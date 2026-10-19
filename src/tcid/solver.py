"""Steady anisotropic heat conduction with finite particle-matrix interface resistance.

Governing problem (no heat generation in the adhesive)::

    q = -K grad T,      div(K grad T) = 0

Material interfaces carry an area-specific resistance R (m^2 K/W): the normal flux
is continuous and the temperature jumps by ``[T] = R q_n``. This is represented by
duplicating the degrees of freedom on the filler side of every particle-matrix
interface and adding the interface-conductance term to the weak form::

    a(T, v) = sum_regions int K grad T . grad v dA + int_Gamma (1/R) [T][v] ds

Perfect contact (R = 0) is handled by tying the degrees of freedom (no duplication),
never by an arbitrarily large conductance.

Bulk stiffness is assembled with scikit-fem P1 elements; the interface and Robin
contact couplings are assembled directly (exact 1D P1 edge mass matrices).

Boundary problem: isothermal hot (T = 1 K) and cold (T = 0) faces, insulated
sidewalls. Optionally each face may instead couple to an isothermal reservoir
through an outer contact resistance (Robin condition).
"""

from __future__ import annotations

import time
from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from skfem import Basis, BilinearForm, ElementTriP1, MeshTri, asm

from .materials import rotated_tensor
from .mesh import FILLER, MATRIX, Mesh2D

SOLVER_VERSION = "solver-1.0"


@BilinearForm
def _conduction(u, v, w):
    gu, gv = u.grad, v.grad
    return (w["kxx"] * gu[0] + w["kxy"] * gu[1]) * gv[0] + (w["kxy"] * gu[0] + w["kyy"] * gu[1]) * gv[1]


@dataclass
class MaterialParams:
    k_matrix: float
    k_filler_par: float
    k_filler_perp: float
    r_interface: float  # m^2 K/W; 0 = perfect contact
    r_hot: float = 0.0  # outer contact resistance at the hot face (0 = isothermal face)
    r_cold: float = 0.0


@dataclass
class Discretization:
    mesh: Mesh2D
    split_interfaces: bool
    dof_of_elem: np.ndarray  # (m, 3)
    dof_node: np.ndarray  # (n_dof,) original mesh node of each dof
    iface_dofs: np.ndarray  # (k, 4): matrix-side (a, b), filler-side (a, b)
    iface_len: np.ndarray  # (k,)
    face_dofs: dict  # "hot"/"cold" -> (k, 2) dofs of boundary edges
    face_len: dict
    basis: Basis
    iface_unit: sp.csr_matrix  # interface coupling for unit conductance

    @property
    def n_dof(self) -> int:
        return int(self.dof_node.size)

    @property
    def interface_length(self) -> float:
        return float(self.iface_len.sum())


def _edge_mass(dofs: np.ndarray, lengths: np.ndarray, n: int) -> sp.csr_matrix:
    """Assemble int phi_i phi_j ds over P1 edges."""
    m = lengths[:, None, None] / 6.0 * np.array([[2.0, 1.0], [1.0, 2.0]])
    rows = np.repeat(dofs, 2, axis=1)
    cols = np.tile(dofs, (1, 2))
    return sp.coo_matrix((m.reshape(-1), (rows.reshape(-1), cols.reshape(-1))), shape=(n, n)).tocsr()


def discretize(mesh: Mesh2D, split_interfaces: bool) -> Discretization:
    """Build the (possibly broken) P1 space over ``mesh``."""
    tri = mesh.triangles
    n_nodes = mesh.points.shape[0]
    kind = mesh.region_kind[mesh.region]

    # all element edges with the owning element
    e_local = np.array([[0, 1], [1, 2], [2, 0]])
    edges = np.sort(tri[:, e_local].reshape(-1, 2), axis=1)
    owner = np.repeat(np.arange(tri.shape[0]), 3)
    key = edges[:, 0].astype(np.int64) * n_nodes + edges[:, 1]
    order = np.argsort(key, kind="stable")
    ks = key[order]
    same = ks[1:] == ks[:-1]
    first, second = order[:-1][same], order[1:][same]  # interior edges: two owners
    e1, e2 = owner[first], owner[second]
    diff = kind[e1] != kind[e2]
    if_edges = edges[first[diff]]
    if_m = np.where(kind[e1[diff]] == MATRIX, e1[diff], e2[diff])
    if_f = np.where(kind[e1[diff]] == MATRIX, e2[diff], e1[diff])

    dof_of_elem = tri.copy()
    dof_node = np.arange(n_nodes)
    if split_interfaces and if_edges.size:
        iface_nodes = np.zeros(n_nodes, bool)
        iface_nodes[if_edges.ravel()] = True
        fil = kind == FILLER
        mask = fil[:, None] & iface_nodes[tri]
        el_idx, loc = np.nonzero(mask)
        pair = tri[el_idx, loc].astype(np.int64) * (mesh.region.max() + 1) + mesh.region[el_idx]
        uniq, inv = np.unique(pair, return_inverse=True)
        new = n_nodes + inv
        dof_of_elem[el_idx, loc] = new
        dof_node = np.concatenate([dof_node, (uniq // (mesh.region.max() + 1)).astype(int)])

    n_dof = dof_node.size

    def edge_dofs(el: np.ndarray, nodes: np.ndarray) -> np.ndarray:
        # map the two nodes of each edge to the dofs used by element ``el``
        out = np.empty_like(nodes)
        for j in range(2):
            match = tri[el] == nodes[:, j : j + 1]
            out[:, j] = dof_of_elem[el][match]
        return out

    pts = mesh.points
    if split_interfaces and if_edges.size:
        dm = edge_dofs(if_m, if_edges)
        df = edge_dofs(if_f, if_edges)
        iface_dofs = np.hstack([dm, df])
        iface_len = np.linalg.norm(pts[if_edges[:, 0]] - pts[if_edges[:, 1]], axis=1)
        jump = iface_dofs
        # jump [T] = T_m - T_f gives the block [[M, -M], [-M, M]] with the edge mass M
        m = iface_len[:, None, None] / 6.0 * np.array([[2.0, 1.0], [1.0, 2.0]])
        mm = np.zeros((iface_len.size, 4, 4))
        mm[:, :2, :2] = m
        mm[:, 2:, 2:] = m
        mm[:, :2, 2:] = -m
        mm[:, 2:, :2] = -m
        rows = np.repeat(jump, 4, axis=1)
        cols = np.tile(jump, (1, 4))
        iface_unit = sp.coo_matrix(
            (mm.reshape(-1), (rows.reshape(-1), cols.reshape(-1))), shape=(n_dof, n_dof)
        ).tocsr()
    else:
        iface_dofs = np.zeros((0, 4), int)
        iface_len = np.linalg.norm(pts[if_edges[:, 0]] - pts[if_edges[:, 1]], axis=1) if if_edges.size else np.zeros(0)
        iface_unit = sp.csr_matrix((n_dof, n_dof))

    # boundary edges -> owning element -> dofs
    face_dofs, face_len = {}, {}
    lookup = dict(zip(ks.tolist(), owner[order].tolist()))
    for tag in ("hot", "cold"):
        be = np.sort(mesh.boundary[tag], axis=1)
        bkey = be[:, 0].astype(np.int64) * n_nodes + be[:, 1]
        el = np.array([lookup[k] for k in bkey.tolist()], int)
        face_dofs[tag] = edge_dofs(el, be)
        face_len[tag] = np.linalg.norm(pts[be[:, 0]] - pts[be[:, 1]], axis=1)

    smesh = MeshTri(pts[dof_node].T.copy(), dof_of_elem.T.copy(), sort_t=False)
    basis = Basis(smesh, ElementTriP1())
    return Discretization(
        mesh=mesh,
        split_interfaces=bool(split_interfaces and if_edges.size),
        dof_of_elem=dof_of_elem,
        dof_node=dof_node,
        iface_dofs=iface_dofs,
        iface_len=iface_len,
        face_dofs=face_dofs,
        face_len=face_len,
        basis=basis,
        iface_unit=iface_unit,
    )


def element_tensors(mesh: Mesh2D, mp: MaterialParams) -> np.ndarray:
    """(m, 2, 2) conductivity tensor of each element in global coordinates."""
    kind = mesh.region_kind[mesh.region]
    theta = mesh.region_theta[mesh.region]
    k = np.zeros((mesh.n_elements, 2, 2))
    k[kind == MATRIX] = np.eye(2) * mp.k_matrix
    f = kind == FILLER
    k[f] = rotated_tensor(mp.k_filler_par, mp.k_filler_perp, theta[f])
    return k


@dataclass
class Solution:
    T: np.ndarray  # dof values
    q_hot: float  # heat entering through the hot face, W per metre depth
    q_cold: float  # heat leaving through the cold face
    q_volume: float  # (1/L) int q . e dA: independent flux estimate
    delta_t_faces: float  # mean hot-face minus mean cold-face temperature (1 for isothermal faces)
    r_layer: float  # area-specific layer resistance, m^2 K/W
    r_assembly: float  # including outer contact resistances (series with Robin faces)
    residual: float
    heat_balance_error: float
    flux_consistency_error: float
    n_dof: int
    solve_time_s: float

    def summary(self) -> dict:
        return {k: v for k, v in self.__dict__.items() if k != "T"}


def solve(disc: Discretization, mp: MaterialParams, t_hot: float = 1.0, t_cold: float = 0.0) -> Solution:
    t0 = time.perf_counter()
    mesh = disc.mesh
    if mp.r_interface < 0:
        raise ValueError("interface resistance must be non-negative")
    if mp.r_interface > 0 and not disc.split_interfaces and disc.iface_len.size:
        raise ValueError("finite interface resistance requires a split discretization")
    if mp.r_interface == 0 and disc.split_interfaces:
        raise ValueError("perfect contact must use a tied (unsplit) discretization")

    kt = element_tensors(mesh, mp)
    nqp = disc.basis.X.shape[1]
    w = {c: np.repeat(kt[:, i, j][:, None], nqp, axis=1) for c, (i, j) in {"kxx": (0, 0), "kxy": (0, 1), "kyy": (1, 1)}.items()}
    A = asm(_conduction, disc.basis, **w).tocsr()
    if disc.split_interfaces:
        A = A + disc.iface_unit * (1.0 / mp.r_interface)
    n = disc.n_dof
    b = np.zeros(n)

    robin = {"hot": mp.r_hot > 0, "cold": mp.r_cold > 0}
    res_t = {"hot": t_hot, "cold": t_cold}
    res_r = {"hot": mp.r_hot, "cold": mp.r_cold}
    A_full = A.copy()  # conduction only (for reactions)
    for tag in ("hot", "cold"):
        if robin[tag]:
            h = 1.0 / res_r[tag]
            A = A + h * _edge_mass(disc.face_dofs[tag], disc.face_len[tag], n)
            np.add.at(b, disc.face_dofs[tag].ravel(), np.repeat(h * res_t[tag] * disc.face_len[tag] / 2.0, 2))

    fixed = np.zeros(n, bool)
    u = np.zeros(n)
    for tag in ("hot", "cold"):
        if not robin[tag]:
            d = np.unique(disc.face_dofs[tag])
            fixed[d] = True
            u[d] = res_t[tag]
    free = ~fixed
    A = A.tocsr()
    rhs = b[free] - A[free][:, fixed] @ u[fixed]
    Aff = A[free][:, free].tocsc()
    u[free] = spla.spsolve(Aff, rhs, permc_spec="COLAMD")
    res = float(np.linalg.norm(Aff @ u[free] - rhs) / max(np.linalg.norm(rhs), 1e-300))

    # heat through each face
    q = {}
    for tag, sign in (("hot", 1.0), ("cold", -1.0)):
        if robin[tag]:
            fd, fl = disc.face_dofs[tag], disc.face_len[tag]
            t_mean_edge = u[fd].mean(axis=1)
            q[tag] = sign * float(np.sum((res_t[tag] - t_mean_edge) / res_r[tag] * fl))
        else:
            d = np.unique(disc.face_dofs[tag])
            r = A_full @ u  # reaction = net conduction out of the dof
            q[tag] = sign * float(r[d].sum())

    # face-average temperatures
    def face_mean(tag):
        fd, fl = disc.face_dofs[tag], disc.face_len[tag]
        return float(np.sum(u[fd].mean(axis=1) * fl) / fl.sum())

    dT = face_mean("hot") - face_mean("cold")

    # independent flux estimate via the volume average of q . e
    qv = volume_flux(disc, u, kt)

    W = mesh.face_length
    q_mean = 0.5 * (q["hot"] + q["cold"])
    r_layer = dT * W / q_mean
    r_assembly = (t_hot - t_cold) * W / q_mean
    return Solution(
        T=u,
        q_hot=q["hot"],
        q_cold=q["cold"],
        q_volume=qv,
        delta_t_faces=dT,
        r_layer=float(r_layer),
        r_assembly=float(r_assembly),
        residual=res,
        heat_balance_error=float(abs(q["hot"] - q["cold"]) / abs(q_mean)),
        flux_consistency_error=float(abs(qv - q_mean) / abs(q_mean)),
        n_dof=n,
        solve_time_s=time.perf_counter() - t0,
    )


def element_gradients(disc: Discretization, u: np.ndarray) -> np.ndarray:
    p = disc.basis.mesh.p.T
    t = disc.dof_of_elem
    x = p[t]
    v = u[t]
    d1 = x[:, 1] - x[:, 0]
    d2 = x[:, 2] - x[:, 0]
    det = d1[:, 0] * d2[:, 1] - d1[:, 1] * d2[:, 0]
    du1 = v[:, 1] - v[:, 0]
    du2 = v[:, 2] - v[:, 0]
    gx = (du1 * d2[:, 1] - du2 * d1[:, 1]) / det
    gy = (du2 * d1[:, 0] - du1 * d2[:, 0]) / det
    return np.column_stack([gx, gy]), 0.5 * np.abs(det)


def element_flux(disc: Discretization, u: np.ndarray, kt: np.ndarray) -> np.ndarray:
    g, _ = element_gradients(disc, u)
    return -np.einsum("eij,ej->ei", kt, g)


def volume_flux(disc: Discretization, u: np.ndarray, kt: np.ndarray) -> float:
    """(1/L) int q . e dA. Equals the face heat flow for any div-free field with insulated sides."""
    g, area = element_gradients(disc, u)
    q = -np.einsum("eij,ej->ei", kt, g)
    e = disc.mesh.flow_dir
    return float(np.sum((q @ e) * area) / disc.mesh.flow_length)


class LayerModel:
    """A meshed layer with cached discretizations, solved for many material scenarios."""

    def __init__(self, mesh: Mesh2D):
        self.mesh = mesh
        self._disc: dict[bool, Discretization] = {}

    def disc(self, split: bool) -> Discretization:
        if split not in self._disc:
            self._disc[split] = discretize(self.mesh, split)
        return self._disc[split]

    def solve(self, mp: MaterialParams, **kw) -> Solution:
        return solve(self.disc(mp.r_interface > 0), mp, **kw)
