"""Conforming triangular meshes with tagged material regions and interfaces (Gmsh).

Particles are meshed as polygons whose vertices lie on the ellipse boundary.
The polygon is the geometry that is actually simulated: its area defines the
recorded loading and its separations define the recorded minimum gap. Because
every particle boundary curve is shared by the particle surface and the matrix
surface, the mesh is conforming along each particle-matrix interface; the solver
then duplicates interface degrees of freedom to represent interface resistance.

No rasterisation is involved, so no artificial perfect contacts can appear.
"""

from __future__ import annotations

import math
import time
from dataclasses import dataclass, field

import gmsh
import numpy as np
import shapely

from .config import MeshSettings
from .geometry import Geometry, ellipse_perimeter, ellipse_vertices

MATRIX, FILLER = 0, 1
MESHER_VERSION = "mesh-1.0"


@dataclass
class Mesh2D:
    points: np.ndarray  # (n, 2)
    triangles: np.ndarray  # (m, 3)
    region: np.ndarray  # (m,) region id per element
    region_kind: np.ndarray  # (n_regions,) MATRIX or FILLER
    region_theta: np.ndarray  # (n_regions,) orientation of the material frame (rad)
    boundary: dict  # tag -> (k, 2) node pairs; tags "hot", "cold", "side"
    flow_dir: np.ndarray  # unit vector from hot face to cold face
    flow_length: float  # distance between hot and cold faces (layer thickness)
    face_length: float  # length of the hot face (layer width)
    info: dict = field(default_factory=dict)

    @property
    def n_elements(self) -> int:
        return int(self.triangles.shape[0])

    def element_areas(self) -> np.ndarray:
        p = self.points[self.triangles]
        return 0.5 * np.abs(
            (p[:, 1, 0] - p[:, 0, 0]) * (p[:, 2, 1] - p[:, 0, 1])
            - (p[:, 2, 0] - p[:, 0, 0]) * (p[:, 1, 1] - p[:, 0, 1])
        )

    def filler_fraction(self) -> float:
        ar = self.element_areas()
        return float(ar[self.region_kind[self.region] == FILLER].sum() / ar.sum())

    def rotated(self, angle: float, origin=(0.0, 0.0)) -> "Mesh2D":
        """Rigidly rotate the mesh, its material frames and its boundary problem."""
        c, s = math.cos(angle), math.sin(angle)
        q = np.array([[c, -s], [s, c]])
        o = np.asarray(origin)
        pts = (self.points - o) @ q.T + o
        return Mesh2D(
            points=pts,
            triangles=self.triangles.copy(),
            region=self.region.copy(),
            region_kind=self.region_kind.copy(),
            region_theta=self.region_theta + angle,
            boundary={k: v.copy() for k, v in self.boundary.items()},
            flow_dir=q @ self.flow_dir,
            flow_length=self.flow_length,
            face_length=self.face_length,
            info={**self.info, "rotated_by": float(angle)},
        )

    def to_npz_dict(self) -> dict:
        d = {
            "points": self.points,
            "triangles": self.triangles,
            "region": self.region,
            "region_kind": self.region_kind,
            "region_theta": self.region_theta,
            "flow_dir": self.flow_dir,
            "flow_length": np.array(self.flow_length),
            "face_length": np.array(self.face_length),
        }
        for k, v in self.boundary.items():
            d[f"boundary_{k}"] = v
        return d

    @classmethod
    def from_npz(cls, z) -> "Mesh2D":
        boundary = {k[len("boundary_") :]: z[k] for k in z.files if k.startswith("boundary_")}
        return cls(
            points=z["points"],
            triangles=z["triangles"],
            region=z["region"],
            region_kind=z["region_kind"],
            region_theta=z["region_theta"],
            boundary=boundary,
            flow_dir=z["flow_dir"],
            flow_length=float(z["flow_length"]),
            face_length=float(z["face_length"]),
        )


class MeshError(RuntimeError):
    pass


# ----------------------------------------------------------------------------- gmsh helpers


def _init():
    if not gmsh.isInitialized():
        gmsh.initialize(interruptible=False)
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setNumber("General.Verbosity", 1)
    gmsh.option.setNumber("General.NumThreads", 1)
    gmsh.clear()
    gmsh.model.add("layer")


def _common_options(h_max: float):
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 1)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 1)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeMax", h_max)
    gmsh.option.setNumber("Mesh.Algorithm", 6)  # Frontal-Delaunay
    gmsh.option.setNumber("Mesh.ElementOrder", 1)
    gmsh.option.setNumber("Mesh.RandomSeed", 1)


def _polyline_loop(pts: np.ndarray, size):
    size = np.broadcast_to(np.asarray(size, float), (len(pts),))
    tags = [gmsh.model.geo.addPoint(float(x), float(y), 0.0, float(h)) for (x, y), h in zip(pts, size)]
    lines = [gmsh.model.geo.addLine(tags[i], tags[(i + 1) % len(tags)]) for i in range(len(tags))]
    return gmsh.model.geo.addCurveLoop(lines), lines


def _extract(surfaces: dict[int, int], lines_by_tag: dict[str, list[int]]):
    """Nodes, triangles with region ids, and tagged boundary edges from the current gmsh model."""
    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    order = np.argsort(node_tags)
    node_tags = node_tags[order]
    xy = coords.reshape(-1, 3)[order, :2]
    index = lambda t: np.searchsorted(node_tags, t)  # noqa: E731

    tris, regs, etags = [], [], []
    for surf, reg in surfaces.items():
        types, el_tags, el_nodes = gmsh.model.mesh.getElements(2, surf)
        for ty, et, en in zip(types, el_tags, el_nodes):
            if ty != 2:
                raise MeshError(f"unexpected element type {ty}")
            tris.append(index(en.reshape(-1, 3)))
            regs.append(np.full(len(et), reg))
            etags.append(et)
    tri = np.concatenate(tris)
    reg = np.concatenate(regs)
    etag = np.concatenate(etags)
    q = np.asarray(gmsh.model.mesh.getElementQualities(etag, "minSICN"))

    boundary = {}
    for tag, lines in lines_by_tag.items():
        edges = []
        for ln in lines:
            types, _, el_nodes = gmsh.model.mesh.getElements(1, ln)
            for ty, en in zip(types, el_nodes):
                edges.append(index(en.reshape(-1, 2)))
        boundary[tag] = np.concatenate(edges) if edges else np.zeros((0, 2), int)

    # drop nodes not used by any triangle (none expected), keep indices compact
    used = np.unique(tri)
    if used.size != xy.shape[0]:
        remap = -np.ones(xy.shape[0], int)
        remap[used] = np.arange(used.size)
        xy = xy[used]
        tri = remap[tri]
        boundary = {k: remap[v] for k, v in boundary.items()}
    # consistent counter-clockwise orientation
    p = xy[tri]
    det = (p[:, 1, 0] - p[:, 0, 0]) * (p[:, 2, 1] - p[:, 0, 1]) - (p[:, 2, 0] - p[:, 0, 0]) * (p[:, 1, 1] - p[:, 0, 1])
    neg = det < 0
    tri[neg] = tri[neg][:, [0, 2, 1]]
    return xy, tri, reg, q, boundary


# ----------------------------------------------------------------------------- particle layers


def _limit_gradation(h: np.ndarray, s: np.ndarray, total: float, grade: float) -> np.ndarray:
    """Limit h growth along a closed curve: h_k <= h_j + grade * dist(s_k, s_j)."""
    n = h.size
    h = h.copy()
    for _ in range(2):
        for k in range(1, 2 * n):
            i, j = k % n, (k - 1) % n
            ds = (s[i] - s[j]) % total
            h[i] = min(h[i], h[j] + grade * ds)
        for k in range(2 * n - 2, -1, -1):
            i, j = k % n, (k + 1) % n
            ds = (s[j] - s[i]) % total
            h[i] = min(h[i], h[j] + grade * ds)
    return h


def particle_boundaries(g: Geometry, ms: MeshSettings, lx: float, ly: float, grade: float = 0.25):
    """Adaptive boundary vertices of every particle.

    The local element size along a particle boundary is limited by
      * the local distance to the nearest other particle (``n_across_gap`` elements across
        every narrow gap),
      * the local distance to the layer boundary (``face_elements_in_margin`` across the margin),
      * geometric resolution (at least ``n_perimeter`` vertices, so the simulated polygon
        loses < 0.3% of the ellipse area),
    and is graded smoothly along the boundary. Vertices lie on the exact ellipse, so the
    simulated polygon is inscribed and every simulated gap is >= the generated gap.
    Gap and margin resolution are never coarsened by ``size_factor`` > 1.
    """
    sf = ms.size_factor
    fine = min(sf, 1.0)
    n = g.n
    fine_polys = [shapely.Polygon(ellipse_vertices(g.cx[i], g.cy[i], g.a[i], g.b[i], g.theta[i], 96)) for i in range(n)]
    out_pts, out_h = [], []
    for i in range(n):
        a, b = g.a[i], g.b[i]
        perim = float(ellipse_perimeter(a, b))
        h_geo = min(perim / ms.n_perimeter, ms.h_bulk) * fine
        m_ref = max(256, int(8 * perim / h_geo))
        t = np.linspace(0.0, 2 * math.pi, m_ref, endpoint=False)
        P = ellipse_vertices(g.cx[i], g.cy[i], a, b, g.theta[i], m_ref)
        seg = np.linalg.norm(np.roll(P, -1, axis=0) - P, axis=1)
        sarc = np.concatenate([[0.0], np.cumsum(seg)[:-1]])
        total = float(seg.sum())
        cutoff = ms.n_across_gap * h_geo * 1.5
        dist = np.hypot(g.cx - g.cx[i], g.cy - g.cy[i])
        nb = np.nonzero((dist < g.a + a + cutoff) & (np.arange(n) != i))[0]
        d_nb = np.full(m_ref, np.inf)
        if nb.size:
            union = shapely.union_all([fine_polys[j] for j in nb])
            d_nb = shapely.distance(shapely.points(P), union)
        d_face = np.minimum.reduce([P[:, 0], lx - P[:, 0], P[:, 1], ly - P[:, 1]])
        h = np.minimum.reduce(
            [
                d_nb / ms.n_across_gap * fine,
                d_face / ms.face_elements_in_margin * fine,
                np.full(m_ref, h_geo),
            ]
        )
        h = np.maximum(h, ms.h_min * fine)
        h = _limit_gradation(h, sarc, total, grade)
        # place vertices equidistant in the metric ds / h
        phi = np.concatenate([[0.0], np.cumsum(seg / (0.5 * (h + np.roll(h, -1))))])
        nv = max(ms.n_perimeter, int(math.ceil(phi[-1])))
        targets = np.linspace(0.0, phi[-1], nv, endpoint=False)
        s_all = np.concatenate([sarc, [total]])
        sv = np.interp(targets, phi, s_all)
        tv = np.interp(sv, s_all, np.concatenate([t, [2 * math.pi]]))
        c, s_ = math.cos(g.theta[i]), math.sin(g.theta[i])
        x, y = a * np.cos(tv), b * np.sin(tv)
        out_pts.append(np.column_stack([g.cx[i] + c * x - s_ * y, g.cy[i] + s_ * x + c * y]))
        out_h.append(np.interp(sv, s_all, np.concatenate([h, [h[0]]])))
    return out_pts, out_h


def mesh_particles(g: Geometry, ms: MeshSettings, flow_axis: str = "y") -> Mesh2D:
    """Mesh a particle-in-matrix layer. ``flow_axis='x'`` is used for rotated checks."""
    if not g.ok:
        raise MeshError(f"cannot mesh a rejected geometry: {g.reason}")
    t0 = time.perf_counter()
    lay = g.layer
    if flow_axis == "y":
        lx, ly = lay["width"], lay["thickness"]
    else:
        lx, ly = lay["thickness"], lay["width"]
    sf = ms.size_factor
    h_face = min(ms.h_bulk * sf, lay["margin"] / ms.face_elements_in_margin * min(sf, 1.0))
    polys, sizes = particle_boundaries(g, ms, lx, ly)

    _init()
    try:
        _common_options(ms.h_bulk * sf)
        corners = [(0, 0), (lx, 0), (lx, ly), (0, ly)]
        pt = [gmsh.model.geo.addPoint(x, y, 0.0, h_face) for x, y in corners]
        ln = [gmsh.model.geo.addLine(pt[i], pt[(i + 1) % 4]) for i in range(4)]
        outer = gmsh.model.geo.addCurveLoop(ln)
        loops, part_lines = [], []
        for i, poly in enumerate(polys):
            loop, lines = _polyline_loop(poly, sizes[i])
            loops.append(loop)
            part_lines.append(lines)
        matrix_surf = gmsh.model.geo.addPlaneSurface([outer] + loops)
        surfaces = {matrix_surf: 0}
        for i, loop in enumerate(loops):
            surfaces[gmsh.model.geo.addPlaneSurface([loop])] = i + 1
        gmsh.model.geo.synchronize()
        if flow_axis == "y":
            tags = {"hot": [ln[0]], "side": [ln[1], ln[3]], "cold": [ln[2]]}
        else:
            tags = {"hot": [ln[3]], "side": [ln[0], ln[2]], "cold": [ln[1]]}
        gmsh.model.mesh.generate(2)
        xy, tri, reg, q, boundary = _extract(surfaces, tags)
        # nodes on each particle boundary (for resolution checks)
        part_node_sets = []
        node_tags_all, _, _ = gmsh.model.mesh.getNodes()
        sorted_tags = np.sort(node_tags_all)
        for lines in part_lines:
            s = set()
            for lt in lines:
                nt, _, _ = gmsh.model.mesh.getNodes(1, lt, includeBoundary=True)
                s.update(np.searchsorted(sorted_tags, nt).tolist())
            part_node_sets.append(s)
    except Exception as exc:  # gmsh raises plain Exceptions
        raise MeshError(f"gmsh failed: {exc}") from exc
    finally:
        gmsh.clear()

    region_kind = np.array([MATRIX] + [FILLER] * g.n)
    region_theta = np.concatenate([[0.0], g.theta])
    flow_dir = np.array([0.0, 1.0]) if flow_axis == "y" else np.array([1.0, 0.0])
    m = Mesh2D(
        points=xy,
        triangles=tri,
        region=reg,
        region_kind=region_kind,
        region_theta=region_theta,
        boundary=boundary,
        flow_dir=flow_dir,
        flow_length=lay["thickness"],
        face_length=lay["width"],
    )
    m.info = {
        "mesher_version": MESHER_VERSION,
        "n_nodes": int(xy.shape[0]),
        "n_elements": int(tri.shape[0]),
        "h_face": float(h_face),
        "h_particle_min": float(min(h.min() for h in sizes)) if g.n else float("nan"),
        "h_particle_median": float(np.median(np.concatenate(sizes))) if g.n else float("nan"),
        "n_boundary_vertices": int(sum(len(p) for p in polys)),
        "size_factor": float(sf),
        "min_quality": float(q.min()),
        "mesh_time_s": time.perf_counter() - t0,
    }
    m.info.update(resolution_checks(m, polys, part_node_sets, ms))
    return m


def resolution_checks(m: Mesh2D, polys, part_node_sets, ms: MeshSettings) -> dict:
    """Detect unresolved shortcuts and record the simulated geometry.

    * ``bridge_elements``: matrix triangles whose vertices all lie on the boundaries of
      two or more different particles, i.e. a gap spanned by a single element.
    * ``margin_bridges``: matrix triangles touching both a particle and a hot/cold face.
    """
    node_part = -np.ones(m.points.shape[0], int)
    for i, s in enumerate(part_node_sets):
        node_part[list(s)] = i
    face_nodes = np.zeros(m.points.shape[0], bool)
    for tag in ("hot", "cold"):
        face_nodes[m.boundary[tag].ravel()] = True
    mat = m.region_kind[m.region] == MATRIX
    tp = node_part[m.triangles[mat]]
    on_part = tp >= 0
    all_on = on_part.all(axis=1)
    distinct = all_on & (tp.min(axis=1) != tp.max(axis=1))
    bridges = int(distinct.sum())
    tf = face_nodes[m.triangles[mat]]
    margin_bridges = int((on_part.any(axis=1) & tf.any(axis=1)).sum())

    sp = [shapely.Polygon(p) for p in polys]
    area = float(sum(p.area for p in sp))
    domain = m.flow_length * m.face_length
    min_gap = float("inf")
    if len(sp) > 1:
        tree = shapely.STRtree(sp)
        for i, p in enumerate(sp):
            idx = tree.query(p.buffer(10 * ms.h_bulk))
            idx = idx[idx != i]
            if idx.size:
                min_gap = min(min_gap, float(np.min(shapely.distance(p, [sp[j] for j in idx]))))
    ok = bridges == 0 and margin_bridges == 0 and m.info.get("min_quality", 1.0) >= ms.min_quality
    reason = []
    if bridges:
        reason.append(f"{bridges} single-element gap bridges")
    if margin_bridges:
        reason.append(f"{margin_bridges} unresolved margin elements")
    if m.info.get("min_quality", 1.0) < ms.min_quality:
        reason.append(f"element quality {m.info['min_quality']:.3f} < {ms.min_quality}")
    return {
        "bridge_elements": bridges,
        "margin_bridges": margin_bridges,
        "loading_polygon": area / domain,
        "min_gap_meshed": min_gap,
        "mesh_ok": bool(ok),
        "mesh_reason": "; ".join(reason),
    }


# ----------------------------------------------------------------------------- layered verification meshes


def mesh_strips(
    levels: list[float],
    kinds: list[int],
    width: float,
    h: float,
    orientation: str = "series",
    thetas: list[float] | None = None,
) -> Mesh2D:
    """Stack of strips for analytic verification.

    ``orientation='series'``: strips are horizontal (normal to the heat flow) at y levels
    ``levels`` and span the full ``width``; total thickness is ``levels[-1]``.
    ``orientation='parallel'``: strips are vertical at x levels ``levels`` spanning the
    thickness ``width`` (here meaning the through-thickness extent).
    """
    n = len(kinds)
    assert len(levels) == n + 1
    thetas = thetas or [0.0] * n
    _init()
    try:
        _common_options(h)
        if orientation == "series":
            left = [gmsh.model.geo.addPoint(0.0, y, 0.0, h) for y in levels]
            right = [gmsh.model.geo.addPoint(width, y, 0.0, h) for y in levels]
            horiz = [gmsh.model.geo.addLine(left[k], right[k]) for k in range(n + 1)]
            rside = [gmsh.model.geo.addLine(right[k], right[k + 1]) for k in range(n)]
            lside = [gmsh.model.geo.addLine(left[k + 1], left[k]) for k in range(n)]
            surfaces = {}
            for k in range(n):
                loop = gmsh.model.geo.addCurveLoop([horiz[k], rside[k], -horiz[k + 1], lside[k]])
                surfaces[gmsh.model.geo.addPlaneSurface([loop])] = k
            tags = {"hot": [horiz[0]], "cold": [horiz[n]], "side": rside + lside}
            flow_length, face_length, flow_dir = levels[-1], width, np.array([0.0, 1.0])
        else:
            thick = width
            bot = [gmsh.model.geo.addPoint(x, 0.0, 0.0, h) for x in levels]
            top = [gmsh.model.geo.addPoint(x, thick, 0.0, h) for x in levels]
            vert = [gmsh.model.geo.addLine(bot[k], top[k]) for k in range(n + 1)]
            bseg = [gmsh.model.geo.addLine(bot[k], bot[k + 1]) for k in range(n)]
            tseg = [gmsh.model.geo.addLine(top[k + 1], top[k]) for k in range(n)]
            surfaces = {}
            for k in range(n):
                loop = gmsh.model.geo.addCurveLoop([bseg[k], vert[k + 1], tseg[k], -vert[k]])
                surfaces[gmsh.model.geo.addPlaneSurface([loop])] = k
            tags = {"hot": bseg, "cold": tseg, "side": [vert[0], vert[n]]}
            flow_length, face_length, flow_dir = thick, levels[-1], np.array([0.0, 1.0])
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.generate(2)
        xy, tri, reg, q, boundary = _extract(surfaces, tags)
    except Exception as exc:
        raise MeshError(f"gmsh failed: {exc}") from exc
    finally:
        gmsh.clear()
    m = Mesh2D(
        points=xy,
        triangles=tri,
        region=reg,
        region_kind=np.asarray(kinds),
        region_theta=np.asarray(thetas, float),
        boundary=boundary,
        flow_dir=flow_dir,
        flow_length=float(flow_length),
        face_length=float(face_length),
    )
    m.info = {"n_nodes": int(xy.shape[0]), "n_elements": int(tri.shape[0]), "min_quality": float(q.min())}
    return m
