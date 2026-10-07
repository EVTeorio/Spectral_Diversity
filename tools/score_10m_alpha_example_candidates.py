from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull, Delaunay

ROOT = Path("reports/tables/figure_sources/10m_alpha_examples")


def tri_area(p):
    return abs(np.cross(p[1] - p[0], p[2] - p[0])) / 2


def circumradius(p):
    a = np.linalg.norm(p[1] - p[0])
    b = np.linalg.norm(p[2] - p[1])
    c = np.linalg.norm(p[0] - p[2])
    area2 = abs(np.cross(p[1] - p[0], p[2] - p[0]))
    if area2 <= 1e-12:
        return np.inf
    return a * b * c / (2 * area2)


def poly_area(p):
    x = p[:, 0]
    y = p[:, 1]
    return abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))) / 2


def score_file(path):
    pts = pd.read_csv(path)[["PC1", "PC2"]].to_numpy(float)
    tri = Delaunay(pts)
    radii = np.array([circumradius(pts[s]) for s in tri.simplices])
    areas = np.array([tri_area(pts[s]) for s in tri.simplices])
    finite = radii[np.isfinite(radii)]
    convex_area = poly_area(pts[ConvexHull(pts).vertices])

    best = None
    for q in np.linspace(0.30, 0.98, 69):
        threshold = np.quantile(finite, q)
        included = np.isfinite(radii) & (radii <= threshold)
        if not included.any():
            continue
        used = np.unique(tri.simplices[included].ravel())
        coverage = len(used) / len(pts)
        ratio = areas[included].sum() / convex_area
        if coverage >= 0.98 and ratio < 0.82:
            best = (ratio, coverage, q, included.sum(), len(pts))
            break
    return best


def main():
    rows = []
    for path in ROOT.glob("*.csv"):
        if path.name == "candidate_summary.csv":
            continue
        try:
            best = score_file(path)
        except Exception:
            continue
        if best:
            ratio, coverage, q, n_tri, n_pts = best
            rows.append((ratio, coverage, q, path.stem, n_tri, n_pts))
    for row in sorted(rows)[:30]:
        print(row)
    print("nresults", len(rows))


if __name__ == "__main__":
    main()
