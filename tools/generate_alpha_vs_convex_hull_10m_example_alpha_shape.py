from pathlib import Path
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.path import Path as MplPath
from matplotlib.patches import Polygon, Rectangle
from scipy.spatial import ConvexHull, Delaunay


PROJECT_DIR = Path(
    r"C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity"
)
SCORE_CSV = PROJECT_DIR / "reports/tables/figure_sources/20_alpha_hull_vs_convex_hull_10m_example_scores.csv"
OUT_PATH = PROJECT_DIR / "Documents/Tables and Figures/20_alpha_hull_vs_convex_hull_10m_example.png"


def triangle_circumradius(tri_points):
    a = np.linalg.norm(tri_points[1] - tri_points[0])
    b = np.linalg.norm(tri_points[2] - tri_points[1])
    c = np.linalg.norm(tri_points[0] - tri_points[2])
    area2 = abs(np.cross(tri_points[1] - tri_points[0], tri_points[2] - tri_points[0]))
    if area2 <= 1e-12:
        return np.inf
    return (a * b * c) / (2.0 * area2)


def boundary_loops(boundary_edges):
    adjacency = defaultdict(list)
    for a, b in boundary_edges:
        adjacency[a].append(b)
        adjacency[b].append(a)

    loops = []
    unused = {tuple(sorted(edge)) for edge in boundary_edges}
    while unused:
        start_edge = unused.pop()
        start, current = start_edge
        loop = [start, current]
        previous = start

        while current != start:
            candidates = [node for node in adjacency[current] if node != previous]
            if not candidates:
                break
            next_node = candidates[0]
            edge_key = tuple(sorted((current, next_node)))
            if edge_key in unused:
                unused.remove(edge_key)
            previous, current = current, next_node
            if current != start:
                loop.append(current)
            if len(loop) > len(boundary_edges) + 5:
                break

        if len(loop) >= 3:
            loops.append(loop)
    return loops


def polygon_area(poly):
    x = poly[:, 0]
    y = poly[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def alpha_shape(points):
    tri = Delaunay(points)
    simplices = tri.simplices
    radii = np.array([triangle_circumradius(points[s]) for s in simplices])
    finite = radii[np.isfinite(radii)]

    best_poly = None
    best_area = np.inf
    convex_area = polygon_area(points[ConvexHull(points).vertices])

    for q in np.linspace(0.58, 0.99, 84):
        threshold = np.quantile(finite, q)
        edge_counts = defaultdict(int)
        for simplex, radius in zip(simplices, radii):
            if radius <= threshold:
                edges = (
                    tuple(sorted((simplex[0], simplex[1]))),
                    tuple(sorted((simplex[1], simplex[2]))),
                    tuple(sorted((simplex[2], simplex[0]))),
                )
                for edge in edges:
                    edge_counts[edge] += 1

        boundary_edges = [edge for edge, count in edge_counts.items() if count == 1]
        loops = boundary_loops(boundary_edges)
        if not loops:
            continue

        for loop in loops:
            poly = points[np.array(loop)]
            area = polygon_area(poly)
            if area <= 0:
                continue
            inside = MplPath(poly, closed=True).contains_points(points, radius=1e-9)
            on_boundary = np.zeros(len(points), dtype=bool)
            on_boundary[np.array(loop)] = True
            if np.all(inside | on_boundary):
                # Prefer the tightest enclosing loop, but avoid a visually indistinguishable convex hull.
                if area < best_area and area < convex_area * 0.98:
                    best_poly = poly
                    best_area = area

    if best_poly is None:
        best_poly = points[ConvexHull(points).vertices]
    return best_poly


def setup_panel(ax, points, extra_right=0.0):
    xmin, ymin = points.min(axis=0)
    xmax, ymax = points.max(axis=0)
    dx = xmax - xmin
    dy = ymax - ymin
    ax.set_xlim(xmin - 0.18 * dx, xmax + (0.18 + extra_right) * dx)
    ax.set_ylim(ymin - 0.18 * dy, ymax + 0.18 * dy)
    ax.axis("off")
    axis_color = "#7a858d"
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    ax.plot([x0 + 0.02 * (x1 - x0), x1 - 0.04 * (x1 - x0)], [y0 + 0.08 * (y1 - y0)] * 2, color=axis_color, lw=1.7)
    ax.plot([x0 + 0.04 * (x1 - x0)] * 2, [y0 + 0.06 * (y1 - y0), y1 - 0.05 * (y1 - y0)], color=axis_color, lw=1.7)
    ax.text(0.52, -0.08, "PC1", transform=ax.transAxes, ha="center", va="top", color=axis_color, fontsize=11)
    ax.text(-0.08, 0.52, "PC2", transform=ax.transAxes, ha="center", va="center", rotation=90, color=axis_color, fontsize=11)


def draw_points(ax, points):
    ax.scatter(points[:, 0], points[:, 1], s=7, facecolor="#244b61", edgecolor="white", linewidth=0.28, alpha=0.72, zorder=5)


def draw_legend(ax):
    handles = [
        plt.Line2D([0], [0], marker="o", color="none", markerfacecolor="#244b61", markeredgecolor="white", markersize=5, label="Pixel scores"),
        Rectangle((0, 0), 1, 1, facecolor="#cfd6d9", edgecolor="#4f5c63", label="Convex-hull area"),
        Rectangle((0, 0), 1, 1, facecolor="#5fa7b0", edgecolor="#075866", label="Alpha-shape area"),
    ]
    ax.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.98, 0.86), frameon=False, fontsize=10, handlelength=1.4)


def main():
    df = pd.read_csv(SCORE_CSV)
    points = df[["PC1", "PC2"]].to_numpy(float)
    convex = points[ConvexHull(points).vertices]
    alpha = alpha_shape(points)

    plt.rcParams.update({"font.family": "Arial", "savefig.transparent": True})
    fig, axes = plt.subplots(1, 2, figsize=(12, 7), dpi=200)
    fig.patch.set_alpha(0)

    for ax in axes:
        ax.set_facecolor((1, 1, 1, 0))

    setup_panel(axes[0], points)
    axes[0].add_patch(Polygon(convex, closed=True, facecolor="#cfd6d9", edgecolor="#4f5c63", linewidth=1.8, alpha=0.86, zorder=1))
    axes[0].add_patch(Polygon(alpha, closed=True, facecolor="#5fa7b0", edgecolor="#075866", linewidth=1.4, alpha=0.74, zorder=2))
    draw_points(axes[0], points)
    axes[0].set_title("A. Convex hull\nOuter envelope around 10 m pixel scores", fontsize=13, fontweight="bold", color="#111820", pad=16)

    setup_panel(axes[1], points, extra_right=0.45)
    axes[1].add_patch(Polygon(alpha, closed=True, facecolor="#5fa7b0", edgecolor="#075866", linewidth=1.8, alpha=0.78, zorder=2))
    draw_points(axes[1], points)
    draw_legend(axes[1])
    axes[1].set_title("B. Alpha-shape area\nBoundary fit to the occupied PC1-PC2 space", fontsize=13, fontweight="bold", color="#111820", pad=16)

    fig.suptitle("Observed 10 m quadrat example in vector-normalized PCA spectral space", x=0.02, y=0.97, ha="left", fontsize=15, fontweight="bold", color="#111820")
    fig.text(0.025, 0.91, "Example tile: 700_a; retained illuminated pixels: 4,745; plotted sample: 1,200 pixels", fontsize=10.5, color="#5d6870")
    fig.subplots_adjust(left=0.06, right=0.97, top=0.80, bottom=0.11, wspace=0.18)

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PATH)
    print(f"Created {OUT_PATH}")


if __name__ == "__main__":
    main()
