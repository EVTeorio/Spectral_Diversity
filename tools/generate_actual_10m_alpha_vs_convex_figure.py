from collections import defaultdict
from pathlib import Path
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Polygon, Rectangle
from scipy.spatial import ConvexHull, Delaunay

from score_10m_alpha_example_candidates import circumradius


PROJECT_DIR = Path(
    r"C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity"
)
TILE_ID = sys.argv[1] if len(sys.argv) > 1 else "112_a"
OUT_NAME = sys.argv[2] if len(sys.argv) > 2 else "20_alpha_hull_vs_convex_hull_10m_example.png"
ALPHA_QUANTILE = float(sys.argv[3]) if len(sys.argv) > 3 else 0.95
SCORE_CSV = PROJECT_DIR / f"reports/tables/figure_sources/10m_alpha_examples/{TILE_ID}.csv"
OUT_PATH = PROJECT_DIR / "Documents/Tables and Figures" / OUT_NAME

COL = {
    "text": "#111820",
    "subtext": "#5d6870",
    "axis": "#7a858d",
    "point": "#244b61",
    "convex_fill": "#cfd6d9",
    "convex_edge": "#4f5c63",
    "alpha_fill": "#5fa7b0",
    "alpha_edge": "#075866",
}


def alpha_complex(points, q=ALPHA_QUANTILE):
    tri = Delaunay(points)
    radii = np.array([circumradius(points[s]) for s in tri.simplices])
    finite = radii[np.isfinite(radii)]
    threshold = np.quantile(finite, q)
    included = np.isfinite(radii) & (radii <= threshold)
    used = np.unique(tri.simplices[included].ravel())
    return tri.simplices[included], used


def boundary_edges(triangles):
    counts = defaultdict(int)
    for tri in triangles:
        for edge in (
            tuple(sorted((tri[0], tri[1]))),
            tuple(sorted((tri[1], tri[2]))),
            tuple(sorted((tri[2], tri[0]))),
        ):
            counts[edge] += 1
    return [edge for edge, count in counts.items() if count == 1]


def polygon_area(points):
    x = points[:, 0]
    y = points[:, 1]
    return abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))) / 2


def boundary_loops(edges):
    adjacency = defaultdict(list)
    for a, b in edges:
        adjacency[a].append(b)
        adjacency[b].append(a)

    unused = {tuple(sorted(edge)) for edge in edges}
    loops = []
    while unused:
        start, current = unused.pop()
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
            if len(loop) > len(edges) + 5:
                break

        if len(loop) >= 3:
            loops.append(loop)
    return loops


def largest_connected_triangle_component(triangles):
    edge_to_triangles = defaultdict(list)
    for idx, tri in enumerate(triangles):
        for edge in (
            tuple(sorted((tri[0], tri[1]))),
            tuple(sorted((tri[1], tri[2]))),
            tuple(sorted((tri[2], tri[0]))),
        ):
            edge_to_triangles[edge].append(idx)

    neighbors = defaultdict(set)
    for triangle_ids in edge_to_triangles.values():
        if len(triangle_ids) < 2:
            continue
        for i in triangle_ids:
            for j in triangle_ids:
                if i != j:
                    neighbors[i].add(j)

    seen = set()
    components = []
    for start in range(len(triangles)):
        if start in seen:
            continue
        stack = [start]
        seen.add(start)
        component = []
        while stack:
            current = stack.pop()
            component.append(current)
            for nxt in neighbors[current]:
                if nxt not in seen:
                    seen.add(nxt)
                    stack.append(nxt)
        components.append(component)

    return triangles[max(components, key=len)]


def setup_ax(ax, points):
    xmin, ymin = points.min(axis=0)
    xmax, ymax = points.max(axis=0)
    dx = xmax - xmin
    dy = ymax - ymin
    ax.set_xlim(xmin - 0.12 * dx, xmax + 0.14 * dx)
    ax.set_ylim(ymin - 0.14 * dy, ymax + 0.12 * dy)
    ax.axis("off")
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    ax.plot([x0 + 0.02 * (x1 - x0), x1 - 0.04 * (x1 - x0)], [y0 + 0.07 * (y1 - y0)] * 2, color=COL["axis"], lw=1.8)
    ax.plot([x0 + 0.04 * (x1 - x0)] * 2, [y0 + 0.04 * (y1 - y0), y1 - 0.05 * (y1 - y0)], color=COL["axis"], lw=1.8)
    ax.text(0.52, -0.08, "PC1", transform=ax.transAxes, ha="center", va="top", color=COL["axis"], fontsize=11)
    ax.text(-0.07, 0.52, "PC2", transform=ax.transAxes, ha="center", va="center", rotation=90, color=COL["axis"], fontsize=11)


def draw_alpha_complex(ax, points, triangles):
    edges = boundary_edges(triangles)
    loops = boundary_loops(edges)
    if loops:
        outer = max(loops, key=lambda loop: polygon_area(points[np.array(loop)]))
        outer_points = points[np.array(outer)]
        ax.add_patch(
            Polygon(
                outer_points,
                closed=True,
                facecolor=COL["alpha_fill"],
                edgecolor=COL["alpha_edge"],
                linewidth=1.6,
                alpha=0.78,
                zorder=2,
            )
        )


def draw_points(ax, points):
    ax.scatter(points[:, 0], points[:, 1], s=8, facecolor=COL["point"], edgecolor="white", linewidth=0.25, alpha=0.78, zorder=5)


def draw_legend(ax):
    handles = [
        plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=COL["point"], markeredgecolor="white", markersize=5, label="Pixel scores"),
        Rectangle((0, 0), 1, 1, facecolor=COL["convex_fill"], edgecolor=COL["convex_edge"], label="Convex-hull area"),
        Rectangle((0, 0), 1, 1, facecolor=COL["alpha_fill"], edgecolor=COL["alpha_edge"], label="Alpha-hull area"),
    ]
    return handles


def main():
    raw_points = pd.read_csv(SCORE_CSV)[["PC1", "PC2"]].to_numpy(float)
    original_triangles, used = alpha_complex(raw_points)
    points = raw_points

    triangles = original_triangles
    convex = points[ConvexHull(points).vertices]

    plt.rcParams.update({"font.family": "Arial", "savefig.transparent": True})
    fig, axes = plt.subplots(1, 2, figsize=(12, 6.4), dpi=220)
    fig.patch.set_alpha(0)
    for ax in axes:
        ax.set_facecolor((1, 1, 1, 0))

    setup_ax(axes[0], points)
    axes[0].add_patch(Polygon(convex, closed=True, facecolor=COL["convex_fill"], edgecolor=COL["convex_edge"], linewidth=2.0, alpha=0.86, zorder=1))
    draw_points(axes[0], points)
    axes[0].set_title("A. Convex hull\nOuter envelope around actual 10 m pixel scores", fontsize=13, fontweight="bold", color=COL["text"], pad=14)

    setup_ax(axes[1], points)
    draw_alpha_complex(axes[1], points, triangles)
    draw_points(axes[1], points)
    axes[1].set_title("B. Alpha-hull area\nBoundary fit around actual occupied PC1-PC2 space", fontsize=13, fontweight="bold", color=COL["text"], pad=14)

    fig.suptitle("Observed 10 m quadrat example in vector-normalized PCA spectral space", x=0.02, y=0.98, ha="left", fontsize=15, fontweight="bold", color=COL["text"])
    fig.text(
        0.025,
        0.91,
        f"Example tile: {TILE_ID}; plotted source sample: {len(raw_points):,} pixels; continuous alpha boundary contains all plotted points",
        fontsize=10.5,
        color=COL["subtext"],
    )
    handles = draw_legend(axes[1])
    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.52, 0.045),
        ncol=3,
        frameon=False,
        fontsize=10,
        handlelength=1.4,
    )
    fig.subplots_adjust(left=0.06, right=0.97, top=0.78, bottom=0.18, wspace=0.16)
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = OUT_PATH.with_name(OUT_PATH.stem + "_tmp.png")
    fig.savefig(tmp_path)
    os.replace(tmp_path, OUT_PATH)
    print(f"Created {OUT_PATH}")
    print(f"Tile: {TILE_ID}; displayed pixels: {len(points)}; source sample: {len(raw_points)}")


if __name__ == "__main__":
    main()
