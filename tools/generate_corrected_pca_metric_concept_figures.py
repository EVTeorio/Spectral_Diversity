from pathlib import Path
from itertools import combinations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon, Rectangle
from scipy.spatial import ConvexHull


PROJECT_DIR = Path(
    r"C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity"
)
FIG_DIR = PROJECT_DIR / "Documents/Tables and Figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

OUT_20 = FIG_DIR / "20_alpha_hull_vs_convex_hull_10m_example.png"
OUT_21 = FIG_DIR / "21_pca_spectral_metric_concept_diagram.png"

COL = {
    "text": "#111820",
    "subtext": "#5d6870",
    "axis": "#7a858d",
    "point": "#244b61",
    "convex_fill": "#cfd6d9",
    "convex_edge": "#4f5c63",
    "alpha_fill": "#5fa7b0",
    "alpha_edge": "#075866",
    "pair": "#9aa2a8",
    "centroid": "#c91f1f",
}

# Conceptual retained-pixel scores in PCA space. The cloud is shaped so that the
# alpha boundary can follow an indentation while still enclosing every point.
POINTS = np.array(
    [
        [-1.25, 0.48],
        [-1.05, 0.18],
        [-0.92, -0.28],
        [-0.66, -0.55],
        [-0.34, -0.34],
        [-0.12, -0.08],
        [0.20, -0.46],
        [0.58, -0.30],
        [0.90, 0.00],
        [0.62, 0.34],
        [0.26, 0.56],
        [-0.08, 0.30],
        [-0.42, 0.14],
        [-0.70, 0.34],
    ],
    dtype=float,
)

ALPHA_POLY = np.array(
    [
        [-1.34, 0.54],
        [-0.72, 0.44],
        [-0.42, 0.22],
        [0.24, 0.64],
        [0.68, 0.40],
        [1.00, 0.03],
        [0.64, -0.38],
        [0.20, -0.56],
        [-0.13, -0.18],
        [-0.60, -0.66],
        [-1.00, -0.36],
        [-1.16, 0.12],
    ],
    dtype=float,
)


def convex_poly(points=POINTS):
    return points[ConvexHull(points).vertices]


def setup_ax(ax, extra_right=0.0):
    ax.set_xlim(-1.55, 1.20 + extra_right)
    ax.set_ylim(-0.85, 0.82)
    ax.axis("off")
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    ax.plot([x0 + 0.06, x1 - 0.05], [y0 + 0.08, y0 + 0.08], color=COL["axis"], lw=1.8)
    ax.plot([x0 + 0.10, x0 + 0.10], [y0 + 0.04, y1 - 0.06], color=COL["axis"], lw=1.8)
    ax.text(0.52, -0.08, "PC1", transform=ax.transAxes, ha="center", va="top", color=COL["axis"], fontsize=11)
    ax.text(-0.07, 0.52, "PC2", transform=ax.transAxes, ha="center", va="center", rotation=90, color=COL["axis"], fontsize=11)


def draw_points(ax, points=POINTS, size=42):
    ax.scatter(points[:, 0], points[:, 1], s=size, facecolor=COL["point"], edgecolor="white", linewidth=0.9, alpha=0.92, zorder=6)


def draw_legend(ax, alpha_label="Alpha-hull area"):
    handles = [
        plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=COL["point"], markeredgecolor="white", markersize=7, label="Pixel scores"),
        Rectangle((0, 0), 1, 1, facecolor=COL["convex_fill"], edgecolor=COL["convex_edge"], label="Convex-hull area"),
        Rectangle((0, 0), 1, 1, facecolor=COL["alpha_fill"], edgecolor=COL["alpha_edge"], label=alpha_label),
    ]
    ax.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.98, 0.93), frameon=False, fontsize=10, handlelength=1.4)


def save_alpha_vs_convex():
    fig, axes = plt.subplots(1, 2, figsize=(12, 6.4), dpi=220)
    fig.patch.set_alpha(0)
    for ax in axes:
        ax.set_facecolor((1, 1, 1, 0))
        setup_ax(ax)

    conv = convex_poly()

    axes[0].add_patch(Polygon(conv, closed=True, facecolor=COL["convex_fill"], edgecolor=COL["convex_edge"], linewidth=2.0, alpha=0.88, zorder=1))
    axes[0].add_patch(Polygon(ALPHA_POLY, closed=True, facecolor=COL["alpha_fill"], edgecolor=COL["alpha_edge"], linewidth=1.6, alpha=0.75, zorder=2))
    draw_points(axes[0])
    axes[0].set_title("A. Convex hull\nOuter envelope around the same pixel scores", fontsize=13, fontweight="bold", color=COL["text"], pad=14)

    axes[1].cla()
    setup_ax(axes[1], extra_right=0.65)
    axes[1].add_patch(Polygon(ALPHA_POLY, closed=True, facecolor=COL["alpha_fill"], edgecolor=COL["alpha_edge"], linewidth=2.0, alpha=0.80, zorder=2))
    draw_points(axes[1])
    draw_legend(axes[1])
    axes[1].set_title("B. Alpha-hull area\nBoundary follows occupied PC1-PC2 space", fontsize=13, fontweight="bold", color=COL["text"], pad=14)

    fig.suptitle("Conceptual 10 m quadrat example in vector-normalized PCA spectral space", x=0.02, y=0.98, ha="left", fontsize=15, fontweight="bold", color=COL["text"])
    fig.text(0.025, 0.91, "The same pixel scores are shown in both panels; the alpha hull follows the occupied spectral region more tightly than the convex hull.", fontsize=10.5, color=COL["subtext"])
    fig.subplots_adjust(left=0.06, right=0.97, top=0.78, bottom=0.12, wspace=0.16)
    fig.savefig(OUT_20, transparent=True)
    plt.close(fig)


def save_three_metric_concept():
    fig, axes = plt.subplots(1, 3, figsize=(14.2, 5.7), dpi=220)
    fig.patch.set_alpha(0)
    for ax in axes:
        ax.set_facecolor((1, 1, 1, 0))
        setup_ax(ax)

    # Rao's Q: pairwise distances among all pixel scores.
    for i, j in combinations(range(len(POINTS)), 2):
        axes[0].plot([POINTS[i, 0], POINTS[j, 0]], [POINTS[i, 1], POINTS[j, 1]], color=COL["pair"], lw=0.75, alpha=0.65, zorder=1)
    draw_points(axes[0])
    axes[0].set_title("A. Spectral Rao's Q\nPairwise dissimilarity among pixels", fontsize=12, fontweight="bold", color=COL["text"], pad=12)

    # Mean Euclidean distance: distance from points to centroid.
    centroid = POINTS.mean(axis=0)
    for p in POINTS:
        axes[1].plot([centroid[0], p[0]], [centroid[1], p[1]], color=COL["centroid"], lw=1.0, alpha=0.55, zorder=1)
    draw_points(axes[1])
    axes[1].scatter([centroid[0]], [centroid[1]], s=95, facecolor=COL["centroid"], edgecolor="white", linewidth=1.0, zorder=7)
    axes[1].text(centroid[0] + 0.08, centroid[1] + 0.06, "centroid", color=COL["centroid"], fontsize=8.5, fontweight="bold")
    axes[1].set_title("B. Mean Euclidean distance\nDistance from pixels to the centroid", fontsize=12, fontweight="bold", color=COL["text"], pad=12)

    # Alpha hull: occupied area in the first two PCA axes.
    axes[2].add_patch(Polygon(ALPHA_POLY, closed=True, facecolor=COL["alpha_fill"], edgecolor=COL["alpha_edge"], linewidth=2.0, alpha=0.80, zorder=2))
    draw_points(axes[2])
    axes[2].set_title("C. Alpha-hull area\nOccupied area in PC1-PC2 spectral space", fontsize=12, fontweight="bold", color=COL["text"], pad=12)

    fig.suptitle("Conceptual representation of PCA-based spectral heterogeneity metrics", x=0.02, y=0.98, ha="left", fontsize=15, fontweight="bold", color=COL["text"])
    fig.text(0.025, 0.91, "Each point represents a retained illuminated pixel projected into spectral PCA space.", fontsize=10.5, color=COL["subtext"])
    fig.subplots_adjust(left=0.055, right=0.98, top=0.78, bottom=0.12, wspace=0.20)
    fig.savefig(OUT_21, transparent=True)
    plt.close(fig)


if __name__ == "__main__":
    save_alpha_vs_convex()
    save_three_metric_concept()
    print(f"Created {OUT_20}")
    print(f"Created {OUT_21}")
