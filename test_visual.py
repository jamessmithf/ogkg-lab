"""Безголова візуалізація прикладів: зберігає PNG із знайденим шляхом."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon

import importlib.util
spec = importlib.util.spec_from_file_location("main_mod", "/home/viktor/Downloads/ogkg/main.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def render(name, polygons, S, T, out_path):
    path, length, all_points, edges = mod.find_shortest_path(S, T, polygons)
    fig, ax = plt.subplots(figsize=(11, 7))
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.set_title(f"{name}: L = {length:.4f}" if path else f"{name}: path NOT found")

    for poly in polygons:
        ax.add_patch(MplPolygon(poly, closed=True, facecolor='#9aa0a6',
                                edgecolor='black', alpha=0.75, linewidth=1.2))

    # граф видимості (легкий)
    for (i, j) in edges:
        p, q = all_points[i], all_points[j]
        ax.plot([p[0], q[0]], [p[1], q[1]], color='#dadada',
                linewidth=0.4, zorder=1)

    if path:
        xs = [p[0] for p in path]
        ys = [p[1] for p in path]
        ax.plot(xs, ys, color='#0066cc', linewidth=2.5, zorder=4)
        if len(xs) > 2:
            ax.plot(xs[1:-1], ys[1:-1], 'o', color='#0066cc',
                    markersize=5, zorder=4)

    ax.plot(S[0], S[1], 'o', color='#22aa22', markersize=12, zorder=5)
    ax.annotate('S', S, textcoords='offset points', xytext=(8, 8),
                fontsize=14, color='#22aa22', fontweight='bold')
    ax.plot(T[0], T[1], 'o', color='#cc2222', markersize=12, zorder=5)
    ax.annotate('T', T, textcoords='offset points', xytext=(8, 8),
                fontsize=14, color='#cc2222', fontweight='bold')

    # автопідбір меж
    all_x = [p[0] for poly in polygons for p in poly] + [S[0], T[0]]
    all_y = [p[1] for poly in polygons for p in poly] + [S[1], T[1]]
    pad = 1.0
    ax.set_xlim(min(all_x) - pad, max(all_x) + pad)
    ax.set_ylim(min(all_y) - pad, max(all_y) + pad)

    fig.savefig(out_path, dpi=120, bbox_inches='tight')
    print(f"{name}: довжина = {length:.4f}, вершин шляху = {len(path) if path else 0}, "
          f"ребер видимості = {len(edges)} -> {out_path}")


if __name__ == "__main__":
    # Приклад 1: опуклі перешкоди
    polygons1 = [
        [(3.5, 6.0), (5.0, 8.5), (6.5, 6.5), (5.5, 5.0)],
        [(8.0, 3.0), (10.5, 4.5), (9.5, 7.5)],
        [(12.0, 6.0), (14.0, 8.5), (15.5, 6.0), (14.0, 4.0)],
        [(7.0, 9.5), (9.5, 11.0), (8.0, 9.0)],
    ]
    render("Example 1 (convex)", polygons1, (1.0, 6.0), (18.5, 6.0),
           "/home/viktor/Downloads/ogkg/example1.png")

    # Приклад 2: неопуклі перешкоди
    polygons2 = [
        [(4.0, 2.5), (4.0, 9.0), (5.5, 9.0), (5.5, 4.5),
         (8.0, 4.5), (8.0, 9.0), (9.5, 9.0), (9.5, 2.5)],
        [(12.5, 4.0), (16.0, 4.0), (16.0, 8.0), (12.5, 8.0)],
        [(13.5, 9.5), (15.0, 11.0), (16.5, 9.5)],
    ]
    render("Example 2 (non-convex)", polygons2, (1.5, 6.0), (18.5, 6.0),
           "/home/viktor/Downloads/ogkg/example2.png")

    # Приклад 3: ситуація як на ілюстрації задачі (S зліва, T справа, кілька форм)
    polygons3 = [
        [(2.5, 5.5), (4.5, 7.0), (6.0, 5.8), (4.5, 6.2), (4.0, 5.5)],  # стріла
        [(3.5, 3.0), (4.5, 4.5), (5.5, 3.0), (4.5, 1.5)],  # ромб
        [(8.0, 4.0), (10.0, 6.0), (10.5, 3.0)],  # трикутник
    ]
    render("Example 3 (mixed)", polygons3, (1.0, 5.0), (12.0, 5.0),
           "/home/viktor/Downloads/ogkg/example3.png")
