"""
Найкоротші шляхи на множині перешкод у 2D.

Постановка: На заданій множині h-перешкод (полігони, сумарна кількість вершин n)
для кожного розташування запитних точок S i T знайти найкоротший шлях між ними.

Алгоритм:
    1. Будуємо граф видимості G(V,E):
         V = {S, T} ∪ (усі вершини перешкод)
         (u,v) ∈ E ⇔ відрізок uv не проходить крізь жодну перешкоду.
    2. Запускаємо алгоритм Дейкстри від S до T по графу видимості.
       (Відомий результат: оптимальний шлях у середовищі полігональних перешкод
        складається лише з ребер графа видимості.)

Складність наївної реалізації: O(N^2) на побудову графа (N = n+2 вершин,
перевірка одного ребра — O(n)) і O(N^2) на Дейкстру → O(N^2) = O((n+h)^2).

GUI на matplotlib: клацанням мишею задаються вершини перешкод (натисніть
"Завершити перешкоду" аби замкнути контур), точки S і T; натискання
"Обрахувати шлях" будує граф видимості та найкоротший шлях.

Запуск:  python3 main.py
"""

import math
import heapq

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, RadioButtons
from matplotlib.patches import Polygon as MplPolygon

EPS = 1e-9


# =============================================================================
# Геометричні примітиви
# =============================================================================

def cross_v(o, a, b):
    """Знаковий векторний добуток (a-o) × (b-o)."""
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def segments_cross_properly(p1, p2, p3, p4):
    """True, якщо відкриті відрізки (p1,p2) і (p3,p4) перетинаються трансверсально
    (тобто власне у внутрішніх точках обох). Дотик у спільній вершині не вважається."""
    d1 = cross_v(p3, p4, p1)
    d2 = cross_v(p3, p4, p2)
    d3 = cross_v(p1, p2, p3)
    d4 = cross_v(p1, p2, p4)
    if (((d1 > EPS and d2 < -EPS) or (d1 < -EPS and d2 > EPS)) and
            ((d3 > EPS and d4 < -EPS) or (d3 < -EPS and d4 > EPS))):
        return True
    return False


def point_on_segment(p, a, b):
    """True, якщо точка p лежить на замкненому відрізку ab."""
    if abs(cross_v(a, b, p)) > EPS:
        return False
    return (min(a[0], b[0]) - EPS <= p[0] <= max(a[0], b[0]) + EPS and
            min(a[1], b[1]) - EPS <= p[1] <= max(a[1], b[1]) + EPS)


def point_on_polygon_boundary(p, polygon):
    m = len(polygon)
    for i in range(m):
        if point_on_segment(p, polygon[i], polygon[(i + 1) % m]):
            return True
    return False


def point_in_polygon_strict(p, polygon):
    """Ray casting. True, якщо точка p лежить *строго* всередині многокутника
    (не на межі)."""
    if point_on_polygon_boundary(p, polygon):
        return False
    x, y = p
    n = len(polygon)
    if n < 3:
        return False
    inside = False
    j = n - 1
    for i in range(n):
        xi, yi = polygon[i]
        xj, yj = polygon[j]
        if (yi > y) != (yj > y):
            denom = yj - yi
            if abs(denom) > EPS:
                x_int = (xj - xi) * (y - yi) / denom + xi
                if x < x_int:
                    inside = not inside
        j = i
    return inside


def visible(u, v, polygons):
    """Перевірка видимості точки v з точки u серед множини перешкод-полігонів.

    1. Відрізок uv не повинен трансверсально перетинати жодне ребро жодного полігона.
    2. Середина відрізка не має лежати строго всередині якоїсь перешкоди
       (закриває випадок, коли uv — діагональ опуклого полігона: проходить
        крізь внутрішність, але не перетинає ребер трансверсально).
    """
    for poly in polygons:
        m = len(poly)
        for i in range(m):
            a = poly[i]
            b = poly[(i + 1) % m]
            if segments_cross_properly(u, v, a, b):
                return False
    midpoint = ((u[0] + v[0]) * 0.5, (u[1] + v[1]) * 0.5)
    for poly in polygons:
        if point_in_polygon_strict(midpoint, poly):
            return False
    return True


def euclidean(p, q):
    return math.hypot(p[0] - q[0], p[1] - q[1])


# =============================================================================
# Граф видимості та Дейкстра
# =============================================================================

def build_visibility_graph(points, polygons):
    """Будує граф видимості над списком точок. Повертає (adj, list_of_edges).

    points[0] = S, points[1] = T, далі — вершини перешкод."""
    n = len(points)
    adj = [[] for _ in range(n)]
    edges = []
    for i in range(n):
        for j in range(i + 1, n):
            if visible(points[i], points[j], polygons):
                d = euclidean(points[i], points[j])
                adj[i].append((j, d))
                adj[j].append((i, d))
                edges.append((i, j))
    return adj, edges


def dijkstra(adj, source, target):
    """Знаходить найкоротший шлях у зваженому графі. Повертає (path, length)."""
    n = len(adj)
    dist = [math.inf] * n
    parent = [-1] * n
    dist[source] = 0.0
    heap = [(0.0, source)]
    while heap:
        d, u = heapq.heappop(heap)
        if d > dist[u]:
            continue
        if u == target:
            break
        for v, w in adj[u]:
            nd = d + w
            if nd < dist[v]:
                dist[v] = nd
                parent[v] = u
                heapq.heappush(heap, (nd, v))
    if math.isinf(dist[target]):
        return None, math.inf
    path = []
    cur = target
    while cur != -1:
        path.append(cur)
        cur = parent[cur]
    path.reverse()
    return path, dist[target]


def find_shortest_path(S, T, polygons):
    """Зручна обгортка: повертає (path_points, length, all_points, vg_edges)."""
    all_points = [S, T]
    for poly in polygons:
        for p in poly:
            all_points.append(p)
    adj, edges = build_visibility_graph(all_points, polygons)
    path_idx, length = dijkstra(adj, 0, 1)
    if path_idx is None:
        return None, math.inf, all_points, edges
    path_points = [all_points[i] for i in path_idx]
    return path_points, length, all_points, edges


# =============================================================================
# GUI
# =============================================================================

class ShortestPathApp:
    XLIM = (0.0, 20.0)
    YLIM = (0.0, 12.0)

    MODE_LABELS = {
        'Вершина перешкоди': 'polygon',
        'Точка S': 'set_s',
        'Точка T': 'set_t',
    }

    def __init__(self):
        self.polygons = []
        self.current_polygon = []
        self.S = None
        self.T = None
        self.mode = 'polygon'
        self.show_vg = False
        self.vg_edges = []
        self.all_points = []
        self.path_points = None
        self.path_len = None

        self.fig = plt.figure(figsize=(13.5, 8.5))
        try:
            self.fig.canvas.manager.set_window_title(
                "Найкоротший шлях на множині перешкод")
        except Exception:
            pass

        # Основна полотно
        self.ax = self.fig.add_axes([0.05, 0.10, 0.66, 0.85])
        self._setup_axes()

        # Радіо: режим вводу
        ax_radio = self.fig.add_axes([0.74, 0.78, 0.24, 0.17])
        ax_radio.set_title("Режим вводу", fontsize=10)
        self.radio = RadioButtons(
            ax_radio, tuple(self.MODE_LABELS.keys()), active=0)
        self.radio.on_clicked(self._set_mode)

        # Кнопки
        bx, bw, bh, gap = 0.74, 0.24, 0.045, 0.012
        by = 0.74

        def add_btn(label, callback):
            nonlocal by
            by -= (bh + gap)
            axb = self.fig.add_axes([bx, by, bw, bh])
            btn = Button(axb, label)
            btn.on_clicked(callback)
            return btn

        self.btn_close = add_btn('Завершити перешкоду', self._close_polygon)
        self.btn_compute = add_btn('Обрахувати шлях', self._compute)
        self.btn_undo = add_btn('Скасувати останню', self._undo)
        self.btn_clear = add_btn('Очистити все', self._clear_all)
        self.btn_vg = add_btn('Показати/сховати граф видимості',
                              self._toggle_vg)
        self.btn_load1 = add_btn('Приклад 1 (опуклі)', self._load_example1)
        self.btn_load2 = add_btn('Приклад 2 (неопуклі)', self._load_example2)

        # Рядок стану
        self.status = self.fig.text(
            0.05, 0.02,
            "Клацайте на полі, щоб додати вершини перешкоди. "
            "Натисніть 'Завершити перешкоду', щоб замкнути контур.",
            fontsize=10)

        self.fig.canvas.mpl_connect('button_press_event', self._on_click)
        self._redraw()

    # ---- допоміжне ----

    def _setup_axes(self):
        self.ax.set_xlim(*self.XLIM)
        self.ax.set_ylim(*self.YLIM)
        self.ax.set_aspect('equal')
        self.ax.grid(True, alpha=0.3)
        self.ax.set_title("Найкоротший шлях на множині перешкод у 2D",
                          fontsize=12)
        self.ax.set_xlabel("x")
        self.ax.set_ylabel("y")

    def _set_status(self, text):
        self.status.set_text(text)
        self.fig.canvas.draw_idle()

    # ---- callbacks ----

    def _set_mode(self, label):
        self.mode = self.MODE_LABELS[label]
        instr = {
            'polygon': "Клацайте, щоб додати вершину поточної перешкоди. "
                       "Потім натисніть 'Завершити перешкоду'.",
            'set_s':   "Клацніть, щоб встановити стартову точку S.",
            'set_t':   "Клацніть, щоб встановити кінцеву точку T.",
        }
        self._set_status(instr[self.mode])

    def _on_click(self, event):
        if event.inaxes is not self.ax:
            return
        if event.button != 1:
            return
        if event.xdata is None or event.ydata is None:
            return
        x, y = event.xdata, event.ydata
        if self.mode == 'polygon':
            self.current_polygon.append((x, y))
        elif self.mode == 'set_s':
            self.S = (x, y)
        elif self.mode == 'set_t':
            self.T = (x, y)
        # Будь-яка зміна вхідних даних робить попередній шлях/граф застарілим
        self.path_points = None
        self.path_len = None
        self.vg_edges = []
        self._redraw()

    def _close_polygon(self, event):
        if len(self.current_polygon) >= 3:
            self.polygons.append(self.current_polygon[:])
            self._set_status(
                f"Додано перешкоду з {len(self.current_polygon)} вершин. "
                f"Загалом перешкод: {len(self.polygons)}, "
                f"сумарно вершин: {sum(len(p) for p in self.polygons)}.")
        else:
            self._set_status("Перешкода повинна мати щонайменше 3 вершини.")
        self.current_polygon = []
        self.path_points = None
        self.path_len = None
        self.vg_edges = []
        self._redraw()

    def _undo(self, event):
        if self.current_polygon:
            self.current_polygon.pop()
            self._set_status("Скасовано останню вершину поточної перешкоди.")
        elif self.polygons:
            removed = self.polygons.pop()
            self._set_status(f"Видалено останню перешкоду ({len(removed)} вершин).")
        else:
            self._set_status("Нічого скасовувати.")
        self.path_points = None
        self.path_len = None
        self.vg_edges = []
        self._redraw()

    def _clear_all(self, event):
        self.polygons = []
        self.current_polygon = []
        self.S = None
        self.T = None
        self.path_points = None
        self.path_len = None
        self.vg_edges = []
        self.show_vg = False
        self._set_status("Очищено все.")
        self._redraw()

    def _toggle_vg(self, event):
        self.show_vg = not self.show_vg
        self._redraw()

    def _load_example1(self, event):
        self.polygons = [
            [(3.5, 6.0), (5.0, 8.5), (6.5, 6.5), (5.5, 5.0)],
            [(8.0, 3.0), (10.5, 4.5), (9.5, 7.5)],
            [(12.0, 6.0), (14.0, 8.5), (15.5, 6.0), (14.0, 4.0)],
            [(7.0, 9.5), (9.5, 11.0), (8.0, 9.0)],
        ]
        self.S = (1.0, 6.0)
        self.T = (18.5, 6.0)
        self.current_polygon = []
        self.path_points = None
        self.path_len = None
        self.vg_edges = []
        self._set_status("Завантажено приклад 1. Натисніть 'Обрахувати шлях'.")
        self._redraw()

    def _load_example2(self, event):
        # Неопуклий "U"-полігон + прямокутник
        self.polygons = [
            [(4.0, 2.5), (4.0, 9.0), (5.5, 9.0), (5.5, 4.5),
             (8.0, 4.5), (8.0, 9.0), (9.5, 9.0), (9.5, 2.5)],
            [(12.5, 4.0), (16.0, 4.0), (16.0, 8.0), (12.5, 8.0)],
            [(13.5, 9.5), (15.0, 11.0), (16.5, 9.5)],
        ]
        self.S = (1.5, 6.0)
        self.T = (18.5, 6.0)
        self.current_polygon = []
        self.path_points = None
        self.path_len = None
        self.vg_edges = []
        self._set_status("Завантажено приклад 2 (неопуклі). "
                         "Натисніть 'Обрахувати шлях'.")
        self._redraw()

    def _compute(self, event):
        if self.S is None or self.T is None:
            self._set_status("Помилка: задайте обидві точки S і T.")
            return
        for i, poly in enumerate(self.polygons):
            if point_in_polygon_strict(self.S, poly):
                self._set_status(f"Помилка: точка S всередині перешкоди #{i + 1}.")
                return
            if point_in_polygon_strict(self.T, poly):
                self._set_status(f"Помилка: точка T всередині перешкоди #{i + 1}.")
                return

        path_points, length, all_points, edges = find_shortest_path(
            self.S, self.T, self.polygons)
        self.all_points = all_points
        self.vg_edges = edges
        if path_points is None:
            self.path_points = None
            self.path_len = None
            self._set_status("Шлях не знайдено (S і T у різних компонентах зв'язності).")
        else:
            self.path_points = path_points
            self.path_len = length
            n_total = sum(len(p) for p in self.polygons)
            self._set_status(
                f"Знайдено! Довжина шляху = {length:.4f}.  "
                f"h = {len(self.polygons)} перешкод, n = {n_total} вершин.  "
                f"|V| = {len(all_points)}, |E| = {len(edges)},  "
                f"шлях має {len(path_points)} точок.")
        self._redraw()

    # ---- малювання ----

    def _redraw(self):
        self.ax.clear()
        self._setup_axes()

        # Перешкоди
        for poly in self.polygons:
            patch = MplPolygon(poly, closed=True,
                               facecolor='#9aa0a6', edgecolor='black',
                               alpha=0.75, zorder=2, linewidth=1.2)
            self.ax.add_patch(patch)

        # Поточний (незавершений) полігон
        if self.current_polygon:
            xs = [p[0] for p in self.current_polygon]
            ys = [p[1] for p in self.current_polygon]
            self.ax.plot(xs, ys, '--', color='#1f77b4', linewidth=1, zorder=3)
            self.ax.plot(xs, ys, 'o', color='#1f77b4', markersize=4, zorder=3)

        # Граф видимості
        if self.show_vg and self.vg_edges:
            for (i, j) in self.vg_edges:
                p, q = self.all_points[i], self.all_points[j]
                self.ax.plot([p[0], q[0]], [p[1], q[1]],
                             color='#d0d0d0', linewidth=0.5, zorder=1)

        # Шлях
        if self.path_points and len(self.path_points) >= 2:
            xs = [p[0] for p in self.path_points]
            ys = [p[1] for p in self.path_points]
            self.ax.plot(xs, ys, color='#0066cc', linewidth=2.6, zorder=4,
                         label=f'Найкоротший шлях, L = {self.path_len:.3f}')
            if len(xs) > 2:
                self.ax.plot(xs[1:-1], ys[1:-1], 'o', color='#0066cc',
                             markersize=5, zorder=4)
            self.ax.legend(loc='upper right', fontsize=10)

        # S
        if self.S is not None:
            self.ax.plot(self.S[0], self.S[1], 'o',
                         color='#22aa22', markersize=12, zorder=5)
            self.ax.annotate('S', self.S, textcoords="offset points",
                             xytext=(8, 8), fontsize=14,
                             color='#22aa22', fontweight='bold')
        # T
        if self.T is not None:
            self.ax.plot(self.T[0], self.T[1], 'o',
                         color='#cc2222', markersize=12, zorder=5)
            self.ax.annotate('T', self.T, textcoords="offset points",
                             xytext=(8, 8), fontsize=14,
                             color='#cc2222', fontweight='bold')

        self.fig.canvas.draw_idle()

    def run(self):
        plt.show()


# =============================================================================
# Точка входу
# =============================================================================

def main():
    app = ShortestPathApp()
    app.run()


if __name__ == "__main__":
    main()
