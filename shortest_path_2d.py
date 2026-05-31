#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Найкоротші шляхи на множині перешкод у 2D.

Лабораторна робота з курсу «Обчислювальна геометрія та комп'ютерна графіка».

Задача: на заданій множині h перешкод (опуклих/простих полігонів, сумарна
кількість вершин яких дорівнює n) для кожного розташування запитних точок S і T
знайти найкоротший (евклідів) шлях між ними, що не перетинає жодну перешкоду.

Метод: граф видимості + пошук A* з лінивим (on-demand) обчисленням сусідів та
просторовою сіткою для прискорення перевірок видимості. Два режими: ТОЧНИЙ
(гарантований оптимум) і ШВИДКИЙ (лише ~k найближчих кандидатів — для дуже
щільних сцен 1000+ перешкод; шлях коректний, на <1% довший за оптимум, але
у десятки-сотні разів швидший).

Запуск GUI потребує робочого Tk. На цьому Mac (macOS 15.6) Tk 8.5 системного
/usr/bin/python3 аварійно завершується, тому використовуйте Python 3.13 з Tk 9.0:
    /opt/homebrew/bin/python3.13 shortest_path_2d.py         # графічний інтерфейс
(один раз встановити Tk:  brew install python-tk@3.13)

Режими без GUI працюють під будь-яким Python 3.x (зокрема /usr/bin/python3):
    python3 shortest_path_2d.py --selftest      # перевірка коректності алгоритму
    python3 shortest_path_2d.py --bench 1000     # заміри ефективності

Залежності: лише стандартна бібліотека Python 3.x (tkinter, math, heapq, random, time).
"""

import math
import heapq
import random
import time
import sys

EPS = 1e-9
INF = float("inf")


# =============================================================================
#  1. Геометричні примітиви
# =============================================================================

def orient(ax, ay, bx, by, cx, cy):
    """Знак векторного добутку (b-a) x (c-a): >0 — ліворуч, <0 — праворуч, 0 — колінеарні."""
    v = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
    if v > EPS:
        return 1
    if v < -EPS:
        return -1
    return 0


def segments_properly_intersect(ax, ay, bx, by, cx, cy, dx, dy):
    """True, якщо відрізки ab та cd перетинаються у внутрішній (власній) точці обох."""
    d1 = orient(cx, cy, dx, dy, ax, ay)
    d2 = orient(cx, cy, dx, dy, bx, by)
    d3 = orient(ax, ay, bx, by, cx, cy)
    d4 = orient(ax, ay, bx, by, dx, dy)
    return d1 * d2 < 0 and d3 * d4 < 0


def point_in_polygon(px, py, poly):
    """Тест 'точка строго всередині простого полігона' методом кидання променя."""
    inside = False
    n = len(poly)
    j = n - 1
    for i in range(n):
        xi, yi = poly[i]
        xj, yj = poly[j]
        if (yi > py) != (yj > py):
            xint = (xj - xi) * (py - yi) / (yj - yi) + xi
            if px < xint:
                inside = not inside
        j = i
    return inside


def dist(ax, ay, bx, by):
    return math.hypot(ax - bx, ay - by)


# =============================================================================
#  2. Просторова сітка (рівномірна) для прискорення перевірок видимості
# =============================================================================

class SpatialGrid:
    """Рівномірна сітка: ребра перешкод розкладаються по комірках, які вони
    перетинають. Для перевірки відрізка переглядаються тільки ребра у комірках
    уздовж нього (а не всі n ребер сцени)."""

    def __init__(self, minx, miny, maxx, maxy, cell):
        self.minx = minx
        self.miny = miny
        self.cell = cell
        self.ncols = max(1, int((maxx - minx) / cell) + 1)
        self.nrows = max(1, int((maxy - miny) / cell) + 1)
        self.cells = {}            # (cx, cy) -> list[edge_id]
        self.edges = []            # edge_id -> (ax, ay, bx, by, vid_a, vid_b)
        self._epoch = []           # позначки «вже перевірено» (дедуплікація)
        self._cur = 0

    def _cell_of(self, x, y):
        cx = int((x - self.minx) / self.cell)
        cy = int((y - self.miny) / self.cell)
        return cx, cy

    def _cells_on_segment(self, ax, ay, bx, by):
        """Комірки, які перетинає відрізок ab (обхід комірок Amanatides–Woo)."""
        cx, cy = self._cell_of(ax, ay)
        cx1, cy1 = self._cell_of(bx, by)
        out = [(cx, cy)]
        dx = bx - ax
        dy = by - ay
        stepx = 1 if dx > 0 else -1
        stepy = 1 if dy > 0 else -1

        if abs(dx) > EPS:
            nb = self.minx + (cx + (1 if dx > 0 else 0)) * self.cell
            t_max_x = (nb - ax) / dx
            t_delta_x = self.cell / abs(dx)
        else:
            t_max_x = INF
            t_delta_x = INF
        if abs(dy) > EPS:
            nb = self.miny + (cy + (1 if dy > 0 else 0)) * self.cell
            t_max_y = (nb - ay) / dy
            t_delta_y = self.cell / abs(dy)
        else:
            t_max_y = INF
            t_delta_y = INF

        guard = self.ncols + self.nrows + 4
        while (cx, cy) != (cx1, cy1) and guard > 0:
            if t_max_x < t_max_y:
                t_max_x += t_delta_x
                cx += stepx
            else:
                t_max_y += t_delta_y
                cy += stepy
            out.append((cx, cy))
            guard -= 1
        return out

    def add_edge(self, ax, ay, bx, by, vid_a, vid_b):
        eid = len(self.edges)
        self.edges.append((ax, ay, bx, by, vid_a, vid_b))
        for key in self._cells_on_segment(ax, ay, bx, by):
            self.cells.setdefault(key, []).append(eid)
        return eid

    def seg_blocked(self, ax, ay, bx, by, skip_a, skip_b):
        """True, якщо відрізок ab власне перетинає якесь ребро перешкоди.

        Комірки обходяться у порядку від a до b, тому перешкода виявляється
        рано і перевірка завершується достроково (ранній вихід). Ребра,
        інцидентні до вершин skip_a/skip_b, ігноруються."""
        if len(self._epoch) != len(self.edges):
            self._epoch = [0] * len(self.edges)
        self._cur += 1
        ep = self._cur
        epoch = self._epoch
        edges = self.edges
        for key in self._cells_on_segment(ax, ay, bx, by):
            bucket = self.cells.get(key)
            if not bucket:
                continue
            for eid in bucket:
                if epoch[eid] == ep:
                    continue
                epoch[eid] = ep
                ex0, ey0, ex1, ey1, ea, eb = edges[eid]
                if ea == skip_a or ea == skip_b or eb == skip_a or eb == skip_b:
                    continue
                if segments_properly_intersect(ax, ay, bx, by,
                                                ex0, ey0, ex1, ey1):
                    return True
        return False


# =============================================================================
#  3. Сцена: перешкоди, вершини, побудова графу видимості та пошук A*
# =============================================================================

class Vertex:
    __slots__ = ("x", "y", "poly", "idx", "size")

    def __init__(self, x, y, poly, idx, size):
        self.x = x          # координати
        self.y = y
        self.poly = poly    # id перешкоди
        self.idx = idx      # локальний індекс у полігоні
        self.size = size    # кількість вершин полігона


class Scene:
    """Геометрична сцена та алгоритм найкоротшого шляху.

    Відокремлена від GUI, тому її можна тестувати/заміряти без екрана."""

    def __init__(self):
        self.polygons = []     # list[list[(x, y)]]
        self.vertices = []     # list[Vertex]  (глобальна нумерація 0..V-1)
        self.grid = None
        self.vis_cache = {}    # vid -> list[(vid, w)]  видимість вершина↔вершина (точний режим)
        self.fast_cache = {}   # vid -> list[(vid, w)]  лише найближчі (швидкий режим)
        self.vgrid = None      # сітка вершин (для пошуку найближчих кандидатів)
        self.vcell = None
        self.S = None
        self.T = None
        self.S_ID = -1
        self.T_ID = -2

    # ---- побудова сцени -----------------------------------------------------

    def clear(self):
        self.polygons.clear()
        self.vertices.clear()
        self.grid = None
        self.vis_cache.clear()
        self.fast_cache.clear()
        self.vgrid = None
        self.S = None
        self.T = None

    def add_polygon(self, pts):
        """Додає перешкоду (список вершин у порядку обходу). Дублікати ігноруються."""
        if len(pts) < 3:
            return
        pid = len(self.polygons)
        self.polygons.append(list(pts))
        size = len(pts)
        for i, (x, y) in enumerate(pts):
            self.vertices.append(Vertex(x, y, pid, i, size))
        self.grid = None
        self.vis_cache.clear()
        self.fast_cache.clear()
        self.vgrid = None

    def build_grid(self):
        """Передобробка: побудова просторової сітки з усіх ребер перешкод."""
        if not self.polygons:
            self.grid = SpatialGrid(0, 0, 1, 1, 1)
            return
        xs = [p[0] for poly in self.polygons for p in poly]
        ys = [p[1] for poly in self.polygons for p in poly]
        if self.S:
            xs.append(self.S[0]); ys.append(self.S[1])
        if self.T:
            xs.append(self.T[0]); ys.append(self.T[1])
        minx, maxx = min(xs) - 10, max(xs) + 10
        miny, maxy = min(ys) - 10, max(ys) + 10
        # розмір комірки ~ середня довжина ребра (баланс пам'ять/швидкість)
        total_edges = sum(len(p) for p in self.polygons)
        span = max(maxx - minx, maxy - miny)
        cell = max(8.0, span / max(8, int(math.sqrt(total_edges)) + 1))
        self.grid = SpatialGrid(minx, miny, maxx, maxy, cell)
        # глобальний індекс першої вершини кожного полігона
        base = 0
        for poly in self.polygons:
            m = len(poly)
            for i in range(m):
                a = poly[i]
                b = poly[(i + 1) % m]
                self.grid.add_edge(a[0], a[1], b[0], b[1], base + i, base + (i + 1) % m)
            base += m
        # сітка вершин для пошуку найближчих кандидатів (швидкий режим)
        self.vcell = cell * 1.5
        self.vgrid = {}
        for i, v in enumerate(self.vertices):
            key = (int((v.x - minx) / self.vcell), int((v.y - miny) / self.vcell))
            self.vgrid.setdefault(key, []).append(i)
        self._vminx, self._vminy = minx, miny

    # ---- видимість ----------------------------------------------------------

    def _pos(self, vid):
        if vid == self.S_ID:
            return self.S
        if vid == self.T_ID:
            return self.T
        v = self.vertices[vid]
        return (v.x, v.y)

    def _adjacent(self, i, j):
        """Чи є вершини i, j сусідніми у спільному полігоні (ребром перешкоди)."""
        vi, vj = self.vertices[i], self.vertices[j]
        if vi.poly != vj.poly:
            return False
        d = abs(vi.idx - vj.idx)
        return d == 1 or d == vi.size - 1

    def visible(self, uid, vid):
        """True, якщо відкритий відрізок між вершинами uid та vid не перетинає
        жодну перешкоду (тобто вони взаємно видимі)."""
        ax, ay = self._pos(uid)
        bx, by = self._pos(vid)

        # випадок двох вершин одного полігона: хорда всередині перешкоди => невидимі
        if uid >= 0 and vid >= 0:
            vu, vv = self.vertices[uid], self.vertices[vid]
            if vu.poly == vv.poly and not self._adjacent(uid, vid):
                mx, my = (ax + bx) / 2.0, (ay + by) / 2.0
                if point_in_polygon(mx, my, self.polygons[vu.poly]):
                    return False

        # перетин з ребрами перешкод (ребра, інцидентні до uid/vid, пропускаємо)
        return not self.grid.seg_blocked(ax, ay, bx, by, uid, vid)

    def vertex_neighbors(self, vid):
        """Список видимих вершин-перешкод (vid, вага) з кешуванням. Незалежний
        від S/T, тому обчислюється один раз і використовується у всіх запитах."""
        cached = self.vis_cache.get(vid)
        if cached is not None:
            return cached
        v = self.vertices[vid]
        res = []
        for j in range(len(self.vertices)):
            if j == vid:
                continue
            if self.visible(vid, j):
                w = self.vertices[j]
                res.append((j, dist(v.x, v.y, w.x, w.y)))
        self.vis_cache[vid] = res
        return res

    # ---- швидкий (наближений) режим: лише найближчі кандидати ----------------

    def _nearest_candidates(self, px, py, k):
        """~k найближчих вершин до (px, py) — кільцевий пошук по сітці вершин.
        Адаптивний: у щільних зонах кандидати близькі, у відкритих їх мало,
        тож захоплюються й далекі (тому шлях лишається майже оптимальним)."""
        cell = self.vcell
        cx = int((px - self._vminx) / cell)
        cy = int((py - self._vminy) / cell)
        found = []
        r = 0
        rmax = max(self.grid.ncols, self.grid.nrows) + 2
        while True:
            for gx in range(cx - r, cx + r + 1):
                for gy in range(cy - r, cy + r + 1):
                    if max(abs(gx - cx), abs(gy - cy)) != r:
                        continue          # лише «кільце» радіуса r
                    bucket = self.vgrid.get((gx, gy))
                    if bucket:
                        found.extend(bucket)
            if (len(found) >= k and r >= 2) or r > rmax:
                break
            r += 1
        return found

    def vertex_neighbors_fast(self, vid, k):
        """Видимі сусіди серед ~k найближчих вершин (з кешуванням). Швидкий
        режим: O(k) перевірок замість O(n) — для дуже щільних сцен 1000+."""
        cached = self.fast_cache.get(vid)
        if cached is not None:
            return cached
        v = self.vertices[vid]
        res = []
        for j in self._nearest_candidates(v.x, v.y, k):
            if j == vid:
                continue
            if self.visible(vid, j):
                w = self.vertices[j]
                res.append((j, dist(v.x, v.y, w.x, w.y)))
        self.fast_cache[vid] = res
        return res

    def build_full_visibility(self):
        """Повна передобробка: граф видимості для всіх вершин (для візуалізації
        та точних замірів). Складність O(n²) — придатна лише для малих сцен."""
        if self.grid is None:
            self.build_grid()
        for i in range(len(self.vertices)):
            self.vertex_neighbors(i)

    def build_knn_visibility(self, k=90):
        """k-найближчий граф видимості — той самий, що бачить ШВИДКИЙ режим A*.
        Складність O(n·k) замість O(n²), тож будується у десятки разів швидше
        і дає набагато менше ребер. Використовується для візуалізації великих
        сцен, де повний граф був би і повільним, і нечитабельним «клубком»."""
        if self.grid is None:
            self.build_grid()
        for i in range(len(self.vertices)):
            self.vertex_neighbors_fast(i, k)

    # ---- пошук найкоротшого шляху A* ----------------------------------------

    def _point_inside_any(self, p):
        for poly in self.polygons:
            if point_in_polygon(p[0], p[1], poly):
                return True
        return False

    def shortest_path(self, S, T, fast=False, k=90):
        """Найкоротший шлях S->T методом A* з лінивим обчисленням видимості.

        fast=False — ТОЧНИЙ режим (гарантований оптимум).
        fast=True  — ШВИДКИЙ режим: на кожному кроці розглядаються лише ~k
                     найближчих вершин, тож для дуже щільних сцен (1000+ перешкод)
                     час падає у десятки-сотні разів. Шлях завжди коректний (без
                     перетину перешкод), але може бути на ~1% довшим за оптимум.

        Повертає словник зі шляхом, довжиною, часом і статистикою."""
        self.S, self.T = S, T
        if self.grid is None:
            self.build_grid()
        t0 = time.perf_counter()

        if self._point_inside_any(S) or self._point_inside_any(T):
            return {"path": None, "length": 0.0, "time": 0.0,
                    "expanded": 0, "status": "S або T всередині перешкоди"}

        Tx, Ty = T

        def h(vid):                       # евристика — пряма відстань до T
            x, y = self._pos(vid)
            return dist(x, y, Tx, Ty)

        def s_visible_vertices():
            if fast:
                cand = self._nearest_candidates(S[0], S[1], k)
            else:
                cand = range(len(self.vertices))
            out = []
            for j in cand:
                if self.visible(self.S_ID, j):
                    w = self.vertices[j]
                    out.append((j, dist(S[0], S[1], w.x, w.y)))
            return out

        def neighbors(uid):
            if uid == self.S_ID:
                out = s_visible_vertices()
                if self.visible(self.S_ID, self.T_ID):
                    out.append((self.T_ID, dist(S[0], S[1], Tx, Ty)))
                return out
            if fast:
                out = list(self.vertex_neighbors_fast(uid, k))
            else:
                out = list(self.vertex_neighbors(uid))
            if self.visible(uid, self.T_ID):
                x, y = self._pos(uid)
                out.append((self.T_ID, dist(x, y, Tx, Ty)))
            return out

        g = {self.S_ID: 0.0}
        came = {}
        closed = set()
        openh = [(h(self.S_ID), 0.0, self.S_ID)]
        expanded = 0

        while openh:
            f, gu, u = heapq.heappop(openh)
            if u in closed:
                continue
            closed.add(u)
            expanded += 1
            if u == self.T_ID:
                break
            for v, w in neighbors(u):
                ng = gu + w
                if ng < g.get(v, INF) - EPS:
                    g[v] = ng
                    came[v] = u
                    heapq.heappush(openh, (ng + h(v), ng, v))

        elapsed = time.perf_counter() - t0
        if self.T_ID not in came and self.T_ID not in closed:
            if fast:
                # у швидкому режимі обмеження кандидатів (~k найближчих) могло
                # розірвати граф і дати хибне «шляху немає» — відкат до точного
                return self.shortest_path(S, T, fast=False)
            return {"path": None, "length": 0.0, "time": elapsed,
                    "expanded": expanded, "status": "Шлях не існує"}

        # відновлення шляху
        path_ids = [self.T_ID]
        cur = self.T_ID
        while cur != self.S_ID:
            cur = came[cur]
            path_ids.append(cur)
        path_ids.reverse()
        path = [self._pos(vid) for vid in path_ids]
        length = g[self.T_ID]
        return {"path": path, "length": length, "time": elapsed,
                "expanded": expanded, "status": "OK",
                "mode": "швидкий" if fast else "точний"}


# =============================================================================
#  4. Генератор перешкод (для авто-режиму та замірів ефективності)
# =============================================================================

def random_convex_polygon(cx, cy, r, k, rng):
    """Випадковий опуклий полігон із k вершин навколо центру (cx, cy)."""
    angles = sorted(rng.uniform(0, 2 * math.pi) for _ in range(k))
    pts = []
    for a in angles:
        rad = r * rng.uniform(0.6, 1.0)
        pts.append((cx + rad * math.cos(a), cy + rad * math.sin(a)))
    return pts


def generate_scene(scene, n_obstacles, width, height, rng,
                   vmin=3, vmax=6, rmin=7, rmax=17):
    """Заповнює сцену n_obstacles випадковими опуклими перешкодами."""
    scene.clear()
    for _ in range(n_obstacles):
        cx = rng.uniform(rmax, width - rmax)
        cy = rng.uniform(rmax, height - rmax)
        r = rng.uniform(rmin, rmax)
        k = rng.randint(vmin, vmax)
        scene.add_polygon(random_convex_polygon(cx, cy, r, k, rng))
    # S, T — у вільних точках
    def free_point():
        for _ in range(2000):
            p = (rng.uniform(5, width - 5), rng.uniform(5, height - 5))
            if not scene._point_inside_any(p):
                return p
        return (5, 5)
    scene.S = free_point()
    scene.T = free_point()
    scene.build_grid()


def generate_worst_case(scene, n_walls, width, height):
    """Найгірший ввід — «слалом». Кожна колонка-стіна сягає обох країв (верхнього
    й нижнього) і має єдиний прохід, чия висота почергово то вгорі, то внизу.
    Тому пройти крізь сцену неможливо інакше, ніж зигзагом крізь усі проходи:
    шлях має Θ(n_walls) ланок, а A* розкриває усі вершини — це й демонструє
    найгіршу (квадратичну) поведінку алгоритму."""
    scene.clear()
    ft = 60.0                                  # товщина рамки
    rect = lambda x0, y0, x1, y1: scene.add_polygon(
        [(x0, y0), (x1, y0), (x1, y1), (x0, y1)])
    # замкнена рамка: S і T опиняються «в пастці», шлях не може її обійти
    rect(0, 0, width, ft)                      # верх
    rect(0, height - ft, width, height)        # низ
    rect(0, 0, ft, height)                     # ліво
    rect(width - ft, 0, width, height)         # право

    yc0, yc1 = ft, height - ft                 # межі коридору
    corr = yc1 - yc0
    upper = yc0 + 0.30 * corr                  # центр верхнього проходу
    lower = yc0 + 0.70 * corr                  # центр нижнього проходу
    gap_half = 0.12 * corr
    step = (width - 2 * ft) / (n_walls + 1)
    # ширина стіни — завжди частка кроку, тож між колонами лишається проміжок
    # (інакше при великій к-сті стін колони злипаються й коридор замуровується)
    wall_w = step * 0.4
    for i in range(n_walls):
        x = ft + step * (i + 1)
        c = upper if i % 2 == 0 else lower
        xL, xR = x - wall_w / 2, x + wall_w / 2
        # колона з двох частин (перекриває рамку, аби не було щілин для ковзання)
        rect(xL, ft - 10, xR, c - gap_half)            # верхня частина
        rect(xL, c + gap_half, xR, height - ft + 10)   # нижня частина
    scene.S = (ft + step * 0.5, upper)
    scene.T = (width - ft - step * 0.5,
               upper if (n_walls - 1) % 2 == 0 else lower)
    scene.build_grid()


# =============================================================================
#  5. Графічний інтерфейс (Tkinter)
# =============================================================================

def run_gui():
    import tkinter as tk
    from tkinter import ttk

    W, H = 980, 720

    scene = Scene()
    rng = random.Random(12345)

    state = {
        "mode": "draw",          # draw | set_s | set_t
        "current": [],           # вершини полігона, що малюється
        "show_graph": False,
        "graph_knn": False,      # True → малюємо k-найближчий граф (великі сцени)
        "last_result": None,
    }

    root = tk.Tk()
    root.title("Найкоротші шляхи на множині перешкод у 2D")

    left = ttk.Frame(root, padding=8)
    left.pack(side="left", fill="y")
    canvas = tk.Canvas(root, width=W, height=H, bg="white", highlightthickness=1,
                       highlightbackground="#888")
    canvas.pack(side="right", fill="both", expand=True)

    # ---- елементи керування -------------------------------------------------
    ttk.Label(left, text="Режим вводу", font=("TkDefaultFont", 10, "bold")).pack(anchor="w")
    mode_var = tk.StringVar(value="draw")

    def set_mode():
        state["mode"] = mode_var.get()
        status(f"Режим: {mode_var.get()}")

    ttk.Radiobutton(left, text="Малювати перешкоду", variable=mode_var,
                    value="draw", command=set_mode).pack(anchor="w")
    ttk.Radiobutton(left, text="Поставити S (старт)", variable=mode_var,
                    value="set_s", command=set_mode).pack(anchor="w")
    ttk.Radiobutton(left, text="Поставити T (фініш)", variable=mode_var,
                    value="set_t", command=set_mode).pack(anchor="w")

    ttk.Separator(left, orient="horizontal").pack(fill="x", pady=6)
    ttk.Label(left, text="Ручний ввід (мишкою):").pack(anchor="w")
    ttk.Label(left, text="ЛКМ — додати вершину\nПКМ / кнопка — завершити перешкоду",
              foreground="#555").pack(anchor="w")
    ttk.Button(left, text="Завершити перешкоду",
               command=lambda: finish_polygon()).pack(fill="x", pady=2)

    ttk.Separator(left, orient="horizontal").pack(fill="x", pady=6)
    ttk.Label(left, text="Авто-генерація", font=("TkDefaultFont", 10, "bold")).pack(anchor="w")
    ttk.Label(left, text="К-сть перешкод:").pack(anchor="w")
    n_var = tk.StringVar(value="300")
    ttk.Entry(left, textvariable=n_var, width=10).pack(anchor="w")
    ttk.Button(left, text="Згенерувати (випадково)",
               command=lambda: do_generate()).pack(fill="x", pady=2)
    ttk.Button(left, text="Найгірший випадок (слалом)",
               command=lambda: do_worst()).pack(fill="x", pady=2)

    ttk.Separator(left, orient="horizontal").pack(fill="x", pady=6)
    ttk.Button(left, text="▶ Знайти найкоротший шлях",
               command=lambda: do_run()).pack(fill="x", pady=2)
    ttk.Label(left, text="(швидкий режим вмикається\nавтоматично при >1000 вершин)",
              foreground="#777").pack(anchor="w")
    show_var = tk.BooleanVar(value=False)
    ttk.Checkbutton(left, text="Показати граф видимості",
                    variable=show_var,
                    command=lambda: toggle_graph()).pack(anchor="w")
    ttk.Button(left, text="Очистити робочу область",
               command=lambda: do_clear()).pack(fill="x", pady=2)

    ttk.Separator(left, orient="horizontal").pack(fill="x", pady=6)
    info = tk.Text(left, width=30, height=12, font=("TkFixedFont", 9),
                   relief="flat", bg="#f4f4f4")
    info.pack(fill="x")

    status_var = tk.StringVar(value="Готово.")
    ttk.Label(left, textvariable=status_var, foreground="#06c",
              wraplength=210).pack(anchor="w", pady=(6, 0))

    def status(msg):
        status_var.set(msg)

    def set_info(lines):
        info.delete("1.0", "end")
        info.insert("1.0", "\n".join(lines))

    # ---- малювання ----------------------------------------------------------
    def redraw():
        canvas.delete("all")
        # перешкоди
        for poly in scene.polygons:
            flat = [c for p in poly for c in p]
            if len(flat) >= 6:
                canvas.create_polygon(*flat, fill="#c9ccd3", outline="#5b5f66")
        # граф видимості: повний (малі сцени) або k-найближчий (великі)
        graph_src = scene.fast_cache if state["graph_knn"] else scene.vis_cache
        if state["show_graph"] and graph_src:
            for i, lst in graph_src.items():
                ax, ay = scene._pos(i)
                for j, _ in lst:
                    if j > i:
                        bx, by = scene._pos(j)
                        canvas.create_line(ax, ay, bx, by, fill="#e7c7c7")
        # поточний полігон, що малюється
        cur = state["current"]
        for (x, y) in cur:
            canvas.create_oval(x - 3, y - 3, x + 3, y + 3, fill="#333")
        for i in range(len(cur) - 1):
            canvas.create_line(*cur[i], *cur[i + 1], fill="#333", dash=(3, 2))
        # шлях
        res = state["last_result"]
        if res and res.get("path"):
            pts = res["path"]
            flat = [c for p in pts for c in p]
            canvas.create_line(*flat, fill="#1565ff", width=3)
        # S, T
        if scene.S:
            x, y = scene.S
            canvas.create_oval(x - 6, y - 6, x + 6, y + 6, fill="#16a34a", outline="")
            canvas.create_text(x, y - 12, text="S", fill="#16a34a",
                               font=("TkDefaultFont", 11, "bold"))
        if scene.T:
            x, y = scene.T
            canvas.create_oval(x - 6, y - 6, x + 6, y + 6, fill="#dc2626", outline="")
            canvas.create_text(x, y - 12, text="T", fill="#dc2626",
                               font=("TkDefaultFont", 11, "bold"))

    def refresh_info():
        nv = len(scene.vertices)
        nh = len(scene.polygons)
        lines = [f"Перешкод (h):   {nh}",
                 f"Вершин (n):     {nv}"]
        res = state["last_result"]
        if res:
            lines.append("")
            lines.append(f"Статус: {res['status']}")
            if res.get("mode"):
                lines.append(f"Режим:  {res['mode']}")
            if res.get("path"):
                lines.append(f"Довжина шляху:  {res['length']:.1f}")
                lines.append(f"Ланок у шляху:  {len(res['path']) - 1}")
            lines.append(f"Вершин розкрито: {res['expanded']}")
            lines.append(f"Час запиту:     {res['time'] * 1000:.1f} мс")
        set_info(lines)

    # ---- обробники подій ----------------------------------------------------
    def on_left(ev):
        x, y = ev.x, ev.y
        m = state["mode"]
        if m == "draw":
            state["current"].append((x, y))
        elif m == "set_s":
            scene.S = (x, y)
            scene.grid = None
            status("S встановлено")
        elif m == "set_t":
            scene.T = (x, y)
            scene.grid = None
            status("T встановлено")
        redraw()

    def finish_polygon():
        cur = state["current"]
        if len(cur) >= 3:
            scene.add_polygon(cur)
            status(f"Додано перешкоду з {len(cur)} вершин")
        else:
            status("Потрібно щонайменше 3 вершини")
        state["current"] = []
        state["last_result"] = None
        redraw()
        refresh_info()

    def on_right(ev):
        finish_polygon()

    def do_generate():
        try:
            n = int(n_var.get())
        except ValueError:
            status("Некоректне число")
            return
        n = max(1, min(n, 6000))
        t0 = time.perf_counter()
        generate_scene(scene, n, canvas.winfo_width() or W,
                       canvas.winfo_height() or H, rng)
        dt = time.perf_counter() - t0
        state["current"] = []
        state["last_result"] = None
        status(f"Згенеровано {n} перешкод ({len(scene.vertices)} вершин) за {dt*1000:.0f} мс")
        redraw()
        refresh_info()

    def do_worst():
        try:
            nw = int(n_var.get())
        except ValueError:
            status("Некоректне число")
            return
        nw = max(2, min(nw, 400))
        generate_worst_case(scene, nw, canvas.winfo_width() or W,
                            canvas.winfo_height() or H)
        state["current"] = []
        state["last_result"] = None
        status(f"Слалом: {nw} стін ({len(scene.vertices)} вершин) — найгірший ввід")
        redraw()
        refresh_info()

    def do_run():
        if scene.S is None or scene.T is None:
            status("Спочатку задайте S і T")
            return
        if scene.grid is None:
            scene.build_grid()
        # швидкий режим вмикається автоматично для великих сцен (>1000 вершин)
        auto_fast = len(scene.vertices) > 1000
        res = scene.shortest_path(scene.S, scene.T, fast=auto_fast)
        state["last_result"] = res
        mode = res.get("mode", "")
        status(f"{res['status']} ({mode}) | {res['time']*1000:.1f} мс | розкрито {res['expanded']}")
        redraw()
        refresh_info()

    # пороги візуалізації графу видимості
    GRAPH_FULL_MAX = 400      # ≤ — повний граф O(n²) (гарантований, але важкий)
    GRAPH_KNN_MAX = 3000      # ≤ — k-найближчий граф O(n·k) (швидкий, читабельний)

    def toggle_graph():
        state["show_graph"] = show_var.get()
        if not state["show_graph"]:
            redraw()
            return
        n = len(scene.vertices)
        if scene.grid is None:
            scene.build_grid()
        if n <= GRAPH_FULL_MAX:
            # малі сцени — повний граф видимості (точний)
            state["graph_knn"] = False
            scene.build_full_visibility()
            status(f"Повний граф видимості ({n} вершин)")
        elif n <= GRAPH_KNN_MAX:
            # великі сцени — k-найближчий граф (той, що бачить швидкий режим A*)
            state["graph_knn"] = True
            t0 = time.perf_counter()
            scene.build_knn_visibility(90)
            dt = (time.perf_counter() - t0) * 1000
            edges = sum(len(l) for l in scene.fast_cache.values()) // 2
            status(f"k-найближчий граф ({n} вершин, {edges} ребер) за {dt:.0f} мс")
        else:
            # надто щільно навіть для k-NN — малювати марно (нечитабельний клубок)
            status(f"Граф видимості вимкнено: {n} вершин > {GRAPH_KNN_MAX} "
                   f"(був би нечитабельний клубок)")
            state["show_graph"] = False
            show_var.set(False)
        redraw()

    def do_clear():
        scene.clear()
        state["current"] = []
        state["last_result"] = None
        status("Очищено")
        redraw()
        refresh_info()

    canvas.bind("<Button-1>", on_left)
    canvas.bind("<Button-3>", on_right)      # ПКМ
    canvas.bind("<Button-2>", on_right)      # середня кнопка / трекпад
    root.bind("<Return>", lambda e: do_run())
    root.bind("<Escape>", lambda e: do_clear())

    refresh_info()
    redraw()
    root.mainloop()


# =============================================================================
#  6. Перевірка коректності та заміри ефективності (без GUI)
# =============================================================================

def selftest():
    print("=== Перевірка коректності ===")
    sc = Scene()
    # один квадрат-перешкода у центрі; S зліва, T справа
    sc.add_polygon([(40, 40), (60, 40), (60, 60), (40, 60)])
    sc.build_grid()
    res = sc.shortest_path((10, 50), (90, 50))
    assert res["status"] == "OK", res
    straight = 80.0
    print(f"  Квадрат: довжина={res['length']:.3f} (пряма={straight}), шлях обходить перешкоду: "
          f"{res['length'] > straight}")
    assert res["length"] > straight, "шлях має бути довшим за пряму крізь перешкоду"
    # точки на одній лінії без перешкод
    sc2 = Scene(); sc2.build_grid()
    r2 = sc2.shortest_path((0, 0), (100, 0))
    print(f"  Без перешкод: довжина={r2['length']:.3f} (очікувано 100.0)")
    assert abs(r2["length"] - 100.0) < 1e-6
    # S всередині перешкоди
    sc3 = Scene()
    sc3.add_polygon([(0, 0), (100, 0), (100, 100), (0, 100)])
    sc3.build_grid()
    r3 = sc3.shortest_path((50, 50), (200, 200))
    print(f"  S всередині: статус='{r3['status']}'")
    assert r3["path"] is None
    print("  Усі перевірки пройдено ✔")


def bench(n_obstacles):
    print(f"=== Заміри ефективності: {n_obstacles} перешкод ===")
    sc = Scene()
    rng = random.Random(7)
    # світ зростає ∝ sqrt(h) — заміри при сталій густині перешкод
    side = int(70 * math.sqrt(n_obstacles))
    generate_scene(sc, n_obstacles, side, side, rng)
    n = len(sc.vertices)
    t0 = time.perf_counter()
    sc.build_grid()
    t_grid = time.perf_counter() - t0
    # перший запит (наповнює кеш видимості — роль передобробки)
    res = sc.shortest_path(sc.S, sc.T)
    # серія повторних запитів для різних S, T (модель задачі)
    times = []
    for _ in range(8):
        S = sc.S
        for _ in range(200):
            p = (rng.uniform(0, side), rng.uniform(0, side))
            if not sc._point_inside_any(p):
                T = p
                break
        r = sc.shortest_path(S, T)
        times.append(r["time"])
    times.sort()
    median = times[len(times) // 2]
    print(f"  h={n_obstacles} перешкод, n={n} вершин")
    print(f"  Сітка (передобробка):   {t_grid*1000:8.1f} мс")
    print(f"  1-й запит (з прогрівом): {res['time']*1000:8.1f} мс  (розкрито {res['expanded']})")
    print(f"  повторні запити, медіана:{median*1000:8.1f} мс  "
          f"(діапазон {times[0]*1000:.0f}..{times[-1]*1000:.0f} мс)")
    if res.get("path"):
        print(f"  Довжина шляху:    {res['length']:.1f}, ланок {len(res['path'])-1}")


def bench_worst(n_walls):
    print(f"=== Найгірший випадок (слалом): {n_walls} стін ===")
    sc = Scene()
    side = max(800, n_walls * 24)
    generate_worst_case(sc, n_walls, side, side)
    n = len(sc.vertices)
    res = sc.shortest_path(sc.S, sc.T)
    print(f"  стін={n_walls}, n={n} вершин")
    print(f"  запит: {res['time']*1000:8.1f} мс  (розкрито {res['expanded']}, "
          f"ланок {len(res['path'])-1 if res.get('path') else 0}, статус {res['status']})")
    if res.get("path"):
        print(f"  довжина шляху: {res['length']:.1f}")


def bench_fast(n_obstacles):
    """Порівняння точного і швидкого режимів на щільній сцені (фіксоване полотно)."""
    print(f"=== Точний vs Швидкий режим: {n_obstacles} перешкод ===")
    sc = Scene()
    rng = random.Random(0)
    generate_scene(sc, n_obstacles, 980, 720, rng)
    n = len(sc.vertices)
    sc.vis_cache.clear(); sc.fast_cache.clear()
    ex = sc.shortest_path(sc.S, sc.T, fast=False)
    fa = sc.shortest_path(sc.S, sc.T, fast=True)
    print(f"  h={n_obstacles}, n={n} вершин")
    print(f"  точний:  {ex['length']:.1f}  за {ex['time']*1000:8.1f} мс")
    print(f"  швидкий: {fa['length']:.1f}  за {fa['time']*1000:8.1f} мс")
    if ex.get("path") and fa.get("path"):
        err = (fa['length'] - ex['length']) / ex['length'] * 100
        print(f"  похибка швидкого: {err:+.3f}% | прискорення: "
              f"{ex['time']/max(fa['time'],1e-9):.0f}x")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--selftest":
        selftest()
    elif len(sys.argv) > 1 and sys.argv[1] == "--bench":
        n = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
        bench(n)
    elif len(sys.argv) > 1 and sys.argv[1] == "--bench-worst":
        n = int(sys.argv[2]) if len(sys.argv) > 2 else 30
        bench_worst(n)
    elif len(sys.argv) > 1 and sys.argv[1] == "--bench-fast":
        n = int(sys.argv[2]) if len(sys.argv) > 2 else 1500
        bench_fast(n)
    else:
        run_gui()
