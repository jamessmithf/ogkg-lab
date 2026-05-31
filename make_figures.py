#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Генерація ілюстрацій для звіту (PIL). Запуск: python3.14 make_figures.py
Сцени та шляхи обчислюються тим самим алгоритмом, що й у програмі."""

import random
from PIL import Image, ImageDraw, ImageFont
import shortest_path_2d as m

OUT = "report/img"

def font(sz, bold=False):
    for p in ("/System/Library/Fonts/Helvetica.ttc",
              "/System/Library/Fonts/Supplemental/Arial.ttf",
              "/Library/Fonts/Arial.ttf"):
        try:
            return ImageFont.truetype(p, sz)
        except Exception:
            pass
    return ImageFont.load_default()

GRAY = (201, 204, 211)
EDGE = (91, 95, 102)
PATHC = (21, 101, 255)
SC = (22, 163, 74)
TC = (220, 38, 38)

def draw_scene(d, scene, res=None, show_graph=False, off=(0, 0)):
    ox, oy = off
    if show_graph and scene.vis_cache:
        for i, lst in scene.vis_cache.items():
            ax, ay = scene._pos(i)
            for j, _ in lst:
                if j > i:
                    bx, by = scene._pos(j)
                    d.line([ox+ax, oy+ay, ox+bx, oy+by], fill=(231, 199, 199))
    for poly in scene.polygons:
        d.polygon([(ox+x, oy+y) for x, y in poly], fill=GRAY, outline=EDGE)
    if res and res.get("path"):
        pts = [(ox+x, oy+y) for x, y in res["path"]]
        d.line(pts, fill=PATHC, width=3, joint="curve")
    if scene.S:
        x, y = scene.S
        d.ellipse([ox+x-6, oy+y-6, ox+x+6, oy+y+6], fill=SC)
        d.text((ox+x-4, oy+y-22), "S", fill=SC, font=font(15))
    if scene.T:
        x, y = scene.T
        d.ellipse([ox+x-6, oy+y-6, ox+x+6, oy+y+6], fill=TC)
        d.text((ox+x-4, oy+y-22), "T", fill=TC, font=font(15))


def fig_visgraph():
    """Малий приклад: граф видимості + найкоротший шлях."""
    sc = m.Scene()
    sc.add_polygon([(150, 90), (230, 70), (250, 150), (180, 170)])
    sc.add_polygon([(330, 200), (420, 190), (430, 290), (350, 300)])
    sc.add_polygon([(180, 280), (260, 270), (250, 360), (170, 350)])
    sc.S = (60, 120); sc.T = (470, 350)
    sc.build_grid(); sc.build_full_visibility()
    res = sc.shortest_path(sc.S, sc.T)
    img = Image.new("RGB", (540, 420), "white")
    d = ImageDraw.Draw(img)
    draw_scene(d, sc, res, show_graph=True)
    img.save(f"{OUT}/fig_visgraph.png")
    print("fig_visgraph:", res["status"], round(res["length"], 1))


def fig_tangent():
    """Шлях огинає перешкоду, торкаючись опуклих вершин (пунктир — пряма SТ)."""
    sc = m.Scene()
    sc.add_polygon([(200, 80), (330, 120), (300, 300), (170, 260)])
    sc.S = (60, 320); sc.T = (470, 90)
    sc.build_grid()
    res = sc.shortest_path(sc.S, sc.T)
    img = Image.new("RGB", (540, 380), "white")
    d = ImageDraw.Draw(img)
    # пряма S-T (заблокована)
    for t in range(0, 100, 4):
        x1 = sc.S[0] + (sc.T[0]-sc.S[0])*t/100
        y1 = sc.S[1] + (sc.T[1]-sc.S[1])*t/100
        x2 = sc.S[0] + (sc.T[0]-sc.S[0])*(t+2)/100
        y2 = sc.S[1] + (sc.T[1]-sc.S[1])*(t+2)/100
        d.line([x1, y1, x2, y2], fill=(170, 170, 170), width=1)
    draw_scene(d, sc, res)
    img.save(f"{OUT}/fig_tangent.png")
    print("fig_tangent:", res["status"], round(res["length"], 1))


def fig_interface():
    """Макет інтерфейсу: ліва панель + робоча область зі шляхом."""
    rng = random.Random(3)
    sc = m.Scene()
    m.generate_scene(sc, 40, 640, 470, rng, rmin=14, rmax=34)
    res = sc.shortest_path(sc.S, sc.T)
    W, H = 870, 500
    img = Image.new("RGB", (W, H), (236, 236, 236))
    d = ImageDraw.Draw(img)
    # ліва панель
    d.rectangle([0, 0, 210, H], fill=(236, 236, 236))
    fb = font(13, True); f = font(12)
    y = 14
    d.text((10, y), "Режим вводу", fill=(20, 20, 20), font=fb); y += 26
    for t in ("(*) Малювати перешкоду", "( ) Поставити S (старт)",
              "( ) Поставити T (фініш)"):
        d.text((14, y), t, fill=(40, 40, 40), font=f); y += 22
    y += 8; d.line([8, y, 202, y], fill=(180, 180, 180)); y += 10
    d.text((10, y), "Авто-генерація", fill=(20, 20, 20), font=fb); y += 24
    d.text((14, y), "К-сть перешкод: [ 300 ]", fill=(40, 40, 40), font=f); y += 26
    for t in ("[ Згенерувати ]", "[ > Знайти шлях ]",
              "[ ] Граф видимості", "[ Очистити ]"):
        d.rectangle([14, y-2, 196, y+18], outline=(170, 170, 170))
        d.text((20, y), t, fill=(40, 40, 40), font=f); y += 26
    y += 6; d.line([8, y, 202, y], fill=(180, 180, 180)); y += 8
    for t in (f"Перешкод (h):   {len(sc.polygons)}",
              f"Вершин (n):     {len(sc.vertices)}",
              f"Статус: {res['status']}",
              f"Довжина: {res['length']:.0f}",
              f"Розкрито: {res['expanded']}",
              f"Час: {res['time']*1000:.0f} мс"):
        d.text((12, y), t, fill=(60, 60, 60), font=font(11)); y += 18
    # робоча область
    d.rectangle([216, 8, W-8, H-8], fill="white", outline=(140, 140, 140))
    draw_scene(d, sc, res, off=(216, 8))
    img.save(f"{OUT}/fig_interface.png")
    print("fig_interface:", res["status"], round(res["length"], 1))


def fig_large():
    """Велика сцена (>1000 вершин) з обчисленим шляхом."""
    rng = random.Random(11)
    sc = m.Scene()
    W, Hh = 760, 560
    m.generate_scene(sc, 320, W, Hh, rng, rmin=8, rmax=18)

    def free_near(tx, ty):
        best, bd = None, 1e18
        for _ in range(4000):
            p = (rng.uniform(10, W - 10), rng.uniform(10, Hh - 10))
            if sc._point_inside_any(p):
                continue
            dd = (p[0] - tx) ** 2 + (p[1] - ty) ** 2
            if dd < bd:
                best, bd = p, dd
        return best
    sc.S = free_near(30, Hh - 30)        # лівий нижній кут
    sc.T = free_near(W - 30, 30)         # правий верхній кут
    sc.grid = None
    res = sc.shortest_path(sc.S, sc.T)
    img = Image.new("RGB", (760, 560), "white")
    d = ImageDraw.Draw(img)
    draw_scene(d, sc, res)
    d.text((10, 10), f"h={len(sc.polygons)} перешкод, n={len(sc.vertices)} вершин, "
                     f"шлях={res['length']:.0f}, {res['time']*1000:.0f} мс",
           fill=(20, 20, 20), font=font(15))
    img.save(f"{OUT}/fig_large.png")
    print("fig_large:", res["status"], len(sc.vertices), "вершин",
          round(res["time"]*1000), "мс")


def fig_worst():
    """Найгірший ввід — слалом у замкненій рамці: шлях змушений зигзагом
    проходити крізь усі проходи."""
    sc = m.Scene()
    m.generate_worst_case(sc, 14, 640, 460)
    res = sc.shortest_path(sc.S, sc.T)
    img = Image.new("RGB", (640, 460), "white")
    d = ImageDraw.Draw(img)
    draw_scene(d, sc, res)
    img.save(f"{OUT}/fig_worst.png")
    print("fig_worst:", res["status"], "links", len(res["path"]) - 1)


if __name__ == "__main__":
    import os
    os.makedirs(OUT, exist_ok=True)
    fig_visgraph()
    fig_tangent()
    fig_interface()
    fig_large()
    fig_worst()
    print("Готово.")
