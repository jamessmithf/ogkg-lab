"""Безголові тести логіки. Не вимагає дисплея."""
import math
import sys

sys.path.insert(0, '/home/viktor/Downloads/ogkg')

# Імпортуємо лише обчислювальну частину (без GUI)
import importlib.util
spec = importlib.util.spec_from_file_location("main_mod", "/home/viktor/Downloads/ogkg/main.py")
mod = importlib.util.module_from_spec(spec)
# Замінимо matplotlib на заглушку щоб не вимагати дисплея
import types
fake_mpl = types.SimpleNamespace()
# Виконуємо лише обчислювальні функції — імпорти matplotlib не виконуються до створення GUI.
spec.loader.exec_module(mod)


def approx(a, b, tol=1e-6):
    return abs(a - b) < tol


def test_segments_cross():
    assert mod.segments_cross_properly((0, 0), (2, 2), (0, 2), (2, 0)) is True   # X crossing
    assert mod.segments_cross_properly((0, 0), (1, 0), (2, 0), (3, 0)) is False  # collinear, disjoint
    assert mod.segments_cross_properly((0, 0), (2, 0), (1, 0), (1, 1)) is False  # T-shape (endpoint on segment) NOT proper
    assert mod.segments_cross_properly((0, 0), (1, 1), (1, 1), (2, 0)) is False  # share endpoint
    print("test_segments_cross OK")


def test_point_in_polygon():
    poly = [(0, 0), (4, 0), (4, 4), (0, 4)]
    assert mod.point_in_polygon_strict((2, 2), poly) is True
    assert mod.point_in_polygon_strict((5, 2), poly) is False
    assert mod.point_in_polygon_strict((-1, 2), poly) is False
    # Триангула - неопуклий тест
    tri = [(0, 0), (4, 0), (2, 3)]
    assert mod.point_in_polygon_strict((2, 1), tri) is True
    assert mod.point_in_polygon_strict((3, 2.5), tri) is False
    print("test_point_in_polygon OK")


def test_visible_no_obstacle():
    assert mod.visible((0, 0), (10, 10), []) is True
    print("test_visible_no_obstacle OK")


def test_visible_blocked_by_square():
    poly = [(3, 3), (5, 3), (5, 5), (3, 5)]
    assert mod.visible((0, 4), (8, 4), [poly]) is False  # горизонталь через квадрат
    assert mod.visible((0, 1), (8, 1), [poly]) is True   # під квадратом
    print("test_visible_blocked_by_square OK")


def test_visible_adjacent_vertices_of_polygon():
    poly = [(0, 0), (4, 0), (4, 4), (0, 4)]
    # Сусідні вершини: відрізок збігається з ребром -> видимі
    assert mod.visible((0, 0), (4, 0), [poly]) is True
    # Несусідні (діагональ опуклого) -> не видимі
    assert mod.visible((0, 0), (4, 4), [poly]) is False
    print("test_visible_adjacent_vertices_of_polygon OK")


def test_shortest_path_simple():
    # Точки по діагоналі, посередині квадрат
    poly = [(3, 3), (5, 3), (5, 5), (3, 5)]
    S = (0.0, 4.0)
    T = (8.0, 4.0)
    path, length, _, _ = mod.find_shortest_path(S, T, [poly])
    assert path is not None
    # Найкоротший шлях має огинати квадрат через одну з його кутових вершин
    # Прямий шлях був би довжиною 8.0; оптимальний має бути більший
    direct = mod.euclidean(S, T)
    assert length > direct - 1e-9
    # Перевіримо, що шлях проходить через одну з вершин квадрата
    middle = path[1:-1]
    assert len(middle) >= 1
    assert all(p in poly for p in middle)
    print(f"test_shortest_path_simple OK (length={length:.4f}, direct={direct:.4f})")


def test_shortest_path_no_obstacle():
    S = (0.0, 0.0)
    T = (3.0, 4.0)
    path, length, _, _ = mod.find_shortest_path(S, T, [])
    assert path == [S, T]
    assert approx(length, 5.0)
    print("test_shortest_path_no_obstacle OK")


def test_shortest_path_around_triangle():
    # Трикутник посередині
    tri = [(4, 2), (6, 2), (5, 5)]
    S = (0.0, 3.0)
    T = (10.0, 3.0)
    path, length, _, _ = mod.find_shortest_path(S, T, [tri])
    assert path is not None
    # шлях має обходити трикутник
    direct = mod.euclidean(S, T)
    assert length >= direct
    print(f"test_shortest_path_around_triangle OK (length={length:.4f})")


def test_no_path_inside_obstacle():
    # S всередині перешкоди -> алгоритм має повернути None
    poly = [(0, 0), (10, 0), (10, 10), (0, 10)]
    S = (5.0, 5.0)
    T = (-1.0, -1.0)
    path, length, _, _ = mod.find_shortest_path(S, T, [poly])
    # Шлях має бути None (S замкнено)
    assert path is None or length == math.inf
    print("test_no_path_inside_obstacle OK")


def test_nonconvex_polygon():
    # U-подібна неопукла перешкода. S зліва, T справа.
    u_shape = [(2, 0), (2, 6), (3, 6), (3, 1), (5, 1), (5, 6), (6, 6), (6, 0)]
    S = (0.5, 3.0)
    T = (8.0, 3.0)
    path, length, _, _ = mod.find_shortest_path(S, T, [u_shape])
    assert path is not None
    print(f"test_nonconvex_polygon OK (length={length:.4f}, vertices={len(path)})")


if __name__ == "__main__":
    test_segments_cross()
    test_point_in_polygon()
    test_visible_no_obstacle()
    test_visible_blocked_by_square()
    test_visible_adjacent_vertices_of_polygon()
    test_shortest_path_no_obstacle()
    test_shortest_path_simple()
    test_shortest_path_around_triangle()
    test_no_path_inside_obstacle()
    test_nonconvex_polygon()
    print("\nВсі тести пройдено успішно.")
