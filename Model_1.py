from py_linq import Enumerable
from scipy import constants
from scipy.optimize import fsolve
import math
import matplotlib.pyplot as plt
import numpy as np


def distance(x1, x2, y1, y2):
    res = math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
    return res


def intersect(x1, x2, y1, y2, eps):
    if len(x1) != len(y1) or len(x2) != len(y2) or len(x2) != len(x1):
        assert ("len x and len y should be equal")

    points = []
    for i in range(len(x1)):
        dist = distance(x1[i], x2[i], y1[i], y2[i])
        if dist <= eps:
            points.append((x1[i], y1[i], dist))

    points = Enumerable(points)
    min_dist = points.select(lambda x: x[2]).min()
    min_point = points.where(lambda x: x[2] == min_dist).first()
    if min_point is None:
        assert ('No intersection found')

    return (min_point[0], min_point[1])


def print_points_table(points):
    print('\begin{tabular}{| c | c | c |}')
    print('\hline\nn  & x  & y \\')
    for i in range(len(points)):
        print(f'\hline\n{i + 1} & {points[i][0]:.3f}  & {points[i][1]:.3f} \\')
    print('\hline\n\end{tabular}')


def gorgeous_scientific_print(x):
    for i in x:
        print(f'{i:.3e}', end=" ")


a = 2e-10
U = 15.168e-18
omega = 1
m = constants.electron_mass
eps = int(1e7)

def V(x: float) -> float:
  if(abs(x) <= a):
    return -U
  else:
    return 0

def Vo(x: float ) -> float:
  return 1/2 * m * omega ** 2 * x ** 2

xMin = -a - a
xMax = a + a
x = np.arange(xMin, xMax, 1e-12)

y = []

for i in range(len(x)):
    y.append(V(x[i]))

plt.plot(x, y)
plt.show()

k_2 = math.sqrt(2 * m * U) / constants.hbar

plt.figure(figsize=(10, 7))

print(f"max k_2: {k_2}")

lefty = np.arange(0, k_2, eps, dtype=float) * a
leftx = np.arange(0, k_2, eps)
n = math.ceil(k_2 * a / constants.pi)

plt.plot(leftx, lefty)

res = []

for i in range(1, n + 1):
    rightx = []
    righty = []
    for j in range(0, int(k_2), eps):
        rightx.append(j)
        righty.append(
            constants.pi * i - 2 * math.asin((constants.hbar * j) / math.sqrt(2 * constants.electron_mass * U)))

    plt.plot(rightx, righty, color='orange')
    point = intersect(leftx, rightx, lefty, righty, eps)
    plt.plot(point[0], point[1], 'bo', color='purple')
    res.append((point[0], point[1]))

plt.show()
print_points_table(res)

plt.figure(figsize=(10, 7))
plt.plot(x, y)
energy = [0] * len(res)


def psi(x: float, n: int, scale: float = 1e23) -> float:
    if abs(x) > a:
        return 0
    return math.sqrt(2 / a) * math.sin(constants.pi * (n + 1) * x / (2 * a) + (n + 1) * constants.pi / 2) / scale


for i in range(len(res)):
    energy[i] = constants.hbar ** 2 * res[i][0] ** 2 / (2 * m)
    plt.axhline(y=energy[i] - U, linestyle="--", linewidth=0.5)

    y_psi = []
    x_psi = []
    for j in range(len(x)):
        tmp = psi(x[j], i)
        if tmp != 0:
            y_psi.append(tmp + (energy[i] - U))
            x_psi.append(x[j])

    plt.plot(x_psi, y_psi, color='orange')

plt.show()
gorgeous_scientific_print(energy)


def hermite(xi: float, n: int) -> float:
    array = [0] * (n) + [1]
    hermite = np.polynomial.hermite.Hermite(array)
    return hermite(xi)


def psi_osc(x, E: float, n: int):
    points = []
    for i in range(len(x)):
        xi = x[i] * math.sqrt(constants.electron_mass * omega / constants.hbar)
        coef = 1 / (math.sqrt(2 ** n * math.factorial(n))) * (
                    constants.electron_mass * omega / (constants.hbar * constants.pi)) ** 0.25
        point = coef * constants.electron_mass / 15e4 * math.exp(-0.5 * xi ** 2) * hermite(xi, n) + E
        points.append(point)

    return points


plt.figure(figsize=(10, 7))
eps = 1e-4
x2 = np.arange(-3e-2, 3e-2, eps)
y2 = []

for i in range(len(x2)):
    y2.append(Vo(x2[i]))

plt.plot(x2, y2)

res = []

n = 0
while True:
    E = (n + 0.5) * constants.hbar * omega

    if (E > y2[0]):
        break

    res.append(E)

    plt.plot(x2, [E] * len(x2), linestyle="--", color='pink', linewidth=0.5)
    plt.plot(x2, psi_osc(x2, E, n), color='orange')
    n += 1

plt.show()
gorgeous_scientific_print(res)