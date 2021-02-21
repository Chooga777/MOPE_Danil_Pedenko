import numpy as np
import random
from prettytable import PrettyTable

list_a = np.array([17, 19, 1, 18])
matr_x = [[random.randint(0, 20) for j in range(3)] for i in range(8)]
matr_y = [list_a[0] + list_a[1] * matr_x[i][0] + list_a[2] * matr_x[i][1] + list_a[3] * matr_x[i][2] for i in
          range(len(matr_x))]
TX = np.array(matr_x).transpose()
matr_x0 = [float((min(TX[i]) + max(TX[i])) / 2) for i in range(len(TX))]
matr_x0.append((min(matr_y) + max(matr_y)) / 2)
yet = list_a[0] + list_a[1] * matr_x0[0] + list_a[2] * matr_x0[1] + list_a[3] * matr_x0[2]
matr_dx = [max(TX[i]) - matr_x0[i] for i in range(len(TX))]
matr_xn = [[(matr_x[i][j] - matr_x0[j]) / matr_dx[j] for j in range(3)] for i in range(8)]
a = []
c = 0
for i in range(len(matr_y)):
    temp = matr_y[i] - matr_x0[3]
    if temp < 0:
        a.append([])
        a[len(a) - 1].append(temp)
        a[len(a) - 1].append(i)
maximum = a[0][0]
koef = 0
for i in range(len(a)):
    for j in a[i]:
        if j < 0:
            if maximum < j:
                maximum = j
                koef = i
ta = PrettyTable()
ta.field_names = ["#", "X1", "X2", "X3", "Y", " ", "XN1", "XN2", "XN3"]
ta.add_rows(
    [
        ["1", matr_x[0][0], matr_x[0][1], matr_x[0][2], matr_y[0], " ", matr_xn[0][0], matr_xn[0][1], matr_xn[0][2]],
        ["2", matr_x[1][0], matr_x[1][1], matr_x[1][2], matr_y[1], " ", matr_xn[1][0], matr_xn[1][1], matr_xn[1][2]],
        ["3", matr_x[2][0], matr_x[2][1], matr_x[2][2], matr_y[2], " ", matr_xn[2][0], matr_xn[2][1], matr_xn[2][2]],
        ["4", matr_x[3][0], matr_x[3][1], matr_x[3][2], matr_y[3], " ", matr_xn[3][0], matr_xn[3][1], matr_xn[3][2]],
        ["5", matr_x[4][0], matr_x[4][1], matr_x[4][2], matr_y[4], " ", matr_xn[4][0], matr_xn[4][1], matr_xn[4][2]],
        ["6", matr_x[5][0], matr_x[5][1], matr_x[5][2], matr_y[5], " ", matr_xn[5][0], matr_xn[5][1], matr_xn[5][2]],
        ["7", matr_x[6][0], matr_x[6][1], matr_x[6][2], matr_y[6], " ", matr_xn[6][0], matr_xn[6][1], matr_xn[6][2]],
        ["8", matr_x[7][0], matr_x[7][1], matr_x[7][2], matr_y[7], " ", matr_xn[7][0], matr_xn[7][1], matr_xn[7][2]],
        ["X0", matr_x0[0], matr_x0[1], matr_x0[2], matr_x0[3], " ", " ", " ", " "],
        ["dx", matr_dx[0], matr_dx[1], matr_dx[2], " ", " ", " ", " ", " "]
    ]
)
print(ta)
print(
    "Початкові коефіцієнти: a0 = {0}, a1 = {1}, a2 = {2}, a3 = {3}".format(list_a[0], list_a[1], list_a[2], list_a[3]))

print("Сгенерована матриця X\n", np.array(matr_x))
print("Відповідні значення Y\n", np.array(matr_y))
print("Отриманні значення X0\n", np.array(matr_x0))
print("Y еталонне = ", yet)
print("Інтервал зміни фактора dx\n", np.array(matr_dx))
print("Нормовані значення XN\n", np.array(matr_xn))
print("Точка плану, що задовольняє критерію вибору\n Y = ", matr_y[a[koef][1]])
