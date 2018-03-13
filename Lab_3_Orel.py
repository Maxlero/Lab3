"""
    Автор: Орел Максим
    Группа: КБ-161
    Вариант: 11
    Дата создания: 28/02/2018
    Python Version: 3.6
"""
import warnings

import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy.misc import derivative

# Constants
accuracy = 0.00001
START_X = -10
END_X = 10
START_Y = -30
END_Y = 100

# IMPORTANT! Look at your function
a0, a1, a2, a3 = 1, 3, -24, 10


def f(x):
    global a0, a1, a2, a3
    return a0 * x ** 3 + a1 * x ** 2 + a2 * x + a3


def diff_f(point):
    return derivative(f, point, dx=1e-5)


def diff_2_f(point):
    return derivative(f, point, dx=1e-5, n=2)


def is_nan(arg):
    return arg != arg


# Находит первую точку в которой существует функция (слева и справа) на заданом промежутке
def find_func(a, b, precision):
    x = np.arange(a, b, precision)

    ret_x_start = None
    ret_x_end = None

    for xn in x:
        if not is_nan(f(xn)):
            ret_x_start = xn
            break

    for xn in reversed(x):
        if not is_nan(f(xn)):
            ret_x_end = xn
            break

    if ret_x_start is None or is_nan(ret_x_start) or ret_x_end is None or is_nan(ret_x_end):
        raise Exception('На данном промежутке функция не существует')

    return ret_x_start, ret_x_end


def build_function(x_from, x_to, dx, func, x1=START_X, x2=END_X, y1=START_Y, y2=END_Y):
    x = np.arange(x_from, x_to, dx)

    plt.plot(x, func(x))
    plt.axis([x1, x2, y1, y2])
    plt.grid(True)
    plt.axhline(y=0, color='k')
    plt.axvline(x=0, color='k')


def show_plot():
    plt.show()


def enter_values():
    while True:
        try:
            a = float(input("Введите левый промежуток: "))
            b = float(input("Введите правый промежуток: "))

            # Сделаем проверку, если ли корень на промежутке
            if f(a) * f(b) < 0:
                print("Кажется на промежутке [" + str(a) + ", " + str(b) + "] есть один или несколько корней.")
                break
            else:
                print("Кажется на введенном промежутке нет корня. Попробуйте еще раз..")

        except ValueError:
            print("Кажется с промежутком что-то не так. Попробуйте другой..")

    return a, b


def lobachevsky_method(a0, a1, a2, a3, iterations):
    if iterations > 0:

        temp_x1 = (a1 / a0) ** (1 / 2 ** iterations)
        temp_x2 = (a2 / a1) ** (1 / 2 ** iterations)
        temp_x3 = (a3 / a2) ** (1 / 2 ** iterations)

        if (abs(f(temp_x1)) <= accuracy or (f(-temp_x1)) <= accuracy) and (
                abs(f(temp_x2)) <= accuracy or (f(-temp_x2)) <= accuracy) and (
                abs(f(temp_x3)) <= accuracy or (f(-temp_x3)) <= accuracy):

            if abs(f(temp_x1)) >= accuracy:
                temp_x1 = -temp_x1

            if abs(f(temp_x2)) >= accuracy:
                temp_x2 = -temp_x2

            if abs(f(temp_x3)) >= accuracy:
                temp_x3 = -temp_x3

            return [[temp_x1, temp_x2, temp_x3], iterations, [f(temp_x1), f(temp_x2), f(temp_x3)]]

    b0 = a0 ** 2
    b1 = a1 ** 2 - 2 * a0 * a2
    b2 = a2 ** 2 - 2 * a1 * a3
    b3 = a3 ** 2

    return lobachevsky_method(b0, b1, b2, b3, iterations + 1)


def lobachevsky_method_interface():
    global a0, a1, a2, a3, accuracy

    print("___________________________________________________________________________________________________________")
    print("Метод Лобачевского:")

    iterations = 0
    res = lobachevsky_method(a0, a1, a2, a3, iterations)

    temp_x1_e = str("%.2e" % res[0][0])
    temp_x1_f = str("%.6f" % res[0][0])
    temp_x2_e = str("%.2e" % res[0][1])
    temp_x2_f = str("%.6f" % res[0][1])
    temp_x3_e = str("%.2e" % res[0][2])
    temp_x3_f = str("%.6f" % res[0][2])

    temp_iterations = str(res[1])

    temp_f_x1 = str("%.5e" % res[2][0])
    temp_f_x2 = str("%.5e" % res[2][1])
    temp_f_x3 = str("%.5e" % res[2][2])

    print(
        ("Корени уравнения x1 = {0} ({1}), x2 = {2} ({3}), x3 = {4} ({5}) \n" +
         "Найдены за {6} итерации(ий) \n" +
         "Функции в этой точке f(x1) = {7}, f(x2) = {8}, f(x3) = {9} <= {10} (точность)").format(
            temp_x1_e, temp_x1_f,
            temp_x2_e, temp_x2_f,
            temp_x3_e, temp_x3_f,
            temp_iterations,
            temp_f_x1, temp_f_x2, temp_f_x3,
            accuracy
        )
    )


if __name__ == "__main__":
    # Отключение вывода некоторых уведомлений
    if not sys.warnoptions:
        warnings.simplefilter("ignore")

    # Устанавливает максимальную глубину рекурсии на 2000
    sys.setrecursionlimit(2000)

    # нарисуем функцию и её производной на промежутке заданном впараметрах проагрммы сверху
    print("Выведен график функции и график производной")
    build_function(START_X, END_X, 0.01, f)
    build_function(START_X, END_X, 0.01, diff_f)
    show_plot()

    try:
        lobachevsky_method_interface()
    except Exception as e:
        print(e)

    print()

    try:
        lobachevsky_method_interface()
    except Exception as e:
        print(e)
