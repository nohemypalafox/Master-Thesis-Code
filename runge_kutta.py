import matplotlib.pyplot as plt
import numpy as np

def f(x):
    return - 0.5 * x ** 2

def runge_kutta(f, x_0, h, n_max):
    sol = np.zeros(n_max)
    sol[0] = x_0

    for j in np.arange(n_max - 1):
        x_j = sol[j]
        k_1 = f(x_j)
        k_2 = f(x_j + 0.5 * h * k_1)
        k_3 = f(x_j + 0.5 * h * k_2)
        k_4 = f(x_j + h * k_3)

        sol[j + 1] = x_j + h / 6.0 * (k_1 + 2 * k_2 + 2 * k_3 + k_4)

    return sol

t_0 = 0
t_f = 1
n_max = 1000

t = np.linspace(t_0,t_f, n_max)
h = t[1] - t[0]
x_0 = 1

x = runge_kutta(f, x_0, h, n_max)

plt.plot(t,x)
plt.show()
