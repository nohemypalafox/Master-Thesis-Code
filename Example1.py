import numpy as np
import matplotlib.pyplot as plt

t_0 = 0
t_f = 1
n_max = 1000

t = np.linspace(t_0,t_f, n_max)
h = t[1] - t[0]
x_0 = 1

x = np.zeros(t.shape[0])
x[0] = x_0

for j in np.arange(n_max -1):

    k_1 = - 0.5 * x[j] ** 2
    k_2 = - 0.5 * (x[j] + 0.5 * h * k_1) ** 2
    k_3 = - 0.5 * (x[j] + 0.5 * h * k_2) ** 2
    k_4 = - 0.5 * (x[j] + h * k_3) ** 2

    x[j+1] = x[j] + h / 6.0 * (k_1 + 2 * k_2 + 2 * k_3 + k_4)

plt.plot(t,x)
plt.show()
