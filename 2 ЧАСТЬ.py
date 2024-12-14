import numpy as np
import matplotlib.pyplot as plt

mu = 398600  #гравитационный параметр
e = 0.736  #эксцентриситет
a = 43200          #большая полуось, км
T = 26121 #период, с

# Временные точки
time = np.linspace(0, T, 1000)

# Вычисление аномалий
M = 2 * np.pi * time / T
E = M
for _ in range(10):
    E = M + e * np.sin(E)
theta = 2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(E / 2))

# Характеристики орбиты
r = a * (1 - e**2) / (1 + e * np.cos(theta))
vr = np.sqrt(mu / a) * e * np.sin(theta) / np.sqrt(1 - e**2)
vt = np.sqrt(mu / a) * (1 + e * np.cos(theta)) / np.sqrt(1 - e**2)
v = np.sqrt(vr**2 + vt**2)

# Построение графиков
plt.figure(figsize=(10, 8))

plt.subplot(2, 2, 1)
plt.plot(time, r)
plt.title("Радиус-вектор от времени")
plt.xlabel("Время, с")
plt.ylabel("Радиус-вектор, км")

plt.subplot(2, 2, 2)
plt.plot(time, vr)
plt.title("Радиальная скорость от времени")
plt.xlabel("Время, с")
plt.ylabel("Радиальная скорость, км/с")

plt.subplot(2, 2, 3)
plt.plot(time, vt)
plt.title("Трансверсальная скорость от времени")
plt.xlabel("Время, с")
plt.ylabel("Трансверсальная скорость, км/с")

plt.subplot(2, 2, 4)
plt.plot(time, v)
plt.title("Модуль скорости от времени")
plt.xlabel("Время, с")
plt.ylabel("Модуль скорости, км/с")

plt.tight_layout()
plt.show()
