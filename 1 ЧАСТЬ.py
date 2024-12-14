import numpy as np
import matplotlib.pyplot as plt
import time
import math

a = 43200   
e = 0.736
T = 26121  

t = np.linspace(0, T, 1000)
M = t / T * 2 * np.pi

def newton_method(M, e):

    E = M
    for _ in range(100):
        E_new = E + (M - (E - e * np.sin(E))) / (1 - e * np.cos(E))
        if abs(E_new - E) < 1e-6:
            break
        E = E_new
    return E

def golden_sel_metod(M, e, tol=1e-6):
    def f(E):
        return E - e * np.sin(E) - M
    a, b = 0, 2 * np.pi
    phi = (1 + np.sqrt(5)) / 2
    while (b - a) > tol:
        c = a + (b - a) / phi
        if f(c) == 0:
            return c
        elif f(a) * f(c) < 0:
            b = c
        else:
            a = c
    return (a + b) / 2

def bisect_method(M, e, tol=1e-6):
    def f(E):
        return E - e * np.sin(E) - M
    a, b = 0, 2 * np.pi
    while (b - a) > tol:
        c = (a + b) / 2
        if f(c) == 0:
            return c
        elif f(a) * f(c) < 0:
            b = c
        else:
            a = c
    return (a + b) / 2

def iterations_method(M, e, tol=1e-6):
    E = M
    E_new = e * np.sin(E) + M
    while E_new - E > tol:
        E = E_new
        E_new = e * np.sin(E) + M
    return E_new

methods = {
    "Ньютона": newton_method,
    "Золотого сечения": golden_sel_metod,
    "Половинного деления": bisect_method,
    "Итераций": iterations_method,
}

results = {}

for name, method in methods.items():
    start_time = time.time()
    E = np.array([method(m, e) for m in M])
    end_time = time.time()
    results[name] = E

for name, E in results.items():
    print(f"Метод {name}: E = {E[-1]:.8f}")

def true_anomaly(E, e):
    return 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))

nu_newton = np.array([true_anomaly(newton_method(m, e), e) for m in M])  
nu_gs = np.array([true_anomaly(golden_sel_metod(m, e), e) for m in M])  
nu_bisection = np.array([true_anomaly(bisect_method(m, e), e) for m in M]) 
nu_iterations = np.array([true_anomaly(iterations_method(m, e), e) for m in M])  

fig, axs = plt.subplots(2, 2, figsize=(15, 12))

# График для метода Ньютона
axs[0, 0].plot(t, M, label='Средняя аномалия M(t)', color='blue')
axs[0, 0].plot(t, results["Ньютона"], label='Эксцентрическая аномалия E(t)', color='red')
axs[0, 0].plot(t, nu_newton, label='Истинная аномалия ν(t)', color='green')
axs[0, 0].set_title('Метод Ньютона')
axs[0, 0].set_xlabel('Время (с)')
axs[0, 0].set_ylabel('Аномалия')
axs[0, 0].legend()
axs[0, 0].grid()

# График для метода золотого сечения
axs[0, 1].plot(t, M, label='Средняя аномалия M(t)', color='blue')
axs[0, 1].plot(t, results["Золотого сечения"], label='Эксцентрическая аномалия E(t)', color='red')
axs[0, 1].plot(t, nu_gs, label='Истинная аномалия ν(t)', color='green')
axs[0, 1].set_title('Метод золотого сечения')
axs[0, 1].set_xlabel('Время (с)')
axs[0, 1].set_ylabel('Аномалия')
axs[0, 1].legend()
axs[0, 1].grid()

# График для метода половинного деления
axs[1, 0].plot(t, M, label='Средняя аномалия M(t)', color='blue')
axs[1, 0].plot(t, results["Половинного деления"], label='Эксцентрическая аномалия E(t)', color='red')
axs[1, 0].plot(t, nu_bisection, label='Истинная аномалия ν(t)', color='green')
axs[1, 0].set_title('Метод половинного деления')
axs[1, 0].set_xlabel('Время (с)')
axs[1, 0].set_ylabel('Аномалия')
axs[1, 0].legend()
axs[1, 0].grid()

# График для метода итераций
axs[1, 1].plot(t, M, label='Средняя аномалия M(t)', color='blue')
axs[1, 1].plot(t, results["Итераций"], label='Эксцентрическая аномалия E(t)', color='red')
axs[1, 1].plot(t, nu_iterations, label='Истинная аномалия ν(t)', color='green')
axs[1, 1].set_title('Метод итераций')
axs[1, 1].set_xlabel('Время (с)')
axs[1, 1].set_ylabel('Аномалия')
axs[1, 1].legend()
axs[1, 1].grid()

plt.tight_layout()
plt.show()
