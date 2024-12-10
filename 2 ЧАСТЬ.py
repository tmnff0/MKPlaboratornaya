import numpy as np
import matplotlib.pyplot as plt

a = 3696  
e = 0.0088  
T = 6720  
mu = 42800 

def anomaly(t, T):

    n = 2 * np.pi / T
    M = n * t
    return M

def newton_method(M, e):

    E = M
    for _ in range(100):
        E_new = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
        if abs(E_new - E) < 1e-6:
            break
        E = E_new
    return E

def ecc_anomaly(M, e):

    E = np.array([newton_method(m, e) for m in M])
    return E

def true_anomaly(E, e):

    beta = np.sqrt((1 + e) / (1 - e))
    nu = 2 * np.arctan(beta * np.tan(E / 2))
    return nu

def rad_vector(nu, a, e):

    p = a * (1 - e**2)
    r = p / (1 + e * np.cos(nu))
    return r

def radial_vel(nu, a, e, mu):

    p = a * (1 - e**2)
    Vr = np.sqrt(mu / p) * e * np.sin(nu)
    return Vr

def transvers_vel(nu, a, e, mu):

    p = a * (1 - e**2)
    Vn = np.sqrt(mu / p) * (1 + e * np.cos(nu))
    return Vn

def speed(Vn, Vr):

    V = np.sqrt(Vn**2 + Vr**2)
    return V

def graphs(t, r, Vr, Vn, V):

    plt.figure(figsize=(10, 10))

    # График радиус-вектора
    plt.subplot(2, 2, 1)
    plt.plot(t, r, label='r (радиус-вектор)', color='blue')
    plt.title('Радиус-вектор')
    plt.xlabel('Время t (в секундах)')
    plt.ylabel('r (м)')
    plt.grid()
    plt.legend()

    # График радиальной скорости
    plt.subplot(2, 2, 2)
    plt.plot(t, Vr, label='Vr(t) (радиальная скорость)', color='red')
    plt.title('Радиальная скорость')
    plt.xlabel('Время t (в секундах)')
    plt.ylabel('Vr (м/с)')
    plt.grid()
    plt.legend()

    # График трансверсальной скорости
    plt.subplot(2, 2, 3)
    plt.plot(t, Vn, label='Vn(t) (трансверсальная скорость)', color='orange')
    plt.title('Трансверсальная скорость')
    plt.xlabel('Время t (в секундах)')
    plt.ylabel('Vn (м/с)')
    plt.grid()
    plt.legend()

    # График модуля скорости
    plt.subplot(2, 2, 4)
    plt.plot(t, V, label='V (модуль скорости)', color='pink')
    plt.title('Модуль скорости')
    plt.xlabel('Время t (в секундах)')
    plt.ylabel('V (м/с)')
    plt.grid()
    plt.legend()

    plt.tight_layout()
    plt.show()

def main():

    t = np.linspace(0, T, 1000)
    M = anomaly(t, T)
    E = ecc_anomaly(M, e)
    nu = true_anomaly(E, e)
    r = rad_vector(nu, a, e)
    Vr = radial_vel(nu, a, e, mu)
    Vn = transvers_vel(nu, a, e, mu)
    V = speed(Vn, Vr)
    graphs(t / (24 * 3600), r, Vr, Vn, V)  

if __name__ == "__main__":
    main()

