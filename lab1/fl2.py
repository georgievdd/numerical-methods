import numpy as np
import matplotlib.pyplot as plt

def iter_up(A, b, x, t, eps, iter_max):
    n = len(b)
    t1 = 1 - t
    it = 0
    norm_old = 0

    for k in range(iter_max):
        it += 1
        norm = 0

        for i in range(n):
            buf = t1 * x[i] + t * b[i] / A[i, i]
            for j in range(n):
                if i != j:
                    buf -= t * A[i, j] * x[j] / A[i, i]
            
            norm = max(norm, abs(buf - x[i]))
            x[i] = buf

        if norm <= eps:
            return x, it, norm

        if norm > norm_old:
            if it > 3:  # Метод расходится
                return None, it, norm
        norm_old = norm

    return x, it, norm

def var_M(B):
    n = B.shape[0]
    A = np.dot(B.T, B)
    return A

def main():
    # Пример системы
    n = 5
    B = np.random.rand(n, n + 1) - 0.5  # Генерация случайной системы
    A = var_M(B[:, :-1])
    b = B[:, -1]

    # Начальные параметры
    iter_max = 100
    eps = 1e-5
    dt = 0.05

    # Массивы для значений t и числа итераций
    t_values = np.arange(0.05, 2.05, dt)
    iter_values = []

    best_t = None
    best_iter = iter_max

    # Поиск оптимального t
    for t in t_values:
        x = np.zeros(n)
        _, it, _ = iter_up(A, b, x, t, eps, iter_max)
        iter_values.append(it)

        if it < best_iter:
            best_iter = it
            best_t = t

    # Подготовка текста для вывода в консоль
    unit_solution = np.zeros(n)
    print(f"Матрица A:\n{A}\n")
    print(f"Вектор b:\n{b}\n")
    print(f"Параметры:\ncode: 0\neps: {eps}\nt: {best_t}\n")
    print(f"Единичное решение:\n{unit_solution}\n")
    print(f"Оптимальные значения: t = {best_t:.2f}, итераций = {best_iter}\n")
    print(f"Массив значений t: {', '.join([f'{t:.2f}' for t in t_values])}\n")
    print(f"Массив значений числа итераций: {', '.join([str(it) for it in iter_values])}\n")

    # Построение графика
    plt.figure(figsize=(10, 6))
    plt.plot(t_values, iter_values, label="Зависимость итераций от t", color='b')
    
    # Отметим оптимальное значение
    plt.axvline(best_t, color='r', linestyle='--', label=f'Оптимальное t = {best_t:.2f}')
    plt.axhline(best_iter, color='g', linestyle='--', label=f'Оптимальное число итераций = {best_iter}')

    # Подписи
    plt.xlabel("Параметр t")
    plt.ylabel("Число итераций")
    plt.title("Зависимость числа итераций от параметра t")
    plt.legend()

    # Показать график
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
