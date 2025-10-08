import math

import sympy as smp
from pick import pick
from sympy.parsing.sympy_parser import implicit_multiplication, standard_transformations, convert_xor


class IteracoesExcedidas(Exception):
    pass


class SemRaizNoIntervalo(Exception):
    pass


def verificar_intervalo(f, a, b):
    if f(a) * f(b) > 0:
        raise SemRaizNoIntervalo


def secante(f, r0, r1, tol, nMax):
    f0 = f(r0)
    f1 = f(r1)
    if math.isclose(f0, 0, abs_tol=tol):
        return r0

    for k in range(0, nMax):
        if f1 == f0:
            raise ZeroDivisionError

        r = r1 - f1 * ((r1 - r0) / (f1 - f0))

        if math.isclose(f(r), 0, abs_tol=tol) or abs(r1 - r0) < tol:
            return r

        r0, r1 = r1, r
        f0, f1 = f1, f(r)

    raise IteracoesExcedidas


def newton_raphson(f, df, r0, tol, nMax):
    if f(r0) == 0.0:
        return r0

    for i in range(0, nMax):
        if df(r0) == 0.0:
            raise ZeroDivisionError
        r1 = r0 - (f(r0) / df(r0))
        if math.isclose(f(r1), 0, abs_tol=tol) or abs(r1 - r0) < tol:
            return r1

        r0 = r1

    raise IteracoesExcedidas


def biseccao(f, a, b, tol, nMax):
    if f(a) * f(b) > 0:
        raise SemRaizNoIntervalo
    if (b - a) / 2 < tol:
        return a + (b - a) / 2

    for i in range(0, nMax):
        r = a + (b - a) / 2
        fr = f(r)

        if math.isclose(fr, 0, abs_tol=tol) or ((b - a) / 2) < tol:
            return r

        if f(a) * fr < 0:
            b = r
        else:
            a = r
    raise IteracoesExcedidas


def main():
    title = 'Selecione um método para encontrar a raiz da função: '
    options = ('Biseccao', 'Newton-Raphson', 'Secante')
    option = pick(options, title)

    try:
        x = smp.symbols('x', real=True)  # Define x como variável real
        f = input("Insira a funcao: ")
        transformations = (standard_transformations + (implicit_multiplication,) + (
            convert_xor,))  # Permite que o parse aceite '^' = ** e 'ax' = a * x
        f = smp.parse_expr(f, local_dict={'x': x},
                           transformations=transformations)  # Converte o input em uma expressão sympy
        df = smp.diff(f, x)  # Deriva f em função de x

        print("Insira o intervalo [a; b] onde existe única uma raiz de f(x): ")
        a = float(input("a: "))
        b = float(input("b: "))

        verificar_intervalo(f, a, b)  # Caso não exista uma única raiz, retorna erro

        raiz = secante(smp.lambdify(x, f), 0, 1, 0.0005, 50)

        print(f"Raiz encontrada: {raiz}")
        print(f"f({raiz}) = {f.subs(x, raiz)}")

    except SemRaizNoIntervalo as e:
        print("Não é possível garantir a existência de uma raiz da função no intervalo dado.")
    except IteracoesExcedidas as e:
        print("O número de iterações excedeu o limite dado.")
    except ZeroDivisionError as e:
        print("Ocorreu um erro de divisao por zero.")


main()
