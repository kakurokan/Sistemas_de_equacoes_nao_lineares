import math

import sympy as smp
from sympy.parsing.sympy_parser import implicit_multiplication, standard_transformations, convert_xor


class IteracoesExcedidas(Exception):
    pass


class SemRaizNoIntervalo(Exception):
    pass


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
    try:
        x = smp.symbols('x', real=True)
        f = input("Insira a funcao: ")
        transformations = (standard_transformations + (implicit_multiplication,) + (convert_xor,))
        f = smp.parse_expr(f, local_dict={'x': x}, transformations=transformations)
        df = smp.diff(f, x)

        raiz = secante(smp.lambdify(x, f), 0, 1, 0.0005, 50)

        print(f"Raiz encontrada: {raiz}")
        print(f"f({raiz:}) = {f.subs(x, raiz)}")
    except SemRaizNoIntervalo as e:
        print("Não existe raiz no intervalo.")
    except IteracoesExcedidas as e:
        print("Numero de iteracoes excedido.")
    except ZeroDivisionError as e:
        print("Divisao por zero.")


main()
