import math
from logging import exception

import sympy as smp
from pick import pick
from sympy.parsing.sympy_parser import convert_xor, standard_transformations, implicit_multiplication_application


class IteracoesExcedidas(Exception):
    pass


class SemRaizNoIntervalo(Exception):
    pass


class MaisDeUmRaiz(Exception):
    pass


def verificar_intervalo(f, df, df2, a, b):
    if f(a) * f(b) > 0:
        raise SemRaizNoIntervalo

    n = 200  # Quantidade de pontos para testar

    crescente = False
    decrescente = False

    # Calcula os valores em df2() para cada ponto
    for i in range(n):
        valor = a + i * (b - a) / (n - 1)
        p = df2(valor)
        if p == 0.0 or df(valor) == 0.0:
            raise MaisDeUmRaiz
        if p > 0 and not decrescente:
            crescente = True
        elif p < 0 and not crescente:
            decrescente = True
        else:
            raise MaisDeUmRaiz


def secante(f, r0, r1, tol, nMax):
    f0 = f(r0)
    f1 = f(r1)
    if math.isclose(f0, 0, rel_tol=tol):
        return r0

    for k in range(0, nMax):
        if f1 == f0:
            raise ZeroDivisionError

        r = r1 - f1 * ((r1 - r0) / (f1 - f0))

        if abs(r1 - r0) < tol or math.isclose(f(r), 0, rel_tol=tol):
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
        if abs(r1 - r0) < tol or math.isclose(f(r1), 0, rel_tol=tol):
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

        if ((b - a) / 2) < tol or math.isclose(fr, 0, rel_tol=tol):
            return r

        if f(a) * fr < 0:
            b = r
        else:
            a = r
    raise IteracoesExcedidas


def main():
    is_running = True

    try:
        while is_running:
            title = 'Selecione um método para encontrar a raiz da função: '
            options = ('Biseccao', 'Newton-Raphson', 'Secante')
            option = pick(options, title)

            x = smp.symbols('x', real=True)  # Define x como variável real
            f = input("Insira a funcao: ")
            transformations = (standard_transformations +
                               (implicit_multiplication_application,
                                convert_xor))  # Permite que o parse aceite '^' = **, 'ax' = a * x, e^a = exp(a)
            locals_vars = {'x': x, 'e': smp.E}  # Permite que o perse reconhece e = exp() e x como variável
            f = smp.parse_expr(f, local_dict=locals_vars,
                               transformations=transformations)  # Converte o input em uma expressão sympy

            print("Insira o intervalo [a; b] onde existe única uma raiz de f(x).")
            a = float(input("a: "))
            b = float(input("b: "))

            df = smp.diff(f, x)  # Deriva f em função de x
            df2 = smp.diff(df, x)  # Segunda derivada de f

            f = smp.lambdify(x, f)  # Converte a função em método
            df = smp.lambdify(x, df)  # Converte a derivada em método
            df2 = smp.lambdify(x, df2)  # Converte a segunda derivada em método

            if not option[0] == 'Biseccao':
                verificar_intervalo(f, df, df2, a, b)  # Caso não exista uma única raiz, retorna erro

            n_max = int(input("Insira o número máximo de iterações: "))
            tol = float(input("Insira a tolerância (absoluta ou relativa): "))

            if option[0] == 'Biseccao':
                raiz = biseccao(f, a, b, tol, n_max)
            elif option[0] == 'Newton-Raphson':
                r0 = a if (df2(a) * f(
                    a) > 0) else b  # Verifica quais dos extremos tem o mesmo sinal da segunda derivada
                raiz = newton_raphson(f, df, r0, tol, n_max)
            else:
                raiz = secante(f, a, b, tol, n_max)

            print(f"Raiz encontrada: {raiz}")
            print(f"f({raiz}) = {f(raiz)}")

            print("\n")
            is_running = input("Deseja continuar? (s/n)").lower() == "s"


    except SemRaizNoIntervalo as e:
        print("Não é possível garantir a existência de uma raiz da função no intervalo dado.")
    except IteracoesExcedidas as e:
        print("O número de iterações excedeu o limite dado.")
    except ZeroDivisionError as e:
        print("Ocorreu um erro de divisão por zero.")
    except MaisDeUmRaiz as e:
        print("Não é possível garantir a unicidade no intervalo")
    except exception as e:
        print(f"Erro inesperado: {e}")


main()
