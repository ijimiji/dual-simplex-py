from fractions import Fraction
from copy import deepcopy
import math

# Use if limitation is absent
infinity = Fraction(10**8)

def transpose(A):
    n = len(A)
    B = [[Fraction() for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            B[i][j] = A[j][i]
    return B

def sign(x):
    if x < 0:
        return -1
    if x > 0:
        return 1
    return 0

def col(A, i):
    return [A[j][i] for j in range(len(A))]

def dot(x, y):
    return Fraction(sum([a*b for (a, b) in zip(x, y)]))

def solve(A, b):
    a, f, n, x = deepcopy(A), deepcopy(b), len(A), [Fraction()] * len(A)

    for k in range(0, n-1):
        for i in range(k+1, n):
            factor = a[i][k] / a[k][k]
            for j in range(k, n):
                a[i][j] = a[i][j] - factor * a[k][j]
            f[i] = f[i] - factor * f[k]
    x[n-1] = f[n-1] / a[n-1][n-1]
    for k in range(n-2, -1, -1):
        sums = f[k]
        for j in range(k+1, n):
            sums = sums - a[k][j] * x[j]
        x[k] = sums / a[k][k]
    return x

def simplex(A, b, c, d_lower, d_upper, basis):
    x = []
    while x == []:
        n = len(c)
        nonbasis = [i for i in range(n) if i not in basis]
        A_base = [[A[i][j] for j in basis] for i in range(len(A))]
        c_base = [c[i] for i in basis]

        u = solve(transpose(A_base), c_base)

        δ = [Fraction() for _ in range(n)]
        for i in nonbasis:
            δ[i] = c[i] - dot(col(A, i), u)

        kappa = [Fraction() for _ in range(n)]
        for i in nonbasis:
            kappa[i] = d_lower[i] if δ[i] <= 0 else d_upper[i]
        
        kappa_base = solve(A_base, [b[i] - sum([kappa[j]*A[i][j] for j in nonbasis]) for i in range(len(A))])
        for (i,value) in zip(basis, kappa_base):
            kappa[i] = value
        
        j_star = -1
        for i in range(n):
            if not(kappa[i] <= d_upper[i] and kappa[i] >= d_lower[i]):
                j_star = i
                break

        if j_star == -1:
            x = kappa
            break

        l = [Fraction() for _ in range(n)]
        l_u = solve(
            [col(A, j) for j in basis], 
            [Fraction() if j != j_star else -sign(kappa[j_star] - (d_lower[j_star] if abs(kappa[j_star] - d_lower[j_star]) < abs(kappa[j_star] - d_upper[j_star]) else d_upper[j_star])) for j in basis]
        )

        for i in nonbasis:
            l[i] = - dot(col(A, i), l_u)

        for (i, val) in zip(basis, l_u):
            l[i] = val

        σ = [-δ[j] / l[j] if l[j]*δ[j] < 0 else infinity for j in nonbasis]

        j_zero = -1
        j_zero = nonbasis[0]
        min_σ = σ[0]
        for (i, val) in zip(nonbasis, σ):
            if min_σ > val:
                min_σ = val
                j_zero = i

        basis = sorted([i for i in basis if i != j_star] + [j_zero])
        y = sum([kappa[i] * c[i] for i in range(len(c))])

    y = sum([x[i] * c[i] for i in range(len(c))])
    return y, x, basis

"""
Example usage

A = [
    [Fraction(1), Fraction(1), Fraction(-1), Fraction()],
    [Fraction(1), Fraction(), Fraction(), Fraction(1)],
]
c = [Fraction(-3), Fraction(1), Fraction(), Fraction()]
b = [Fraction(2), Fraction(4)]
d_lower = [Fraction(-5, 2), Fraction(), Fraction(), Fraction()]
d_upper = [Fraction(infinity), Fraction(16, 3), infinity, infinity]
y, x, basis = simplex(A, b, c, d_lower, d_upper, [0, 1])
xs = ["0" if frac.numerator == 0 else f"{frac.numerator}/{frac.denominator}" for frac in x]
print(f"{y = }; {xs = }; {basis = }")
"""

"""
Get solution in Integers

for i in range(len(x)):
    lower = deepcopy(d_lower)
    upper = deepcopy(d_upper)
    ceil = math.ceil(x[i])
    floor = math.floor(x[i])
    upper[i] = floor
    y_floor = -infinity
    if d_lower[i] < upper[i]:
        y_floor, x, basis = simplex(A, b, c, d_lower, upper, basis)
    lower[i] = ceil
    y_ceil = -infinity
    if lower[i] < d_upper[i]:
        y_ceil, x, basis = simplex(A, b, c, lower, d_upper, basis)
    if y_floor > y_ceil:
        d_upper = upper
    else:
        d_lower = lower
y, x, basis = simplex(A, b, c, d_lower, d_upper, [0, 1])
xs = ["0" if frac.numerator == 0 else f"{frac.numerator}/{frac.denominator}" for frac in x]
print(f"{y = }; {xs = }; {basis = }")
"""