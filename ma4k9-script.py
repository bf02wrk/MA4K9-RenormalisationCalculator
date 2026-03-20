"""A simple Python script that calculates the renormalised amplitudes of a family of simple $\phi^4$ Feynman diagrams.

This is a simple Python program I wrote to handle some of the calculations for my MA4K9 dissertation on the Connes-Kreimer theory of renormalisation. Its purpose is to calculate the renormalised amplitudes of massless $\phi^4$ Feynman graphs in the family >O<, >OO<, >OOO<, etc. 

The function calculate_renormalised_amplitude(N) does this by taking the loop number $N$ of the graph as an input, calculating the appropriate number of terms of the series expansions of $\Gamma(z)$ and $(1-2z)^{-1}$, multiplying them out, discarding unneeded terms and then applying the main theorem of Connes-Kreimer theorem to find an expression for the renormalised part of the Feynman integral in terms of $z$, $\lambda$, $L=\log\left( \frac{p^2}{4\pi \mu^2} \right)$ and the constants $\gamma$ and $\zeta(n)$, $n\geq 2$.

The function evaluate(expression, log, coupling_constant, loop_num) substitutes into the input expression the values $z=0$, $L=$log_val, $\lambda=$coupling_constant and uses $N=$loop_num to replace the zeta values with their numerical values. This gives us a single complex value for each loop number.

Dependencies: sympy, itertools, scipy, time, numpy
"""
import sympy
from itertools import chain, combinations
import scipy
import time
from numpy import euler_gamma, log

lam, gamma, z, L = sympy.symbols('\lambda \gamma z L')

def coeffs(N):   
    """Calculates the coefficients of gamma(z) and 1/gamma(z) following the recursive formula of [Zwillinger & Jeffrey 2000]"""
    S = [0, sympy.symbols('\gamma')] + [sympy.symbols(f'\zeta({i})') for i in range(2, N + 1)]
    C = [1]
    D = [1]
    for n in range(0, N):
        C += [0]
        D += [0]
        for k in range(0, n + 1):
            C[-1] += sympy.Rational(1, n + 1) * (-1)**(k + 1) * S[k + 1] * C[n - k]  # formulae for laurent coefficients
            D[-1] += sympy.Rational(1, n + 1) * (-1)**(k) * S[k + 1] * D[n - k]
    return C, D

def Gamma(x, C):
    """Laurent polynomial for gamma(x) where C is the list of coefficients"""
    total = 0
    for n in range(len(C)):
        total += C[n] * x**(n - 1)
    return total

def invGamma(x, D):
    """Laurent polynomial for 1/gamma(x) where D is the list of coefficients"""
    total = 0
    for n in range(len(D)):
        total += D[n] * x**(n + 1)
    return total

def T(p):
    """Renormalisation scheme - returns the pole part of p"""
    q = sympy.polys.polytools.reduced(p, [z])[1]
    q = q - sympy.polys.polytools.reduced(q, [1/z])[1]
    return q

def evaluate(expression, log_val, coupling_constant, loop_num): 
    """Substitutes values into the expression

    Keyword arguments:
    expression -- the SymPy expression into which we substitute the values
    log_val -- value of L
    coupling_constant -- value of lambda
    loop_num -- value of N, number of loops
    """
    new = expression.subs({lam:coupling_constant, gamma:euler_gamma,  sympy.pi:scipy.pi, L:log_val})
    S = [sympy.symbols(f'\zeta({i})') for i in range(2, loop_num + 1)]
    new = new.subs([(S[i - 2], scipy.special.zeta(i)) for i in range(2, loop_num + 1)])
    return sympy.simplify(new)

def powerset(iterable): 
    """Subsequences of the iterable from shortest to longest. Copied from itertools documentation.
    E.g. powerset([1,2,3]) → () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)
    """
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s) + 1)) 


def calculate_renormalised_amplitude(N):  
    """Does the whole calculation for a graph with N loops"""
    C, D = coeffs(N)
    gamma1 =  sympy.simplify(sympy.expand(-sympy.Rational(1, 2) *z *sum([(2 * z)**i for i in range(N + 1)]) * Gamma(z, C) * Gamma(-z, C)**2 * invGamma(-2 * z, D)))
    gammas = [gamma1]
    for i in range(N - 1):
        nextgamma = sympy.simplify(sympy.polys.polytools.reduced(gammas[-1] * gamma1, [z**(N - 1 - i)])[1])
        gammas += [nextgamma]   # iteratively calculates the gamma part of the integral and adds it to a list for efficiency
    
    
    def phi(n):     # the feynman rule
        if n == 0:
            return 1
        else:
            answer = -sympy.I * lam**(n + 1) / (2**n * (4 * sympy.pi)**(2 * n))     # constant part
            answer *= sum([(-n * z * L)**i / sympy.factorial(i) for i in range(N + 1)])     # exponent part
            answer *= gammas[n - 1]   # gamma part
            answer = sympy.polys.polytools.reduced(answer, [z**(N + 1 - n)])[1]     # truncates
            return  sympy.simplify(answer)

    phi_values = [phi(i) for i in range(1, N + 1)]     # unrenormalised feynman integrals for each order up to n - making a list for efficiency
    phi_minus_values = [-T(phi_values[0])]
    phi_plus_values = [phi_values[0] + phi_minus_values[0]]
    X = [sympy.symbols(f'X_{i}') for i in range(N)]     # formal symbols for phi_-(\Gamma_i) so we can simplify before substituting
    Y = [sympy.symbols(f'Y_{i}') for i in range(N)]     # formal symbols for phi(\Gamma_i)
    for n in range(2, N + 1):
        loops = [i for i in range(n)]   # list of loops in the graph labelled 0 to n-1
        div = list(powerset(loops))     # set of divergent subgraphs
        div.remove(())  # remove 1 and \Gamma
        div.remove(tuple(loops))
        total_phi = 0   # our final answer
        for subset in div:
            complement = [i for i in loops if i not in subset]
            phi_term = 1    # contribution from this subgraph
            counter = 0     # counting the number of loops in each connected component of the subgraph
            for j in subset:    # counts the size of connected components
                if j + 1 in subset:
                    counter += 1
                else:
                    phi_term *= X[counter]
                    counter = 0
            counter = 0
            for j in complement:
                if j + 1 in complement:
                    counter += 1
                else:
                    phi_term *= Y[counter]
                    counter = 0
            total_phi += phi_term
        total_phi = total_phi.subs([(X[i], phi_minus_values[i]) for i in range(n - 1)])    # substituting our phi values into the polynomial
        total_phi = total_phi.subs([(Y[i], phi_values[i]) for i in range(n - 1)])
        phi_minus_values += [-T(sympy.simplify(phi_values[n - 1] + total_phi))]     # adding phi- to the list
        phi_plus_values += [sympy.simplify(phi_values[n - 1] + phi_minus_values[n - 1] + total_phi)]
    return sympy.polys.polytools.reduced(phi_plus_values[N - 1], [z])[1]

if __name__ == "__main__":
    for i in range(1,4):
        t0 = time.time()
        expression = calculate_renormalised_amplitude(i)
        log_val = -log(4*scipy.pi) + (0+1j)*scipy.pi
        print(sympy.latex(evaluate(expression, log_val, 0.1, i)))
        print(time.time() - t0)
        
