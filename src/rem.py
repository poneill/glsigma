import random
from utils import *
from math import *
from matplotlib import pyplot as plt
from scipy.stats import norm

def test_lemma(n):
    return max(random.gauss(0,1) for i in xrange(2**n))/sqrt(n)

def test_inequality():
    def f(u):
        return 1/u*exp(-u**2/2.0) * (1-2/u**2)
    def g(u):
        return sqrt(2*pi)*(1-norm.cdf(u))
    def h(u):
        return 1/u*exp(-u**2/2.0)
    plt.plot(*pl(f,myrange(0.1,10,.1)))
    plt.plot(*pl(g,myrange(0.1,10,.1)))
    plt.plot(*pl(h,myrange(0.1,10,.1)))
    plt.show()
