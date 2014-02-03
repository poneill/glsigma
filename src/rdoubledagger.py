"""
The purpose of this script is to explore the new statistic
R_doubledagger (= (rdagger_high + rdagger_lo)/2) in the singleton
regulatory scenario.
"""
from utils import *
G = 10000
n = 16
alpha = 0.5
alpha_high = 0.25
w = 10
epsilon = 2
beta = 1.61
const = 1
Zb = G*((exp(-beta*(epsilon + const))+exp(-beta*const))/2.0)**w
Z = Zb/(1-alpha)
e1 = (log(alpha_high) + log(Z))/-beta
ei = (log(alpha-alpha_high) - log(n-1) + log(Z))/-beta
mismatches_1 = e1/(epsilon)
mismatches_i = ei/(epsilon)

p_1 = exp(-beta*e1)/Z
p_i = exp(-beta*ei)/Z

def sample_motif():
    return ([[int(random.random() > mismatches_1/w) for i in range(w)]] +
            [[int(random.random() > mismatches_i/w) for i in range(w)]
             for j in range(n-1)])

print "loaded"
