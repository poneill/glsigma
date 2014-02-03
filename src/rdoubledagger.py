"""
The purpose of this script is to explore the new statistic
R_doubledagger (= (rdagger_high + rdagger_lo)/2) in the singleton
regulatory scenario: one site at alpha_high, n-1 sites at (alpha-alpha_high)/(n-1).
"""

from utils import *
beta = 1.61
alphabet = 4
def rseq(G=10000,n=16,alpha=0.5,alpha_high=0.25,w=10,ep=2):
    Zb = G*(((alphabet - 1)*exp(-beta*(ep))+1)/alphabet)**w
    Z = Zb/(1-alpha)
    e1 = (log(alpha_high) + log(Z))/-beta
    ei = (log(alpha-alpha_high) - log(n-1) + log(Z))/-beta
    # check to see that probabilities are satisfied
    p_1 = exp(-beta*e1)/Z
    p_i = exp(-beta*ei)/Z
    print "alpha_high:",alpha_high,p_1
    print "alpha-alpha_high:",alpha-alpha_high,(n-1)*p_i

    mismatches_1 = e1/(ep)
    mismatches_i = ei/(ep)
    print "mismatches:",mismatches_1,mismatches_i
    mismatch_prob1 = mismatches_1/w
    mismatch_probi = mismatches_i/w
    mismatch_prob = (mismatch_prob1 + (n-1)*mismatch_probi)/float(n)
    match_prob = 1 - mismatch_prob
    column_entropy = -(match_prob*log2(match_prob) +
                      ((alphabet-1)*(mismatch_prob/(alphabet-1))*
                       log2(mismatch_prob/(alphabet-1))))
    prior_entropy = log2(alphabet)
    return w*(prior_entropy-column_entropy)

def rdd(G=10000,n=16,alpha=0.5,alpha_high=0.25,w=10,ep=2):
    alpha_low = alpha - alpha_high
    rd_low = log2(g) - log2(15) + log2(alpha_low)
    rd_high = log2(g) - log2(1) + log2(alpha_high)
    return (rd_low + rd_high)/2


def sample_Z():
    """Check to confirm Zb: it works"""
    return sum(exp(-beta*(ep*sum(random.random() > 1/float(alphabet)
                                 for j in range(w))))
               for i in range(G))
def sample_motif():
    return ([[int(random.random() < mismatches_1/w) for i in range(w)]] +
            [[int(random.random() < mismatches_i/w) for i in range(w)]
             for j in range(n-1)])

print "loaded"

def make_plot():
    alpha_start = 0.1
    alpha_stop= 0.7
    alpha_step = 0.1
    alpha_high_start = 0.01
    alpha_high_step = 0.001
    for a in myrange(alpha_start,alpha_stop,alpha_step):
 	plt.plot(*pl(lambda ah:rseq(alpha=a,alpha_high=ah),
                     myrange(alpha_high_start,a,alpha_high_step)))

