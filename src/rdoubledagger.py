"""
The purpose of this script is to explore the new statistic
R_doubledagger (= (rdagger_high + rdagger_lo)/2) in the singleton
regulatory scenario: one site at alpha_high, n-1 sites at (alpha-alpha_high)/(n-1).
"""
from math import *
from utils import *
beta = 1.61
alphabet = 4
epsilon = 2
w = 16

def find_Zb(G,w,epsilon):
    return G*(((alphabet - 1)*1+exp(-beta*(-epsilon)))/alphabet)**w
def rseq(G=10000,n=16,alpha=0.5,alpha_high=0.25,w=10,epsilon=2):
    Zb = find_Zb(G,w,epsilon)
    Z = Zb/(1-alpha)
    e1 = (log(alpha_high) + log(Z))/-beta
    ei = (log(alpha-alpha_high) - log(n-1) + log(Z))/-beta
    print "e1:",e1,"ei:",ei
    # check to see that probabilities are satisfied
    p_1 = exp(-beta*e1)/Z
    p_i = exp(-beta*ei)/Z
    print "alpha_high:",alpha_high,p_1
    print "alpha-alpha_high:",alpha-alpha_high,(n-1)*p_i

    matches_1 = -e1/(epsilon)
    matches_i = -ei/(epsilon)
    print "matches:",matches_1,matches_i
    match_prob1 = matches_1/w
    match_probi = matches_i/w
    match_prob = (match_prob1 + (n-1)*match_probi)/float(n)
    mismatch_prob = 1 - match_prob
    column_entropy = -(match_prob*log2(match_prob) +
                      ((alphabet-1)*(mismatch_prob/(alphabet-1))*
                       log2(mismatch_prob/(alphabet-1))))
    prior_entropy = log2(alphabet)
    return w*(prior_entropy-column_entropy)

def rdd(G=10000,n=16,alpha=0.5,alpha_high=0.25,w=10,epsilon=2):
    """compute Rdd by levels"""
    alpha_low = alpha - alpha_high
    rd_low = log2(G) - log2(n-1) + log2(alpha_low)
    rd_high = log2(G) - log2(1) + log2(alpha_high)
    return (rd_low + rd_high)/2

def rddd(G=10000,n=16,alpha=0.5,alpha_high=0.25,w=10,epsilon=2):
    """compute Rdd by sites"""
    alpha_low = alpha - alpha_high
    rd_low = log2(G) - log2(n-1) + log2(alpha_low)
    rd_high = log2(G) - log2(1) + log2(alpha_high)
    print "rd_low,rd_high:",rd_low,rd_high
    return ((n-1)*rd_low + rd_high)/n + 1.6

def new_stat(alphas,G=10000,w=16,epsilon=2):
    print "in new stat:"
    print "alphas:",alphas
    print "G:",G
    print "w:",w
    n = len(alphas)
    print "n:",n
    # Step 1: find ep_i for each alpha_i
    Zb = G*(((alphabet - 1)*1 + exp(-beta*-epsilon))/alphabet)**w
    print "Zb:",Zb
    alpha = sum(alphas)
    Z = Zb*(1/(1-alpha))
    eps = [(log(a) + log(Z))/-beta for a in alphas]
    print "eps:",eps
    recovered_alphas = [exp(-beta*ep)/Z for ep in eps]
    print recovered_alphas
    sse = sum([(a-re)**2 for a,re in zip(alphas,recovered_alphas)])
    print "SSE:",sse
    assert sse < 1e-16
    # Step 2: compute Rseq from ep_is
    ps = [-ep/(epsilon*w) for ep in eps] # match probabilities
    print "ps:",ps
    p = mean(ps) # average mismatch prob
    pcheck = 1/(epsilon*beta)*(1.0/(n*w)*sum(log(a) for a in alphas) +
                                (1.0/w)*(log(G)-log(1-alpha)) +
                                log((3+exp(beta*epsilon))/4))
    
    
    q = 1 - p # match prob
    
    print "match prob p:",p,"pcheck:",pcheck
    col_entropy = -(p*log2(p) + 3*(q/3.0)*log2(q/3.0))
    print "col_entropy:",col_entropy
    Hprior = log2(alphabet)
    rseq = w*(Hprior - col_entropy)
    return rseq,mean(eps),p

def test_col_ent(m,w):
    q = m/float(w)
    p = 1 - q
    return -(p*log2(p) + 3*(q/3.0)*log2(q/3.0))

def alpha_sabot(alpha,alpha_high,n):
    site_alpha_low = (alpha - alpha_high)/(n-1)
    return [alpha_high] + [site_alpha_low]*(n-1)

def rd(G=10000,n=16,alpha=0.5,alpha_high=0.25,w=10,epsilon=2):
    return log2(G) - log2(n) + log2(alpha)

def sample_Zb(w=16):
    """Check to confirm Zb: it works"""
    return sum(exp(-beta*(epsilon*sum(random.random() > 1/float(alphabet)
                                 for j in range(w))))
               for i in range(G))

def sample_motif(qs,w):
    return ([[int(random.random() < q) for i in range(w)]
             for q in qs])

print "loaded"

def make_plot():
    alpha_start = 0.1
    alpha_stop= 0.7
    alpha_step = 0.1
    alpha_high_start = 0.01
    alpha_high_step = 0.001
    n = 16
    for a in myrange(alpha_start,alpha_stop,alpha_step):
 	plt.plot(*pl(lambda ah:rseq(alpha=a,alpha_high=ah),
                     myrange(alpha_high_start,a,alpha_high_step)))
        # plt.plot(*pl(lambda ah:rddd(alpha=a,alpha_high=ah),
        #              myrange(alpha_high_start,a,alpha_high_step)),linestyle='--')
        plt.plot(*pl(lambda ah:new_stat(alpha_sabot(alpha=a,alpha_high=ah,n=n,w=w)),
                     myrange(alpha_high_start,a,alpha_high_step)),linestyle='--')
        plt.scatter(a/n,rseq(alpha=a,alpha_high=a/n))
    plt.show()

def make_all_alphas(trials,vary_n=False):
    all_alphas = []
    m = n
    for i in xrange(trials):
        if vary_n:
            m = random.choice(range(4,16))
        pre_rs = [0] + sorted([random.random() for _ in range(m-1)]) + [1]
        rs = [y-x for (x,y) in pairs(pre_rs)]
        alpha_factor = random.random()
        all_alphas.append([r*alpha_factor for r in rs])
    return all_alphas

def random_alphas_plot():
    trials = 1000
    #statistic = lambda alphas:tan(pi*sum(alphas)-pi/2)
    #statistic = lambda alphas:(lambda x:log2(x) + log2(1/(1-x) + log2(G)))(sum(alphas))
    statistic = sum
    statistic2 = lambda x:log(x) + log(1/(1-x)) + log(G)
    all_alphas = filter(lambda alphas:0.01 < sum(alphas) < 0.99,
                        make_all_alphas(trials,vary_n=False))
    plt.scatter(*transpose([(show(statistic(alphas)),new_stat(alphas))
                            for alphas in all_alphas]))
    plt.plot(*pl(statistic2,myrange(0.01,1,0.01)),color='r')
    plt.show()
