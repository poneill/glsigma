"""
Reduce the glsigma model to simplest possible setting: binary alphabet, fixed epsilon.
"""
import random
from utils import *

G = 10000
w = 10
n = 16
epsilon = 2
beta = 1.61

def score(site):
    return sum(site)*epsilon

def make_genome(G):
    return make_site(G)

def score_genome(genome):
    return [score(genome[i:i+w]) for i in range(G-w + 1)]

def make_site(w):
    return [int(random.random() < 0.5) for _ in range(w)]
    

def compute_props(genome):
    return [exp(-beta*ep) for ep in score_genome(genome)]

def compute_Z(genome):
    return sum(compute_props(genome))

def predict_Z():
    return G*((exp(-beta*epsilon) + 1)/2.0)**w

def compute_epsilon_i(alpha):
    """Estimate average number of mismatches if P(s_i) = alpha/n"""
    Zb = predict_Z()
    return (1/beta)*log((n*(1/float(alpha)-1))/Zb)

def compute_rseq(alpha,n=n,w=w):
    """Assuming mismatches are distributed independently, compute
    expected Rseq for motif satisfying alpha"""
    epsilon_i = compute_epsilon_i(alpha)
    p = epsilon_i/float(w) # prob of mismatch in column, assuming independence
    safelog2 = lambda x: log2(x) if x > 0 else 0
    col_entropy = -(p*safelog2(p) + (1-p)*safelog2(1-p))
    rseq = (1-col_entropy)*w
    return rseq

def fw_approx(mu,sigma_sq,n):
    """Suppose Z = sum(X), where X_i~LN(mu,sigma_sq),1<=i<=n.  Assume
    Z is also log normal and return mu_z,sigma_sq_z.  See Cobb and
    Rumi, "Approximating the Distribution of a Sum of Log-normal
    Random Variables".
    """
    sigma_sq_z = log((exp(sigma_sq)-1)/float(n) + 1)
    mu_z = log(n*exp(mu)) + (sigma_sq-sigma_sq_z)/2.0
    return mu_z,sigma_sq_z

def fw_test(mu,sigma_sq,n,trials=100):
    """Verifies that fw_approx basically works"""
    Zs = [sum([exp(random.gauss(mu,sqrt(sigma_sq))) for i in range(n)])
          for t in range(trials)]
    mu_z,sigma_sq_z = fw_approx(mu,sigma_sq,n)
    pred_mean_Z = ln_mean(mu_z,sigma_sq_z)
    pred_var_Z = ln_var(mu_z,sigma_sq_z)
    print "pred mu:",pred_mean_Z,"obs mu:",mean(Zs)
    print "pred sigma_sq:",pred_var_Z,"obs sigma_sq:",variance(Zs)
    
def params_of_Z(G=G):
    """Z = sum_{i=1}^G exp(-beta*epsilon).  Return mu_z,sigma_sq_z"""
    mu = -beta*epsilon*w/2.0
    sigma_sq = beta**2 * epsilon**2*w/4.0
    print "mu:",mu,"sigma_sq:",sigma_sq
    mu_z,sigma_sq_z = fw_approx(mu,sigma_sq,G)
    return mu_z,sigma_sq_z

def approx_Z_test(trials=100,G=G):
    Zs = [compute_Z(make_genome(G)) for i in verbose_gen(range(trials))]
    logZs = map(log,Zs)
    print "M-W:",mannwhitneyu(logZs,normal_model(logZs))
    mu_z,sigma_sq_z = params_of_Z(G)
    pred_mean_Z = ln_mean(mu_z,sigma_sq_z)
    pred_var_Z = ln_var(mu_z,sigma_sq_z)
    pred_log_Z = log(G * ((exp(-beta*epsilon)+1)/2.0)**w)
    print "pred mu:",pred_mean_Z,"obs mu:",mean(Zs)
    print "pred sigma_sq:",pred_var_Z,"obs sigma_sq:",variance(Zs)
    print "pred <logZ>:",pred_log_Z,"obs <log Z>",mean(logZs)
    
def ln_mean(mu,sigma_sq):
    """Compute mean of log normal rv"""
    return exp(mu + sigma_sq/2.0)

def ln_var(mu,sigma_sq):
    """Compute variance of log normal rv"""
    return (exp(sigma_sq)-1)*exp(2*mu+sigma_sq)

def compute_rfreq(alpha,G=G,n=n):
    fg_term = alpha*(log(alpha)-log(n))
    h_in_nats = beta*epsilon + log(G) + w*log((exp(-beta*epsilon)+1)/2.0) + fg_term
    h_in_bits = h_in_nats * log2(exp(1))
    return log2(G+n) - h_in_bits

def sample_rfreq(alpha,G=G):
    ei = compute_epsilon_i(alpha)
    genome = make_genome(G)
    site_props = [exp(-beta*ei) for i in range(n)]
    props = compute_props(genome)
    Z = sum(site_props + props)
    ps = [pr/Z for pr in site_props + props]
    prior_ent = log2(G+n)
    post_ent = -sum(p*log2(p) for p in ps)
    rfreq = prior_ent - post_ent
    print alpha,sum(p/Z for p in site_props),ei,prior_ent,post_ent,rfreq
    return rfreq
    
    
