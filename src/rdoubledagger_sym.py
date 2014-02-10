from sympy import *

p, q, n, alpha_i,alpha_o,alpha, ep_i,ep_o,Z, Zb = var("p, q, n, alpha_i,alpha_o, alpha, ep_i,ep_o,Z, Zb",real=True)

G,w,epsilon,beta,A,logomega = var("G,w,epsilon,beta,A,logomega",positive=True)
omega = (A-1 + exp(-beta*-epsilon))/A
#Zb = G*omega**w
Zb = G*exp(logomega)**w
Z = Zb/(1-alpha)
alpha = alpha_i + alpha_o
ep_i = (log(alpha_i) + log(Z))/-beta
ep_o = (log(alpha_o) + log(Z))/-beta
m_i = ep_i/(-epsilon)
p_i = m_i/w
m_o = ep_o/(-epsilon)
p_o = m_o/w
p = (p_i + (n-1)*p_o)/n
q = (1-p)/(A-1)
H = -(p*log(p,2) + q*log(q/3,2))
