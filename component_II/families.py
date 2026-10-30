"""Fixed-family definitions shared by builder and checker.
A family: witness value  f = n*j + n^2*alpha + n*kappa + K(n,w) - mu,  (alpha,kappa) in [0,1]^2, mu in [mu0,mu1].
'modes' are pairs (c1,c0): the box succeeds if max(f - c1*mu) < c0 (each mode is a proved valid target).
'tau_modes' are extra modes valid only on boxes inside 'nbhd' (list of (alo,ahi,klo,khi,mulo)) -- used only
for the (2,2) resonant families, whose complement inside nbhd is covered by a Farkas corner tree."""
from fractions import Fraction as F
from itertools import product

def single(t,e,c,Z=12):
    """one-variable resonant family: K = n e (z+c) - t (z+c)^2"""
    return lambda n: [((z,), n*e*(z+c)-t*(z+c)**2) for z in range(-Z,Z+1)]

def diag_block(terms,W):
    """sum over coordinates of n*b_i*(w_i+c_i) - q_i*(w_i+c_i)^2, w in W^k"""
    def K(n):
        out=[]
        for w in product(W,repeat=len(terms)):
            out.append((w, sum(n*b*(x+c)-q*(x+c)**2 for x,(q,b,c) in zip(w,terms))))
        return out
    return K

def gram_block(Mq,b,c,W):
    """large-a: K = n b.w - phi(w+c), phi(x)=x^T Mq x"""
    def K(n):
        out=[]
        for w in product(W,repeat=3):
            x=[w[i]+c[i] for i in range(3)]
            phi=sum(Mq[i][k]*x[i]*x[k] for i in range(3) for k in range(3))
            out.append((w, n*sum(b[i]*w[i] for i in range(3))-phi))
        return out
    return K
