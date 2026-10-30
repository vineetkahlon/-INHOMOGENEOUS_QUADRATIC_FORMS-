"""Numerical discovery helpers, never imported by the verifier."""
import numpy as np
from scipy.optimize import linprog
from fractions import Fraction as F
from exact_linear import solve_exact, verify_leaf

def farkas(rows,dim):
    mat=np.array([[*a,b] for a,b in rows],dtype=float).T
    r=linprog(np.ones(len(rows)),A_eq=mat,b_eq=[0]*dim+[-1],bounds=(0,None),method='highs')
    if not r.success:raise ValueError('No numerical Farkas multiplier')
    support=np.flatnonzero(r.x>1e-9).tolist()
    weights=solve_exact([[*rows[i][0],rows[i][1]] for i in support],[F(0)]*dim+[F(-1)])
    cert=[[i,str(w)] for i,w in zip(support,weights) if w]
    verify_leaf(rows,cert,dim);return cert
