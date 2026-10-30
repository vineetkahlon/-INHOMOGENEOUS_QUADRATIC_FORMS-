"""Exact linear algebra; Python standard library only."""
from fractions import Fraction as F

def solve_exact(columns,target):
    """Solve a small rational linear system, setting free variables to zero."""
    n=len(columns)
    M=[[columns[j][i] for j in range(n)]+[target[i]] for i in range(len(target))]
    pivot=[];r=0
    for col in range(n):
        p=next((i for i in range(r,len(M)) if M[i][col]),None)
        if p is None:continue
        M[r],M[p]=M[p],M[r]
        scale=M[r][col];M[r]=[v/scale for v in M[r]]
        for i in range(len(M)):
            if i!=r and M[i][col]:
                q=M[i][col];M[i]=[a-q*b for a,b in zip(M[i],M[r])]
        pivot.append(col);r+=1
        if r==len(M):break
    if any(not any(row[:-1]) and row[-1] for row in M):
        raise ValueError('Inconsistent exact recovery')
    out=[F(0)]*n
    for i,p in enumerate(pivot):out[p]=M[i][-1]
    return out

def verify_leaf(rows,certificate,dim):
    coeff=[F(0)]*dim;constant=F(0)
    for index,raw in certificate:
        if type(index) is not int or not 0<=index<len(rows):raise ValueError('Bad row index')
        w=F(raw)
        if w<0:raise ValueError('Negative multiplier')
        a,b=rows[index]
        coeff=[c+w*x for c,x in zip(coeff,a)];constant+=w*b
    if any(coeff) or constant>=0:raise ValueError('Invalid exact Farkas identity')
