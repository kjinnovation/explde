"""magnus_ref.py -- reference implementation for centered Magnus quadrature in Legendre moments.

Companion code for
  [1] K. C. Johnson, "The Centered Magnus Exponent and Magnus Quadrature Formulas of Orders 2, 4, 6
      and 8," Version 09-Nov-2026, magnus_order_2-8.pdf, https://ai.vixra.org/abs/2611.0005
  [2] K. C. Johnson, "Magnus quadrature formulas of order 10 on a centered interval,"
      Version 09-Nov-2026, magnus_order_10.pdf, https://ai.vixra.org/abs/2611.0004;
      coefficients in magnus10coef.m.
  Introductory tutorial: K. C. Johnson, "Quadrature approximation to the Magnus exponent for
      numerical solution of non-autonomous coupled linear differential equations,"
      Version 09-Nov-2026, explde.pdf, https://vixra.org/abs/2604.0065
Repository: https://github.com/kjinnovation/explde
Written with the assistance of Claude Opus 5.5 (Anthropic).

What it does (pure Python, exact rational arithmetic with fractions.Fraction):
  1. centered_series(D): the Magnus exponent Omega(t) through t^D in centered Legendre moments g_k,
     derived from the centered Magnus theorem ([1], Eq. 6).
  2. lie_decompose(): rewrite a Lie polynomial in the Lyndon basis (as printed in [1], Eq. 9).
  3. evaluate(): expand a nested-commutator formula in the free algebra of the moments, so that it can
     be compared word by word with the series (the exact verification used in [1] and [2]).
  4. quad_moments(): express quadrature moments in exact moments, to verify a quadrature rule.

Conventions (centered interval [-t, t]):
  g_k(t) = int_{-t}^{t} P_k(tau/t) A(tau) dtau,  weight(g_k) = k+1,  g_k = O(t^(k+1)).
  Omega(t) is the Magnus exponent of the propagator from -t to t; it is odd in t.
  Centered Magnus theorem: Omega' = sum_{k even} B_k/k! ad_Omega^k (A(t)+A(-t)) - 1/2 [Omega, A(t)-A(-t)].
  A word is a tuple of letter indices (k = moment g_k); a polynomial is a dict word -> Fraction.
  A formula is a dict with keys 'Xi.gk', 'Xi.Cj', 'Yi.gk', 'Yi.Cj', 'ci':  C_i = [X_i, Y_i],
  Omega = g0 + sum_i c_i C_i.
Run "python3 magnus_ref.py" for a self-test (series through weight 7, orders 4, 6, 8)."""
from fractions import Fraction as F
from collections import defaultdict
from math import factorial, comb

# ---------- basic algebra ----------
def bernoulli(n):                      # B_0..B_n with B_1 = -1/2
    B=[F(1)]
    for m in range(1,n+1): B.append(-sum(comb(m+1,k)*B[k] for k in range(m))/(m+1))
    return B
def wt(w): return sum(k+1 for k in w)
def add(p,q,c=1):
    r=defaultdict(F,p)
    for k,v in q.items(): r[k]+=c*v
    return {k:v for k,v in r.items() if v}
def mul(p,q,W):
    r=defaultdict(F)
    for a,x in p.items():
        for b,y in q.items():
            if wt(a)+wt(b)<=W: r[a+b]+=x*y
    return {k:v for k,v in r.items() if v}
def br(p,q,W): return add(mul(p,q,W),mul(q,p,W),-1)
def scale(p,c): return {k:c*v for k,v in p.items() if c*v}
def legendre(k):                       # coefficients of P_k (Fractions, ascending powers)
    c0,c1=[F(1)],[F(0),F(1)]
    if k==0: return c0
    for n in range(1,k):
        c2=[F(0)]*(n+2)
        for i,v in enumerate(c1): c2[i+1]+=F(2*n+1,n+1)*v
        for i,v in enumerate(c0): c2[i]-=F(n,n+1)*v
        c0,c1=c1,c2
    return c1
def mat_inv(M):
    n=len(M); A=[list(r)+[F(int(i==j)) for j in range(n)] for i,r in enumerate(M)]
    for c in range(n):
        p=next(r for r in range(c,n) if A[r][c]); A[c],A[p]=A[p],A[c]; A[c]=[v/A[c][c] for v in A[c]]
        for r in range(n):
            if r!=c and A[r][c]: f=A[r][c]; A[r]=[a-f*b for a,b in zip(A[r],A[c])]
    return [r[n:] for r in A]

# ---------- 1. centered series in Legendre moments ----------
def centered_series(D):
    """Omega through t^D, as a dict word(g-letters) -> Fraction (weights <= D)."""
    B=bernoulli(D+1)
    # Taylor form: letter j = A_j = A^{(j)}(0); key (word, degree in t)
    def tadd(p,q,c=1):
        r=defaultdict(F,p)
        for k,v in q.items(): r[k]+=c*v
        return {k:v for k,v in r.items() if v}
    def tmul(p,q):
        r=defaultdict(F)
        for (w1,d1),v1 in p.items():
            for (w2,d2),v2 in q.items():
                if d1+d2<=D: r[(w1+w2,d1+d2)]+=v1*v2
        return {k:v for k,v in r.items() if v}
    tbr=lambda p,q: tadd(tmul(p,q),tmul(q,p),-1)
    Aev={((j,),j):F(2,factorial(j)) for j in range(0,D+1,2)}   # A(t)+A(-t)
    Aod={((j,),j):F(2,factorial(j)) for j in range(1,D+1,2)}   # A(t)-A(-t)
    Om={}
    for _ in range(D+1):                                        # fixed-point iteration, one order per pass
        rhs=dict(Aev); ad=Aev
        for k in range(1,D+1):
            ad=tbr(Om,ad)
            if not ad: break
            if k%2==0: rhs=tadd(rhs,ad,B[k]/factorial(k))
        rhs=tadd(rhs,tbr(Om,Aod),F(-1,2))
        Om={(w,d+1):v/(d+1) for (w,d),v in rhs.items() if d+1<=D}
    # X_j = t^(j+1) A_j in terms of centered moments: g_k = sum_j X_j/j! * int_{-1}^1 P_k u^j du
    im=lambda p: F(0) if p%2 else F(2,p+1)
    N=D+1
    Cm=[[sum(c*im(i+j) for i,c in enumerate(legendre(k)))/factorial(j) for j in range(N)] for k in range(N)]
    Ci=mat_inv(Cm)                                              # X_j = sum_k Ci[j][k] g_k
    out=defaultdict(F)
    for (w,d),v in Om.items():
        cur={():F(1)}
        for j in w:
            nxt=defaultdict(F)
            for a,ca in cur.items():
                for k in range(N):
                    if Ci[j][k] and wt(a+(k,))<=D: nxt[a+(k,)]+=ca*Ci[j][k]
            cur=nxt
        for b,c in cur.items(): out[b]+=v*c
    return {k:v for k,v in out.items() if v}

# ---------- 2. Lyndon basis (for printing Lie polynomials) ----------
def is_lyndon(w): return all(w<w[i:]+w[:i] for i in range(1,len(w)))
def std_fact(w):
    for i in range(1,len(w)):
        if is_lyndon(w[i:]): return w[:i],w[i:]
def expand(w,W=99):
    if len(w)==1: return {w:F(1)}
    u,v=std_fact(w); return br(expand(u),expand(v),W)
def bracket(w):
    if len(w)==1: return 'g%d'%w[0]
    u,v=std_fact(w); return '[%s,%s]'%(bracket(u),bracket(v))
def lie_decompose(poly):
    """Coefficients of a Lie polynomial in the Lyndon (standard-bracketing) basis; exact."""
    rest=dict(poly); out=[]
    while rest:
        w=min(rest,key=lambda w:(wt(w),len(w),w))             # smallest word is Lyndon with leading coeff 1
        assert is_lyndon(w), w
        c=rest[w]; out.append((w,c)); rest=add(rest,expand(w),-c)
    return out

# ---------- 3. commutator formulas ----------
def evaluate(formula,G,W):
    """formula: dict 'Xi.name'/'Yi.name'/'ci' -> Fraction, names g0..,C1..; G: list of polys (moments)."""
    pool={'g%d'%k:G[k] for k in range(len(G))}
    n=max(int(k[1:].split('.')[0]) for k in formula if k[0] in 'XY')
    for i in range(1,n+1):
        X={}; Y={}
        for k,v in formula.items():
            if k.startswith('X%d.'%i): X=add(X,pool[k.split('.')[1]],v)
            if k.startswith('Y%d.'%i): Y=add(Y,pool[k.split('.')[1]],v)
        pool['C%d'%i]=br(X,Y,W)
    out=dict(G[0])
    for k,v in formula.items():
        if k[0]=='c': out=add(out,pool['C'+k[1:]],v)
    return out
def letter(k): return {(k,):F(1)}
def quad_moments(nodes,W_k,M):
    """Moments of a rule g_k ~ 2t*sum_i W_k[i] A(u_i t) (as in [1], Eq. 10, where W_k[i] = w_i P_k(u_i)
    with sum_i w_i = 1), expressed in exact moments g_0..g_{M-1} (exact when the u_i, W_k[i] are rational)."""
    G=[]
    for Wk in W_k:
        g={}
        for m in range(M):
            q=(2*m+1)*sum(w*sum(c*u**e for e,c in enumerate(legendre(m))) for u,w in zip(nodes,Wk))
            if q: g=add(g,letter(m),q)
        G.append(g)
    return G

ORDER4={'X1.g0':F(1),'Y1.g1':F(1),'c1':F(-1,2)}
ORDER6={'X1.g0':F(1),'Y1.g1':F(1),'X2.g0':F(1),'Y2.g2':F(1),'Y2.C1':F(1,10),
        'X3.g0':F(1),'X3.g2':F(-1),'X3.C1':F(-3,10),'Y3.g1':F(1),'Y3.C2':F(-1,6),'c3':F(-1,2)}

ORDER8={'X1.g0':F(-1,2),'X1.g2':F(5,7),'Y1.g1':F(1),'Y1.g3':F(-1),'X2.g0':F(1),'X2.g2':F(-5,2),'Y2.g2':F(5),'Y2.C1':F(-1),
 'X3.g0':F(1,2),'X3.g2':F(-5,7),'X3.C1':F(3,14),'Y3.g1':F(-2),'Y3.g3':F(2),'Y3.C2':F(1,21),'X4.C1':F(3,14),'Y4.g1':F(-1),
 'Y4.g3':F(7,2),'X5.g0':F(2),'X5.g2':F(-5),'X5.C1':F(15,14),'Y5.g2':F(10),'Y5.C3':F(-1),'Y5.C4':F(1),'X6.g0':F(1),
 'X6.C1':F(-1),'X6.C3':F(1),'Y6.g1':F(3),'Y6.g3':F(7),'Y6.C2':F(1,2),'Y6.C5':F(-1,6),'c3':F(7,20),'c6':F(-1,20)}

if __name__=='__main__':
    W=9
    S=centered_series(W)
    for w,c in lie_decompose({k:v for k,v in S.items() if wt(k)<=7}): print('%8s  %s'%(c,bracket(w)))
    for name,form,order in (('order 4',ORDER4,4),('order 6',ORDER6,6),('order 8',ORDER8,8)):
        G=[letter(k) for k in range(order//2)]
        R=add(evaluate(form,G,order+1),{k:v for k,v in S.items() if wt(k)<=order+1},-1)
        print(name,'exact moments: mismatches through weight',order,'=',sum(1 for k in R if wt(k)<=order))
    # Newton-Cotes, order 6: 5 equally spaced nodes, interpolatory weights (1/2) int P_k L_i
    u=[F(-1),F(-1,2),F(0),F(1,2),F(1)]
    im=lambda p: F(0) if p%2 else F(2,p+1)
    V=mat_inv([[x**j for j in range(5)] for x in u])            # L_i(x) = sum_j V[j][i] x^j
    Wk=[[sum(V[j][i]*sum(c*im(j+e) for e,c in enumerate(legendre(k))) for j in range(5))/2 for i in range(5)] for k in range(3)]
    G=quad_moments(u,Wk,9)
    R=add(evaluate(ORDER6,G,7),{k:v for k,v in S.items() if wt(k)<=7},-1)
    print('order 6 Newton-Cotes moments: mismatches through weight 6 =',sum(1 for k in R if wt(k)<=6),
          '; weight-7 residual terms =',sum(1 for k in R if wt(k)==7))
