<!--
  Magnus_quadrature_research_guide.md
  Research notes and reproduction guide for centered Magnus quadrature formulas of orders 2-10.

  Companion documents (cite these for the results; this file records the method):
    [1] K. C. Johnson, "The centered Magnus exponent and Magnus quadrature formulas of orders 2, 4, 6 and 8,"
        Version 02-Nov-2026, magnus_order_2-8.pdf
    [2] K. C. Johnson, "Magnus quadrature formulas of order 10 on a centered interval,"
        Version 02-Nov-2026, magnus_order_10.pdf
    Repository:       https://github.com/kjinnovation/explde
    Coefficient file: magnus10coef.m  ([x,y] = magnus10coef(); Appendices A and B of [2])
    Reference code:   magnus_ref.py   (reproduced in full in the Appendix below)

  Author: Kenneth C. Johnson, ken@kjinnovation.com
  Prepared with the assistance of Claude Opus 5.5 (Anthropic). Version 02-Nov-2026.
-->

# Centered Magnus quadrature in Legendre moments: research notes and reproduction guide

This file describes how the Magnus quadrature formulas in [1] and [2] were derived and verified, and what was learned along the way, so that other researchers can reproduce the work and extend it, for example in sessions with Claude or another AI assistant. The papers give the results. This file gives the procedure, the checks, the dead ends, and open questions.

## 1. Problem and notation

The linear system $F'(t)=A(t)\,F(t)$ is advanced over a centered interval $[-t,t]$ (step $h=2t$) by $F(t)=\exp(\Omega(t))\,F(-t)$. The propagator $U(t)=V(t)V(-t)^{-1}$ satisfies $U'=A(t)U+UA(-t)$ and $U(-t)=U(t)^{-1}$, so $\Omega$ is odd in $t$. (Here $V$ is the conventional propagator, $V'=AV$, $V(0)=I$. Do not write $F(t)F(-t)^{-1}$: $F$ is a vector.)

**Centered Magnus theorem** ([1], Eq. 34):

$$\Omega'=\sum_{k\ \mathrm{even}}\frac{B_k}{k!}\,\mathrm{ad}_\Omega^k\big(A(t)+A(-t)\big)-\tfrac12\big[\Omega,\;A(t)-A(-t)\big].$$

**Centered Legendre moments** ([1], Eq. 35): $g_k(t)=\int_{-t}^{t}P_k(\tau/t)\,A(\tau)\,d\tau$. $g_k=O(t^{k+1})$, *weight* $k+1$; $g_0,g_2,g_4$ are odd in $t$, $g_1,g_3$ even. They equal the Legendre moments over $[0,h]$ (with $P_k(2\tau/h-1)$) for $h=2t$, so formulas derived on $[0,h]$ about the midpoint carry over unchanged.

**Series** ([1], Eq. 41; [2], Eq. 2). The only term outside commutators is $g_0$; only odd weights occur; $g_k$ first appears at weight $2k+1$. An order-$2n$ formula (residual $O(t^{2n+1})$) needs $g_0,\dots,g_{n-1}$ and must match the series at odd weights $\le 2n-1$.

| Order | Moments | Order conditions (weights 3, 5, 7, 9) | Total | Commutators used |
|---|---|---|---|---|
| 2 | $g_0$ | - | 0 | 0 |
| 4 | $g_0,g_1$ | 1 | 1 | 1 |
| 6 | $g_0..g_2$ | 1, 5 | 6 | 3 |
| 8 | $g_0..g_3$ | 1, 5, 14 | 20 | 6 |
| 10 | $g_0..g_4$ | 1, 5, 16, 48 | 70 | 9 |

(Counts are dimensions of the free Lie algebra on the moments at each weight, not counting single letters.)

## 2. Results in brief

* **Order 4:** $\Omega=g_0-\tfrac12[g_0,g_1]$.
* **Order 6** ([1], Eq. 47): $C_1=[g_0,g_1]$, $C_2=[g_0,\,g_2+\tfrac1{10}C_1]$, $C_3=[g_0-g_2-\tfrac3{10}C_1,\;g_1-\tfrac16C_2]$, $\Omega=g_0-\tfrac12C_3$. Its leading-error norm is 11% below that of Blanes-Casas-Ros (BCR 2002).
* **Order 8** ([1], Eq. 48): the BCR 2002 six-commutator scheme rewritten in Legendre moments (`ORDER8` in the Appendix code). A search over all six-commutator parity patterns found six solution families (5-12 free parameters); minimizing the leading error gave a norm 1.8 times smaller but no improvement on tests, so BCR's scheme is effectively optimal.
* **Order 10** ([2]): nine commutators, parity pattern OE, OE, OO, OE, OO, OE, OO, OE, OE; 81 coefficients $x_k$ (25 rational, 56 irrational algebraic numbers). No eight-commutator formula was found. Coefficients: Appendix A of [2] and `magnus10coef.m`.
* **Embedded pair** ([2]): on 7 Gauss-Kronrod nodes, an order-6 estimate from the 3 Gauss nodes (3 commutators, rational) and an order-10 estimate with 7 more commutators (96 coefficients $y_k$; 35 rational, 61 irrational). $\|M_{10}-M_6\|$ estimates the error.

**Quadrature criterion** ([1]): a rule of polynomial degree $2n-1$ (exact when $P_k A$ is a polynomial of degree $2n-1$) gives moments accurate enough for order $2n$, because the moment errors raise weights by at least $2n-2k$ and $g_k$ appears only at weight $\ge 2k+1$.

| Rule | Order 2 | 4 | 6 | 8 | 10 | New samples/step | Notes |
|---|---|---|---|---|---|---|---|
| Gauss-Legendre | 1 | 2 | 3 | 4 | 5 | $n$ | weights on products $P_kA$ |
| Newton-Cotes | 2 | 3 | 5 | 7 | - | nodes − 1 | interpolatory or product weights both work; split step reuses all nodes |
| Gauss-Lobatto | 2 | 3 | 4 | 5 | 6 | $n$ | **product weights required** (interpolatory fails for $g_1$ when $n\ge3$) |
| Nine-node "modified Newton-Cotes" | | | | | 9 | 8 | $u=0,\pm\frac13,\pm\frac23,\pm\frac56,\pm1$ (sixths grid; avoids Runge / alternating NC-9 weights) |

**Gauss-Kronrod pairs** (error estimation, $2n+1$ new samples per step): 3 nodes give orders 2 and 6; 5 nodes ($0,\pm1/\sqrt3,\pm\sqrt{6/7}$; weights $14/45, 27/110, 49/495$, normalized to sum to 1) give orders 4 and 8; 7 nodes give orders 6 and 10 (the embedded pair of [2]).

## 3. Pipeline, step by step

Each step lists what was done and the check that confirmed it. The Appendix has working code for steps 1, 2, 3 and 7 (series, Lyndon basis, formula evaluation, quadrature substitution).

**Step 1 - Series.** Expand $A(\tau)=\sum_j A_j\tau^j/j!$, iterate the centered Magnus theorem as a fixed point in exact rationals (one order of $t$ per pass), then substitute $X_j=t^{j+1}A_j$ in terms of $g_k$ by inverting $g_k=\sum_j X_j\,m_{kj}/j!$, $m_{kj}=\int_{-1}^1P_k(u)u^j\,du$. *Check:* the result agrees term by term with the conventional midpoint series over $[0,h]$ through weight 9, and has only odd weights. Takes under a second through $t^9$ in pure Python; $t^{11}$ is feasible.

**Step 2 - Lyndon basis.** Rewrite the series as a Lie polynomial in the Lyndon basis (standard bracketing) by repeatedly peeling off the smallest remaining word. The coefficient of a Lyndon word $\ell$ in the word expansion of a Lie polynomial determines it uniquely, so the order conditions can be imposed as "coefficient of each odd-weight Lyndon word equals the series coefficient".

**Step 3 - Ansatz.** $C_i=[X_i,Y_i]$, where $X_i,Y_i$ are linear combinations of the moments and earlier commutators, and $\Omega=g_0+\sum c_iC_i$. Require each operand to have definite parity (odd: $g_0,g_2,g_4$ and odd commutators; even: $g_1,g_3$ and even commutators). Then every even-weight term vanishes automatically and only odd-weight conditions remain. OO and EE brackets are even and appear only as operands ($c_i=0$). Each bracket has two rescaling ("gauge") freedoms ($X_i\to aX_i$, $Y_i\to bY_i$, compensated downstream), which the solver should fix or exploit.

**Step 4 - Minimum number of commutators.** For each parity pattern, compute the generic rank of the Jacobian of the order conditions with respect to the free coefficients at random points (a necessary condition for solvability is rank = number of conditions, but solutions may lie on special sets where the rank drops). Then search for exact solutions with Levenberg-Marquardt (analytic Jacobian, many random restarts, continuation from larger solutions by driving one bracket to zero). Order-10 statistics ([2]): 8 commutators, max generic rank 64 of 70 over all 2,187 patterns, 322 restarts + 9 continuations, best residual $9\times10^{-4}$; 9 commutators, rank 69 over 6,561 patterns, exact solutions in 3 families (Jacobian rank 55-56 at the solutions, 16-20 free parameters beyond 24 gauge freedoms). The failure at 8 is numerical evidence, not a proof.

**Step 5 - Choosing a member of a family.** Minimize the leading-error coefficients (the residual at weight $2n+1$, in the word basis) with each moment scaled by its leading size ($g_k\approx c_{kk}\,t^{k+1}A^{(k)}(0)/k!$, $c_{kk}=2^{k+1}(k!)^2/(2k+1)!$ on $[-t,t]$). Part of the leading error is a fixed floor (terms with $g_n$ and quadrature error) that no coefficient choice removes. Then fix coefficients one at a time to simple rationals (0, ±1, small fractions), re-solving after each fix, while keeping the weighted error within ~25% of the minimum; finally apply gauge rescalings to simplify. At order 8 this recovers fully rational formulas; at order 10 some coefficients always ended up determined at irrational values, in every fixing order tried.

**Step 6 - Root polishing and conditioning.** With the rationals fixed, the remaining coefficients satisfy a polynomial system with rational coefficients (order 10: 70 equations in 56 unknowns, degree ≤ 13, 4,951 terms). Polish by Gauss-Newton, $x\leftarrow x-(J^TJ)^{-1}J^TE(x)$, evaluating $E$ and $J$ from exact rational coefficients in high-precision decimal arithmetic (320 digits). Convergence is quadratic (from 17 digits: $10^{-15}\to10^{-31}\to10^{-53}\to10^{-107}\to10^{-206}$). Full column rank of $J$ shows the solution is isolated, so the coefficients are real algebraic numbers. Order 10: $\kappa(J)=2.0\times10^8$ ($2.3\times10^7$ column-scaled); Kronrod pair: $78\times61$, $\kappa=3.2\times10^{10}$ ($1.9\times10^8$ scaled). The condition number limits recovering the coefficients from residuals, not the use of the formula: rounding the coefficients to double changes the residuals by about $10^{-15}$. Integer-relation (LLL) tests excluded minimal polynomials of degree ≤ 12 with coefficients below about $10^{10}$.

**Step 7 - Quadrature rules.** For a rule $g_k\approx 2t\sum_iW_{ki}A(u_it)$ (the form of [1], Eq. 42), express each quadrature moment in exact moments, $\tilde g_k=\sum_m Q_{km}g_m$ with $Q_{km}=(2m+1)\sum_iW_{ki}P_m(u_i)$, substitute into the formula, and compare every word through weight $2n$ (order 10: letters $g_0..g_9$, weight ≤ 10) with the series. Product weights are $W_{ki}=w_iP_k(u_i)$, with the weights normalized so that $\sum_iw_i=1$ as in [1] and [2]; interpolatory weights are $W_{ki}=\frac12\int_{-1}^1P_kL_i$. A deliberately weak rule (e.g. 5 equally spaced nodes at order 10) should fail this test, which is a useful control.

**Step 8 - Numerical tests** ([1], Results). $4\times4$ matrix polynomials $A(\tau)=\sum_{k=0}^{12}C_k\tau^k$, $C_k=s\,d^k N_k$, $N_k$ standard normal from `numpy.random.default_rng(seed)`, with (seed, $s$, $d$) = (42, 1, 0.8) random, (9, 2, 0.5) commutator-heavy, (7, 0.3, 1.5) rapidly varying. Reference: $\Omega=\log(V(t)V(-t)^{-1})$ with $V$ from its Taylor series (200 terms) and `scipy.linalg.logm`. Use $t=0.4,0.2,0.1,0.05$; error ratios per halving should approach $2^{2n+1}$. At order 10 the rapidly varying test reaches rounding level (~$2\times10^{-15}$) at $t=0.05$.

**Step 9 - Documents.** The papers were generated from Python content files into LaTeX and LyX (format 544) so that equations are editable in LyX; coefficients and tables were inserted from the verified data, never typed by hand.

## 4. Lessons and pitfalls

* **Verify every formula exactly**, on all words through the target weight, with exact moments *and* with each quadrature rule. Numerical order tests alone can miss a wrong coefficient that only affects one problem class.
* **Gauss-Lobatto weights must be applied to the products** $P_kA$. Interpolating through the $n+1$ samples is exact only for $A$ of degree $n$, which fails for $g_1$ at orders ≥ 6. For Newton-Cotes both product and interpolatory weights work (accuracy within a few percent); an earlier claim that product weights lose accuracy was wrong.
* **The Kronrod pair needs the Kronrod nodes**: a degree-9 seven-node set with endpoints suffices for the single order-10 formula but not for the embedded pair, which uses $g_4$ at weight 7.
* **One matrix product per commutator** ($YX=(XY)^\dagger$) works only when the operands stay in a Lie algebra closed under the relevant symmetry: skew-Hermitian, real skew-symmetric, or Hamiltonian ($YX=-J(XY)^TJ$). It does **not** work for Hermitian or symmetric $A$ (their commutators are skew).
* **Exact moments can be much better** than quadrature moments when $A$ varies rapidly relative to its non-commutativity (3-7 times in the tests), in agreement with BCR 2002.
* **Step splitting:** Newton-Cotes (and the nine-node sixths-grid rule) reuse all parent samples, halving the new evaluations; Gauss-Lobatto reuses endpoints (and the midpoint for odd node counts); Gauss-Legendre reuses nothing.
* **Literature check:** the centered representation of the step itself appears in Blanes, Casas, Oteo and Ros (2009), Eqs. 63-66, and Newton-Cotes moment formulas appear in BCR 2000 and the 2009 review; read the primary sources before claiming novelty.
* **Report numbers carefully:** two statements in an earlier internal summary were wrong (a singular value described as relative when it was absolute, so $\kappa$ was understated; and a condition count of 84 that is actually 78). Re-derive key numbers from data before publishing them.
* **Computing environment:** the sandbox used had no `mpmath` or `sympy` and no network access to PyPI, so everything was done with `fractions.Fraction` and `decimal.Decimal` (note `Decimal(0)**0` raises an error). Truncated free-algebra products are the main cost; precomputed index tables for word concatenation make Decimal evaluation fast.

## 5. Working with Claude on this kind of research

What worked well:

* **Start each session with a handoff file** (like this one): notation, established results, the exact task, deliverables, and the verification required. Ask Claude to produce such a file at the end of a long session.
* **Ask for exact verification and controls**: "verify on all words through weight $2n$ in exact arithmetic, and show a rule that should fail". Ask for numbers to be generated from data rather than retyped.
* **Ask for the conditioning and convergence history** of any numerical solve, not just the final residual.
* **Supply primary sources** (PDFs of the papers you compare against) and ask Claude to check equation numbers and claims against them.
* **Keep coefficient files machine-readable** (`.m`, JSON, or the polynomial-system text file) and regenerate documents from them.
* **Review the prose for overclaiming**: AI drafts tend to state novelty or accuracy too strongly; challenge each comparative statement.

Example prompts:

* "Using magnus_ref.py, derive the centered series through $t^{11}$ and list the number of order conditions at order 12."
* "Search for an order-12 formula: for each parity pattern with $N$ commutators, compute the generic Jacobian rank of the order conditions; report the smallest $N$ that reaches full rank, then run Levenberg-Marquardt from random starts."
* "Polish these coefficients by Gauss-Newton in 300-digit arithmetic; report the residual per iteration and the singular values of the Jacobian."
* "Verify this quadrature rule against the order-8 formula of [1] by substituting its moments in exact moments and comparing all words through weight 8; include a control rule that should fail."

## 6. Open questions

* Prove (or refute) that eight commutators cannot give order 10 (e.g. by Gröbner bases or a dimension argument on the special stratum).
* Find closed forms or minimal polynomials for the irrational order-10 coefficients, or a family member with all coefficients rational using more commutators.
* Order 12 and beyond: count conditions, find minimum commutator numbers, and check whether irrational coefficients become unavoidable earlier.
* Embedded pairs at order 8 (e.g. on 5 Kronrod nodes) that reuse the low-order commutators, and pairs on Lobatto or sixths-grid nodes suited to step splitting.
* Adaptive step-size control strategies using these pairs and sample reuse.
* Structure-preserving variants (one product per commutator) and their accuracy for Hermitian or non-normal problems.

## 7. References

1. K. C. Johnson, "The centered Magnus exponent and Magnus quadrature formulas of orders 2, 4, 6 and 8," Version 02-Nov-2026 (magnus_order_2-8.pdf).
2. K. C. Johnson, "Magnus quadrature formulas of order 10 on a centered interval," Version 02-Nov-2026 (magnus_order_10.pdf).
3. S. Blanes, F. Casas, J. Ros, "Improved high order integrators based on the Magnus expansion," *BIT* 40 (2000) 434-450.
4. S. Blanes, F. Casas, J. Ros, "High order optimized geometric integrators for linear differential equations," *BIT* 42 (2002) 262-284.
5. S. Blanes, F. Casas, J. A. Oteo, J. Ros, "The Magnus expansion and some of its applications," *Physics Reports* 470 (2009) 151-238.
6. A. S. Kronrod, *Nodes and Weights of Quadrature Formulas: Sixteen-Place Tables*, Consultants Bureau, New York, 1965.

## Appendix: reference code (magnus_ref.py)

Self-contained (Python 3, standard library only). `python3 magnus_ref.py` prints the series through weight 7 in the Lyndon basis (it should match [1], Eq. 41) and verifies the order-4, 6 and 8 formulas with exact moments and the order-6 formula with five-node Newton-Cotes moments (expected output: 0 mismatches each).

```python
"""magnus_ref.py -- reference implementation for centered Magnus quadrature in Legendre moments.

Companion code for
  [1] K. C. Johnson, "The centered Magnus exponent and Magnus quadrature formulas of orders 2, 4, 6
      and 8," Version 02-Nov-2026 (magnus_order_2-8.pdf).
  [2] K. C. Johnson, "Magnus quadrature formulas of order 10 on a centered interval,"
      Version 02-Nov-2026 (magnus_order_10.pdf); coefficients in magnus10coef.m.
Repository: https://github.com/kjinnovation/explde
Written with the assistance of Claude Opus 5.5 (Anthropic).

What it does (pure Python, exact rational arithmetic with fractions.Fraction):
  1. centered_series(D): the Magnus exponent Omega(t) through t^D in centered Legendre moments g_k,
     derived from the centered Magnus theorem ([1], Eq. 34).
  2. lie_decompose(): rewrite a Lie polynomial in the Lyndon basis (as printed in [1], Eq. 41).
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
    """Moments of a rule g_k ~ 2t*sum_i W_k[i] A(u_i t) (as in [1], Eq. 42, where W_k[i] = w_i P_k(u_i)
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
```
