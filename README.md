# explde

explde numerically solves the linear matrix differential equation F'(t) = A(t)\*F(t), with constant or non-constant coefficients, using the Magnus exponent method with Gauss-Legendre quadrature of accuracy order 2, 4, 6, 8, or 10 in the integration step size. Paged operation is supported for parallel ensemble computations.

Author: Kenneth C. Johnson (KJ Innovation), ken@kjinnovation.com

Also on MATLAB File Exchange: https://www.mathworks.com/matlabcentral/fileexchange/183695-explde

## Files

| File | Contents |
|---|---|
| `explde.m` | Solver: `F = explde(A,t,ord,F0)` for constant A (array) or variable A (function handle) |
| `magnus_exp.m` | Magnus exponent over one step, orders 2, 4, 6, 8, 10 (used by explde) |
| `magnus10coef.m` | Coefficients of the order-10 formulas (`[x,y] = magnus10coef()`) |
| `demo_explde.m` | Examples: constant A, and a variable-A test case with an exact solution |
| `docs/magnus_order_2-8.pdf` | The centered Magnus exponent and Magnus quadrature formulas of orders 2, 4, 6 and 8 |
| `docs/magnus_order_10.pdf` | Magnus quadrature formulas of order 10 on a centered interval |
| `docs/*.lyx`, `docs/*.tex` | LyX and LaTeX sources of the two documents |
| `research/Magnus_quadrature_research_guide.md` | How the formulas were derived and verified, and how to continue the research |
| `research/magnus_ref.py` | Python reference code (exact rational arithmetic) for the series and formula verification |

## Requirements

- MATLAB R2020b or later (uses `pagemtimes`).
- `pageexpm` from MATLAB File Exchange: https://www.mathworks.com/matlabcentral/fileexchange/107959-pageexpm

Put `explde.m`, `magnus_exp.m`, `magnus10coef.m` and `pageexpm.m` on the MATLAB path.

## Quick start

```matlab
% Constant A: explde(A,[0,t]) equals expm(A*t)
F = explde([1,2;2,1],[0,1]);

% Variable A: F'(t) = A(t)*F(t), F(0) = F0, ten steps of 0.1, order 8
J = [0,-1;1,0];
B = [1,1;1,1];
R = @(t) [cos(t),-sin(t);sin(t),cos(t)];
A = @(t) J+R(t)*B*R(-t);              % exact solution: F(t) = R(t)*expm(t*B)*F0
F0 = [1;0];
t = linspace(0,1,11);
[F,F_] = explde(A,t,8,F0);            % F at t(end); F_{j} at t(j)
err = F-R(1)*expm(B)*F0
```

See `demo_explde.m` for error plots at all orders and for paged (ensemble) operation, and `help explde` for the full syntax.

## Method

Each step from t0 to t1 is advanced by F(t1) = expm(Omega)\*F(t0). Omega is computed from Legendre moments of A over the step, which are evaluated by Gauss-Legendre quadrature and combined in nested commutators:

| Order | Samples of A per step | Commutators | Local error |
|---|---|---|---|
| 2 | 1 | 0 | O(dt^3) |
| 4 | 2 | 1 | O(dt^5) |
| 6 | 3 | 3 | O(dt^7) |
| 8 | 4 | 6 | O(dt^9) |
| 10 | 5 | 9 | O(dt^11) |

The documents in `docs/` give the derivation, the formulas, the alternative quadrature rules (Newton-Cotes, Gauss-Lobatto, Gauss-Kronrod error-estimating pairs), and test results.

## Citation

K. C. Johnson, "The centered Magnus exponent and Magnus quadrature formulas of orders 2, 4, 6 and 8," and "Magnus quadrature formulas of order 10 on a centered interval," https://github.com/kjinnovation/explde

## License

BSD 3-Clause; see `LICENSE`.

## Acknowledgment

This work was carried out with the assistance of Claude Opus 5.5 (Anthropic), a large language model. Under the author's direction, Claude derived the Magnus series in Legendre moments, searched for minimal-commutator formulas and optimized their leading error, and drafted portions of the text.
