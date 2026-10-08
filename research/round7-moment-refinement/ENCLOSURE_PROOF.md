# Independent enclosure of the frozen bit moment

This is arithmetic analysis of an unchanged histogram, not a new network or
proof that an actual integer multiplication algorithm realizes that histogram.

The pinned original `moment.certify` numerically estimates a root, then checks
a rational upper bound on a denominator10^9 grid. Its decisive inequality is

`exp(a log(m/t)) <= 1/(1−a L_t)`,

where`L_t` is its30-term rational upper bound for`log(m/t)` and`a L_t<1`.
The numerical search chooses a starting grid point only; the final rational
upper-bound check establishes feasibility. It does not establish the exact
root, the largest feasible grid point for the true moment, or optimality of the
histogram. The reciprocal relaxation itself has second-order slack.

The independent checker imports no upstream arithmetic code. For rational
`x>=1`, repeatedly divide by3/2 until`1<=y<3/2`, obtaining

`log x = k log(3/2) + 2 atanh((y−1)/(y+1))`.

Both atanh arguments lie in`[0,1/5]`. After48 terms,

`2 sum_(j=0)^47 z^(2j+1)/(2j+1) <= 2 atanh z`

and the omitted positive tail is at most
`2 z^97/[97(1−z²)]`. Thus two rational endpoints enclose each logarithm.
Rounding the lower endpoint down and upper endpoint up to denominator10^64
preserves enclosure.

For`0<=x_lo<=x_hi<1`, the degree12 Taylor polynomial at`x_lo` is a lower
bound for`exp(x)`. At the upper endpoint, the omitted tail is bounded by

`(x_hi^13/13!)/(1−x_hi/14)`,

since every subsequent term ratio is at most`x_hi/14`. Each positive weighted
term of

`F(a) = (1/W) sum_t c_t (t/m) exp(a log(m/t))`

is separately rounded outward to denominator10^64 before summation. All
decisive calculations use integers and rational numbers; no floating operation
appears in the checker.

The checker also encloses the mathematical reciprocal surrogate using the
same independent log interval. A rejected reciprocal point concerns that
surrogate, not the true exponential moment; these are distinct claims.

`independent_assembly.py` additionally encloses the literal original30-term
log-majorant reciprocal. If the original base2 scaling uses`k` reductions,
its log overestimate is at most

`2(k+1)(1/3)^61/[61(1−1/9)]`.

Adding this explicit allowance to the independent upper log endpoint covers
the original majorant without importing its implementation. The independent
lower log endpoint is also a lower bound for that majorant. Directed
reciprocal bounds then certify its accepted and rejected adjacent points.

The assembly audit reproduces every original floor/ceiling operation with
integers. In particular, independent log endpoints force the same integer
`ceil(576 log119453132304)`, hence`C1=14694`; no floating guard calculation
is trusted. Original epsilon andκ grids remain10^−16 and10^−17. The audit
checks all strict rational assembly conditions and cost margins, and compares
every reported variant's selected parameters with the baseline certificate.

Strict monotonicity of`F(a)` follows because every actual child has
`0<t<m`, positive multiplicity and hence`log(m/t)>0`. Therefore an accepted
grid point and rejected adjacent grid point isolate the maximal feasible
grid point for this **fixed histogram**, provided both strict inequalities
are certified. They do not establish a globally best construction or assembly.
