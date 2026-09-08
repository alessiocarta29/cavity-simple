# cavity-simple

A finite-volume solver for the steady incompressible Navier–Stokes equations,
applied to the lid-driven cavity. Uniform Cartesian grid, collocated variables,
SIMPLE coupling.

Written in C++ while working through Ferziger, Perić & Street, *Computational
Methods for Fluid Dynamics*. The aim is not speed. The aim is a solver where
every choice can be explained, and where the claims made below are backed by
something in the repository.

## Problem

A square cavity of side `L`. The lid slides at constant speed `U`; the other
three walls are fixed. No-slip everywhere. Steady laminar flow at
`Re = ρUL/μ`.

## Numerical method

| Component | Choice |
|---|---|
| Discretisation | Finite volume, cell-centred |
| Grid | Uniform Cartesian, **collocated** |
| Convection | Upwind in the matrix, deferred correction towards CDS (blending factor `gamma`) |
| Diffusion | Central differences (2nd order) |
| Pressure–velocity coupling | SIMPLE |
| Checkerboard suppression | Rhie–Chow momentum interpolation |
| Linear solver | Stone's SIP, α = 0.92 |
| Streamfunction | Poisson equation `lap(psi) = -omega`, same SIP |
| Walls | Diffusive coefficient `2μ` (the wall sits `Δx/2` from the cell centre), no mass flux |

Four choices deserve an explanation, because they are where a collocated SIMPLE
solver is usually either right or wrong.

### Why Rhie–Chow is needed

SIMPLE needs a mass flux on each face, and the natural way to build one is to
average the velocities of the two cells that share the face:
`u_e = ½(u_P + u_E)`. That average turns out to be unusable, and the reason is
not obvious until the algebra is written out.

Each cell velocity comes from its own momentum equation, and carries the
pressure gradient evaluated at that cell centre. On a uniform grid that gradient
is `(p_E − p_W)/(2Δx)` at `P`, and `(p_EE − p_P)/(2Δx)` at `E`. Averaging the
two velocities averages the two gradients, and what comes out is

    (p_EE + p_E − p_P − p_W) / (4Δx)

Look at which pressures appear. The face between `P` and `E` is straddled by a
difference that reaches two cells out on either side and never uses `p_E − p_P`
on its own. The mass flux through a face is blind to the pressure difference
across that very face.

The consequence is easiest to see with a pressure field that alternates cell by
cell: `+a, −a, +a, −a`. Substituting into the expression above gives exactly
zero. That field produces no face pressure gradient anywhere, so it drives no
mass flux, so the continuity equation cannot see it. Pressure is determined only
up to an arbitrary amount of it, and in practice the solution grows a
checkerboard pattern that nothing removes. This is the odd–even decoupling that
staggered grids avoid by construction and that a collocated grid has to deal
with explicitly.

Rhie–Chow interpolation fixes it by subtracting the averaged gradient from the
face velocity and adding back the compact two-point one:

    u_e = ½(u_P + u_E) − (V/A_P)_e · [ (p_E − p_P)/Δx − ½( (dp/dx)_P + (dp/dx)_E ) ]

For the alternating field the bracket is now `2a/Δx` rather than zero, so the
oscillation finally produces a mass flux, continuity objects, and the mode is
damped. On a smooth field the compact gradient and the averaged one nearly
agree, the bracket is small, and almost nothing changes.

Two things are worth being explicit about.

The bracket is a difference between two approximations of the same derivative,
so it does not vanish as the grid is refined in the sense of being exactly zero —
it is a deliberate inconsistency of order `Δx²` added to the discretisation. That
is the price of using a collocated grid, and it is paid knowingly.

The coefficient `(V/A_P)_f` is not chosen for convenience. The same coefficient
appears in the pressure-correction equation, because the velocity correction
follows from the same momentum equation. The two are consistent by construction,
and that consistency is what makes the pressure correction actually drive the
mass imbalance to zero rather than merely reduce it. It also means the
Rhie–Chow damping inherits whatever is done to `A_P`, including under-relaxation —
see limitation 1 below.

### Why the convection scheme is deferred, not selected

The earlier version of this solver put central-difference coefficients directly
into the matrix. That works until it does not: on a 32 × 32 grid at Re = 1000 it
diverges outright, while the finer grids converge. The cell Peclet number there
is above 30, and §5.8 explains what goes wrong — of the two schemes only upwind
satisfies `A_P ≥ Σ|A_l|`, the sufficient condition for iterative solvers to
converge. Past a cell Peclet number of 2 the central matrix loses that property,
and SIP is left with a system it has no guarantee of handling.

Deferred correction keeps the good matrix and gets the good scheme anyway. The
face value we want is a blend, and it can be rewritten without approximation as

    φ_f = γ·φ_f^CDS + (1 − γ)·φ_f^UDS  =  φ_f^UDS + γ·( φ_f^CDS − φ_f^UDS )

The first term goes into the matrix, so the matrix is always upwind whatever `γ`
is. The second is evaluated with the previous iterate and moved to the source
term. At convergence the outer iteration stops changing the solution, so the
lagged term equals what it would have been implicitly: the explicit upwind
contribution cancels the implicit one, and what is actually solved is the
blended scheme. Lagging a term changes how fast the iteration converges, not
what it converges to.

So the change buys three things. The solver now converges where it previously
crashed. One code path covers CDS, UDS and everything between, instead of two
implementations that can drift apart. And it is the mechanism production codes
use, OpenFOAM included, so the work transfers.

The cost is that the converged solution is no longer identical to what the
central-coefficient matrix produced. Changing the matrix changes `A_P`, and
`A_P` feeds both the Rhie–Chow damping and the pressure-correction equation. At
convergence the mass imbalance vanishes, which leaves `A_P → 4μ + A_wall` for
central differences but `A_P → Σ max(m,0) + 4μ + A_wall` for upwind — much
larger where convection dominates. These are two different discretisations
heading for the same continuum solution. The results section quantifies the gap
and shows it shrinking with refinement, as it must.

### The pressure system is singular

Pressure is only defined up to a constant, so the pressure-correction matrix has
no unique solution. Cell 0 is pinned with `p' = 0`. The compatibility condition
holds because every boundary face is a wall: the mass imbalances of the
individual cells cancel when summed over the domain.

### Under-relaxation

Applied implicitly, `A_P ← A_P/α_u` and `Q ← Q + (1−α_u)·A_P·φ_old`, not by
mixing successive solutions.

## Build and run

```bash
g++ -O2 -o cavity cavity_simple.cpp
./cavity
```

No external dependencies. Parameters are given as `key=value` on the command
line; anything left out keeps its default.

```bash
./cavity N=128 Re=1000 gamma=1 nswUV=3 nswP=10 tol=1e-4
./cavity --spazzate            # inner sweep study
```

| Key | Default | Meaning |
|---|---|---|
| `N` | 128 | cells per direction |
| `Re` | 1000 | Reynolds number |
| `U`, `rho`, `L` | 5, 1000, 1 | lid speed, density, cavity side |
| `au`, `ap` | 0.8, 0.2 | under-relaxation, velocity and pressure |
| `gamma` | 1 | 1 = CDS, 0 = UDS, in between a blend |
| `nswUV`, `nswP` | 3, 20 | SIP sweeps per outer iteration |
| `tol` | 1e-4 | required *relative* residual reduction |
| `maxIter` | 5000 | iteration cap |

Output files: `campo.csv` (one row per cell: `x, y, u, v, p, psi`),
`parametri.csv` (settings, iteration count, residual, wall-clock time, vortex
strengths) and `spazzate.csv` (one row per run of the sweep study).

The figures below are produced by `grafici.py`, which reads those files. It
takes the grid size from the data rather than having it written in, so the
figures stay correct when the grid changes.

## Convergence criterion

For each equation, the residual is the L1 norm of `q − A·φ` taken before the
inner sweeps. It measures how far the current field is from satisfying the
discrete equation. Each residual is divided by its own value at the first outer
iteration, and the largest of the three ratios is compared against `tol`.

Normalising is not cosmetic. A plain sum over cells makes the tolerance depend
on the grid: the same number would mean something different at N = 64 and at
N = 256. Any study comparing iteration counts between configurations would then
be partly measuring its own stopping rule. §5.7 gives the targets: inner
iterations can stop after a reduction of one or two orders of magnitude, outer
iterations should not stop before three to five.



## Results

Reference values from §8.4.1 at Re = 1000: `psi_min = -0.11893`,
`psi_max = 0.00173`.

### What the deferred correction is for

An earlier version of this solver put central-difference coefficients directly
into the matrix. It is faster where it works — and it stops working sooner than
one might expect.

Grid fixed at 64 × 64, Reynolds number increasing, everything else identical
(`nswP = 20`, `nswUV = 3`, `α_u = 0.8`, `α_p = 0.2`, four orders of residual
reduction):

| Re | cell Peclet | CDS in the matrix | deferred correction, `gamma = 1` |
|---|---|---|---|
| 1000 | ~16 | 682 iter | 890 iter |
| 2000 | ~31 | diverges | 1533 iter |
| 3200 | ~50 | diverges | 2244 iter |
| 5000 | ~78 | diverges | 3165 iter |
| 7500 | ~117 | diverges | diverges |

Doubling the Reynolds number on the same grid is enough to lose the central
matrix entirely. The deferred-correction form runs to Re = 5000, an operating
range five times wider.

The governing parameter is the cell Peclet number, not the grid size. At
Re = 2000 on 64 × 64 the cell Peclet number is about 31 — the same value reached
at Re = 1000 on 32 × 32, where the central matrix also diverges. Both failures
happen at the same place on the same axis. §5.8 gives the reason: of the two
schemes only upwind satisfies `A_P ≥ Σ|A_l|`, the sufficient condition for
iterative solvers to converge, and the central matrix loses that property above
a cell Peclet number of 2. That the solver survives to 16 is luck; the condition
is sufficient, not necessary, so exceeding it guarantees nothing either way.

The blending factor extends the range further. At Re = 7500 nothing converges
with `gamma = 1`, but `gamma = 0` converges in 580 iterations. Trading accuracy
for stability is available as a dial precisely because the matrix does not
change when the dial moves. With central coefficients in the matrix there is no
dial.

(At Re = 10000 on this grid nothing converges, including `gamma = 0`. Whether a
steady solution is reachable there at this resolution, or the under-relaxation
is simply too aggressive, has not been established.)

Two things this does **not** buy. It does not improve accuracy: at `gamma = 1`
the scheme solved is the same central scheme, and the agreement with reference
data is unchanged wherever both formulations converge. And it is not free —
lagging a term costs about 30% more outer iterations at Re = 1000, measured
below. The trade is a slower solver over a much wider operating range.

### Iteration counts against the reference

Table 12.1 of the book reports outer iteration counts for this problem on a
single grid with `α_u = 0.8`, `α_p = 0.2`, non-uniform grid, zero initial field,
four orders of residual reduction: 250 iterations at 32², 433 at 64², 1352 at
128².

This solver takes 890 at 64² under nominally matching settings — a factor of
2.1. Part of that is accounted for and part is not.

**Measured: the deferred correction costs a factor of 1.30.** Running the
earlier central-matrix version with the same relative stopping criterion and
identical settings gives 682 iterations against 890. This isolates the cost of
lagging the correction term, since nothing else differs between the two runs.

**Unaccounted: a factor of 1.58 remains** between 682 and 433. Three candidates,
none yet demonstrated:

- *Inner solver effort.* The book stops the inner iterations on a relative
  residual criterion rather than a fixed sweep count. The sweep study in this
  README shows outer iterations ranging from 2996 to 804 as `nswP` alone varies
  from 2 to 50, a factor of 3.7 — larger than the entire gap being explained. If
  the reference solve on the pressure equation is effectively harder than
  `nswP = 20`, this alone could account for most of it. Testable now.
- *Grid.* The reference uses a non-uniform grid, this solver a uniform one. The
  expected sign is unclear: stretched grids have worse aspect ratios and usually
  converge more slowly, not faster. Not testable until the non-uniform grid is
  implemented.
- *Implementation details* not visible in the published description.

The gap that remains unexplained is larger than the part accounted for, so no
conclusion is drawn here. It is recorded as an open question rather than
attributed to a plausible cause.

### The two formulations are not identical, and should not be

The gap between the two should grow with convection and shrink with grid
refinement. It does. Largest difference in `u`, as a fraction of the peak
velocity, at N = 64:

| Re | 10 | 100 | 1000 |
|---|---|---|---|
| difference | 0.05% | 0.50% | 4.81% |

At Re = 1000 the difference drops from 4.81% at N = 64 to 2.48% at N = 128.
These are two different discretisations heading for the same continuum solution,
not two implementations of the same one.

### How much the upwind scheme smears the solution

Same grid (64 × 64), same Reynolds number (1000), both converged to a residual
reduction of nine orders.

![Centreline profiles, CDS against UDS](figures/schemi_profili.png)

The two profiles share the same shape, but the upwind curves are flatter
everywhere. The peaks are cut down and the near-wall gradients are softened.
Nothing here is unstable or noisy — the solution just looks like a lower
Reynolds number than the one requested. That is numerical diffusion: the upwind
scheme adds an artificial viscosity proportional to the cell size, and on this
grid it is not small.

![Streamlines, CDS against UDS](figures/linee_corrente.png)

The streamlines make the cost easier to see. The primary vortex is weaker and
its centre has drifted. The two corner vortices, drawn with the same contour
levels in both panels, are visibly smaller on the right.

| | `gamma = 1` (CDS) | `gamma = 0` (UDS) |
|---|---|---|
| `psi_min` | -0.11351 | -0.09011 |
| `psi_max` | 0.00183 | 0.00104 |

Upwind underestimates the primary vortex by about 24% here, and the corner
vortices by considerably more. Those corner vortices are the demanding
quantity, which is why the book reports both.

### How many inner sweeps to spend

`./cavity --spazzate` varies `nswP` at fixed `nswUV`, then `nswUV` at the best
`nswP`, on two grids.

§12.2.2.1 sets out the trade-off. Solving the inner systems accurately is wasted
effort, because the coefficients will be rebuilt many times before the
non-linear problem is anywhere near solved. But stopping too early pushes the
work onto the outer loop, and the total cost goes up again. Where the optimum
sits depends on the problem, so it has to be measured.

![Outer iterations and wall-clock time against nswP](figures/spazzate.png)

The two curves tell different stories, and that is the whole point. Outer
iterations (blue) fall the whole way: more sweeps on the pressure always means
fewer outer iterations. Wall-clock time (red) has a clear minimum around
`nswP = 10`, then climbs again as each iteration gets more expensive than it is
worth. Anyone measuring only the iteration count would conclude that more sweeps
are always better, which is the wrong answer.

At N = 64:

| `nswP` | outer iterations | wall-clock |
|---|---|---|
| 2 | 2996 | 2.55 s |
| 5 | 1432 | 1.37 s |
| 10 | 1034 | **1.16 s** |
| 20 | 890 | 1.31 s |
| 30 | 841 | 1.59 s |
| 50 | 804 | 2.05 s |

The minimum is at `nswP = 10` on both grids tested. Varying `nswUV` instead, the
outer count stops improving after two or three sweeps while the cost keeps
rising, so there is nothing to gain beyond `nswUV = 3`.

The asymmetry between the two systems is the one described in §7.1.7. The
momentum equations carry convection and have Dirichlet conditions on all four
walls, and under-relaxation adds further to their diagonal — they are easy. The
pressure correction is a pure Poisson problem with no Dirichlet condition
anywhere, singular up to a constant and anchored only by the reference cell —
it is not.

## Validation

Two independent reference solutions are used, and the choice is deliberate.

**Ghia, Ghia & Shin (1982)** is the conventional benchmark: second-order finite
differences on a 129 × 129 uniform grid, coupled strongly implicit multigrid,
tabulated centreline velocities. It is the reference a reader expects to see.

**Erturk, Corke & Gökçöl (2005)** is used because it fails differently. They
solve the streamfunction–vorticity form on a uniform 601 × 601 grid, driving the
maximum absolute residual of the governing equations below 1e-10.

That formulation contains no pressure. There is no pressure–velocity coupling to
get wrong, no Rhie–Chow interpolation, no splitting error from SIMPLE, no
pressure boundary condition. Continuity is satisfied identically by
construction. None of the failure modes of this solver exist in theirs.

### Centreline velocities at Re = 1000

Mean absolute deviation over the tabulated stations. The computed profile is
interpolated onto the stations of each paper, not the other way round: the
tabulated positions are the data, our grid is an arbitrary choice.

| Grid | vs Ghia, u | vs Ghia, v | vs Erturk, u | vs Erturk, v |
|---|---|---|---|---|
| 32 × 32 | 0.03404 | 0.04321 | 0.05183 | 0.06049 |
| 64 × 64 | 0.01031 | 0.01001 | 0.01658 | 0.01878 |
| 128 × 128 | 0.00170 | 0.00398 | 0.00409 | 0.00453 |

All four columns fall monotonically. Against Erturk the successive ratios are
3.1 and 4.1 for u, 3.2 and 4.1 for v — the factor of four expected of a
second-order method.

### Why the agreement with Ghia is not the better result

The deviation from Ghia at 128 × 128 is less than half the deviation from
Erturk. The naive reading is that this solver agrees better with Ghia. Station
by station near the lid tells a different story:

| y | this solver (128²) | Ghia | Erturk | vs Ghia | vs Erturk |
|---|---|---|---|---|---|
| 0.9531 | 0.46789 | 0.46604 | 0.47432 | +0.00185 | −0.00643 |
| 0.9688 | 0.57747 | 0.57492 | 0.58192 | +0.00255 | −0.00445 |
| 0.9766 | 0.66151 | 0.65928 | 0.66747 | +0.00223 | −0.00596 |

This solver and Ghia sit on the same side. Erturk sits on the other, and the gap
between Ghia and Erturk is of the same size as ours. Erturk has closely spaced
stations in this region, so this is genuine disagreement between the references
rather than an interpolation artefact — the shear layer under the moving lid,
where the boundary layer is thinnest.

The explanation is the reason for using two references in the first place. Ghia
is second-order finite differences on 129 × 129; this solver is second-order
finite volumes on 128 × 128. They under-resolve the same layer in the same
direction. A 601 × 601 solution does not. The closer agreement with Ghia
measures a shared failure mode, not accuracy.

Erturk et al. make the same point about their own comparison, describing the
results of Ghia et al. as somewhat under-resolved, and noting that at Re = 1000
even a 401 × 401 second-order grid can be considered so.

### Grid convergence

Three systematically refined uniform grids, refinement ratio 2, all converged to
a relative residual reduction of nine orders. Richardson extrapolation and the
Grid Convergence Index of Roache, with the standard safety factor of 1.25.

| | 32² | 64² | 128² | observed p | extrapolated | GCI (fine) |
|---|---|---|---|---|---|---|
| `psi_min` | -0.10233 | -0.11351 | -0.11748 | 1.494 | -0.11967 | 2.33% |
| `psi_max` | 0.00210 | 0.00183 | 0.00176 | 1.948 | 0.001736 | 1.74% |

Reference values, from three independent sources that agree to four significant
figures: Erturk's own Richardson extrapolation gives -0.118942, the Chebyshev
spectral solution of Botella & Peyret gives -0.1189366, and §8.4.1 of Ferziger
gives -0.11893. Table V of Erturk gives 0.0017281 for the larger secondary
vortex.

Against these, the extrapolated values here are off by 0.61% and 0.46%. The best
single computation, 128², is off by 1.24% on `psi_min`. Extrapolation halves the
error using data that was already available.

The GCI on the fine grid is 2.33% against a true error of 1.24% for `psi_min`,
and 1.74% against 1.73% for `psi_max`. The index estimates the error without
knowing the answer, and in both cases it is correct or conservative.

**The observed orders are diagnostic.** `psi_max` converges at 1.95, essentially
second order. `psi_min` converges at 1.49. A second-order scheme should give two
for both. The discrepancy is consistent with limitation 2 below — the
first-order pressure treatment at the walls — degrading the global quantity
while leaving the corner vortex, which sits away from the moving lid, largely
alone. The grid convergence study identified a defect without any reference
value being involved.

### Iterative convergence

Discretisation error can only be measured if the iterative error is much smaller
than it. On the finest grid the default iteration cap stops the solver at a
residual reduction of seven orders rather than the nine requested. Repeating
that run to full convergence moves `psi_min` from -0.11747 to -0.11748, the
observed order from 1.497 to 1.494, and the extrapolated value by 0.02%.

The contamination is therefore negligible here, but it is shown rather than
assumed. The solver prints a warning when it exits on the iteration cap, and
that warning is the reason the check was made.
## Known limitations

Current properties of the code, not open questions.

1. **Rhie–Chow uses the under-relaxed `A_P`.** The damping therefore scales with
   `α_u`, and the converged solution keeps a weak dependence on the relaxation
   factor (Majumdar, 1988). This has to be fixed before the under-relaxation
   study, or the optimum found will not be comparable with the published one.
2. **Pressure at wall faces uses the cell value as the ghost value.** The
   one-sided gradient this produces is effectively taken over `2Δx` instead of
   `Δx`, so it is first-order at the walls. Linear extrapolation would restore
   second order. Expect this to spoil the observed order in a refinement study.
3. **Uniform grid only.** Non-uniform spacing needs face interpolation factors,
   variable cell volumes, and in particular a weighted interpolation of
   `(V/A_P)_f` inside the Rhie–Chow term.
4. **No stagnation or divergence detection** beyond the iteration cap.

## References

- Ferziger, J. H., Perić, M. & Street, R. L., *Computational Methods for Fluid
  Dynamics*, Springer.
- Ghia, U., Ghia, K. N. & Shin, C. T. (1982), "High-Re solutions for
  incompressible flow using the Navier-Stokes equations and a multigrid method",
  *Journal of Computational Physics*, 48(3), 387–411.
- Rhie, C. M. & Chow, W. L. (1983), "Numerical study of the turbulent flow past
  an airfoil with trailing edge separation", *AIAA Journal*, 21(11), 1525–1532.
- Stone, H. L. (1968), "Iterative solution of implicit approximations of
  multidimensional partial differential equations", *SIAM Journal on Numerical
  Analysis*, 5(3), 530–558.
- Majumdar, S. (1988), "Role of underrelaxation in momentum interpolation for
  calculation of flow with nonstaggered grids", *Numerical Heat Transfer*,
  13(1), 125–132.

## Note on tooling

The solver was written by hand from the equations. An AI assistant was used for
code review and refactoring once the numerical method was working. The
discretisation choices, the debugging decisions and the convergence settings are
the author's own, and are documented in the comments.
