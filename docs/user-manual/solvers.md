# ODE Solvers

OASIS provides three families of implicit ODE solvers for time integration. The choice of solver affects accuracy, stability, and computational cost.

## Solver Selection

Set `timeIntMethod` in `dataProblem.dat` (line 18) or YAML `solver.method`:

| Value | Solver | Order | Adaptive | Best for |
|---|---|---|---|---|
| 1 | BDF1 | 2 (internal) | Automatic | — |
| 2 | BDFN | 1–6 (configurable) | Optional | Stiff mooring problems |
| 3 | **ESDIRK46** | 4 (embedded 5th) | Optional | General use (**default**) |

## ESDIRK46

Explicit Singly Diagonal Implicit Runge-Kutta method with 6 stages. This is the **recommended default solver**.

Based on: *"Fourth-order Runge–Kutta schemes for fluid mechanics applications"* (Carpenter, 2005).

**Characteristics:**

- 4th-order main scheme with embedded 5th-order error estimator
- Handles stiff and non-stiff problems
- Built-in adaptive time stepping with error control
- No initialisation required (self-starting)

**Configuration:**

```
3               // timeIntMethod (ESDIRK46)
2               // timeIntOrder (unused for ESDIRK)
1               // timeIntAdaptivity (adaptive)
0               // timeIntJacNumStepsMax
1e-5            // timeIntAbsTol
1e-3            // timeIntRelTol
20              // maxIterStep
```

**Adaptivity options** (from `ESDIRK_adaptive_time_step_data.dat`, optional):

| Parameter | Default | Description |
|---|---|---|
| `adaptivity_type` | 2 | Error estimation method (2 = Noventa 2018) |
| `adaptivity_smooth_flag` | 2 | Smoothing (2 = atan smoothing) |
| `rho_min` | 2/3 | Minimum step multiplier to accept |

## BDFN

Backward Differentiation Formula with configurable order $N$ (1 to 6).

**Characteristics:**

- Implicit, unconditionally A-stable (order 1-2), A(α)-stable (order 3-6)
- Higher order = more accurate but potentially less robust
- Stores $N+1$ previous states
- Self-initialises by running $N$ steps of BDF1

**Order coefficients:**

| Order | Formula | Stability |
|---|---|---|
| 1 | $y_{n+1} = y_n + \Delta t \, f_{n+1}$ | A-stable |
| 2 | $y_{n+1} = \frac{4}{3}y_n - \frac{1}{3}y_{n-1} + \frac{2}{3}\Delta t \, f_{n+1}$ | A-stable |
| 3-6 | Higher-order multi-step | A(α)-stable |

**Configuration:**

```
2               // timeIntMethod (BDFN)
3               // timeIntOrder (order 3)
0               // timeIntAdaptivity (fixed step)
0               // timeIntJacNumStepsMax
1e-5            // timeIntAbsTol
1e-3            // timeIntRelTol
20              // maxIterStep
```

!!! note "Fixed Time Step Constraints"
    When `timeIntAdaptivity = 0`:
    
    - `writeTimeStep` must be a multiple of `maxTimeStep`
    - `simulationTime` must be a multiple of `maxTimeStep`

## BDF1

First-order BDF (backward Euler). Internally implemented via the BDF2 class.

- Unconditionally stable but only first-order accurate
- Use `timeIntMethod = 2` with `timeIntOrder = 1` for standalone BDF1

## Common Parameters

All solvers share these convergence parameters:

| Parameter | Description | Typical value |
|---|---|---|
| `timeIntAbsTol` | Absolute error tolerance | $10^{-5}$ |
| `timeIntRelTol` | Relative error tolerance | $10^{-3}$ |
| `maxIterStep` | Max Newton iterations per step | 20 |
| `timeIntJacNumStepsMax` | Steps between Jacobian recomputation (0 = every step) | 0 |

## Solver Guidelines

**General dynamics (free bodies + waves):**

- Use ESDIRK46 with adaptivity ON — most robust, excellent error control

**Stiff mooring problems:**

- Start with ESDIRK46. If convergence issues, try BDFN order 3–4

**Wave excitation analysis (locked bodies):**

- ESDIRK46 with fixed $\Delta t = 0.1$ s

## Convergence Troubleshooting

If the solver fails to converge:

1. Reduce `maxTimeStep` by 1–2 orders of magnitude
2. Increase `maxIterStep` (e.g., 20 → 50)
3. Tighten tolerances (`timeIntAbsTol` $10^{-5} \to 10^{-6}$)
4. Change solver order (BDFN: try order 2 or 3)
5. Switch solver family (ESDIRK46 → BDFN or vice versa)

## Debug Output

When debug flags are enabled, the solver writes:

| File | Columns | Description |
|---|---|---|
| `time_adaptivity_debug.txt` | $t$, $\Delta t$, EWT, LTE, $\rho$, rejected | Adaptive time stepping diagnostics |

## Related Examples

| Example | Solver |
|---|---|
| `solvers/bdf1` | BDF order 1 |
| `solvers/bdfn_order2` | BDF order 2 |
| `solvers/bdfn_order4` | BDF order 4 |
| `solvers/esdirk46` | ESDIRK 4th order |
