# The Non-Dimensional Mass Flow Function: Sizing, Matching, and Sweeping

This page covers how `turbo-design` uses Mattingly's non-dimensional mass flow function ($\tilde m$) across three related capabilities: sizing a single component's annulus before any blade row has been solved, matching a compressor and turbine on a common shaft, and sweeping a fixed geometry across operating points to build a performance map.

---

## 1. The relation

Continuity ($\dot m = \rho A V$) combined with the isentropic relations for a calorically perfect gas collapses into one function of Mach number and $\gamma$ alone:

$$\tilde m(M,\gamma) = M\left(1+\frac{\gamma-1}{2}M^2\right)^{-\frac{\gamma+1}{2(\gamma-1)}}$$

Implemented as `mass_flow_function(M, gamma)` in `turbodesign/isentropic.py`. It carries no gas constant and no absolute scale — the same curve applies to air, combustion products, or any other working fluid, and to any physical size of machine.

$\tilde m$ **peaks once, at $M=1$, and falls off on both sides** — subsonic and supersonic. That peak,

$$\tilde m_{max}(\gamma) = \left(\frac{\gamma+1}{2}\right)^{-\frac{\gamma+1}{2(\gamma-1)}}$$

(`mass_flow_function_max(gamma)`), is the physical ceiling on how much mass a given area can pass — the sonic/choked limit. It is the identity behind everything else on this page:

$$\tilde m(M,\gamma)\cdot \frac{A}{A^*}(M,\gamma) = \tilde m_{max}(\gamma)$$

which ties `mass_flow_function` to the already-used `A_As` area-ratio relation and gives the choke margin below a physical reading.

### Dimensional form

Multiplying by $\sqrt{\gamma/R}$ recovers Mattingly's dimensional mass flow parameter, and inverting it gives the area a given massflow, total state, and target Mach need:

$$MFP(M,\gamma,R) = \sqrt{\gamma/R}\;\tilde m(M,\gamma) \qquad\qquad A = \frac{\dot m\sqrt{T_0}}{P_0\cdot MFP(M,\gamma,R)}$$

`mass_flow_parameter(M, gamma, R)` and `area_for_massflow(massflow, P0, T0, M, gamma, R, blockage=0.0)` in `turbodesign/isentropic.py`. `min_area_for_massflow(...)` is the same relation at $M=1$ — the sonic/choked area.

### Choke margin

$$\text{choke\_margin}(M,\gamma) = 1-\frac{\tilde m(M,\gamma)}{\tilde m_{max}(\gamma)}$$

`choke_margin(M, gamma)` in `isentropic.py`. Because $\tilde m$ peaks at $M=1$ on **both** sides, this margin is $\ge 0$ everywhere and exactly $0$ only at $M=1$ — it does not by itself distinguish subsonic from supersonic. Every blade row gets this computed automatically each solve (`flow_math.update_choke_diagnostics`, wired into `stator_calc`/`rotor_calc` in both `compressor_math.py` and `turbine_math.py`), using the relative-frame Mach for rotors and the absolute-frame Mach for stators/IGVs — populating `row.mass_flow_function`, `row.choke_margin`, and `row.choke_margin_min`.

### The feasibility guard

Because no exit angle or pressure split can beat $\tilde m_{max}$, a target massflow that exceeds it at a row's current area/$P_0$/$T_0$ can never converge — regardless of angle-matching or pressure-balance mode. `flow_math.assert_flow_capacity` checks this for every row up front, in both `CompressorSpool.balance_pressure`/`_angle_match` and `TurbineSpool._balance_pressure`/`_angle_match`, and raises a `ValueError` naming the row, the required vs. maximum $\tilde m$, and the minimum area that would fix it — before the solver spends any iterations discovering the same thing the hard way. (`flow_math.explain_infeasible_massflow` produces the identical diagnosis at the one place a bounded Mach solve can still fail *inside* `initialize()`, before the guard's own checkpoint runs.)

---

## 2. Sizing a component

`turbo-design` solves stage-by-stage radial equilibrium — it never assumes an overall pressure ratio up front. What it *does* need, before it can take a single step, is a physical annulus to integrate across. $\tilde m$ supplies exactly that at the two stations where no blade row exists yet to hand you an area.

### Station 1/2 — component inlet

Pick a target massflow, total conditions, and inlet Mach; `area_for_massflow` returns the annulus area, and `flow_math.radii_for_area(area, mean_radius)` splits it into hub and shroud radii around a chosen mean radius — the thin-annulus complement, $h = A/(4\pi r_{mean})$. That pair is the first control point of a `Passage`.

```python
area1 = area_for_massflow(massflow_kg_s, P01_Pa, T01_K, M1, gamma=1.4, R=287.0)
rhub1, rshroud1 = radii_for_area(area1, rmean)
```

Everything downstream of that first station — rotor and stator exit areas, pressure rise, exit angles — is a stage design choice (velocity triangles, loading), not a Mach target, and is left to the normal row-by-row design.

### Station 4 — the NGV throat, choked

A first-stage nozzle guide vane is the engine's metering orifice, and it runs choked ($M=1$) across essentially the entire normal operating range — the critical pressure ratio for combustion gas ($\gamma\approx1.33$) is only about $1.85$, so turbine inlet-to-throat ratios exceed it from idle through max power. That makes the sizing relation exact, not an approximation:

$$A_4 = \frac{\dot m_4\sqrt{R\,T_{04}}}{P_{04}\cdot \tilde m_{max}(\gamma)}$$

`ComponentSizingSpec.ngv_choked` (default `True`, in `turbodesign/shaft_match.py`) drives station 4's sizing Mach to exactly $1.0$; setting it `False` falls back to a user-supplied `inlet_mach` for a non-metering station instead.

---

## 3. Component matching

`ShaftMatch` (`turbodesign/shaft_match.py`) sizes and solves a turbine to match an already-solved `CompressorSpool` on a common shaft. Currently only that direction is implemented (compressor known, turbine sized) — the classic gas-generator problem.

**Phase 1 — `size()`**, pure algebra, no solver: shaft power balance and mass continuity with fuel/bleed,

$$\dot m_4 = \dot m_2\cdot(1+\text{far}-\text{bleed}) \qquad\qquad P_{turbine}\cdot\eta_{mech} = P_{compressor}$$

(`turbine_massflow(compressor_massflow, fuel_air_ratio, bleed_fraction)`), combustor transition ($T_{04}$ specified, $P_{04}=P_{03}(1-\text{loss})$), and $\tilde m$-based area sizing at each station — including the choked NGV throat from §2. Self-consistency is checked before returning: every sized station's area, fed back through `Massflow(P0, T0, A, M, gamma, R)`, reproduces the target massflow exactly.

**Phase 2 — `build()`/`match()`**: a user-supplied *reference* turbine's blade angles, loss coefficients, and stage layout are reused verbatim (MFP gives an area, not a blade shape — synthesizing angles from scratch is explicit follow-up work); only its annulus is rescaled to hit the sizing target.

The matching solve is genuinely two-dimensional, not one. For a fixed-angle turbine in pressure-balance mode, `spool.massflow` is only a solver **seed** — the achieved massflow and power are both outputs of the annulus geometry and the boundary pressures, not independently dictated. `match()` therefore searches two knobs together: annulus area (sets flow capacity, i.e. massflow) and exit static pressure (sets the expansion ratio, i.e. power), as nested bounded 1-D searches, until both the massflow and power residuals close.

---

## 4. Sweeping an operating map

`sweep_operating_points` (`turbodesign/operating_map.py`) mutates a fixed, already-designed spool in place across a range of $(\dot m, N)$ pairs and collects the resulting pressure ratio, efficiency, and choke margin at each — the standard compressor/turbine characteristic. It reuses the same choke-margin diagnostic from §1: a point whose target massflow exceeds $\tilde m_{max}$ at the swept condition is recorded as `converged=False` with the reason, and the sweep continues rather than aborting.

```python
points = sweep_operating_points(spool, massflow_points, rpm_points)
```

`spool.set_rpm(rpm)` pushes shaft speed onto every row (`rpm` is only stamped at construction time — `solve()` never re-reads `self.rpm` on repeat calls), and `flow_math.reset_rotor_power` clears each rotor's `power`/`power_mean` before every re-solve, since a stale value from a very different prior operating point can degrade the next solve's initial guess.

**Meaningful for `CompressorSpool`, not (by itself) for `TurbineSpool`.** A compressor's `massflow_points` genuinely drives its row-to-row pressure balance. A turbine in pressure-balance mode has the same "massflow is just a seed" property described in §3 — sweeping `massflow_points` alone can converge to the *same* physical point every time, since nothing forces the achieved massflow to track the request. Sweeping `rpm_points` is meaningful for either machine, since it changes blade speed directly. A genuine turbine massflow/pressure-ratio map would need the exit pressure to vary too, the same two-knob idea `ShaftMatch.match()` already uses — not implemented in the sweep today.

**The loss model determines whether the map is trustworthy away from design.** `sweep_operating_points` has no knowledge of loss, incidence, or Mach — `LossBaseClass.__call__(row, upstream)` already re-evaluates from each row's live solved state every call, so an incidence-aware loss model makes a sweep physically honest with zero sweep-side changes. Today, axial rows generally use `FixedPressureLoss`/`FixedPolytropicEfficiency` — constants, not responsive to off-design incidence — so a swept map's massflow/choke axis is honest immediately (pure continuity), but the efficiency/pressure-ratio axis silently reuses the design-point loss coefficient at every point. A real axial incidence-loss correlation (candidates: port `NASA23B20` from `turbodesign/loss/turbine/otac.py`, or translate the centrifugal module's working `ImpellerIncidenceConrad` in `turbodesign/centrifugal/losses.py`) is the prerequisite for a trustworthy efficiency map — separate work, not part of the sweep itself.

---

## See also

- [`mass_flow_function.md`](mass_flow_function.md) (this page)
- [Entropy-Based Efficiency](entropy_based_efficiency.md) — the efficiency definition `overall_polytropic_efficiency` and its turbine counterpart build on
- [`loss.md`](loss.md) — TD2 loss correlation, one of the pluggable loss models referenced in §4
- Tutorial notebook: [`examples/shaft-match/compressor_turbine_match.ipynb`](https://colab.research.google.com/github/nasa/turbo-design/blob/main/examples/shaft-match/compressor_turbine_match.ipynb) — sections 1-4 of this page, run end to end
