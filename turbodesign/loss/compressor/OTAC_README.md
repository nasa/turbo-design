OTAC compressor loss models
===========================

The classes in `otac.py` mirror the legacy OTAC `*.int` loss models. They are
currently placeholders that return zeros and emit a warning so the loss API can
be wired without failing.

When you translate an OTAC model, these mappings are a good starting point:

- `Fl_IR` → upstream `BladeRow` (inlet to the current row)
- `Fl_OR` → current `BladeRow` (outlet of the current row)
- `V`, `W`, `U`, `Vz`, `Vr`, `Vtheta` → `row.V`, `row.W`, `row.U`, `row.Vx`, `row.Vr`, `row.Vt`
- `alpha`, `beta` (flow angles) → `row.alpha1/alpha2` (absolute), `row.beta1/beta2` (relative)
- `ts`, `Tt`, `T` → static temperature `row.T`, total temperature `row.T0`
- `Ps`, `Pt`, `P` → static pressure `row.P`, total pressure `row.P0`
- `rho`, `rhos`, `rhot` → density `row.rho`
- `mu` → dynamic viscosity `row.mu`
- `ht` (total enthalpy) → `row.Cp * row.T0`
- `MN`, `M` → Mach number `row.M` (absolute) or `row.M_rel` (relative)
- `radius`, `radiusTipInlet`, `radiusExit` → `row.r` entries (`row.r[0]` hub, `row.r[-1]` tip)
- `bwidth`, `pitch`, `chord`, `throat` → `row.pitch`, `row.chord`, `row.throat`

Most OTAC models assume enthalpy loss; set the loss type accordingly in the
class `__init__`. If a model clamps the loss (e.g., to 25–50% of the available
enthalpy rise), keep those guards to avoid runaway penalties.

Recommended workflow:

1. Pick a model, open the matching `otac/*.int` file, and translate `calculate`
   into a vectorized NumPy computation inside the corresponding class in
   `otac.py`.
2. Replace the base class with `LossBaseClass` directly once implemented
   (remove `_OTACStub` inheritance) and drop the warning.
3. Use `row`/`upstream` attributes from `BladeRow` for inputs; add any new
   geometry fields you need to `BladeRow` with sensible defaults.
4. Return an array shaped like `row.r` (use `np.full_like(row.r, value)` when
   the loss is a scalar).

This file is intentionally brief; keep notes here as you refine the mappings.
