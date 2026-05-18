# Solving the Radial Equilibrium Equation

**Appendix A derivation — Equations (15) → (38)**

This document walks through the full derivation in the appendix of
*Turbo Design – Open-Source Radial Equilibrium Turbomachinery Solver: Part I. Turbines*,
showing how the radial equilibrium equation (15) is rearranged into an explicit
first-order ODE for the meridional velocity (38) that the solver integrates from
hub to shroud.

Every step below is the algebra implied by the paper, written out in full so the
chain can be checked term by term.

---

## Nomenclature

| Symbol | Meaning | Units |
|---|---|---|
| $P$ | Static pressure | Pa |
| $P_0$ | Total (stagnation) pressure | Pa |
| $T$ | Static temperature | K |
| $T_0$ | Total (stagnation) temperature | K |
| $\rho$ | Density | kg/m³ |
| $V$ | Absolute velocity magnitude | m/s |
| $V_M$ | Meridional velocity | m/s |
| $V_T$ | Tangential (whirl) velocity | m/s |
| $V_r$ | Radial velocity component | m/s |
| $V_{ax}$ | Axial velocity component | m/s |
| $r$ | Radius | m |
| $r_M$ | Meridional radius of curvature of the streamline | m |
| $M$ | Distance along the meridional streamline | m |
| $\alpha$ | Absolute flow angle (between $V$ and $V_M$) | rad |
| $\phi$ | Streamline inclination (between $V_M$ and the axial direction) | rad |
| $\gamma$ | Ratio of specific heats | – |
| $C_p$ | Specific heat at constant pressure | J/(kg·K) |

---

## 1. Starting point — the radial equilibrium equation (15)

The balance of pressure and inertial forces on a fluid element, written in the
radial ($\hat{e}_r$) frame, is

$$
\frac{1}{\rho}\frac{dP}{dr} \;=\; \frac{V_T^{2}}{r}\;-\;\frac{V_M^{2}}{r_M}\cos\phi\;-\;V_r\frac{dV_M}{dM}
\tag{15}
$$

For convenience, name the right-hand side

$$
\mathcal{R}\;\equiv\;\frac{V_T^{2}}{r}\;-\;\frac{V_M^{2}}{r_M}\cos\phi\;-\;V_r\frac{dV_M}{dM},
\qquad\text{so that}\qquad
\frac{dP}{dr}=\rho\,\mathcal{R}.
$$

The three terms on the right are, in order, the **whirl** (centripetal)
contribution, the **meridional curvature** contribution, and the
**streamwise acceleration** contribution. The last term is the one most
analyses drop for axial machines (where $V_r$ is small) but which must be
retained for radial machines, where hub and shroud radii change appreciably.

---

## 2. The closure problem

Equation (15) mixes a pressure gradient ($dP/dr$) with velocities. To integrate
it radially we need the pressure gradient expressed in terms of quantities the
solver tracks: the meridional velocity $V_M$, the total pressure $P_0$, and the
total temperature $T_0$. The next steps remove $P$ in favour of those.

---

## 3. Pressure in terms of total conditions (30) – (31)

From the definition of total temperature for a calorically perfect gas,

$$
T_0 = T + \frac{V^{2}}{2C_p}
\qquad\Longrightarrow\qquad
\frac{T}{T_0}=1-\frac{V^{2}}{2C_p T_0}
\tag{30}
$$

and the isentropic pressure–temperature relation gives

$$
\frac{P}{P_0}=\left(\frac{T}{T_0}\right)^{\frac{\gamma}{\gamma-1}}
=\left(1-\frac{V^{2}}{2C_p T_0}\right)^{\frac{\gamma}{\gamma-1}}.
\tag{31}
$$

Solving for the static pressure,

$$
P=P_0\left(1-\frac{V^{2}}{2C_p T_0}\right)^{\frac{\gamma}{\gamma-1}}.
$$

---

## 4. Differentiate with respect to radius (32)

Taking the radial derivative of the expression for $P$ gives the form of the
static-pressure gradient that appears in (15):

$$
\frac{dP}{dr}=\frac{d}{dr}\!\left[\,P_0\left(1-\frac{V^{2}}{2C_p T_0}\right)^{\frac{\gamma}{\gamma-1}}\right].
\tag{32}
$$

Both $P_0$ and the bracket vary with radius, so this is a product that will need
the product rule — done in step 6.

---

## 5. The $B$ and $C$ substitutions (33)

Using the velocity decomposition (from the main text, with $\alpha$ measured
between $V$ and the meridional velocity),

$$
V^{2}=V_M^{2}\left(1+\tan^{2}\alpha\right),
$$

define the pressure-ratio factor

$$
B=\left(1-\frac{V_M^{2}\left(1+\tan^{2}\alpha\right)}{2C_p T_0}\right)^{\frac{\gamma}{\gamma-1}}.
\tag{33}
$$

To compress the algebra, also define the dimensionless group

$$
C\;\equiv\;\frac{\left(1+\tan^{2}\alpha\right)}{2C_p}\,\frac{V_M^{2}}{T_0}.
$$

With these definitions,

$$
B=(1-C)^{\frac{\gamma}{\gamma-1}},\qquad P=P_0\,B.
$$

> **Physical meaning of $C$.** $C$ is the local kinetic energy scaled by total
> enthalpy — essentially a scaled Mach number squared. As $C\to 1$ the bracket
> $(1-C)\to 0$ and the formulation becomes singular; this is the choking limit
> discussed in the Solution Stability section of the paper.

---

## 6. Product rule (34)

Because $P=P_0 B$, the static-pressure gradient (32) expands as

$$
\frac{dP}{dr}=\frac{dP_0}{dr}B+P_0\frac{dB}{dr}.
\tag{34}
$$

The two remaining unknown derivatives are $dB/dr$ (step 7) and, hidden inside it,
$dC/dr$ (step 8).

---

## 7. Derivative of $B$ (35)

With $B=(1-C)^{\frac{\gamma}{\gamma-1}}$, the chain rule gives

$$
\frac{dB}{dr}
=\frac{\gamma}{\gamma-1}\,(1-C)^{\frac{\gamma}{\gamma-1}-1}\cdot\left(-\frac{dC}{dr}\right).
$$

The exponent reduces cleanly:

$$
\frac{\gamma}{\gamma-1}-1=\frac{\gamma-(\gamma-1)}{\gamma-1}=\frac{1}{\gamma-1},
$$

so

$$
\frac{dB}{dr}=-\,\frac{\gamma}{\gamma-1}\,(1-C)^{\frac{1}{\gamma-1}}\,\frac{dC}{dr}.
\tag{35}
$$

---

## 8. Derivative of $C$ (36)

Write $C=\kappa\,\dfrac{V_M^{2}}{T_0}$ with $\kappa=\dfrac{1+\tan^{2}\alpha}{2C_p}$
treated as locally constant along the streamline. Applying the quotient rule to
$V_M^{2}/T_0$,

$$
\frac{d}{dr}\!\left[\frac{V_M^{2}}{T_0}\right]
=\frac{2V_M\dfrac{dV_M}{dr}\,T_0-V_M^{2}\dfrac{dT_0}{dr}}{T_0^{2}}
=\frac{2V_M}{T_0}\frac{dV_M}{dr}-\frac{V_M^{2}}{T_0^{2}}\frac{dT_0}{dr},
$$

hence

$$
\frac{dC}{dr}=\frac{\left(1+\tan^{2}\alpha\right)}{2C_p}
\left[\,\frac{2V_M}{T_0}\frac{dV_M}{dr}-\frac{V_M^{2}}{T_0^{2}}\frac{dT_0}{dr}\,\right].
\tag{36}
$$

---

## 9. Assembling the full equation (37)

Substitute (35) and (36) into (34), then set $dP/dr=\rho\,\mathcal{R}$ from (15).
This produces the radial equilibrium equation written entirely in terms of
$P_0$, $V_M$, and $T_0$ and their radial derivatives:

$$
\frac{dP_0}{dr}B
\;-\;P_0\,\frac{\gamma}{\gamma-1}\,(1-C)^{\frac{1}{\gamma-1}}
\frac{\left(1+\tan^{2}\alpha\right)}{2C_p}
\left[\frac{2V_M}{T_0}\frac{dV_M}{dr}-\frac{V_M^{2}}{T_0^{2}}\frac{dT_0}{dr}\right]
\;=\;
\rho\left(\frac{V_T^{2}}{r}-\frac{V_M^{2}}{r_M}\cos\phi-V_r\frac{dV_M}{dM}\right).
\tag{37}
$$

This is the implicit form: $dV_M/dr$ appears inside the bracket on the left and
$\mathcal{R}$ still carries a $dV_M/dM$ term on the right.

---

## 10. The cleaned, solved form (38)

Collect the constant-like coefficient into a single symbol:

$$
A\;\equiv\;-\,P_0\,\frac{\gamma}{\gamma-1}\,(1-C)^{\frac{1}{\gamma-1}}\,\frac{\left(1+\tan^{2}\alpha\right)}{2C_p}.
$$

Substituting $A$ into (37) collapses the left side to

$$
\frac{dP_0}{dr}B+A\left[\frac{2V_M}{T_0}\frac{dV_M}{dr}-\frac{V_M^{2}}{T_0^{2}}\frac{dT_0}{dr}\right]=\rho\,\mathcal{R}.
$$

Expanding and isolating the $dV_M/dr$ term,

$$
A\,\frac{2V_M}{T_0}\frac{dV_M}{dr}
=\rho\,\mathcal{R}-B\frac{dP_0}{dr}+A\,\frac{V_M^{2}}{T_0^{2}}\frac{dT_0}{dr},
$$

and finally solving for the meridional-velocity gradient,

$$
\boxed{\;
\frac{dV_M}{dr}=\frac{T_0}{2V_M A}
\left[\,\rho\left(\frac{V_T^{2}}{r}-\frac{V_M^{2}}{r_M}\cos\phi-V_r\frac{dV_M}{dM}\right)
-B\frac{dP_0}{dr}+A\,\frac{V_M^{2}}{T_0^{2}}\frac{dT_0}{dr}\,\right]\;}
\tag{38}
$$

Equivalently, pulling the last term out of the bracket (since
$\frac{T_0}{2V_M A}\cdot A\frac{V_M^{2}}{T_0^{2}}=\frac{V_M}{2T_0}$):

$$
\frac{dV_M}{dr}=\frac{T_0}{2V_M A}\left[\,\rho\,\mathcal{R}-B\frac{dP_0}{dr}\,\right]
+\frac{V_M}{2T_0}\frac{dT_0}{dr}.
$$

> **Note on the final term — and what the code actually does.**
> The implementation in
> [`turbodesign/radeq.py`](https://github.com/nasa/turbo-design/blob/main/turbodesign/radeq.py)
> uses the **pulled-out form** verbatim:
>
> ```python
> dVm_dr = T0/(2*Vm*A) * (rho*eqn15_rhs - B*dP0_dr) + Vm/(2*T0) * dT0_dr
> ```
>
> i.e. the temperature contribution is the separate additive term
> $\dfrac{V_M}{2T_0}\dfrac{dT_0}{dr}$, exactly as derived above. This confirms
> the algebra: solving (37) for $dV_M/dr$ yields
> $+A\,\dfrac{V_M^{2}}{T_0^{2}}\dfrac{dT_0}{dr}$ inside the bracket, equivalently
> $+\dfrac{V_M}{2T_0}\dfrac{dT_0}{dr}$ outside it. The typeset Equation (38) in
> the manuscript writes that term as $\dfrac{V_M^{2}}{T_0}\dfrac{dT_0}{dr}$
> *inside* the bracket, which drops the $A/T_0$ factor and is not consistent
> with the substitution of $A$ into (37). The code and the derivation agree;
> the manuscript's boxed equation should be reconciled to match (a typesetting
> issue in the paper, not an error in the solver).

---

## 11. Closure — $V_r$ and $V_T$ in terms of $V_M$

Equation (38) still contains $V_T$ and $V_r$. Both close directly through the
velocity decomposition, so the only independent unknown velocity is $V_M$:

$$
V_r=V_M\sin\phi,\qquad
V_T=V_M\tan\alpha,\qquad
V^{2}=V_M^{2}\left(1+\tan^{2}\alpha\right).
$$

With these substitutions, (38) is a first-order ordinary differential equation
for $V_M(r)$ along a radial cut. In the implementation it is integrated as one
member of a coupled system together with $P_0$ and $T_0$ — the ODE state vector
is $y=[P_0,\,T_0,\,V_M]$ and the right-hand side returns
$[\,dP_0/dr,\;dT_0/dr,\;dV_M/dr\,]$, with $dP_0/dr$ and $dT_0/dr$ obtained by
differentiating the total-condition relations (Section 3) using the frozen
static profiles $T(r)$ and $P(r)$ for the current row. The geometric quantities
$r_M$ and $\phi$ come from the passage geometry.

---

## 12. How the solver uses (38)

The radial march is implemented in
[`turbodesign/radeq.py`](https://github.com/nasa/turbo-design/blob/main/turbodesign/radeq.py)
in the function `radeq(row, upstream, downstream)`, which wraps the inner
right-hand-side `ode_radeq_streamtube(r, y)`:

1. At a blade-row station the radial profiles of the static state $T(r)$,
   $P(r)$, blade angle $\alpha(r)$, streamline curvature $r_M(r)$, and
   inclination $\phi(r)$ are available on the row.
2. Initial conditions are taken at **mid-span** (50 % hub-to-shroud):
   $y_0=[P_{0,m},\,T_{0,m},\,V_{M,m}]$.
3. The system is integrated outward with `scipy.integrate.solve_ivp` in two
   sweeps from mid-span — once toward the **shroud** and once toward the
   **hub** — and the two solutions are concatenated into a hub-to-shroud
   profile. (Integrating from the middle outward keeps both ends well-posed.)
4. Inside the right-hand side, the meridional coupling $dV_M/dM$ is estimated
   with a shape-preserving (`PchipInterpolator`) fit of $V_M$ versus meridional
   distance $m$ across the upstream, current, and downstream rows, then
   differentiated.
5. If the curvature is negligible ($|r_M|\le\varepsilon$, i.e. an essentially
   axial passage) the solver falls back to the **simple radial equilibrium**
   right-hand side $\mathcal{R}_{\text{simple}} = V_T^{2}/r$, dropping the
   curvature and streamwise-acceleration terms.
6. The streamline positions and the pressure split between blade rows are then
   adjusted in an outer loop until mass flow is balanced — see
   `adjust_streamlines()` in
   [`turbodesign/solve_radeq.py`](https://github.com/nasa/turbo-design/blob/main/turbodesign/solve_radeq.py)
   and the *Solver Modes* section of the main text.

The choking guard is explicit in the code: if $C>1$ the right-hand side raises
an exception rather than returning `NaN`. The singularity sits in the
$(1-C)^{1/(\gamma-1)}$ factor of $A$ and $B$, so $C$ must stay below unity —
this defines the solver's operating envelope.

---

## 13. Mapping the derivation to the code

Every symbol in the derivation appears, with the same meaning, in
`ode_radeq_streamtube` inside
[`turbodesign/radeq.py`](https://github.com/nasa/turbo-design/blob/main/turbodesign/radeq.py):

| Derivation | Code | Notes |
|---|---|---|
| $V_T=V_M\tan\alpha$ | `Vt = Vm*np.tan(alpha)` | closure, Section 11 |
| $V_r=V_M\sin\phi$ | `Vr = Vm*np.sin(phi)` | closure, Section 11 |
| $C=\dfrac{(1+\tan^{2}\alpha)V_M^{2}}{2C_p T_0}$ | `C = (1 + np.tan(alpha)**2) * Vm**2/(2*Cp*T0)` | Section 5 |
| $B=(1-C)^{\gamma/(\gamma-1)}$ | `B = (1-C)**(gamma/(gamma-1))` | Eq (33) |
| $A=-P_0\dfrac{\gamma}{\gamma-1}(1-C)^{1/(\gamma-1)}\dfrac{1+\tan^{2}\alpha}{2C_p}$ | `A = -P0 * gamma/(gamma-1) * (1-C)**(1/(gamma-1)) * (1 + np.tan(alpha)**2)/(2*Cp)` | Section 10 |
| $\mathcal{R}=\dfrac{V_T^{2}}{r}-\dfrac{V_M^{2}}{r_M}\cos\phi-V_r\dfrac{dV_M}{dM}$ | `eqn15_rhs = Vt**2/r - Vm**2/rm*np.cos(phi) - Vr*dVm_dm` | RHS of Eq (15) |
| $\mathcal{R}_{\text{simple}}=\dfrac{V_T^{2}}{r}$ | `eqn15_rhs_simple = Vt**2/r` | axial fallback |
| $\dfrac{dT_0}{dr}=\dfrac{dT}{dr}+\dfrac{V_M(1+\tan^{2}\alpha)}{C_p}\dfrac{dV_M}{dr}$ | `dT0_dr = dT_dr + Vm/Cp * (1 + np.tan(alpha)**2)*dVm_dr` | from $T_0=T+V^{2}/2C_p$ |
| $\dfrac{dP_0}{dr}$ (chain rule on $P_0=P(T_0/T)^{\gamma/(\gamma-1)}$) | `dP0_dr = dP_dr * (T0/T)**(gamma/(gamma-1)) + P*gamma/(gamma-1) * (T0/T)**(1/(gamma-1)) * (T*dT0_dr-T0*dT_dr)/T**2` | inverse of Eq (31) |
| Eq (38), pulled-out form | `dVm_dr = T0/(2*Vm*A) * (rho*eqn15_rhs - B*dP0_dr) + Vm/(2*T0) * dT0_dr` | see *Note on the final term* |

Supporting modules:

- **`turbodesign/radeq.py` → `radeq()`** — top-level radial-equilibrium march
  (Eqs 15, 30–38); inner `ode_radeq_streamtube()` is the ODE right-hand side.
- **`turbodesign/solve_radeq.py` → `adjust_streamlines()`** — outer loop that
  redistributes streamlines to balance mass flow and recomputes $\phi$ and
  $r_M$ via the passage geometry.
- **`turbodesign/passage.py` → `Passage.streamline_curvature()`** — returns
  $\phi$ (inclination) and $r_M$ (meridional radius of curvature) that feed
  Eq (15).
- **`turbodesign/turbine_math.py`** — `compute_gas_constants` ($C_p$, $\gamma$),
  `compute_quantities`, `compute_power`.
- **`turbodesign/flow_math.py` → `compute_massflow`** — radial mass-flow
  integration used by the streamline balancing loop.

> A stale inline comment in the code labels the `dVm_dr` line `# Eqn 21`; this
> refers to an earlier internal numbering and corresponds to Equation (38) in
> the manuscript.

---

*Generated as supplementary documentation for the Turbo-Design repository.
Equation numbers match Appendix A of the manuscript; code references are to
[`nasa/turbo-design`](https://github.com/nasa/turbo-design) `main`.*
