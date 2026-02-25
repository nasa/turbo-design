# Entropy-Based Efficiency in Turbomachinery: The $T_2 \Delta s$ Approach

## 1. Introduction

Efficiency is the central performance metric in turbomachinery design. It quantifies how much of the available energy is converted to useful work (turbines) or how much work is required relative to the ideal (compressors). The classical isentropic efficiency definitions work well for axial machines with small radius change, but they break down for radial and mixed-flow machines. This document derives the entropy-based efficiency $\eta = w/(w + T_2 \Delta s)$ from first principles, explains why it is superior, and demonstrates that it applies to both compressors and turbines.

---

## 2. Conventional Isentropic Efficiency

### 2.1 Turbine

For an adiabatic turbine, the total-to-total isentropic efficiency is defined as the ratio of actual work to ideal (isentropic) work:

$$\eta_{tt} = \frac{h_{01} - h_{02}}{h_{01} - h_{02s}}$$

For a calorically perfect gas ($h_0 = C_p T_0$):

$$\eta_{tt} = \frac{T_{01} - T_{02}}{T_{01} - T_{02s}}$$

The isentropic exit temperature $T_{02s}$ comes from the isentropic relation:

$$\frac{T_{02s}}{T_{01}} = \left(\frac{P_{02}}{P_{01}}\right)^{(\gamma - 1)/\gamma}$$

Therefore:

$$T_{02s} = T_{01} \left(\frac{P_{02}}{P_{01}}\right)^{(\gamma-1)/\gamma}$$

The denominator $h_{01}-h_{02s}$ depends entirely on the **absolute total pressure ratio** $P_{02}/P_{01}$.

### 2.2 Compressor

For an adiabatic compressor:

$$\eta_{tt} = \frac{h_{02s} - h_{01}}{h_{02} - h_{01}} = \frac{T_{02s} - T_{01}}{T_{02} - T_{01}}$$

Again, $T_{02s}$ is computed from the absolute total pressure ratio:

$$T_{02s} = T_{01} \left(\frac{P_{02}}{P_{01}}\right)^{(\gamma-1)/\gamma}$$

### 2.3 The Hidden Assumption

Both definitions rely on $P_{02}/P_{01}$ — the **absolute-frame** total pressure ratio. For an axial machine where the blade speed $U$ is approximately the same at inlet and exit ($r_1 \approx r_2$), the absolute total pressure captures the aerodynamic losses directly. But what happens when $r_1 \neq r_2$?

---

## 3. The Problem: Frame Change in Radial Machines

### 3.1 The Euler Turbomachinery Equation

The specific work extracted (turbine) or input (compressor) is:

$$w = U_1 V_{\theta 1} - U_2 V_{\theta 2}$$

where $U = \omega r$ is the blade speed and $V_\theta$ is the absolute tangential velocity. For a radial inflow turbine, $U_1 \gg U_2$ (flow moves from large to small radius), so the work is dominated by the change in $U$, not the aerodynamic turning.

### 3.2 Rothalpy Conservation

In the relative frame, the **rothalpy** is conserved across a rotor (for steady, adiabatic flow):

$$I = h_{0,\text{rel}} - \frac{1}{2}U^2 = \text{const}$$

Expanding:

$$h + \frac{1}{2}W^2 - \frac{1}{2}U^2 = \text{const}$$

where $W$ is the relative velocity. This means:

$$h_{0,\text{rel},1} - \frac{1}{2}U_1^2 = h_{0,\text{rel},2} - \frac{1}{2}U_2^2$$

### 3.3 Absolute vs. Relative Total Pressure

The absolute total pressure is related to the relative total pressure by:

$$P_0 = P_{0R} \cdot f(U, W, \text{geometry})$$

More precisely, for an ideal gas:

$$\frac{P_0}{P} = \left(1 + \frac{\gamma - 1}{2}M^2\right)^{\gamma/(\gamma-1)}$$

$$\frac{P_{0R}}{P} = \left(1 + \frac{\gamma - 1}{2}M_R^2\right)^{\gamma/(\gamma-1)}$$

where $M$ and $M_R$ are the absolute and relative Mach numbers. Here $V$ is the absolute velocity and $W$ is the relative velocity (i.e. the velocity seen in the rotating frame of the blade). The difference between $P_0$ and $P_{0R}$ is driven by the difference between $V^2 = V_x^2 + V_r^2 + V_\theta^2$ and $W^2 = W_x^2 + W_r^2 + W_\theta^2$, which is:

$$V^2 - W^2 = 2UV_\theta - U^2$$

For a radial turbine where $U_1$ is large (say $U_1 = 400$ m/s) and $U_2$ is small ($U_2 = 150$ m/s), the change in $\frac{1}{2}U^2$ across the rotor is:

$$\Delta\left(\frac{1}{2}U^2\right) = \frac{1}{2}(400^2 - 150^2) = 68{,}750 \text{ J/kg}$$

This is a massive contribution to the absolute total enthalpy (and hence $P_0$) change that has **nothing to do with aerodynamic loss**. It is purely a geometric consequence of the radius change.

### 3.4 Why $\eta \to 1$

The isentropic denominator $T_{01} - T_{02s}$ is computed from the absolute total pressure ratio. When the frame-change contribution dominates $P_{02}/P_{01}$, the isentropic temperature drop $T_{01} - T_{02s}$ becomes very close to the actual temperature drop $T_{01} - T_{02}$, regardless of how much aerodynamic loss exists. The result:

$$\eta_{tt} = \frac{T_{01} - T_{02}}{T_{01} - T_{02s}} \approx 1$$

even when the total pressure loss coefficient $Y_p$ is significant. **The conventional efficiency is blind to the loss.**

---

## 4. Entropy: A Frame-Independent Measure of Loss

### 4.1 The Second Law

For an adiabatic process, the second law of thermodynamics states:

$$\Delta s = s_2 - s_1 \geq 0$$

Equality holds only for a reversible (isentropic) process. Any irreversibility — viscous dissipation, shock waves, mixing, tip leakage — produces entropy. The greater the loss, the greater $\Delta s$.

### 4.2 Entropy Change for an Ideal Gas

For a calorically perfect ideal gas, the entropy change between two states is:

$$\Delta s = C_p \ln\frac{T_2}{T_1} - R \ln\frac{P_2}{P_1}$$

This is implemented in the turbo-design compressor solver (`compressor_math.py`, line 92):

```python
entropy_rise_local = 0.5 * (row.Cp + upstream.Cp) * np.log(T_local / upstream.T) \
                   - row.R * np.log(P_local / upstream.P)
```

### 4.3 The Key Insight: Entropy is a State Function

Entropy is a **thermodynamic state function**. It depends only on the thermodynamic state ($T$, $P$) of the fluid, not on the frame of reference. Whether you compute $\Delta s$ in the absolute frame or the relative frame, you get the same answer:

$$\Delta s = C_p \ln\frac{T_2}{T_1} - R \ln\frac{P_2}{P_1} = C_p \ln\frac{T_{2R}}{T_{1R}} - R \ln\frac{P_{2R}}{P_{1R}}$$

But note: **static** temperature and pressure are the same in both frames ($T = T_R$, $P = P_R$). Therefore, we can also write:

$$\Delta s = C_p \ln\frac{T_{02R}}{T_{01R}} - R \ln\frac{P_{0R,2}}{P_{0R,1}}$$

### 4.4 Simplification for Adiabatic Rotors

For an adiabatic rotor, rothalpy is conserved. Using the relation $h_{0,\text{rel}} = C_p T_{0R}$ and rothalpy conservation, the relative total temperature does not change across the rotor in the ideal case. For the entropy calculation, the key simplification is:

$$\Delta s = R \ln\frac{P_{0R,1}}{P_{0R,2}}$$

This is because any decrease in relative total pressure (at approximately constant relative total temperature) is due entirely to **irreversible loss**. The centrifugal effect does not appear in $P_{0R}$ — it has been removed by working in the relative frame.

This is the formula used in `flow_math.py` (line 159):

```python
ds = row.R * np.log(np.mean(ref.P0R) / np.mean(row.P0R))
```

---

## 5. The Gouy-Stodola Theorem

### 5.1 Classical Statement

The Gouy-Stodola theorem (also known as the lost work theorem) states that for any irreversible process, the **lost work** — the difference between the maximum possible work and the actual work — is:

$$w_{\text{lost}} = T_0 \cdot \Delta s_{\text{total}}$$

where $T_0$ is a reference temperature (typically the environment temperature) and $\Delta s_{\text{total}}$ is the total entropy generated.

### 5.2 Derivation

Consider a control volume with steady flow. From the first law (adiabatic, no potential energy):

$$w = h_{01} - h_{02}$$

From the second law, the maximum work (reversible process) would give:

$$w_{\text{max}} = h_{01} - h_{02s}$$

where state $2s$ has the same pressure as state 2 but $s_{2s} = s_1$ (isentropic). The lost work is:

$$w_{\text{lost}} = w_{\text{max}} - w = h_{02} - h_{02s}$$

For an ideal gas, we can relate the enthalpy difference to the entropy difference. Starting from:

$$dh = T \, ds + v \, dP$$

At constant pressure ($P_2 = P_{2s}$, since both states share the same exit pressure):

$$h_{02} - h_{02s} \approx T_2 \cdot (s_2 - s_1) = T_2 \cdot \Delta s$$

This approximation becomes exact for small $\Delta s$ and is very good for typical turbomachinery losses. Therefore:

$$\boxed{w_{\text{lost}} = T_2 \cdot \Delta s}$$

### 5.3 Physical Interpretation

The term $T_2 \cdot \Delta s$ represents the **work that was destroyed by irreversibility**. Every unit of entropy generated at temperature $T_2$ corresponds to $T_2$ units of lost work per unit mass. This is why entropy generation is the fundamental measure of loss in turbomachinery — it directly translates to lost performance.

---

## 6. Deriving the Entropy-Based Efficiency

### 6.1 Definition

The entropy-based efficiency is defined as:

$$\eta = \frac{w_{\text{actual}}}{w_{\text{actual}} + w_{\text{lost}}} = \frac{w}{w + T_2 \cdot \Delta s}$$

For a **turbine** (work output):

$$\boxed{\eta = \frac{C_p(T_{01} - T_{02})}{C_p(T_{01} - T_{02}) + T_2 \cdot \Delta s}}$$

where $\Delta s = R \ln(P_{0R,1}/P_{0R,2})$ for a rotor.

For a **compressor** (work input), the definition is analogous. The actual work input is $w = C_p(T_{02} - T_{01})$, and the ideal work would be $w - w_{\text{lost}}$:

$$\boxed{\eta = \frac{C_p(T_{02} - T_{01}) - T_2 \cdot \Delta s}{C_p(T_{02} - T_{01})} = 1 - \frac{T_2 \cdot \Delta s}{w}}$$

Or equivalently, using the "ideal work over actual work" convention:

$$\eta = \frac{w - T_2 \Delta s}{w}$$

Both forms give $\eta = 1$ when $\Delta s = 0$ (no loss) and $\eta < 1$ when $\Delta s > 0$.

### 6.2 Implementation in turbo-design

The turbine implementation (`flow_math.py`, lines 152–163):

```python
# Entropy-based total-total efficiency:  η = w / (w + T_exit·Δs)
# The standard isentropic formula η = ΔT0/(T01−T0_is) uses the
# absolute P0 ratio which, for radial machines with large radius
# change, is dominated by the frame change and barely reflects the
# relative-frame loss — giving η ≈ 1 even with significant Yp.
# The entropy-based definition always isolates the irreversibility.
if np.mean(ref.P0R) > 0 and np.mean(row.P0R) > 0 and deltaT > 0:
    ds = row.R * np.log(np.mean(ref.P0R) / np.mean(row.P0R))
    w_per_mass = row.Cp * deltaT
    row.eta_total = w_per_mass / (w_per_mass + row.T.mean() * max(ds, 0.0))
else:
    row.eta_total = (ref.T0.mean() - row.T0.mean()) / \
                    max(ref.T0.mean() - row.T0_is.mean(), 1e-9)
```

Note: the fallback (line 163) uses the conventional isentropic formula when relative-frame data ($P_{0R}$) is unavailable — for example, across a stator where there is no relative frame.

### 6.3 Limiting Cases

**Zero loss** ($\Delta s = 0$):

$$\eta = \frac{w}{w + 0} = 1$$

Both definitions agree.

**Small radius change** ($U_1 \approx U_2$, axial machine): The absolute total pressure ratio $P_{02}/P_{01}$ is no longer dominated by the frame change. In this limit, the conventional isentropic efficiency and the entropy-based efficiency converge to the same value.

**Large radius change** (radial machine): The entropy-based efficiency correctly captures the loss through $P_{0R}$, while the conventional efficiency is misleadingly close to 1.

---

## 7. Applicability: Compressors, Turbines, and Both

### 7.1 Turbines

The entropy-based approach is critical for:
- **Radial inflow turbines** (automotive turbochargers, small gas turbines) where $r_1/r_2$ can exceed 2:1
- **Mixed-flow turbines** where there is moderate radius change
- **Axial turbines** — works correctly, gives the same answer as conventional when $\Delta r$ is small

### 7.2 Compressors

The entropy rise in compressors is computed from static conditions (`compressor_math.py`, line 92):

$$\Delta s = \bar{C}_p \ln\frac{T_2}{T_1} - R \ln\frac{P_2}{P_1}$$

where $\bar{C}_p = (C_{p,1} + C_{p,2})/2$ accounts for temperature-dependent specific heats. This is used for:
- **Centrifugal compressors** (large radius change from inducer to diffuser)
- **Mixed-flow compressors**
- **Axial compressors** (works correctly, entropy-based loss models like Koch-Smith use $\Delta s / C_p$ directly)

### 7.3 Why It Works for Both

The $T_2 \Delta s$ approach works universally because:

1. **Entropy is frame-independent** — no correction needed for rotating vs. stationary frames
2. **Entropy is additive** — losses from different sources (profile, endwall, tip clearance, shock) all contribute to $\Delta s$, which can be summed
3. **Entropy maps directly to lost work** — via the Gouy-Stodola theorem, every Joule of $T \Delta s$ is a Joule of lost work
4. **No dependence on pressure ratio convention** — unlike $Y_p$ which can be defined with inlet or exit dynamic head, $\Delta s$ has only one definition

---

## 8. Summary

| Property | Conventional $\eta_{tt}$ | Entropy-Based $\eta = w/(w + T_2\Delta s)$ |
|---|---|---|
| Depends on | Absolute $P_{02}/P_{01}$ | Relative $P_{0R}$ or static $T$, $P$ |
| Frame-independent | No | Yes |
| Radial machines | Gives $\eta \approx 1$ (misleading) | Gives realistic values |
| Physical meaning | Temperature ratio | Lost work (Gouy-Stodola) |
| Applies to | Axial machines primarily | All machine types |
| Zero-loss limit | $\eta = 1$ | $\eta = 1$ |
| Axial-machine limit | Standard value | Same as conventional |

The entropy-based efficiency definition is the thermodynamically rigorous way to quantify turbomachinery performance. It reduces to the conventional definition for axial machines (where both work fine) and provides correct results for radial and mixed-flow machines (where the conventional definition fails). This is why turbo-design uses $T_2 \Delta s$ as the primary efficiency metric for turbine rotors.
