# A centrifugal compressor meanline module

`turbodesign.centrifugal` adds meanline design and off-design analysis for **radial
(centrifugal) compressors** — impeller, vaneless space, vaned diffuser, deswirl.

It is a **new, self-contained module**. It imports nothing from `compressor_math.py`,
`radeq.py` or `flow_math.py`, and it changes **no existing file**. The axial turbine paths this
library exists to serve cannot regress as a result of merging it.

---

## Read this before you read the code

**The model is falsified on one of the three machines it was validated against, and that machine
is NASA's own CC3.** We are telling you this first, not in a footnote, because it is the single
most important thing about this contribution.

Validated on **one frozen coefficient set, with no per-machine retuning** — that constraint is the
whole point, and it is what makes the failure legible:

| machine | measurand | model PR_tt | measured PR_tt | signed error |
|---|---|---|---|---|
| Eckardt Rotor O | impeller-only, measured (NASA TM-75232, Table 6.1) | 2.158352 | 2.2190 | **−2.73 %** |
| NASA HECC | stage, measured (Braunscheidel 2014) | 4.780524 | 4.6847 | **+2.05 %** |
| NASA CC3 | impeller-only, digitized (Kulkarni 2014, Fig. 8) | 4.813221 | 4.1667 | **+15.52 %** |

**The errors do not share a sign. There is no bias to calibrate out.** The residual is structural,
and it separates Eckardt from {HECC, CC3} — but machine class (splittered vs splitter-free), Mach
number, loading, and the `PR ≤ 3.5` validity range of the Aungier correlations are **all confounded**
across exactly those three machines. Breaking the confound requires a fourth machine that is either
splitter-free *and* high-pressure-ratio, or splittered *and* low-pressure-ratio. **We do not hold
one.**

A dedicated experiment (2026-07) **refuted slip factor as the cause**: no admissible slip model
improves on the shipped one, and the best spread achievable across the admissible set is the
control's own.

So: this module reproduces HECC to about two percent on a frozen coefficient set, and it does not
reproduce CC3. **We are not offering a solved problem. We are offering a working meanline
centrifugal capability, with its failure documented and its cause bounded.**

---

## Why it exists

The library's own paper states that it *"currently does not support compressor design — both axial
and radial."* There is no centrifugal solver, and `turbodesign/deviation/` contains **no slip
factor** — slip is the centrifugal analogue of deviation, and without it an impeller's work input
is unrecoverable.

There *is* a centrifugal loss suite in `turbodesign/loss/compressor/otac.py` — sixteen models,
correctly cited. **None of them can currently run.** Every one declares `LossType.Enthalpy`, which
the compressor path does not consume, so each contributes `Yp = 0`: a loss-free impeller, with no
error and no warning. This module does not depend on that path; it carries its own loss set.

---

## What's in it

| module | what it does |
|---|---|
| `solver.py` | the meanline solve — impeller through deswirl |
| `state.py` | thermodynamic state; rothalpy-consistent relative-frame totals |
| `geometry.py` | flowpath, stations at any meridional angle φ |
| `slip.py` | Wiesner, Stanitz, Busemann, Qiu |
| `losses.py` | internal and parasitic losses, kept **separate** |
| `diffusion.py` | vaneless space, vaned diffuser, deswirl |
| `coefficients.py` | the frozen coefficient set — one set, all machines |
| `components.py`, `candidates.py` | component assembly; alternative loss sets |

Every empirical coefficient cites its source (author, year, equation) or is explicitly marked
unverified. Correlations guard their fitted range rather than extrapolating silently.

### Two things it does differently, on purpose

**Internal and parasitic losses are not the same physics.** Internal losses (incidence, friction,
clearance, mixing) destroy total pressure and do not change work. Parasitic losses (disc friction,
recirculation, leakage) **add work** and add no pressure. They enter efficiency on opposite sides
of the fraction. Lumping them into a single `Yp` makes the efficiency rollover unreproducible at
any coefficient value — which is why four of the sixteen `otac.py` models cannot simply be "wired
up" to a pressure-loss coefficient.

**Choke is on the meridional Mach number.** The absolute Mach number routinely exceeds 1 at a
centrifugal impeller exit *without choking*. `if M > 1: choked` is a bug in a radial machine.

---

## Tests

The unit and oracle suites run against an independent meanline calculation and against measured
geometry from NASA sources. They do not require the performance-map data, which is not
redistributed here.

```bash
uv run pytest tests/unit tests/oracle -q
```

---

## Data

The geometry shipped with this module comes from **US Government works**: HECC blade angles, vane
angles and flowpath from Medic et al., NASA/CR-2014-218114; Eckardt rotor geometry from NASA
TM-75232.

**The measured performance maps are deliberately not included.** They were digitized from figures
in AIAA and ASME publications, and their redistribution rights have not been cleared. The
validation figures above can be reproduced from those sources by anyone who holds them.
