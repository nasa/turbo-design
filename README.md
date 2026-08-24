# Turbo Design 3
This tool is a streamline turbomachinery design tool solving the radial equilibrium equations. It can be used for designing compressors and turbines. The designs can have counter rotating stages, different working fluids, and cooling. The intent of this tool is to enable added flexibility in which loss models are used. Because it's a python, it can connect with custom machine learning based loss models.

## Table of Contents
- [Key Features](#key-features)
- [Code Structure](#code-structure)
  - [Core Solver Classes](#core-solver-classes)
  - [Mathematics & Solvers](#mathematics--solvers)
  - [Loss Model Architecture](#loss-model-architecture)
- [How Loss Models Work](#how-loss-models-work)
  - [Loss Model Interface](#loss-model-interface)
  - [Available Loss Types](#available-loss-types)
  - [Using Loss Models](#using-loss-models)
  - [Built-in Turbine Loss Models](#built-in-turbine-loss-models)
  - [Built-in Compressor Loss Models](#built-in-compressor-loss-models)
  - [Loss Model Data Files](#loss-model-data-files)
  - [Custom Loss Models](#custom-loss-models)
- [Solver Modes](#solver-modes)
  - [Pressure Balance Mode (Default)](#1-pressure-balance-mode-default)
  - [Angle Matching Mode](#2-angle-matching-mode)
- [Typical Workflow](#typical-workflow)
- [Tutorials](#tutorials)
- [Technical Notes](#technical-notes)
- [Understanding Velocity Triangles](#understanding-velocity-triangles)
- [Documentation](#documentation)
  - [Online Documentation](#online-documentation)
  - [Building Documentation Locally](#building-documentation-locally)
  - [CI/CD for Documentation](#cicd-for-documentation)
- [Contributors](#contributors)
- [License](#license)

## Key Features
- Streamline-based radial equilibrium solver for axial and radial turbomachinery
- Pluggable loss model architecture supporting multiple empirical correlations
- Two solver modes: **Pressure Balance** (fixed blade angles) and **Angle Matching** (target massflow)
- Support for cooling flows, counter-rotating stages, and custom working fluids via Cantera
- Integration with geometry generation tools (AGF format) for 3D CFD pre-processing

## Code Structure

The turbodesign package is organized into the following key components:

### Core Solver Classes
- **[TurbineSpool](turbodesign/turbine_spool.py)** - Main turbine solver coordinating streamline calculations
- **[CompressorSpool](turbodesign/compressor_spool.py)** - Main compressor solver with similar architecture
- **[BladeRow](turbodesign/bladerow.py)** - Dataclass representing a single blade row with 150+ flow properties
- **[Passage](turbodesign/passage.py)** - Meridional passage geometry (hub/shroud curves)
- **[Inlet](turbodesign/inlet.py)** / **[Outlet](turbodesign/outlet.py)** - Boundary condition definitions

### Mathematics & Solvers
- **[turbine_math.py](turbodesign/turbine_math.py)** - Turbine flow calculations (stator_calc, rotor_calc)
- **[compressor_math.py](turbodesign/compressor_math.py)** - Compressor flow calculations
- **[solve_radeq.py](turbodesign/solve_radeq.py)** - Radial equilibrium equation solver
- **[flow_math.py](turbodesign/flow_math.py)** - Massflow, area, and power calculations

### Loss Model Architecture
- **[loss/losstype.py](turbodesign/loss/losstype.py)** - Abstract base class `LossBaseClass` defining the loss model interface
- **[loss/turbine/](turbodesign/loss/turbine/)** - Turbine loss correlations (TD2, Ainley-Mathieson, Kacker-Okapuu, Craig-Cox, Traupel)
- **[loss/compressor/](turbodesign/loss/compressor/)** - Compressor loss correlations (Lieblein, OTAC, Diffusion Factor)
- **[loss/fixedpolytropic.py](turbodesign/loss/fixedpolytropic.py)** - Fixed polytropic efficiency model
- **[loss/fixedpressureloss.py](turbodesign/loss/fixedpressureloss.py)** - Fixed pressure loss model

## How Loss Models Work

Turbo Design 3 uses a **pluggable loss model architecture** that allows users to easily swap between different empirical correlations or implement custom loss functions (including ML-based models).

### Loss Model Interface

All loss models inherit from `LossBaseClass` and implement a simple interface:

```python
from turbodesign.loss.losstype import LossBaseClass
from turbodesign.enums import LossType

class MyCustomLoss(LossBaseClass):
    def __init__(self):
        super().__init__(LossType.Pressure)  # or Enthalpy, Entropy, Polytropic

    def __call__(self, row: BladeRow, upstream: BladeRow) -> np.ndarray:
        """Calculate loss coefficient for each streamline.

        Args:
            row: Current blade row being evaluated
            upstream: Upstream blade row providing inlet conditions

        Returns:
            Array of loss coefficients matching row.r shape
        """
        # Your loss calculation here
        return loss_array
```

### Available Loss Types

The `LossType` enum defines four types of loss representations:

1. **Pressure Loss** (`LossType.Pressure`) - Total pressure loss coefficient (Yp)
2. **Enthalpy Loss** (`LossType.Enthalpy`) - Enthalpy loss coefficient
3. **Entropy Loss** (`LossType.Entropy`) - Entropy generation
4. **Polytropic Efficiency** (`LossType.Polytropic`) - Polytropic efficiency (0-1)

### Using Loss Models

Loss models are assigned to individual blade rows:

```python
from turbodesign.loss.turbine import TD2, AinleyMathieson
from turbodesign import make_rotor_row, make_stator_row

# Create blade rows
rotor = make_rotor_row(row_type=RowType.Rotor, hub_location=0.5, stage_id=1)
stator = make_stator_row(row_type=RowType.Stator, hub_location=0.75, stage_id=1)

# Assign loss models
rotor.loss_function = TD2()
stator.loss_function = AinleyMathieson()
```

### Built-in Turbine Loss Models

| Model | Description | Best For | Reference |
|-------|-------------|----------|-----------|
| **TD2** | NASA legacy correlation from TD2 code | Initial estimates, historical validation | NASA SP-290 |
| **Ainley-Mathieson** | Classic cascade-based correlation | Steam turbines, impulse turbines | ARC R&M 2974 (1951) |
| **Kacker-Okapuu** | Modern update to Ainley-Mathieson | Gas turbines, subsonic to transonic | ASME 81-GT-120 (1982) |
| **Craig-Cox** | Alternative correlation | Steam turbines | IMechE (1971) |
| **Traupel** | European-based correlation | Steam turbines | Thermische Turbomaschinen (1977) |
| **FixedPolytropicEfficiency** | Constant efficiency assumption | Preliminary design, sensitivity studies | - |

### Built-in Compressor Loss Models

| Model | Description | Best For | Reference |
|-------|-------------|----------|-----------|
| **Lieblein** | Diffusion factor based loss | Axial compressors | NACA RM E57A28 (1957) |
| **Diffusion Factor** | Simplified diffusion loss | Preliminary design | - |

### Loss Model Data Files

Some correlations (Ainley-Mathieson, Kacker-Okapuu, etc.) rely on digitized charts stored as `.pkl` files:
- Auto-downloaded from GitHub on first use to `~/.cache/TD3_LossModels/`
- Can be regenerated locally by running `python build_dataset.py` in `references/Turbines/<ModelName>/`
- **Note:** Depending on your Python version, pickle binaries may have compatibility issues. If auto-download fails, regenerate locally.

**Building Loss Model Correlations:**

The turbine loss correlations were estimated using axial steam turbine data. Correlation figures are extracted and surface fitted:

- [Ainley Mathieson](https://colab.research.google.com/github/nasa/turbo-design/blob/main/references/Turbines/AinleyMathieson/ainley_mathieson.ipynb)
- [Craig Cox](https://colab.research.google.com/github/nasa/turbo-design/blob/main/references/Turbines/CraigCox/craig_cox.ipynb)
- [Traupel](https://colab.research.google.com/github/nasa/turbo-design/blob/main/references/Turbines/Traupel/traupel.ipynb)
- [KackerOkapuu](https://colab.research.google.com/github/nasa/turbo-design/blob/main/references/Turbines/KackerOkapuu/kacker_okapuu.ipynb)

### Custom Loss Models

You can implement custom loss models (including ML-based) by subclassing `LossBaseClass`:

```python
import joblib
from turbodesign.loss.losstype import LossBaseClass
from turbodesign.enums import LossType

class MLLossModel(LossBaseClass):
    def __init__(self, model_path: str):
        super().__init__(LossType.Pressure)
        self.model = joblib.load(model_path)  # Load your trained ML model

    def __call__(self, row, upstream):
        # Extract features from row and upstream
        features = self._extract_features(row, upstream)
        # Predict loss using ML model
        loss = self.model.predict(features)
        return loss
```

## Solver Modes

Turbo Design 3 supports two primary solving modes, automatically selected based on outlet boundary conditions:

### 1. Pressure Balance Mode (Default)

**When used:** Outlet initialized with `init_static(P, percent_radii)` (no massflow specified)

**How it works:**
- Blade exit angles (beta2) are **fixed** based on geometry or user input
- Solver adjusts static pressures along streamlines to satisfy radial equilibrium
- Massflow is calculated as a result of the pressure distribution

**Best for:**
- Design problems where blade geometry is already defined
- Matching CFD results with known blade angles
- Sensitivity studies varying inlet conditions

```python
outlet = Outlet(num_streamlines=5)
outlet.init_static(P=100000, percent_radii=[0.5])  # Pressure balance mode
spool.solve()  # Blade angles fixed, pressures adjusted
```

### 2. Angle Matching Mode

**When used:** Outlet initialized with `init_static(P, percent_radii, massflow=target)`

**How it works:**
- Target massflow is **specified** by user
- Solver iteratively adjusts blade exit angles (beta2) to match target massflow
- Radial equilibrium maintained while tweaking angles
- Useful for inverse design problems

**Best for:**
- Meeting specific massflow requirements
- Engine cycle matching (where massflow is constrained)
- Preliminary design before geometry is finalized

```python
outlet = Outlet(num_streamlines=5)
outlet.init_static(P=100000, percent_radii=[0.5], massflow=29.4)  # Angle matching mode
spool.solve()  # Angles adjusted to match massflow=29.4 kg/s
```

## Typical Workflow

```python
from turbodesign import TurbineSpool, Inlet, Outlet, Passage, PassageType
from turbodesign.loss.turbine import TD2
from turbodesign import make_rotor_row, make_stator_row
from cantera import Solution

# 1. Define working fluid
fluid = Solution('air.yaml')

# 2. Define meridional passage geometry
passage = Passage(hub_x, hub_r, shroud_x, shroud_r, passageType=PassageType.Axial)

# 3. Set up inlet boundary conditions
inlet = Inlet(hub_location=0, shroud_location=0, beta=[0])
inlet.init_total(P0=200000, T0=500, M=0.5)

# 4. Set up outlet boundary conditions
outlet = Outlet(num_streamlines=5)
outlet.init_static(P=100000, percent_radii=[0.5])  # Pressure balance mode

# 5. Create blade rows with loss models
rotor1 = make_rotor_row(row_type=RowType.Rotor, hub_location=0.5, stage_id=1)
rotor1.loss_function = TD2()
rotor1.num_blades = 36

stator1 = make_stator_row(row_type=RowType.Stator, hub_location=0.75, stage_id=1)
stator1.loss_function = TD2()
stator1.num_blades = 46

# 6. Create spool and solve
spool = TurbineSpool(
    passage=passage,
    massflow=20,
    inlet=inlet,
    outlet=outlet,
    rows=[rotor1, stator1],
    rpm=10000,
    num_streamlines=5,
    fluid=fluid
)

spool.solve()

# 7. Post-process results
spool.plot()
spool.plot_velocity_triangles()
spool.export_properties('results.json')
print(f"Total power: {spool.total_power()} W")
```

## Tutorials

[Turbine Design: EEE High-Pressure Turbine](https://colab.research.google.com/github/nasa/turbo-design/blob/main/examples/EEE-HPT/eee_hpt.ipynb) - Based on GE's 2-stage HPT design [[3]](#ref-eee-hpt). Full CFD results available on https://data.nasa.gov/dataset/eee-2-stage-hpt-cfd-tecplot-results.


[Compressor Design: EEE High-Pressure Compressor](https://colab.research.google.com/github/nasa/turbo-design/blob/main/examples/EEE-HPC/eee_hpc.ipynb)
> **Note on EEE-HPC Geometry and Results:** The massflow and performance predictions may not exactly match the original 1982 publications [[1]](#ref-eee-hpc). The EEE program's original design work (1970s-1980s) used computational tools and geometry definitions that differ from modern reconstructions. The geometry used in this example is based on publicly released data [[2]](#ref-eee-geom), which was reconstructed from contractor reports and may contain discrepancies in blade radii, twist distributions, and exact airfoil coordinates compared to the test hardware described in the original GE publications. Without access to the original CAD files or manufacturing drawings, exact matching is not possible.
>
> **References:**
> - <a name="ref-eee-hpc"></a>[1] Holloway, P. R., et al. (1984). [*Energy Efficient Engine: High Pressure Compressor Detail Design Report*](https://ntrs.nasa.gov/citations/19850002690). NASA CR-165558, General Electric Company.
> - <a name="ref-eee-geom"></a>[2] Claus, R. W., Beach, T., Turner, M., Siddappaji, K., & Hendricks, E. S. (2015). [*Geometry and Simulation Results for a Gas Turbine Representative of the Energy Efficient Engine (EEE)*](https://ntrs.nasa.gov/citations/20150003286). NASA/TM-2015-218408.
> - <a name="ref-eee-hpt"></a>[3] Halila, E. E., Lenahan, D. T., & Thomas, T. T. (1982). [*Energy Efficient Engine High Pressure Turbine Test Hardware Detailed Design Report*](https://ntrs.nasa.gov/citations/19850002687). NASA CR-167955, General Electric Company.


[Simple 1-stage compressor validation](https://colab.research.google.com/github/nasa/turbo-design/blob/main/examples/3RowSteady-1D/3RowSteady.ipynb)

[Radial inflow turbine design](https://colab.research.google.com/github/nasa/turbo-design/blob/main/examples/radial-turbine/radial_turbine-1D.ipynb)

[Turbine optimization with scipy](https://colab.research.google.com/github/nasa/turbo-design/blob/main/examples/optturb-turbine/optturb.ipynb)

[Multi-stage turbine optimization](https://colab.research.google.com/github/nasa/turbo-design/blob/main/examples/optturb-multistage/optturb-multistage.ipynb)

[Sizing and Matching a Compressor and Turbine with the Mass Flow Function](https://colab.research.google.com/github/nasa/turbo-design/blob/main/examples/shaft-match/compressor_turbine_match.ipynb) — Sizes a compressor inlet and a choked NGV throat directly from the non-dimensional mass flow function, matches a turbine to the compressor's shaft power and mass-flow continuity with `ShaftMatch`, and sweeps the result into a performance map with `sweep_operating_points`.

## Technical Notes

[The Non-Dimensional Mass Flow Function: Sizing, Matching, and Sweeping](mass_flow_function.md) — Derives $\tilde m(M,\gamma)$ and the choke-margin feasibility guard, and covers the three places it's used: sizing a component's inlet/NGV throat before any blade row is solved, matching a compressor and turbine on a shaft with `ShaftMatch`, and sweeping a fixed geometry into a performance map with `sweep_operating_points` — including why a turbine's massflow can't be swept the same way a compressor's can.

[Entropy-Based Efficiency: The T2·delta-s Approach](entropy_based_efficiency.md) — Derives why the conventional isentropic efficiency $\eta = \Delta T_0 / (T_{01} - T_{02s})$ breaks down for radial machines with large radius change, and how the entropy-based definition $\eta = w/(w + T_2 \Delta s)$ resolves this by isolating irreversibility from frame-change effects. Includes full derivations from first and second law, the Gouy-Stodola lost-work theorem, and shows applicability to both compressors and turbines.

## Understanding Velocity Triangles

Turbodesign tracks all flow properties leaving each blade row. Work is computed using `Work = U*(Vt1-Vt2)` [J/kg]; `Power = massflow * Work` [W].

For turbines, maximize tangential velocity exiting the stator and minimize it leaving the rotor to extract maximum work. The diagram below shows how velocity data is stored in each `BladeRow` object:

<img src="references/turbine_velocity_triangles.jpg" alt="Velocity Triangle for a Turbine" style="width:400px;"/>

# Documentation

## Online Documentation

The latest HTML documentation is automatically built and deployed to GitHub Pages on every push to the main branch.

**Documentation URL:** [https://nasa.github.io/turbo-design/](https://nasa.github.io/turbo-design/)

## Building Documentation Locally

To build the documentation locally:

```bash
cd docs
make html
```

Then open `docs/build/html/index.html` in your browser.

**Requirements:**
- Sphinx
- sphinx-rtd-theme
- Python dependencies (numpy, scipy, matplotlib, pandas)

## CI/CD for Documentation

Documentation builds are automated using GitHub Actions. The workflow:
- Builds on every push and pull request
- Deploys to GitHub Pages on push to main/master branch
- Uses Python 3.12 and Sphinx

See [.github/workflows/README.md](.github/workflows/README.md) for setup details.

# Contributors

## Fortran Verson
| Person | Contribution/Role | Dates |
| ------ | ------ | ------ |
| Simon Chen | AXOD | - 2020 |
| Arthur Glassman | TD2 | unknown |
| Paht Juangphanich | Maintainer | 2020-2022 |


## Python Version TD3
| Person | Contribution/Role | Dates | Email |
| ------ | ----------------- | ----- | ------|
| Paht Juangphanich | Turbo Design 3 | - | paht.juangphanich@nasa.gov |
| Andress William | TD2 | Summer 2021 | bill.andress@gmail.com |

# License
[NASA Open Source Agreement](https://opensource.org/licenses/NASA-1.3)

