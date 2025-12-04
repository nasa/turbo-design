"""AGF generation and parsing utilities."""

from dataclasses import dataclass, asdict
from typing import List

import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as plt


@dataclass
class Inlet_bcs:
    ptin: float
    ttin: float
    pspan: float
    machin: float
    alpin: float
    phiin: float


@dataclass
class Outlet_bcs:
    rpm: float
    gamma: float
    psout: float
    twall: float
    molwt: float


@dataclass
class Domain:
    xhup: float  # xhub upstream
    rhup: float  # rhub upstream
    xtup: float  # xtip upstream
    rtup: float  # rtip upstream

    xhdw: float  # xhub downwind
    rhdw: float  # rhub downwind
    xtdw: float  # xtip downwind
    rtdw: float  # rtip downwind


@dataclass
class Settings:
    nprof: int = 1  # number of spanwise points that define pitched average inlet profile
    ifang: int = 10  # 0 - adiabatic wall, 1 - temperature wall
    hbl: float = 1  # Inlet hub boundary layer thickness as percent span
    tbl: float = 1  # Inlet tip boundary layer thickness as percent span

    nblades: int = 8  # number of blades
    npts: int = 300  # number of points per blade
    nspans: int = 3  # number of spans/sections
    ity: int = 5  # What format are the blades in.
    # ITY 5 = x rth r
    # ITY 7 = x1,theta1,r1,x2,theta2,r2
    # ITY 10 = x1,y1,z1,x2,y2,z2
    iym: int = 0  # 0 = do not flip airfoil along x-axis
    tcls: int = 1  # 1 = has tip clearance
    hcls: int = 0  # 0 = no hub clearance
    lete: int = 10  # do not modify leading edge or trailing edge
    isplit: int = 0  # no splitters

    nht: int = 1  # Number of axial points defining the endwall


@dataclass
class Clearance:
    tlecl: float  # tip leading edge clearance
    tmccl: float  # tip mid clearance
    ttecl: float  # tip te clearance
    hlecl: float = 0  # hub le clearance
    hmccl: float = 0  # hub mid clearance
    htecl: float = 0  # hub te clearance


class AGF_Setup:
    def __init__(self, template_file: str = "template.agf", name="radial-turbine"):
        self.agf_template = template_file
        self.name = name
        self.endwall: str = ""
        self.sections: str = ""
        self.clearance: Clearance
        self.settings: Settings
        self.domain: Domain
        self.inlet: Inlet_bcs
        self.outlet: Outlet_bcs

    def add_passage(self, hub: npt.NDArray, shroud: npt.NDArray):
        """Add passage hub/shroud geometry."""
        domain = Domain(
            xhup=hub[0, 0],
            rhup=hub[0, 1],
            xtup=shroud[0, 0],
            rtup=shroud[0, 1],
            xhdw=hub[-1, 0],
            rhdw=hub[-1, 1],
            xtdw=shroud[-1, 0],
            rtdw=shroud[-1, 1],
        )
        endwall = []
        for i in range(hub.shape[0]):
            xl = hub[i, 0]
            rl = hub[i, 1]
            xu = shroud[i, 0]
            ru = shroud[i, 1]
            if i < hub.shape[0] - 1:
                endwall.append(f"{xl:.4f}   {rl:.4f}   {xu:.4f}   {ru:.4f}\n")
            else:
                endwall.append(f"{xl:.4f}   {rl:.4f}   {xu:.4f}   {ru:.4f}")
        self.endwall = "".join(endwall)
        self.domain = domain
        self.settings.nht = hub.shape[0]

    def add_blade(self, ss: npt.NDArray, ps: npt.NDArray, IsDuct: bool = False):
        """Add the blade geometry."""
        sections = []
        if IsDuct:
            self.settings.nspans = 0
        else:
            self.settings.nspans = ss.shape[0]
        section_indx = 1

        for i in range(ss.shape[0]):
            x = np.hstack([ss[i, :, 0], ps[i, 1:-1, 0]])
            y = np.hstack([ss[i, :, 1], ps[i, 1:-1, 1]])
            z = np.hstack([ss[i, :, 2], ps[i, 1:-1, 2]])
            n = len(x)  # number of points
            self.settings.npts = n
            r = np.sqrt(y**2 + z**2)
            th = np.arctan2(y, z)
            rth = r * th
            sections.append(f"*SECTION\t{section_indx}\n")
            sections.append(f"- SECTION - {section_indx}\t{n}\n")
            sections.append(">----RAD------XOFF------YOFF------ROTD----CONEANGLE----\n")
            sections.append("0.0000\t0.0000\t0.0000\t0.0000\t0.0000\n")
            sections.append("x       rth        r\n")
            for j in range(len(x)):
                sections.append(f"{x[j]:0.6f}    {rth[j]:0.6f}    {r[j]:0.6f}\n")
            section_indx += 1

        self.sections = "".join(sections)

    def add_clearance(self, clearance: Clearance):
        self.clearance = clearance

    def add_settings(self, settings: Settings):
        self.settings = settings

    def add_inlet(self, inlet: Inlet_bcs):
        self.inlet = inlet

    def add_outlet(self, outlet: Outlet_bcs):
        self.outlet = outlet

    def build(self, output_filename: str = "stator.agf"):
        with open(self.agf_template, "r") as file:
            file_content = file.read()
            file_content = file_content.replace("[name]", f"{self.name}")

            domain_dict = asdict(self.domain)
            clearance_dict = asdict(self.clearance)
            settings_dict = asdict(self.settings)
            inlet_dict = asdict(self.inlet)
            outlet_dict = asdict(self.outlet)

            with open(output_filename, "w") as f:
                for k, v in domain_dict.items():
                    file_content = file_content.replace(f"[{k}]", f"{v:0.4f}")

                for k, v in clearance_dict.items():
                    file_content = file_content.replace(f"[{k}]", f"{v}")

                for k, v in settings_dict.items():
                    file_content = file_content.replace(f"[{k}]", f"{v}")

                for k, v in inlet_dict.items():
                    file_content = file_content.replace(f"[{k}]", f"{v:0.4f}")

                for k, v in outlet_dict.items():
                    file_content = file_content.replace(f"[{k}]", f"{v:0.4f}")

                if self.sections:
                    file_content = file_content.replace("[sections]", self.sections)

                file_content = file_content.replace("[endwall]", self.endwall)
                f.write(file_content)


def read_agf(file_path: str) -> dict:
    """
    Read an AGF file and return its contents in structured form.

    Returns a dict with keys: 'header', 'sections', 'endwall', 'raw'.
    """
    with open(file_path, "r") as f:
        raw = f.read()

    data = {"header": {}, "sections": [], "endwall": [], "raw": raw}
    lines = raw.splitlines()
    section = None
    for line in lines:
        if line.startswith("*SECTION"):
            if section:
                data["sections"].append(section)
            section = {"header": line, "points": []}
        elif section and line.strip() and line[0].isdigit():
            try:
                parts = [float(p) for p in line.split()[:3]]
                section["points"].append(parts)
            except Exception:
                continue
        elif section is None and "=" in line:
            parts = line.split("=")
            if len(parts) == 2:
                key = parts[0].strip("[] \t")
                val = parts[1].strip()
                data["header"][key] = val

    if section:
        data["sections"].append(section)

    if "[endwall]" in raw:
        after = raw.split("[endwall]")[1]
        for ln in after.splitlines():
            parts = [p for p in ln.split() if p]
            if len(parts) == 4:
                try:
                    data["endwall"].append([float(x) for x in parts])
                except Exception:
                    pass

    return data


def plot_airfoil_inputs(nsections: int, npts: int):
    xthr = np.zeros(shape=(nsections, npts, 3))  # section_num x theta r
    with open("AIRFOIL.INPUTS", "r") as f:
        [f.readline() for _ in range(4)]  # skip first 4 lines
        for i in range(nsections):
            for j in range(npts):
                line = f.readline()
                temp = [float(p) for p in line.split(" ") if p]
                xthr[i, j, 0] = temp[0]


def plot_airfoil_inputs_2D(nsections: int, npts: int):
    xthr = np.zeros(shape=(nsections, npts, 3))
    with open("AIRFOIL.INPUTS", "r") as f:
        [f.readline() for _ in range(4)]
        for i in range(nsections):
            for j in range(npts):
                line = f.readline()
                temp = [float(p) for p in line.split(" ") if p]
                xthr[i, j, 0] = temp[0]
                xthr[i, j, 1] = temp[1]
                xthr[i, j, 2] = temp[2]
            [f.readline() for _ in range(2)]
    for i in range(nsections):
        plt.plot(xthr[i, :, 0], xthr[i, :, 1], "r", label=f"section {i}")
    plt.axis("equal")
    plt.xlabel("x-axial")
    plt.ylabel("y")
    plt.show()
