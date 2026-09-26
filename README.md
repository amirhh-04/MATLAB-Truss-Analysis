<div align="center">

# 🏗️ MATLAB Truss Analysis — 2D & 3D Finite Element Method (FEM)

<p>
<img src="https://img.shields.io/badge/License-MIT-green.svg" alt="License">
<img src="https://img.shields.io/badge/Software-MATLAB-orange?logo=mathworks" alt="MATLAB">
<img src="https://img.shields.io/badge/Method-Direct%20Stiffness-blue" alt="Method">
<img src="https://img.shields.io/badge/Analysis-Linear%20Static-success" alt="Type">
<img src="https://img.shields.io/badge/Toolboxes-None%20required-informational" alt="Toolboxes">
</p>

**Linear static analysis of planar (2D) and space (3D) truss structures using the Finite Element Method and the Direct Stiffness Method — implemented from scratch in pure MATLAB.**

[Features](#-features) · [Theory](#-analysis-method) · [Inputs](#-input-definition) · [Results](#-results) · [Getting Started](#-getting-started) · [Assumptions](#-assumptions--scope)

</div>

---

## 📌 Overview

This repository provides two standalone MATLAB scripts for **truss structural analysis**:

| Solver | Structure | DOFs per Node | Element Type |
|:---:|:---:|:---:|:---:|
| `truss2D.m` | Planar truss | 2 (`UX`, `UY`) | 2-node axial bar |
| `truss3D.m` | Space truss | 3 (`UX`, `UY`, `UZ`) | 2-node axial bar |

Both solvers implement the complete **linear static finite element workflow** — from model definition (nodes, elements, loads, boundary conditions) through **global stiffness matrix assembly**, system solution, and post-processing of **axial forces, stresses, strains, and yield-stress utilization** — without relying on any external structural-analysis solver or toolbox.

The project is intended for **educational, computational mechanics, structural analysis, and engineering programming** applications.

---

## ✨ Features

| Feature | 2D | 3D |
|---|:---:|:---:|
| Direct stiffness method | ✅ | ✅ |
| Element stiffness matrix formulation | ✅ | ✅ |
| Global stiffness matrix assembly | ✅ | ✅ |
| Nodal load definition | ✅ (magnitude + angle) | ✅ (components `Fx Fy Fz`) |
| Flexible boundary-condition definition | ✅ | ✅ (per-DOF flags) |
| Static nodal displacement solution | ✅ | ✅ |
| Element elongation (µm) | ✅ | ✅ |
| Axial force (N) | ✅ | ✅ |
| Normal stress (MPa) | ✅ | ✅ |
| Strain (microstrain, µε) | ✅ | ✅ |
| Yield-stress utilization (%) | ✅ | ✅ |
| Members > 100% utilization highlighted in red | ✅ | ✅ |
| Undeformed + deformed geometry plot | ✅ | ✅ |
| Node & element labeling | ✅ | ✅ |
| Command-window result reporting | ✅ | ✅ |

---

## 📐 Analysis Method

The solver follows the standard linear static finite element procedure for truss structures.

The global governing system is:

$$
[K]\{u\} = \{F\}
$$

where:

| Symbol | Meaning |
|:---:|---|
| `[K]` | Global stiffness matrix |
| `{u}` | Global nodal displacement vector |
| `{F}` | Global nodal force vector |

### 2D Truss Element

With `c = cos(θ)` and `s = sin(θ)` derived from the element orientation:

$$
k_e =
\frac{AE}{L}
\begin{bmatrix}
c^2 & cs & -c^2 & -cs \\
cs & s^2 & -cs & -s^2 \\
-c^2 & -cs & c^2 & cs \\
-cs & -s^2 & cs & s^2
\end{bmatrix}
$$

### 3D Truss Element

The 3D formulation uses the **direction cosines** of the element axis:

$$
c_x = \frac{\Delta x}{L},
\qquad
c_y = \frac{\Delta y}{L},
\qquad
c_z = \frac{\Delta z}{L}
$$

The element stiffness matrix is written as:

$$
k_e =
\frac{AE}{L}
\begin{bmatrix}
\lambda & -\lambda \\
-\lambda & \lambda
\end{bmatrix}
$$

with:

$$
\lambda =
\begin{bmatrix}
c_x^2 & c_xc_y & c_xc_z \\
c_xc_y & c_y^2 & c_yc_z \\
c_xc_z & c_yc_z & c_z^2
\end{bmatrix}
$$

### Solution Strategy

1. Assemble the full global stiffness matrix `[K]` by direct stiffness summation.
2. Apply boundary conditions by **eliminating the rows and columns** of restrained DOFs.
3. Solve the reduced system for free-DOF displacements.
4. Recover the complete displacement vector (zero at restrained DOFs).

---

## 🧾 Input Definition

All model data are defined directly at the top of each script.

### Common Units

| Quantity | Unit |
|---|:---:|
| Coordinates | m |
| Cross-sectional area `A` | m² |
| Young's modulus `E` | Pa |
| Yield stress | Pa |
| Forces | N |
| Force angle (2D only) | degrees |

### 2D Inputs (`truss2D.m`)

**Nodes** — `[x y]` per row:

```matlab
data{1,1}.nodes = [0 0
                   3.8 9.12/7
                   ...];
```

**Elements** — `[node1 node2 A E yieldStress]`:

```matlab
data{1,1}.elements = [1 2 0.0005 200e9 20e7
                      ...];
```

**Forces** — `[nodeNumber F angle]`, with components computed as `Fx = F·cos(θ)`, `Fy = F·sin(θ)` (`θ` in degrees):

```matlab
data{1,1}.forces = [1 5700 -90
                    ...];
```

**Supports** — `[nodeNumber type orien]`:

| `type` | `orien` | Restrained DOFs |
|:---:|:---:|---|
| `1` | `1` | `UX` only |
| `1` | `2` | `UY` only |
| `2` | any | `UX` + `UY` (pin) |

```matlab
data{1,1}.supports = [1 2 0      % pin at node 1
                      8 1 2];    % UY restraint at node 8
```

### 3D Inputs (`truss3D.m`)

**Nodes** — `[x y z]` per row.

**Elements** — `[node1 node2 A E yieldStress]` (same format as 2D).

**Forces** — direct Cartesian components `[nodeNumber Fx Fy Fz]`:

```matlab
data{1,1}.forces = [1 0 0 -25000
                    ...];
```

**Supports** — `[nodeNumber Dx Dy Dz]`, where each translational DOF is defined independently:

| Flag | Meaning |
|:---:|---|
| `0` | Free |
| `1` | Restrained |

For example `[3 1 1 1]` fully fixes node 3 in X, Y, and Z.

---

## ⚙️ Analysis Workflow

1. Define nodal coordinates.
2. Define element connectivity, area, elastic modulus, and yield stress.
3. Compute element lengths and orientation parameters (angle / direction cosines).
4. Build element stiffness matrices.
5. Assemble the global stiffness matrix.
6. Assemble the global load vector.
7. Apply boundary conditions.
8. Solve the reduced linear system.
9. Recover the full displacement vector.
10. Compute element elongation and axial force.
11. Compute element stress and strain.
12. Compute yield-stress utilization.
13. Visualize undeformed and deformed configurations.

---

## 📊 Results

Per-element results are reported in the MATLAB command window:

| Quantity | Symbol | Unit | Sign Convention |
|---|:---:|:---:|---|
| Elongation | `ΔL` | µm | positive = elongation |
| Axial force | `N` | N | positive = tension |
| Normal stress | `σ = N/A` | MPa | positive = tension |
| Strain | `ε` | **µε (microstrain)** | positive = elongation |
| Utilization | — | % | always positive |

Yield-stress utilization is computed as:

$$
\text{Utilization} =
\left| \frac{\sigma}{\sigma_y} \right| \times 100
$$

Members whose utilization exceeds **100%** are drawn in red in the structure plot.

---

## 🎨 Visualization

Each script produces a labeled plot of the structure:

| Item | Representation |
|---|---|
| Undeformed members | Blue (red if utilization > 100%) |
| Deformed configuration | Magenta overlay (indicative visual aid) |
| Node numbers | Labeled at node positions |
| Element numbers | Labeled at member midpoints |
| Supports | Marked at restrained nodes |

3D support marker colors indicate the restrained DOF combination:

| `Dx Dy Dz` | Color |
|:---:|---|
| `1 1 1` | ⚫ Black |
| `1 1 0` | ⚪ White |
| `1 0 0` | 🟢 Green |
| `0 1 0` | 🔵 Cyan |
| `0 0 1` | 🟣 Magenta |
| `0 1 1` | 🫒 Olive |
| `1 0 1` | 🟠 Orange |
| `0 0 0` | 🔵 Blue (free) |

---

## 📁 Project Structure

```text
MATLAB-Truss-Analysis/
├── 2D/
│   └── truss2D.m
├── 3D/
│   └── truss3D.m
├── LICENSE
└── README.md
```

---

## 🚀 Getting Started

### Requirements

- MATLAB (base installation — **no toolboxes required**)

### Installation

```bash
git clone https://github.com/amirhh-04/MATLAB-Truss-Analysis.git
cd MATLAB-Truss-Analysis
```

### Usage

```matlab
% 2D planar truss
run('2D/truss2D.m')

% 3D space truss
run('3D/truss3D.m')
```

To model your own structure, edit the Inputs Section at the top of the corresponding script (nodes, elements, forces, supports) and re-run it.

---

## 📖 Assumptions & Scope

The implementation performs linear static analysis of ideal truss structures.

### Model assumptions

- Members carry **axial force only** (pin-jointed idealization)
- **Linear elastic** material behavior
- Constant cross-sectional properties per element
- **Small displacements**
- Static loads
- Translational DOFs only

### Not included

- Geometric nonlinearity
- Material plasticity
- Buckling analysis
- Dynamic / transient / modal analysis
- Fatigue analysis
- Joint slip modeling
- Design-code compliance checking

> **Notes:**
>
> - Reported strain is expressed in **microstrain (µε)**.
> - The deformed-geometry overlay is provided as a **qualitative visual aid** alongside the exact numerical element results.
> - The nodal force vector `[K]{u}` is computed internally; explicit support **reaction reporting** is planned for a future release.

The project should be considered an **educational and computational-mechanics implementation**, not a certified structural design package.

---

## 🗺️ Roadmap

- [ ] Support reaction reporting
- [ ] Deformation scale factor control
- [ ] Improved input validation
- [ ] Verification & benchmark examples
- [ ] Integration of the GUI version
- [ ] Expanded material and section libraries
- [ ] Structural optimization features

---

## 🤝 Contributing

Issues and pull requests are welcome — feel free to open an issue for bug reports, benchmark models, or feature suggestions.

---

## 📄 License

This project is licensed under the MIT License — see the `LICENSE` file for details.

---

## 👤 Author

**Amirhh**

|  |  |
|---|---|
| 🐙 GitHub | [amirhh-04](https://github.com/amirhh-04) |
| 🌐 Website | [amirhh.ir](https://amirhh.ir) |
| 💼 LinkedIn | [amirhh](https://linkedin.com/in/amirhossein-hassan) |

---

## 🔍 Keywords

MATLAB · Finite Element Method · FEM · FEA · Finite Element Analysis · Truss Analysis · Truss Solver · 2D Truss · 3D Truss · Space Truss · Structural Analysis · Structural Engineering · Direct Stiffness Method · Stiffness Matrix · Computational Mechanics · Mechanical Engineering
