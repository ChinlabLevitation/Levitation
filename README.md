# Code for the Levitation experiment at the Chin Lab at UChicago

## `sim_symmetric.jsx`

A React-based interactive simulator for a rigid body undergoing **thermophoretic levitation** in a harmonic temperature trap.

### Physics Model

Solves a 13-DOF ODE using RK4 integration:

| State variable | Description |
|---|---|
| `x` (3D) | Center-of-mass position |
| `v` (3D) | Translational velocity |
| `q` (quaternion, 4D) | Orientation |
| `ω_b` (3D, body frame) | Angular velocity |

**Equations of motion:**
- **Translation:** `M v̇ = R(−Γ ∇lnT_b − Ξ v_b) − Mg ẑ`
- **Rotation:** `I ω̇_b = −ω_b × (I ω_b) − Λ ∇lnT_b − Cᵣ ω_b`

**Temperature field** is a harmonic trap:
```
∇lnT = (αr·x,  αr·y,  αz·z − Mg/γ)
```
The last term sets the levitation equilibrium height automatically from `Mg/Γ₂₂`.

### The Four Physical Tensors

All are stored as symmetric 6-vectors `[t₁₁, t₁₂, t₁₃, t₂₂, t₂₃, t₃₃]` (24 total DOF), enforcing Onsager reciprocity and energy dissipation symmetry:

| Tensor | Symbol | Physical meaning |
|---|---|---|
| Translational drag | **Ξ** | `F_drag = −Ξ v_body` |
| Thermophoretic force | **Γ** | `F_T = −Γ ∇lnT_body` |
| Thermophoretic torque | **Λ** | `τ_T = −Λ ∇lnT_body` (zero for mirror-symmetric bodies) |
| Rotational drag | **Cᵣ** | `τ_drag = −Cᵣ ω_body` |

### How to Run

The simulator is set up as a Vite/React app at `Theory/levitation-sim/`. To start it:

```bash
cd Theory/levitation-sim
npm run dev
```

Then open the local URL printed by Vite (typically `http://localhost:5173`).

> **Dependencies:** React 19, Recharts 3, Vite 8. Run `npm install` if `node_modules` is missing (it is gitignored).

### UI

Once running, the interface has three panels:

1. **Left panel — Parameters:** Configure mass, inertia `diag(I₁,I₂,I₃)`, trap strengths `αr`/`αz`, and all four tensors (click a tensor to expand its 3×3 symmetric editor with Isotropic/Diag/Zero/Symmetrise shortcuts). Set initial position `x₀`, velocity `v₀`, angular velocity `ω₀`, and initial tilt.

2. **Middle panel — 3D View:** After clicking **▶ Run Simulation**, watch the cylindrical body move and rotate. Use the Elev/Azim sliders to orbit the camera. Play/pause/scrub through the trajectory.

3. **Right panel — Time-series charts:** Position `(x,y,z)`, angular velocity `(ω₁,ω₂,ω₃)`, tilt angle, and speed vs. time, with a red cursor line tracking the current playback frame.

### Default Parameters

| Parameter | Default | Notes |
|---|---|---|
| Mass | 1 kg | |
| Gravity | 9.81 m/s² | |
| Inertia | diag(0.02, 0.02, 0.008) | Cylinder-like |
| Ξ, Γ, Cᵣ | `5·I`, `5·I`, `0.5·I` | Isotropic |
| Λ | 0 | Mirror-symmetric body |
| Initial tilt | 30° | Around x-axis |
| Simulation time | 20 s, 1200 steps | |

