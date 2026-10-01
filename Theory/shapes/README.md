# Shapes series: which particle shape can spin and orbit?

Mathematica notebooks that apply the free-molecule surface stress of the corrected notes (`../ThermophoreticForce_v22.pdf`, Eq. 50) to simple near-spherical shapes in a uniform temperature gradient, with no trap. For each shape they give the Mathematica definition, a diagram, the four response tensors (plus the translation–rotation coupling B), simulated trajectories, and the orbital radius and period where an orbit exists.

| notebook | shape | result |
|---|---|---|
| `Shape1_Ellipsoid.nb` | triaxial ellipsoid, semi-axes 1 + ε e_i | T_th = 0 exactly; orientation neutral; straight-line glide; no orbit |
| `Shape2_EllipsoidOffset.nb` | the same ellipsoid with a CoM offset in any direction | a generic offset is chiral (C1), but T_th = −[r_c]× f_th cannot drive a spin: the particle hangs, then glides or hovers; no orbit |
| `Shape3_Harmonics.nb` | R = 1 + ε s(n), s built from spherical harmonics | chirality first appears at O(ε²) from l = 2 × l = 3. An ellipsoid with an xyz twist spins in place; tilting the ellipsoid against the twist gives a stable orbit (R = 24.5 µm, period 0.436 s for the defaults) |

`ShapeLevitation.wl` is the shared package:

- `symbolicShape`: an exact ε-series of the tensors, the CoM and the inertia, via sphere moments.
- `makeParticle`: SI tensors at finite ε, by quadrature.
- `steadyMotions`: exact steady spin and orbit states, with their linear stability.
- `simulate`: the full nonlinear rigid-body equations.
- Also: `orbitFromRun`, `hoverDrift`, `convexityMargin`, `shapeDiagram`, `bodyMesh` and `drawBody`.

The notebooks load the package from their own directory.

The sources are in `src/*.src`. To rebuild a notebook and its PDF:

```
cd Theory/shapes && SHAPES_DIR=$PWD ../mathematica_tools/build.sh src/Shape3_Harmonics.src Shape3_Harmonics.nb Shape3_Harmonics.pdf
```

In `Shape3_Harmonics.nb`, the interactive explorer (Section 8) starts idle. Tick **compute** to run it.
