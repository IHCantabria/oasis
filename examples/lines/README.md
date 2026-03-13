# Group 3: Mooring Lines Examples

Isolated tests for mooring line dynamics. Each example tests a single aspect of line behavior (material model, BCP type, seabed interaction, friction) while keeping other subsystems minimal. A floating body provides dynamic excitation through heave free-decay unless otherwise noted.

## Examples

| Example | What it tests | BCP types | Key configuration |
|---|---|---|---|
| `single_line` | Basic mooring coupling | 1 Anchor, 1 Fairlead | 1 line, linear elastic, L=120m, tension-only |
| `multi_line` | Multi-point mooring system | 4 Anchors, 4 Fairleads | 4 lines at 90° spacing, L=140m each |
| `joint_connection` | Joint BCP coupling between segments | 2 Anchors, 2 Joints, 2 Fairleads | 4 lines: Anchor→Joint→Fairlead chain |
| `elastic_anchor` | Elastic anchor restoring force | 4 Elastic Anchors, 4 Fairleads | F = −c(1−e^(−kx)), c=50kN, k=0.1/m |
| `prescribed_motion` | Prescribed endpoint dynamics (no body) | 1 Actuator, 1 Anchor | Sinusoidal surge A=2m T=8s, no floating body |
| `viscoelastic` | Viscoelastic material (memory kernel) | 1 Anchor, 1 Fairlead | flag_stiffness=1, polynomial elastic + kernel |
| `tabulated_stiffness` | Tabulated strain-stress curve | 1 Anchor, 1 Fairlead | flag_stiffness=5, 5-point interpolation |
| `tension_symmetric` | Symmetric tension model (allows compression) | 1 Anchor, 1 Fairlead | flag_tension=1 vs default tension-only |
| `seabed_contact` | Ground contact (no friction) | 1 Anchor, 1 Fairlead | L=180m (excess on floor), GK=1e6, GC=0.1 |
| `friction_isotropic` | Isotropic stick-slip friction | 1 Anchor, 1 Fairlead | frictionModel=1, flat seabed |
| `friction_anisotropic` | Anisotropic friction on inclined seabed | 1 Anchor, 1 Fairlead | frictionModel=2, 10° inclined seabed |

## Shared Configuration

All examples use (unless noted):
- **Simulation time:** 30s
- **Integration:** ESDIRK46 (method=3), order=2, dt=0.005s (line dynamics need smaller step)
- **Water depth:** -100m, density=1025 kg/m³, gravity=9.81 m/s²
- **Body:** isobara.ehydb, RAD_DIFF, 6 DOFs, initial heave displacement z=1m, radiation ON
- **Waves:** REG H=0.0 T=10 (no wave excitation)
- **Line type:** Mooring (lineType=1)
- **Seabed:** Flat at -100m (from dataSeaFloor.dat)
- **Mooring IC:** Newton's method (flag=1)

## Running

Run a single example:
```bash
cd <example_name>
run.bat          # Windows
sh run.sl .      # Linux
```

Run all line examples (included in the master test suite):
```bash
cd examples
run_all_examples.bat    # Windows
sh run_all_examples.sh  # Linux
```

## Expected Behavior

- **single_line:** Body oscillates in heave with mooring line restoring force; line tension varies cyclically
- **multi_line:** Symmetric 4-point mooring constrains all DOFs; body decays faster than single_line
- **joint_connection:** Two-segment mooring chain; joint BCP moves freely between upper and lower segments
- **elastic_anchor:** Anchor positions shift under load (non-zero anchor displacement in output)
- **prescribed_motion:** No body output; line tension follows 8s sinusoidal cycle matching actuator motion
- **viscoelastic:** Tension response shows hysteresis (phase lag) compared to elastic baseline
- **tabulated_stiffness:** Non-linear tension-strain relationship; stiffening at higher strains
- **tension_symmetric:** Line can develop negative tension (compression); compare with single_line
- **seabed_contact:** LineIni shows nodes resting on seabed (z ≈ -100m); ground reaction forces present
- **friction_isotropic:** Modified tensions vs no-friction baseline; stick-slip transitions visible
- **friction_anisotropic:** Directional friction effects on inclined seabed; different tangential/normal responses

## Relevant Source Code

- `src/Lines/Lines_Dyn.cpp` — `ReadPropertiesASCII()` for input format; `SEM_computeF()` for force computation
- `src/Lines/Lines_QS.cpp` — Quasi-static catenary solver for initialization
- `src/BCPs/BCPs.cpp` — BCP reading (Actuator/Anchor/Joint/Fairlead/ElasticAnchor)
- `src/SeaFloor/SeaFloor.cpp` — Seabed geometry (flat, inclined, variable bathymetry)
- `src/Simulations/Simulation.cpp` — `ReadBcpsASCII()` for BCP index allocation; `ReadLinesASCII()` for line loading
