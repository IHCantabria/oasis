# Group 4: Spring Examples

Isolated tests for spring connections between bodies and fixed structures. Each example tests a different physical spring configuration while keeping other subsystems minimal. Bodies are excited by regular waves.

## Examples

| Example | What it tests | BCP types | Key configuration |
|---|---|---|---|
| `pile_connector` | Body-to-pile nonlinear spring with clearance gap and friction | 1 Anchor, 1 Fairlead | 1 spring, nonlinear stress-strain (gap+contact), velocity-based friction, line constraint on pile axis |
| `neoprene_connector` | Body-to-body neoprene+wire spring joining two floating dikes | 2 Fairleads | 1 spring, asymmetric nonlinear stiffness (neoprene compression / wire tension), damping, no friction |

## Shared Configuration

All examples use:
- **Simulation time:** 50s
- **Integration:** ESDIRK46 (method=3), order=2, dt=0.02s
- **Water depth:** -8m, density=1025 kg/m³, gravity=9.81 m/s²
- **Body:** dique_cadenas.hydb.h5, RAD_DIFF, 6 DOFs, linear hydrostatics, radiation ON, 1st order excitation with instant position
- **Waves:** REG H=0.5m T=5.0s heading=225°
- **Spring stress model:** Nonlinear (stressModelFlag=2) with stress-strain curves
- **IRF time span:** 40s

## Physical Descriptions

### pile_connector — Body Connected to a Marine Pile

A single floating dike is connected to a vertical marine pile (fixed structure) via a fender guide spring. The spring models:
- **Clearance gap**: The fender guide has ±30mm play in surge (X) and ±40mm play in sway (Y). Within these gaps, no restoring force acts.
- **Contact stiffness**: Beyond the gap, very high stiffness prevents further penetration (simulating steel-on-steel contact).
- **Friction**: Velocity-based friction (μ_d=0.02, μ_s=0.20) resists vertical sliding along the pile.
- **Line constraint**: BCP1 (pile side) uses type=1 (line constraint along Z-axis), allowing the connection point to slide vertically along the pile while constraining horizontal displacement.

### neoprene_connector — Two Dikes Joined by Neoprene + Wire

Two identical floating dikes are placed side by side as port breakwater infrastructure. They are connected at their adjacent edges by a combined neoprene pad and wire/chain connector:
- **Axial (surge) behavior**: Asymmetric — neoprene provides compression stiffness (4.5 MN/m) when dikes push together, wire provides tension stiffness (15 MN/m) when dikes pull apart.
- **Transverse (sway/heave) behavior**: 13-point nonlinear restoring force with cross-coupling to axial tension (transverse displacement causes axial reaction due to wire geometry).
- **Damping**: 500 kN·s/m in all three translational DOFs absorbs energy.
- **No friction**: The neoprene/wire connection is modeled without friction.

## Running

Run a single example:
```bash
cd <example_name>
run.bat          # Windows
sh run.sl .      # Linux
```

Run all spring examples (included in the master test suite):
```bash
cd examples
run_all_examples.bat    # Windows
sh run_all_examples.sh  # Linux
```

## Expected Behavior

- **pile_connector:** Body oscillates under wave excitation within the pile's clearance gap. When displacement exceeds ±30mm (surge) or ±40mm (sway), high contact forces appear. Friction forces resist vertical sliding. Spring forces show flat-zero plateau within gap and sharp ramps at contact.
- **neoprene_connector:** Both dikes respond to waves with slightly different phases due to heading angle. Axial spring force is asymmetric (stiffer in tension than compression). Cross-coupling generates axial loads from transverse relative motion. Damping smooths the dynamic response.

## Relevant Source Code

- `src/Spring/Spring.cpp` — `ReadPropertiesASCII()` for input format; `computeSpringForces()` for force computation
- `src/Spring/Spring.hpp` — Spring class definition, members, and enums
- `src/BCPs/BCPs.cpp` — BCP reading (Actuator/Anchor/Joint/Fairlead/ElasticAnchor)
- `src/Simulations/Simulation.cpp` — `ReadSpringsASCII()` for spring loading; `ReadBcpsASCII()` for BCP allocation
