# Group 10: Winches

Isolated tests for winch-controlled mooring systems. Each example tests a specific aspect of the winch actuation and controller while keeping other subsystems minimal.

## Examples

| Example | What it tests | Lines | Winches | Key configuration |
|---|---|---|---|---|
| `constant_tension` | Winch controller applying constant pretension to all mooring lines | 4 mooring lines (5 nodes each) at ±45° | 4 winches (one per line at fairlead end) | Constant tension = 30000 N per winch, body starts at surge=1m, free decay |
| `horizontal_control` | Winch controller for horizontal plane positioning (surge/sway/yaw) | 4 mooring lines (5 nodes each) at ±45° | 4 winches (one per line at fairlead end) | State-space controller, straight-lines inversor, reference=(0,0,0), body starts at surge=2m |

## Shared Configuration

- **Simulation time**: 30s
- **Integration**: ESDIRK46 (method=3), order=2, dt=0.005s
- **Water depth**: -100m, density=1025 kg/m³, gravity=9.81 m/s²
- **Body**: `one_box.hydb.h5`, RAD_DIFF, 6 DOFs, radiation ON
- **Waves**: REG H=0.0 T=10 (no wave excitation)
- **Lines**: Mooring (lineType=1), tension-only, linear elastic (EA=1e7 N), L=160m, 50 kg/m
- **Winches**: Inertia=50 kg·m², radius=0.25m, drag=0.5 N·s/rad
- **Controller dt**: 0.1s
- **Seabed**: Flat at -100m
- **Mooring IC**: Newton's method (flag=1)

## Physical Layout (horizontal_control)

```
                 Anchor 2 (-80,+80)         Anchor 1 (+80,+80)
                        \                   /
                         \    Body (0,0)   /
                    FL2 (-5,+5)-----FL1 (+5,+5)
                          |    one_box    |
                    FL4 (-5,-5)-----FL3 (+5,-5)
                         /                 \
                        /                   \
                 Anchor 4 (-80,-80)         Anchor 3 (+80,-80)
```

Each fairlead has a winch (winchId 1-4) that controls line length. The state-space controller reads body position (surge, sway, yaw) and commands winch torques via a straight-lines inversor (bounded optimization) to drive the body toward the reference position (0, 0, 0). The controller starts at t=5s. The body begins displaced 2m in surge, allowing free-decay dynamics before active control begins.
