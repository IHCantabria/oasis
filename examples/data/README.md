# examples/data

This folder contains the hydrodynamic databases (`.hydb.h5` files) required to run the example cases included in `examples/`.

| File | Description |
|---|---|
| `one_box.hydb.h5` | Single box — baseline hydrodynamic database |
| `one_box_5p0.hydb.h5` | Single box — 5.0 m draft variant |
| `one_box_7p5.hydb.h5` | Single box — 7.5 m draft variant |
| `one_box_10p0.hydb.h5` | Single box — 10.0 m draft variant |
| `one_box_owc.hydb.h5` | Single box with OWC module |
| `trl_plus.hydb.h5` | TRL+ platform |
| `two_boxes.hydb.h5` | Two-body system |

## How to obtain these files

These files are **not included in the repository source code** due to their size.  
They are distributed as `examples.zip` attached to each [GitHub release](https://github.com/IHCantabria/oasis/releases).

1. Go to the [Releases page](https://github.com/IHCantabria/oasis/releases).
2. Download `examples.zip` from the desired release.
3. Extract the archive and copy the contents of the `data/` folder here.

## Generating your own databases

The `.hydb.h5` files are produced by [**SeaMotions**](https://github.com/IHCantabria/SeaMotions), an open-source BEM solver developed by IHCantabria. SeaMotions computes first-order hydrodynamic coefficients (added mass, radiation damping, wave excitation forces, QTFs) for arbitrary floating bodies, including multi-body layouts and OWCs, and writes the results in the HDF5 format expected by OASIS.
