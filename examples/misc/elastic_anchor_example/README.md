# Elastic Anchor Example

Test case for ElasticAnchorBCP (type 5) - compliant mooring anchors with exponential restoring force.

## Description

Based on generic_example, simplified to focus on elastic anchor behavior:
- **4 elastic anchors** (type 5) replace the original fixed anchors
- **4 mooring lines** connecting directly from elastic anchors to body fairleads (no joints)
- **Linear hydrostatics** (no STL mesh, no sinking)
- **No wind**, no wind turbines, no winches, no springs

## Mooring Configuration

```
  Fairlead 1 (BCP 1)          Fairlead 2 (BCP 2)
     (+10, +2.5, 0)              (-10, +2.5, 0)
         \    BODY    /
          \__________/
          /          \
         /            \
  Fairlead 3 (BCP 3)   Fairlead 4 (BCP 4)
     (+10, -2.5, 0)      (-10, -2.5, 0)
         |                    |
         | Line 1-4           | Line 1-4
         |                    |
  ElasticAnchor 1 (BCP 5)  ElasticAnchor 2 (BCP 6)
     (+20, +25, -12.5)       (-20, +25, -12.5)

  ElasticAnchor 3 (BCP 7)  ElasticAnchor 4 (BCP 8)
     (+20, -25, -12.5)       (-20, -25, -12.5)
```

## Elastic Anchor Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| mass | 500 kg | Anchor mass (dynamic inertia) |
| volume | 0.1 m³ | For drag calculation |
| c_param | 50000 N | Ultimate holding capacity |
| k_param | 0.1 1/m | Rate parameter (63% at 10m displacement) |

## Force Law

```
F_restoring = -c * (1 - exp(-k * x)) * direction
```

Where x is displacement from reference position. Force saturates at c = 50 kN.

## Running

### Windows
```
cd elastic_anchor_example
run.bat
```

### Linux (SLURM)
```
cd elastic_anchor_example
sbatch run.sl
```

## What to Check

1. Simulation completes without errors
2. Elastic anchor positions deviate from reference under mooring tension
3. Restoring forces remain bounded by c_param = 50 kN
4. No NaN or Inf in output
