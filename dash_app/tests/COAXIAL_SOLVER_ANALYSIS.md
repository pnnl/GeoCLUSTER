# Coaxial Solver Numerical Stability Analysis

## Problem Summary
Coaxial SBT2 simulations are failing with invalid temperature values (600-1500°C or negative) even with correct geometry parameters. The solver matrix becomes ill-conditioned (condition number > 5e8), leading to numerical instability.

## Root Causes Identified

### 1. Matrix Conditioning Issues
The linear system `L * Sol = R` becomes ill-conditioned when:
- **Thermal resistance `R_cp` becomes very small**: Terms like `1/R_cp/(A_flow_annulus*rho_f*cp_f)` in matrix `L` become very large
- **Flow areas become very small**: If `A_flow_annulus` or `A_flow_centerpipe` approach zero, velocities become very large, causing instability
- **Hydraulic diameter becomes negative or zero**: If `Dh_annulus = 2*(radius - outerradiuscenterpipe)` ≤ 0, Reynolds numbers become invalid

### 2. Geometry Validation
Flow areas are calculated in `sbt_utils.py::compute_tube_geometry()`:
```python
A_flow_annulus = math.pi*(radius**2 - outerradiuscenterpipe**2)   # [m^2]
A_flow_centerpipe = math.pi*radiuscenterpipe**2                     # [m^2]
Dh_annulus = 2*(radius - outerradiuscenterpipe)                     # [m]
```

**Critical condition**: `outerradiuscenterpipe < radius` must be satisfied, otherwise:
- `A_flow_annulus` ≤ 0 (invalid)
- `Dh_annulus` ≤ 0 (invalid)
- Matrix becomes singular or ill-conditioned

### 3. Current Validation
`clgs.py` has clamping logic (lines 293-312) that prevents `Diameter2 >= Diameter1`, but:
- Clamping happens **before** `set_tube_geometry` divides by 2
- Edge cases where `outerradiuscenterpipe ≈ radius` (very thin annulus) may still cause issues
- No validation of `A_flow_annulus` or `Dh_annulus` after calculation

## Matrix Structure
The coaxial solver builds a 4×N system (4 equations per element):
1. Downflowing fluid heat balance
2. Rock temperature equation
3. SBT algorithm equation
4. Upflowing fluid heat balance

Key matrix terms that can cause instability:
- `L[i, i] = 1/Deltat + u/Deltaz + 1/R_cp/(A_flow*rho_f*cp_f)` - becomes large if `R_cp` is small or `A_flow` is small
- `L[i, j] = -1/R_cp/(A_flow*rho_f*cp_f)` - coupling terms become large if `R_cp` is small

## Current Safeguards
1. **Condition number check** (line 1586): Raises error if `cond_num > 5e8`
2. **Warning threshold** (line 1592): Warns if `cond_num > 1e8`
3. **Solution validation** (lines 1602-1611): Checks for NaN/Inf in solution
4. **Temperature validation** (lines 1626-1635): Ensures temperatures are in reasonable range (-100 to 500°C)

## Recommended Next Steps

### 1. Add Geometry Validation After Calculation
Add checks in `sbt_utils.py::compute_tube_geometry()` after calculating flow areas:
```python
if A_flow_annulus <= 0:
    raise ValueError(f"Invalid annulus flow area: {A_flow_annulus} m². "
                     f"radius={radius} m, outerradiuscenterpipe={outerradiuscenterpipe} m")
if Dh_annulus <= 0:
    raise ValueError(f"Invalid annulus hydraulic diameter: {Dh_annulus} m. "
                     f"radius={radius} m, outerradiuscenterpipe={outerradiuscenterpipe} m")
```

### 2. Add Minimum Flow Area Check
Ensure annulus flow area is not too small (e.g., > 1% of wellbore area):
```python
min_annulus_area = 0.01 * math.pi * radius**2
if A_flow_annulus < min_annulus_area:
    raise ValueError(f"Annulus flow area too small: {A_flow_annulus} m² < {min_annulus_area} m²")
```

### 3. Add Thermal Resistance Validation
Check that `R_cp` is not too small before building matrix:
```python
min_R_cp = 1e-6  # Minimum thermal resistance [K/W]
if np.any(R_cp < min_R_cp):
    raise ValueError(f"Thermal resistance R_cp too small: min={np.min(R_cp):.2e} K/W")
```

### 4. Improve Clamping Logic
Update clamping in `clgs.py` to ensure sufficient annulus clearance:
```python
# Ensure center pipe diameter is at most 80% of wellbore diameter
# This guarantees annulus area > 0.36 * wellbore_area
max_d2_ratio = 0.8
if d2 >= d1 * max_d2_ratio:
    new_d2 = d1 * max_d2_ratio
    # ... clamping logic
```

### 5. Add Matrix Preconditioning (Advanced)
Consider using iterative solvers with preconditioning for ill-conditioned matrices:
- `scipy.sparse.linalg.bicgstab` with ILU preconditioner
- `scipy.sparse.linalg.gmres` with Jacobi preconditioner

### 6. Test Parameter Ranges
Create a test script to systematically test parameter combinations:
- Various `mdot` values (10-200 kg/s)
- Various `radius` and `radiuscenterpipe` ratios
- Identify stable parameter ranges

## Test Results from `test_coaxial_quick.py`
- **Parameters**: `mdot=30`, `Diameter1=0.4572`, `Diameter2=0.254` (correct diameters)
- **Result**: Failed with invalid temperatures (Min=-127.91°C)
- **Conclusion**: Even with correct geometry, solver is unstable for these parameters

## Questions to Investigate
1. What parameter combinations produce stable results for coaxial SBT2?
2. Is there a minimum `mdot` or geometry ratio required for stability?
3. Does the original MATLAB/Python code have different validation or solver settings?
4. Should we use a different solver or numerical method for coaxial geometry?

