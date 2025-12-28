# Coaxial SBT Debugging Notes

## Summary
Coaxial SBT simulations are failing due to numerical instability in the solver, not unit conversion issues.

## Unit Conversion Validation
✅ **Our implementation is CORRECT:**
- `Diameter1` (wellbore) from `radius-vertical-select` slider: Already a diameter (range 0.2159-0.4445 m)
- `Diameter2` (center pipe) from `radius-lateral-select` slider: Radius (range 0.0635-0.174 m) → Converted to diameter (`* 2`) before passing to `run_sbt_final`
- `set_tube_geometry` expects diameters and divides by 2 to get radii: `radius = Diameter1/2`, `radiuscenterpipe = Diameter2/2`

This matches the original `sbt.py` implementation.

## Test Results
- **With correct unit conversions** (`Diameter1=0.4572`, `Diameter2=0.254`):
  - Simulation runs but becomes numerically unstable
  - Condition numbers: ~1.6-3.7e8 (below threshold of 5e8)
  - Eventually produces invalid temperatures (negative values)
  - Fails after many iterations

- **With incorrect units** (as in original tests: `Diameter1=0.2286`, `Diameter2=0.127`):
  - Fails immediately with invalid geometry
  - Produces invalid temperatures right away

## Root Cause
The coaxial solver has inherent numerical stability issues:
1. High condition numbers (~1.6-3.7e8) indicate ill-conditioned matrices
2. Solver becomes unstable over time, even when starting with reasonable values
3. Condition number threshold (5e8) catches extreme cases but not intermediate instabilities

## Current Status
- ✅ Unit conversions are correct
- ✅ Error handling is in place (catches invalid temperatures, ill-conditioned matrices)
- ✅ Summary table handles empty data gracefully
- ⚠️ Coaxial solver has numerical stability issues for certain parameter combinations

## Recommendations
1. **Parameter validation**: Add checks to prevent invalid geometry combinations
2. **Solver improvements**: May need to investigate the coaxial solver implementation in `sbt_v27.py`
3. **User guidance**: Document which parameter combinations work for coaxial geometry
4. **Alternative approach**: Consider if coaxial geometry needs different numerical methods

## Working Parameter Examples
From terminal output, some cases show reasonable temperatures initially:
- `mdot=30`, `Diameter1=0.4572`, `Diameter2=0.254`: Starts with reasonable temps (193-345°C) but becomes unstable
- Need to find parameter combinations that remain stable throughout the simulation

