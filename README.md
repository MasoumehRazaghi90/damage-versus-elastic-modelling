# Damage versus elastic modelling of human back skin

MATLAB/GIBBON and FEBio scripts for inverse analysis of human back skin and simulations of auxetic skin mesh patterns and a Zimmer mesher. The paper compares hyperelastic and damage based Ogden–GOH formulations at applied stretch λ = 1.25 and extends the damage simulations to λ = 1.40.

## Files

| Folder or file | Purpose |
| --- | --- |
| `data/Ni_Annaidh_2012/*_wd.mat` | Experimental parallel and perpendicular stress–stretch data for the predamage fitting range. |
| `data/Ni_Annaidh_2012/HumanBack_stress_stretch_{para,perp}.mat` | Full experimental curves including the damage range. |
| `inverse-analysis/HumanBack_wd_NEO.m` | Predamage neo-Hookean–GOH fit. The matrix uses FEBio's Ogden form with `m1 = 2`. |
| `inverse-analysis/HumanBack_wd_OGDEN.m` | Predamage one-term Ogden–GOH fit. |
| `inverse-analysis/damage_HumanBack_NEO.m` | Unified damage model for the neo-Hookean–GOH mixture. |
| `inverse-analysis/damage_HumanBack_OGDEN.m` | Unified damage model for the Ogden–GOH mixture. |
| `inverse-analysis/twodamage_HumanBack_OGDEN.m` | Initial separate ground matrix and fibre damage model. |
| `inverse-analysis/twodamage_HumanBack_OGDEN_fixed.m` | Intermediate revision addressing the perpendicular fit. |
| `inverse-analysis/twodamage_HumanBack_OGDEN_fixed_Last.m` | Final supplied two-damage Ogden–GOH fit. |
| `mesh-models/damage/auxetic_patterns_unload.m` | Damage based auxetic pattern simulations, including unloading. |
| `mesh-models/damage/zimmer_60_60_3mm_unload_R1.m` | Damage based Zimmer mesher simulation, including unloading. |
| `mesh-models/hyperelastic/auxetic_patterns.m` | Hyperelastic auxetic pattern simulations. |
| `mesh-models/hyperelastic/zimmer_60_60_3mm.m` | Hyperelastic Zimmer mesher simulation. |
| `mesh-models/hyperelastic/allGeometries.m` | Geometry function used by both auxetic scripts. |

`para` and `perp` refer to loading parallel and perpendicular to Langer's lines. `_wd` means without damage in the inverse fit. `_unload` indicates the loading/unloading mesh workflow; those scripts contain damage material definitions.

## Run locally

1. Install MATLAB, GIBBON and FEBio. The scripts construct FEBio 4.0 model specifications.
2. Change the `gibbonFolder` setting in each inverse analysis script to your own GIBBON installation. The packaged inverse scripts point to `data/Ni_Annaidh_2012` in this repository. The attachment filename suffixes `(1)` and `(2)` were removed to match the scripts' `load` calls.
3. Add `mesh-models/hyperelastic` to your MATLAB path before running the damage auxetic script, since it calls `allGeometries`. Configure any remaining machine specific paths and MATLAB working directory for your installation.
4. Run `HumanBack_wd_NEO.m` or `HumanBack_wd_OGDEN.m` for predamage calibration. The `damage_*` and `twodamage_*` scripts use the full curves. Use `twodamage_HumanBack_OGDEN_fixed_Last.m` for the final supplied two-damage calibration, then run the corresponding mesh scripts.

The scripts have not been executed or numerically validated in this packaged layout. The paper draft is not included. Add the final article citation and DOI here when published: **[citation pending]**.

The experimental `.mat` data may have separate source and redistribution terms. Verify permission to redistribute them before making the repository public. The MIT license covers code that the copyright holder has authority to license.

## License

Copyright © 2026 Masoumeh Razaghi. See [LICENSE](LICENSE).
