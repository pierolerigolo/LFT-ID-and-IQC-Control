# LFT identification and IQC robust control (W2 pipeline)

Data-driven robust controller design for a plant whose resonance varies with an
added mass. Run the three scripts in order from the repository root (MATLAB with
the Control System, Robust Control and System Identification toolboxes, plus
YALMIP and MOSEK).

| Step | Script | Input | Output |
|---|---|---|---|
| 1 | `a_MultimodelFRFID.m` | `Data/MIVIDA 261028/*_xsswp.mat` (one sine sweep per mass) | `Data/GFRF_iddata_1028_spa.mat` (set `save_output = true`) |
| 2 | `b_LFT2ndorderID.m` | `Data/GFRF_iddata_1028_spa.mat` | `Data/LFT_GFRF_GdmGLTIW2n<nW2>_1028_spa.mat` |
| 3 | `c_IQCcontrolSynW2.m` | the LFT file of step 2 (`idname`) | `Data/controller/*.mat`, `Data/diaries/*.txt` |
| 3 (alt.) | `c_IQCcontrolSynElli.m` | the measured FRFs and grid of the step-2 file | `Data/controller/ELLI_*.mat`, `Data/diaries/ELLI_*.txt` |

## Model (step 2)

    G(dm, Delta) = G_dm(m0 - dm) * G_LTI * (1 + W2 * Delta)

- `G_dm = s^2/(m s^2 + c s + k)`: grey-box fit on 25-75 Hz, one mass per
  operating point, common `c`, `k`; Tustin-discretised so the LFT in `dm` is exact.
- `m0`, `dm_bound`: centre and half-width of the identified masses,
  `dm` in `[-dm_bound, dm_bound]` (real, centred, not normalised).
- `G_LTI`: mean of the residuals `G_i / G_dm,i`.
- `W2`: `ucover` fit of the relative residual error; `Delta` is a 1x1 LTI block
  with `||Delta|| <= 1`.

## Synthesis (step 3)

Fixed-order discrete state-space controller, frequency-domain LMIs on the LFT
grid (Gupta, Klauser, Karimi 2024; Schuchert, Gupta, Karimi 2024).

1. Robust stability only, constant multiplier at `eta0`.
2. Minimise nominal `||W1 S||_inf` subject to robust stability (stage-1 multiplier).
3. Alternate: largest `eta` with an optimised multiplier `D(w)`, `E(w)` for the
   frozen controller (bisection), then the stage-2 step with that multiplier.

Any step that increases the objective is rejected and rolled back. The saved
`K_final` is certified robustly stable for `|dm| <= eta`; `eta >= dm_bound`
covers the whole identified mass range.

## Code layout

- `+datadriven/+Controller`: controller parametrisation (`SS`, frozen copy `Fixed`).
- `+datadriven/+Objective`: `Hinf`, `IQC_reduced`, multiplier optimisation
  (`IQC_reduced_combinedPI_par_fix`, `..._eta_bisection`).
- `+datadriven/+Constraint/polygonal_chain.m`: nominal stability constraint.
- `+utils`: data loading, `G_dm` grey-box fit, magnitude plots. `Gdm.m` is the
  grey-box model used by `idgrey`.

## Elliptical-set alternative (`c_IQCcontrolSynElli.m`)

Same approach as `exp_controller.m` of the reference code
(gitlab.epfl.ch/ddmac/paper-codes/2024/non-parametric-iqc-for-generalised-plant):
the nominal model and a minimum-volume ellipse per frequency are computed
directly from the measured FRFs (`+datadriven/find_elliptical_uncertainty.m`),
then robust stability only -> nominal H-inf performance subject to robust
stability, with the same settings, stopping rules and OL-stable fallbacks as
`c_IQCcontrolSynW2.m`. Same `W1` and frequency grid as the W2 design.

