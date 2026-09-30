# Classic DMX runner and shared delta-DM overlay

The prenoise runner is now implemented in `classic_dmx.py` and exported as
`pint_pal.BBX.run_classic_dmx`. The processing notebook and running fits are
unchanged. The broader noise/comparison interfaces below remain design proposals.

## Implemented prenoise interface (first draft)

```python
classic = run_classic_dmx(
    model=mo, toas=to, timing_config=tc,
    binning="preserve",               # or "rebuild"
    dmx_kwargs=None,                  # e.g. {"bin_width": 6.5}
    dm_order=2,                       # None: preserve; 0: disable; n: DM1..DMn
    dm_values=None,                   # preserve values; initialize new ones at 0
    solar_wind="replace",             # "preserve", "replace", "remove"
    sw_kwargs={"NE_SW": 6.67, "SWM": 1, "SWP": 2.0},
    sw_free=("SWP",),                 # freeze all other deterministic SW params
    frozen_params=("DM", "DM1", "DM2"),
    unfrozen_params=("DMX_*",),
    fit=True,                         # False: inspect preparation without fitting
    maxiter=None,                     # use YAML iteration count
)
fo_dmx = classic.prenoise
```

Capabilities:

1. Preserve input DMX bins or rebuild with `du.setup_dmx`, using its returned
   TOAs and recording the input/output indices.
2. Keep, add, or truncate DM derivatives through the lite_utils helpers. For
   `dm_order=0`, PINT requires a structural DM1 parameter: it stays zero/frozen,
   while higher derivatives are removed.
3. Start with YAML free parameters, remove explicit frozen names, add explicit
   unfrozen names, and assign the final list once. Exact names/globs are accepted;
   contradictory or unmatched selections raise errors. Intentionally removed
   model parameters are removed from YAML's list and reported in metadata.
4. Preserve, replace or remove solar wind. `sw_free=()` fixes all SW parameters;
   `("NE_SW",)` fits density; `("SWP",)` fits only the SWM=1 exponent; both names
   fit both. `None` follows YAML and explicit overrides. `SWP` cannot be fitted
   with SWM=0; SWM and SWEPOCH remain settings rather than free parameters.
5. Discover all PINT NoiseComponent subclasses and pass their component names
   explicitly to `lu.remove_noise`, including inherited DM/SW GPs.
6. Call `find_empty_masks(freeze=True)` after applying the free list. Reject free
   DM with all populated DMX amplitudes free; an empty fixed bin is not a valid
   reference constraint. Warn about DM derivative/DMX covariance.

Selective DMX example: with YAML `free-dmx: true`, use
`frozen_params=("DM", "DMX_0001")` and `unfrozen_params=()` to freeze that bin
while leaving the other bins free. To override YAML entirely for DMX, list
specific complementary frozen/unfrozen amplitudes. An all-DMX glob cannot also
appear opposite a specific DMX name because overlap is an error.

`ClassicDMXRun.prenoise` is a normal PINT fitter for existing diagnostics.
`fitted` distinguishes preparation from fitting; `metadata` records the actual
policy, removed noise components and selected TOAs. No noise inference or output
writing occurs. Numerical fit failures propagate rather than returning success.

Verification: focused tests use real PINT parameter/component objects, current
repository helpers, and a recording fitter. Numerical fitting and actual
`setup_dmx` selection need validation in the notebook's working environment.

## Scope and module placement

Add `classic_dmx.py` for fitting and `comparison.py` for result adapters and plots
inside `pint_pal.BBX`. Use functions plus small data containers, rather than a new
runner class or changes to the BB segmentation engine. Start with narrowband
validation; do not advertise wideband support until its full path is tested.

The first milestone is a classic prenoise DMX fit and an overlay with an existing
GP reconstruction. Independent vanilla-noise inference is the next increment.

## Proposed runner interface

```python
classic = run_classic_dmx(
    model=classic_model,              # snapshot retaining original DMX
    toas=selected_toas,
    timing_config=tc,
    binning="preserve",               # or "rebuild"
    dmx_kwargs=None,                  # forwarded only when rebuilding
    noise_mode="none",                # "none", "run", or "load"
    noise_model_kwargs=None,          # explicit classic recipe
    sampler_kwargs=None,
    noise_dir=None,                   # exact directory for load/run
    output_dir=None,
    maxiter=None,                     # default tc.get_niter()
)
```

Return `ClassicDMXRun(prenoise, postnoise=None, metadata, paths)`. Each fit is a
`FitResult` with label, method, stage, fitter, fit status, stable TOA identities,
configuration snapshot and provenance. Keep the prenoise fitter independent of
the postnoise fitter. A nonconverged fit remains inspectable and visibly labeled;
it must not silently become a successful validation baseline.

Runner sequence:

1. Deep-copy model, TOAs and configuration. Validate source identity, pulse-number
   state and selection; make no changes to the active GP workflow.
2. Remove inherited stochastic components using an explicit, reviewed component
   policy. The current `lu.remove_noise()` default is not a complete GP cleanup.
   Preserve the classic deterministic DM/SW/FD prescription; do not inject the
   GP branch's DM derivatives or SW assumptions.
3. Preserve existing DMX ranges by default, rejecting missing/invalid DMX instead
   of silently choosing bins. `rebuild` invokes `dmx_utils.setup_dmx()` on the
   copies, retains its returned TOAs and records dropped identities and settings.
4. Construct the fitter through the timing configuration, resolve its free
   parameters against the classic model, and validate the DM/DMX reference
   constraint. Missing configured parameters are an explicit configuration error.
   Do not blindly unfreeze both a common DM offset and every DMX offset.
5. Fit and preserve the prenoise result. `noise_mode="none"` ends here.
6. For `run`, use `noise_utils.model_noise()` with an independently supplied
   vanilla recipe and a new output directory. Initial recipe: timing model,
   white noise and optional red noise; DM GP, SW GP and chromatic GP disabled.
   Deterministic solar wind remains governed by the classic model.
7. For `load`, load existing noise values only after validating source, TOA and
   model/configuration provenance. No implicit inference or fallback rerun.
8. Use the current `get_map_noise_values` / `add_noise_to_model` sequence on a
   copy of the prenoise model, then refit and preserve the postnoise result.
   Select a fitter supporting the attached covariance model.

The initial wrapper targets the current Discovery noise path. The old BBX
`NoiseAnalysis.noise_analysis_from_fitter()` calls an obsolete attachment API and
should not be reused unchanged. Wideband also needs explicit chain-directory
handling: `get_map_noise_values()` currently resolves paths without a wideband
argument. Saving must include configuration and input identity, not just a PAR.

## Shared series contract and adapters

`DMSeries` contains label, source/method, stage, time MJD, DM values with units,
optional intervals, covariance or joint draws, uncertainty kind, selected physical
components, support mask and reference/provenance metadata. Missing uncertainty
is represented as unavailable, never zero. Preserve full draws/covariance for
later projections instead of storing only marginal error bars.

```python
dmx = dm_series_from_fitter(classic.prenoise, quantity="total_dispersive")
gp = dm_series_from_gp(
    payload=post_proc,
    reference_model=gp_reference_model,
    label="GP + BBX",
    stage="postnoise",
    quantity="total_dispersive",
)
fig, ax, comparison = plot_ddm_overlay(
    [dmx, gp],
    reference=DMReference(value=reference_dm),
    mjd_range=None,                   # or (start_mjd, end_mjd)
)
```

These names describe the proposed API. `reference_dm` is an explicit scalar with
DM units; `gp_reference_model` is the exact model used to generate the payload.
Do not substitute a subsequently refitted model, which can double-count timing
parameter corrections. The adapter must validate that provenance.

Physical quantity choices:

- `total_dispersive`: deterministic DM trends + DMX or DM GP + modeled solar
  wind. This is the first comparison target when the two methods partition
  interstellar and solar-wind contributions differently.
- `iism`: deterministic DM trends + DMX or DM GP, excluding explicitly modeled
  SW. Label DMX as a modeled non-SW contribution, not proof that residual solar
  wind contamination is absent.
- Component-only views can follow, explicitly labeled as such. FD, red noise
  and general chromatic noise do not enter the DM sum.

For GP series, reconstruct DM, SW and relevant timing-model contributions using
matched draw indices from the joint conditional samples. Refactor/reuse the
reconstruction and unit conversion behind `plot_utils` rather than scrape plot
artists. Sum contributions before taking quantiles so cross-covariance survives.
At fixed MAP hyperparameters, label the band "conditional uncertainty".

For deterministic DMX, extract uncentered DMX values directly and add the selected
deterministic contributions. Propagate the relevant fitter covariance. Existing
`simple_dmxparse()` independently mean-centers and has a broken exception path;
it is not the shared extraction contract.

## Reference, support and plot behavior

- Subtract one explicit common scalar DM reference from every series by default.
  This preserves fitted offsets. An optional shape-only mode may later remove
  each series' weighted mean using the same reference support and weights; label
  that transformation and propagate its covariance/draws.
- `mjd_range` changes only the display. Zooming must not refit, renormalize or
  change the reference. Clip visible DMX intervals at the view boundary.
- Draw DMX as interval-supported values and GP as a curve/band on available
  evaluation points. Deterministic trends/SW may vary within DMX intervals;
  the total reconstructed signal is therefore not necessarily a flat step.
- Do not linearly interpolate DMX midpoints or extrapolate beyond valid support.
  Mask undefined DMX intervals and visibly break GP curves across declared gaps.
  GP prediction through gaps can be added as an explicitly distinguished view.
- Native-support overlays may show different TOA sets, with counts and a mismatch
  notice. Numerical comparisons require explicit common support and TOA identity
  alignment. MJD alone is not a unique TOA identifier.
- Return the plotted arrays, transformations and metadata with the figure, so
  later tables use exactly the same quantities. Plotting has no fitting or file
  writing side effects; callers choose whether to save.

Later amplitude-difference panels should project both results onto common
intervals using documented weights. Apply that same linear projection to draws
or covariance. Do not divide unmatched uncertainty arrays or interpret a
difference significance as independent: both fits use the same observations.

## Notebook placement and verification

Personal notebook: preserve the classic input; run/load the baseline; adapt the
chosen baseline and GP result; plot full span or zoom. Keep the detailed report
in `BBX_Analysis.ipynb`. Existing BB+DMX fitters use the same fitter adapter.

First implementation checks:

- Runner leaves caller model, TOAs and config unchanged; preserve mode retains
  bins/TOA identities; prenoise and postnoise states cannot alias.
- Noise run/load/none branches cannot accidentally launch inference or overwrite
  an existing run; GP recipe settings cannot leak into the vanilla run.
- Synthetic DMX and GP representations of the same signal agree after conversion
  and referencing, including deterministic trends and SW contributions.
- Zoom does not change values; gaps are not bridged; TOA reordering is detected.
- Joint-draw/covariance transformations retain uncertainty and do not double-count
  timing-model reference values. Validate against one saved real narrowband run.

Full FD/GP covariance reports, residual comparison suites and model selection
statistics remain subsequent increments, using the same result containers.
