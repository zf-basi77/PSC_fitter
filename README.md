# Post-synaptic current (PSC) Waveform Analysis

Python/Jupyter notebook for batch analysis of electrophysiology waveform data exported from PatchView.

The analysis is designed to process PatchView helper files and corresponding waveform files (see examples folder), identify individual events, fit their decay kinetics, and extract several descriptive properties of each event.

*** The fitting function employs a differential evolution algorithm to ensure best fits; parameters can be tuned to optimize fits. Along with safeguards, the algorithm allows for the fitting of all possible PSCs rather than just a 'well-behaved' average; this allows for granular analysis of PSC kinetics across recording length, conditions, etc. Visualization tools are included. Given a fit, different metrics (adj_r^2) are included to assess how well the fit represents the raw input. Additional diagnostic metrics will be included because R^2 alone may not adequately characterize model fit in nonlinear regression, particularly when residuals deviate from the assumptions of the statistical model. ***

## Overview

The main analysis notebook:

- Finds matching PatchView helper and waveform files.
- Extracts cell, series, and date information from filenames.
- Filters waveform data to columns corresponding to recorded events.
- Uses the `InUse` information from the helper file to select events for analysis.
- Extracts individual waveform events from the recorded traces.
- Fits either a single- or double-exponential decay model.
- Optionally uses differential evolution to obtain initial parameters for double-exponential fitting.
- Calculates goodness-of-fit metrics.
- Measures event amplitude and rise kinetics.
- Stores raw and fitted signals alongside the extracted measurements.
- Combines the results into a pandas DataFrame for downstream analysis.

The resulting dataframe contains both event-level measurements and the original waveform data, allowing individual events to be inspected after the batch analysis.

## Models

### Single exponential

The single-exponential model is:

```python
y = y0 + A * exp(-(x - x0) / t)
```

implemented as:

```python
def expFit(x, y0, A, t, x0):
    return y0 + A*np.exp(-(x - x0) / t)
```

The resulting `t` parameter is stored as `tau`.

### Double exponential

The double-exponential model is:

```python
y = y0 + A1 * exp(-(x - x0) / t1) + A2 * exp(-(x - x0) / t2)
```

implemented as:

```python
def dblExpFit(x, y0, A1, x0, A2, t1, t2):
    return y0 + A1*np.exp(-(x - x0) / t1) + A2*np.exp(-(x - x0) / t2)
```

When the double-exponential model is used, the two time constants are retained separately as `fast_tau` and `slow_tau`.

The notebook can use `scipy.optimize.differential_evolution` to estimate initial parameters for the double-exponential fit before refinement with `curve_fit`. This is controlled by the `use_evolution_estimate` flag.

## Analysis workflow

The main analysis proceeds approximately as follows:

```text
PatchView helper files
        │
        ├── identify cell / series / date
        │
        ▼
Find matching waveform file
        │
        ▼
Load waveform + helper data
        │
        ▼
Select events marked InUse
        │
        ▼
Extract individual events
        │
        ▼
Fit exponential decay
        │
        ├── single exponential
        │
        └── double exponential
              │
              └── optional differential-evolution initialization
        │
        ▼
Calculate event properties
        │
        ├── tau / fast tau / slow tau
        ├── R²
        ├── adjusted R²
        ├── RMSE ratio
        ├── amplitude
        ├── rise slope
        └── rise time
        │
        ▼
Store raw and fitted signals
        │
        ▼
Final pandas DataFrame
```

## Extracted measurements

| Column | Description |
|---|---|
| `tau` | Time constant from the single-exponential fit |
| `fast_tau` | Fast component of the double-exponential fit |
| `slow_tau` | Slow component of the double-exponential fit |
| `r2` | R² goodness-of-fit measure |
| `adjR2` | Adjusted R² |
| `rmse_ratio` | RMSE relative to the standard deviation of the data |
| `amplitude(pA)` | Estimated event amplitude in pA |
| `riseSlope` | Linear slope of the rising portion of the event |
| `riseTime_full(ms)` | Time between the event peak and the median pre-event baseline |
| `peak_id` | Identifier for the individual event |
| `cell` | Cell identifier |
| `series` | Series identifier |
| `date` | Recording date |
| `vibes` | Qualitative indicator of fitting success |
| `raw_sig` | Original waveform used for analysis |
| `fit_sig` | Fitted waveform |
| `full_sig` | Full original trace surrounding the event |

Additional downstream processing adds fields such as `dir_check`, `ampi`, `combo`, `orient_flag`, `fast_tau(ms)`, and `slow_tau(ms)`.

## Input data

The notebook expects a directory containing PatchView helper and waveform files.

The input directory is specified with:

```python
master_path = r'M:\HK patchview data\E2'
```

The helper/waveform files are identified automatically by `get_files()`:

- waveform files contain `csv` and `wave` in their filenames
- helper files contain `tsv` and `helper`, or `csv` and `helper`

The notebook then matches helper and waveform files using the recording date, cell, and series identifiers extracted from the filenames.

### Important

The `master_path` is specific to the original analysis environment and should be changed to the location of the PatchView files on your system.

## Main configuration

```python
do_the_double = True
pflag = False
check_min_pos = True

use_evolution_estimate = True
evol_plot = False
```

### `do_the_double`

Controls whether events are fit with a double or single exponential.

`True` uses the double-exponential model and returns `fast_tau` and `slow_tau`.

`False` uses the single-exponential model and returns `tau`.

### `pflag`

Controls diagnostic plotting during the analysis. Plotting is disabled by default because the analysis may process thousands of individual events.

### `check_min_pos`

Controls an additional check on where the minimum of the extracted event occurs.

If the minimum is displaced slightly from the expected starting position, the waveform can be trimmed to begin at that minimum. Events with a substantially displaced minimum are skipped.

### `use_evolution_estimate`

Controls whether differential evolution is used to obtain initial parameters for the double-exponential fit.

This can substantially increase computational time, but provides a more flexible method for finding starting parameters.

### `evol_plot`

Controls diagnostic plots for the differential-evolution fitting procedure. This is intended primarily for inspecting individual fits rather than running large batch analyses.

## Fitting strategy

For ordinary fitting, the notebook uses `scipy.optimize.curve_fit` with the `dogbox` method.

The fitting procedure attempts an initial fit with a limited number of function evaluations. If this fails, it retries with a larger `max_nfev`.

For differential-evolution fitting, the workflow is:

1. Identify the portion of the waveform used for fitting.
2. Define a sum-of-squared-error objective function.
3. Generate parameter bounds.
4. Run `differential_evolution` with a fixed seed.
5. Use the resulting parameters as initial values for `curve_fit`.
6. If necessary, retry the curve fit using the `dogbox` method and a larger iteration limit.

The notebook records a qualitative fitting status in `vibes`:

- `good` — initial fitting attempt succeeded
- `bad` — fitting required the fallback procedure

## Goodness-of-fit metrics

### R²

R² is calculated using `r2_score(yData, y_fit)` and is used as an indicator of fit quality.

### Adjusted R²

Adjusted R² penalizes the six-parameter double-exponential model for its number of parameters.

### RMSE ratio

The notebook defines RMSE relative to the standard deviation of the data as `rmse_ratio`. This provides a measure of model error relative to the variability of the waveform.

## Rise kinetics

The notebook also characterizes the rising portion of each event.

`get_riseSlope_dx()` identifies a baseline region, determines the portion of the rising waveform to analyze, and performs a linear regression over that region.

It returns:

- `riseSlope` — slope of the rising phase
- `dx` — duration of the analyzed rise

The rise time is subsequently converted to milliseconds and stored as `riseTime_full(ms)`.

## Downstream analysis

After the main event-level analysis, the notebook performs several additional checks and transformations, including:

- checking the orientation of fitted signals
- calculating an additional amplitude measure
- generating a combined recording identifier
- converting time constants to milliseconds

The recording identifier is constructed from:

```python
date + '_' + cell + '_' + series
```

and stored as `combo`.

The notebook also includes `plot_sigNfit()`, which can be used to visually compare an individual raw signal with its fitted signal.

## Dependencies

The analysis uses:

```text
numpy
pandas
matplotlib
seaborn
scipy
scikit-learn
tqdm
```

The main scientific functionality comes from:

- `scipy.optimize.curve_fit`
- `scipy.optimize.differential_evolution`
- `scipy.stats.linregress`
- `sklearn.metrics.r2_score`

## Running the analysis

1. Open `postPatchView_analysis.ipynb` in Jupyter Notebook or JupyterLab.
2. Update `master_path` to the directory containing the PatchView files.
3. Set the fitting options at the beginning of the main analysis cell.
4. Run the function-definition cells.
5. Run the main analysis cell.
6. The resulting event-level dataframes are stored in `dfs`.

The individual recording dataframes are concatenated downstream as needed.

## Notes

This notebook was developed for a specific PatchView electrophysiology analysis workflow and therefore contains assumptions about:

- PatchView filename conventions
- waveform sampling rate
- event positioning within extracted traces
- helper-file structure
- the location of the event relative to the extracted waveform
- initial parameters used for exponential fitting

In particular, the analysis currently uses a sampling rate of **20 kHz** when converting waveform indices to time.

The notebook is therefore not intended to be a completely general-purpose PatchView analysis package without adapting these assumptions to the relevant acquisition configuration.

## Repository structure

```text
.
├── README.md
└── postPatchView_analysis.ipynb
```

The notebook contains both the analysis functions and the batch-processing workflow.
