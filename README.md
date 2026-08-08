# rpc-hv-efficiency                    

**Tracklet-based RPC efficiency vs High Voltage (HV) analysis**\
*(CMS/RPC-style workflow for laboratory or QC studies)*

------------------------------------------------------------------------

## 0) What this repository does

This repository provides a **reproducible pipeline** to measure RPC
efficiency as a function of HV using **tracklet reconstruction**. It
supports:

-   **2-chamber** (A → B) efficiency (reference → target)
-   **3-chamber** (A+B → C) efficiency (two references define a
    tracklet; third is target)
-   **4-chamber** (A+B+C → D) efficiency (three references; fourth is
    target)

It produces, for each HV point:

1.  **2D efficiency maps** (ROOT histograms):
    -   `h_total` (denominator: eligible events)
    -   `h_tracklet` (numerator: matched events)
    -   `h_eff` (efficiency map)
2.  **Diagnostic ROOT histograms** for residuals and angles:
    -   `h_dx`, `h_dy`, `h_dt`
    -   `h_dx_dy`
    -   `h_theta_x`, `h_theta_y`, `h_theta`
    -   `h_theta_x_theta_y`
    -   candidate-level residual histograms, when enabled in the producer
3.  **JSON summaries** with ROI and global efficiencies + binomial
    errors
4.  **Optional per-tracklet diagnostic tables**:
    -   CSV output with matched-tracklet quantities
    -   NDJSON output with event-level diagnostic information
5.  **Efficiency vs HV plots** (simple and CMS/RPC-group style)
6.  **Working-point validation plots** near the fitted WP:
    -   beam/local efficiency map
    -   `Δx`, `Δy`, `θx`, `θy` distributions
    -   `θx` vs `θy` 2D angular map
    -   tracklet-level correlation matrix
7.  **3D tracklet visualizations** showing reconstructed trajectories
    through all chambers
8.  **Automatic cut estimation** for `Δx`, `Δy`, `Δt` and optional
    angular consistency, evaluated per working point when desired
9.  **Regime-aware analysis** for `normal` and `beam` running
    conditions
10. **BX-tree-aware analysis**, for example `events`, `events_36BX`,
    `events_24BX`, `events_18BX`, or `events_12BX`, if those TTrees are
    present in the input `data.root` files

------------------------------------------------------------------------

# Download

``` bash
git clone https://github.com/Klaudio-Peqini/rpc-hv-efficiency
cd rpc-hv-efficiency
```

------------------------------------------------------------------------

## 1) Repository layout

    rpc-hv-efficiency/
    ├─ src/
    │  ├─ rpc_tracklet_efficiency.py
    │  ├─ rpc_hv_pipeline.py
    │  ├─ rpc_estimate_tracklet_cuts.py
    │  ├─ rpc_make_eff_summary.py
    │  ├─ rpc_efficiency_vs_hv.py
    │  ├─ rpc_cms_style_eff_plot.py
    │  ├─ rpc_plot_efficiency_maps.py
    │  ├─ rpc_plot_tracklets_3d.py
    │  └─ rpc_wp_validation_plots.py
    ├─ configs/
    │  ├─ beam_roi.json
    │  └─ default_cuts.json
    ├─ examples/
    │  ├─ run_beam_scan.sh
    │  └─ run_normal_scan.sh
    ├─ synthetic/
    │  └─ make_synthetic_rpc_hvscan.py
    ├─ README.md
    ├─ requirements.txt
    ├─ .gitignore
    ├─ CHANGELOG.md
    └─ CITATION.cff

------------------------------------------------------------------------

## 2) Requirements

### Mandatory

-   Python ≥ 3.9
-   ROOT with PyROOT enabled

Test PyROOT:

``` bash
python3 -c "import ROOT; print(ROOT.gROOT.GetVersion())"
```

### Python dependencies

The usual analysis and plotting workflow requires packages such as
`numpy`, `matplotlib`, and `pandas`, in addition to ROOT/PyROOT.
Install the repository dependencies with:

``` bash
python3 -m pip install -r requirements.txt
```

------------------------------------------------------------------------

## 3) Physics definition

Efficiency is defined as:

ε = N_matched / N_eligible

-   Eligible events are defined only by reference chambers.
-   Target chamber never contributes to denominator.
-   Matching requires compatibility within:
    -   `tol-x`
    -   `tol-y`
    -   `tol-t`

In the updated workflow, these tolerances can be used in two different
ways:

-   **fixed-cut mode**: user-provided values are applied identically to
    all HV points
-   **dynamic-cut mode**: a pre-processing step estimates robust cut
    values from the data for each working point, then passes them to the
    tracklet selection script

For beam studies, the same detector can show different apparent behavior
from the natural background case. For this reason the pipeline supports:

-   `--regime normal`
-   `--regime beam`

This is useful both for bookkeeping and for applying distinct cut
estimation, ROI handling, and output naming.

### 3.1 Residuals and angular diagnostics

For a matched anchor-target pair, the residuals are defined as:

-   `dx = xT - x_pred`, after the alignment correction used by the
    producer
-   `dy = yT - y_pred`, after the alignment correction used by the
    producer
-   `dt`, the timing residual used in the matching condition

For a two-chamber diagnostic, the pair angles are:

``` text
θx = atan(dx / Δz)
θy = atan(dy / Δz)
θ  = atan(sqrt(dx² + dy²) / |Δz|)
```

For three- and four-chamber modes, the more physical angular estimate is
obtained from the fitted reference track:

``` text
x(z) = ax z + bx
y(z) = ay z + by
θx = atan(ax)
θy = atan(ay)
θ  = atan(sqrt(ax² + ay²))
```

These quantities are intended as **diagnostic validation variables** near
the working point. They help assess residual alignment, angular spread,
tracklet quality, and possible correlations between geometry, timing,
and reconstructed tracklets.

------------------------------------------------------------------------

## 4) Input TTree selection and BX studies

By default, the analysis reads the TTree named:

``` text
events
```

Some `data.root` files may also contain BX-specific TTrees, for example:

``` text
events_36BX
events_24BX
events_18BX
events_12BX
```

The pipeline and producer support selecting the input tree with:

``` bash
--tree events_36BX
```

Example for a 36BX beam analysis:

``` bash
python3 src/rpc_hv_pipeline.py \
  --chambers Ch201 Ch202 \
  --producer src/rpc_tracklet_efficiency.py \
  --producer-mode optionc \
  --summarizer src/rpc_make_eff_summary.py \
  --hvscript src/rpc_efficiency_vs_hv.py \
  --z 0 100 \
  --tree events_36BX \
  --anchor-index 0 \
  --target-index 1 \
  --hit-mode centroid \
  --require-ncluster-eq1 \
  --regime beam \
  --out-dir out_dynamic_regime_36BX \
  --fit
```

Use a separate output directory for each tree or BX selection to avoid
mixing results:

``` text
out_dynamic_regime/
out_dynamic_regime_36BX/
out_dynamic_regime_24BX/
```

------------------------------------------------------------------------

## 5) 3D Tracklet Visualization

The script:

    src/rpc_plot_tracklets_3d.py

provides a full 3D geometrical visualization of reconstructed tracklets across all chambers.

### 5.1 What the script does

For each event passing the matching criteria:

- Reconstructs a straight-line trajectory using the chamber hit coordinates and z-positions.
- Draws a publication-style 3D tracklet using a central arrow and a light continuation segment.
- Uses color coding based on inclination angle, so near-vertical and more inclined tracklets can be distinguished visually.
- Optionally displays intermediate chamber hits (for 3- or 4-chamber setups).
- Draws RPC chamber contours or semi-transparent planes at their corresponding z-planes.
- Constrains the displayed detector area so that unphysical out-of-range tracks can be rejected or hidden.

This allows:

- Direct visual validation of tracklet geometry
- Inspection of misalignment effects
- Verification of z-position configuration
- Identification of outlier or mismatched tracklets
- Presentation-ready detector visualization

### 5.2 Geometry Requirements

You must provide:

- ```--files```: one ROOT file per chamber (same HV point)
- ```--z```: one z-position (cm) per chamber
- ```--anchor-index```: index of the anchor chamber
- ```--target-index```: index of the target chamber

Also:
- Directional arrows = reconstructed tracklets
- Color coding = track inclination or chosen track category
- Optional intermediate points
- RPC chamber contours or detector planes drawn at each z-plane
- Customizable chamber rectangle via:

``` bash
python3 src/rpc_plot_tracklets_3d.py \
--files path/to/file1/HV6/data.root path/to/file2/HV6/data.root \
--z z1 z2 \
--tree events \
--tag HV6 \
--anchor-index 0 \
--target-index -1 \
--max-tracks 400 \
--show-intermediate \
--pdf \
--chamber-rect 0 54 0 120 \
--chamber-style k \
--chamber-lw 1.2
```

### Example (3 chambers)

``` bash
python3 src/rpc_plot_tracklets_3d.py \
  --files ChA/HV6/data.root ChB/HV6/data.root ChC/HV6/data.root \
  --z 0 50 100 \
  --tree events \
  --anchor-index 0 \
  --target-index 2 \
  --tol-x 7 --tol-y 16 --tol-t 3 \
  --max-tracks 300 \
  --pdf
```

### Example (3D tracklet reconstruction for 2 chambers in a stylized fashion)

```bash
python3 src/rpc_plot_tracklets_3d.py \
  --files /path/to/file1/data.root /path/to/file2/data.root \
  --z 0 10 \
  --tree events \
  --tag <working_point_label> \
  --anchor-index 0 \
  --target-index -1 \
  --max-tracks 1000 \
  --show-intermediate \
  --pdf \
  --chamber-rect 0 60 0 120 \
  --chamber-style k \
  --chamber-lw 1.2
```

------------------------------------------------------------------------

## 6) Running the HV Pipeline

The same pipeline supports 2-, 3-, and 4-chamber studies. The number of
input chamber directories must match the number of z positions given with
`--z`.

### 6.1 Fixed-cut mode

``` bash
python3 src/rpc_hv_pipeline.py \
  --chambers ChA ChB ChC \
  --producer src/rpc_tracklet_efficiency.py \
  --producer-mode optionc \
  --summarizer src/rpc_make_eff_summary.py \
  --hvscript src/rpc_efficiency_vs_hv.py \
  --z 0 50 100 \
  --anchor-index 0 \
  --target-index 2 \
  --tol-x 7.0 \
  --tol-y 20.0 \
  --tol-t 3.0 \
  --hit-mode centroid \
  --require-ncluster-eq1 \
  --out-dir out_hv_3ch \
  --fit
```

### 6.2 Dynamic-cut mode with regime selection

``` bash
python3 src/rpc_hv_pipeline.py \
  --chambers ChA ChB ChC \
  --producer src/rpc_tracklet_efficiency.py \
  --producer-mode optionc \
  --summarizer src/rpc_make_eff_summary.py \
  --hvscript src/rpc_efficiency_vs_hv.py \
  --z 0 50 100 \
  --anchor-index 0 \
  --target-index 2 \
  --hit-mode centroid \
  --require-ncluster-eq1 \
  --regime beam \
  --auto-cuts \
  --auto-cuts-script src/rpc_estimate_tracklet_cuts.py \
  --out-dir out_hv_3ch_beam \
  --fit
```

### 6.3 Using precomputed cut JSON

``` bash
python3 src/rpc_hv_pipeline.py \
  --chambers ChA ChB ChC \
  --producer src/rpc_tracklet_efficiency.py \
  --producer-mode optionc \
  --summarizer src/rpc_make_eff_summary.py \
  --hvscript src/rpc_efficiency_vs_hv.py \
  --z 0 50 100 \
  --anchor-index 0 \
  --target-index 2 \
  --regime normal \
  --cuts-json configs/default_cuts.json \
  --out-dir out_hv_3ch_normal \
  --fit
```

### 6.4 Running with diagnostic CSV/NDJSON output

For approval studies, WP validation, angle plots, and correlation
matrices, rerun the analysis with the diagnostic table switches:

``` bash
--write-csv --write-ndjson
```

These options are important because the final ROOT map file contains
histograms, while the **correlation matrix requires event-by-event or
tracklet-by-tracklet quantities**. The CSV file is therefore the main
input for the numerical correlation matrix.

Example for a 36BX beam scan:

``` bash
python3 src/rpc_hv_pipeline.py \
  --chambers Ch201 Ch202 \
  --producer src/rpc_tracklet_efficiency.py \
  --producer-mode optionc \
  --summarizer src/rpc_make_eff_summary.py \
  --hvscript src/rpc_efficiency_vs_hv.py \
  --z 0 100 \
  --tree events_36BX \
  --anchor-index 0 \
  --target-index 1 \
  --hit-mode centroid \
  --require-ncluster-eq1 \
  --regime beam \
  --out-dir out_dynamic_regime_36BX \
  --fit \
  --write-csv \
  --write-ndjson
```

During each HV point, the producer first writes local diagnostic files:

``` text
tracklets_full.csv
tracklets.ndjson
```

The HV pipeline then moves and renames them into the `maps/` directory:

``` text
out_dynamic_regime_36BX/maps/tracklets_full_HV6_beam.csv
out_dynamic_regime_36BX/maps/tracklets_HV6_beam.ndjson
```

The exact HV tag depends on the processed HV point. For example, `HV6`
may be replaced by `HV5`, `HV7`, etc.

The CSV is the recommended file for:

- checking event-by-event residuals,
- checking corrected and raw `Δx`, `Δy`, and `Δt`,
- plotting or recomputing `θx`, `θy`, and `θ`,
- creating the tracklet-level correlation matrix.

The NDJSON file is useful when one wants a more verbose event-level log
that can be inspected line by line.

### 6.5 Direct CSV generation for one HV point

If you do not want to rerun the full HV pipeline, you can run the
producer directly for a single HV point and ask it to write the CSV and
NDJSON diagnostics:

``` bash
python3 src/rpc_tracklet_efficiency.py \
  --files Ch201/HV6/data.root Ch202/HV6/data.root \
  --z 0 100 \
  --tree events_36BX \
  --tag HV6_beam \
  --anchor-index 0 \
  --target-index 1 \
  --hit-mode centroid \
  --require-ncluster-eq1 \
  --regime beam \
  --tol-x 7.0 \
  --tol-y 20.0 \
  --tol-t 3.0 \
  --write-csv \
  --write-ndjson
```

This direct command writes in the current directory:

``` text
tracklets_efficiency_HV6_beam.root
tracklets_full.csv
tracklets.ndjson
```

For repository-style bookkeeping, move or rename them consistently, for
example:

``` bash
mkdir -p out_dynamic_regime_36BX/maps
mv tracklets_efficiency_HV6_beam.root \
   out_dynamic_regime_36BX/maps/tracklets_efficiency_HV6_beam.root
mv tracklets_full.csv \
   out_dynamic_regime_36BX/maps/tracklets_full_HV6_beam.csv
mv tracklets.ndjson \
   out_dynamic_regime_36BX/maps/tracklets_HV6_beam.ndjson
```

To check that the CSV exists and contains the expected columns:

``` bash
head -n 1 out_dynamic_regime_36BX/maps/tracklets_full_HV6_beam.csv
```

Typical columns include:

``` text
xA, yA, x_pred, y_pred, xT, yT,
dx_raw, dy_raw, dx, dy, dt,
dx0, dy0,
theta_x_deg, theta_y_deg, theta_deg
```

------------------------------------------------------------------------

## 7) Automatic Cut Estimation

The script:

    src/rpc_estimate_tracklet_cuts.py

can be used to determine reasonable cut values before the full HV scan is run.

Typical strategy:

1.  build a loose preselected sample
2.  inspect `Δx`, `Δy`, `Δt` distributions
3.  estimate robust centers and widths
4.  write per-HV values into a JSON file
5.  pass those values to the producer during the main pipeline run

This helps avoid using a single fixed set of cuts when:

- cluster size evolves with HV
- timing broadens at high working points
- beam conditions differ from the natural background regime
- one wants more reproducible and physically motivated selections

------------------------------------------------------------------------

## 8) Beam vs Normal Regimes

The repository supports two analysis regimes:

-   **normal**: intended for natural background / broad detector acceptance studies
-   **beam**: intended for localized beam-region studies and ROI-centered analyses

In practice, this means:

- Output names include the regime label
- Cut estimation can be performed separately for each regime
- Summaries can be produced consistently without overwriting one another
- Beam-region studies can be kept distinct from general detector-efficiency scans

------------------------------------------------------------------------

## 9) Diagnostic ROOT histograms and per-tracklet tables

The updated producer keeps the standard map output unchanged:

``` text
h_total
h_tracklet
h_eff
```

It also writes diagnostic ROOT histograms useful for WP validation:

``` text
h_dx_raw_candidate
h_dy_raw_candidate
h_dx_raw_dy_raw_candidate
h_dx_candidate
h_dy_candidate
h_dt_candidate
h_dx_dy_candidate
h_dx
h_dy
h_dt
h_dx_dy
h_theta_x
h_theta_y
h_theta
h_theta_x_theta_y
```

The matched-tracklet histograms are filled after the matching cuts are
applied. The candidate histograms are filled before the final matching
condition, and are useful for inspecting the residual sample before
accepting the final matched tracklets.

The diagnostic ROOT histograms are sufficient for approval-style plots
such as:

``` text
Δx distribution
Δy distribution
Δt distribution
θx distribution
θy distribution
θ distribution
θx vs θy map
Δx vs Δy map
```

However, the ROOT histograms alone are **not sufficient** for a full
tracklet-level correlation matrix, because a correlation matrix requires
event-by-event values in the same row. For that reason, the analysis
must also write the per-tracklet CSV.

When the pipeline is run with:

``` bash
--write-csv --write-ndjson
```

the per-tracklet CSV can include quantities such as:

``` text
xA, yA, x_pred, y_pred, xT, yT,
dx_raw, dy_raw, dx, dy, dt,
dx0, dy0,
theta_x_deg, theta_y_deg, theta_deg
```

These columns are the preferred input for the numerical
tracklet-correlation matrix. Each row corresponds to one matched
tracklet, so the correlations are computed between variables from the
same reconstructed object.

To check that the diagnostic histograms are present:

``` bash
root -l out_dynamic_regime_36BX/maps/tracklets_efficiency_HV6_beam.root
```

Inside ROOT:

``` cpp
.ls
```

You should see at least:

``` text
h_dx
h_dy
h_theta_x
h_theta_y
h_theta_x_theta_y
```

To check that the diagnostic CSV was created:

``` bash
ls out_dynamic_regime_36BX/maps/tracklets_full_HV*_beam.csv
head -n 1 out_dynamic_regime_36BX/maps/tracklets_full_HV6_beam.csv
```

A quick check of the number of matched tracklets stored in the CSV is:

``` bash
python3 - <<'PY'
import pandas as pd
csv = "out_dynamic_regime_36BX/maps/tracklets_full_HV6_beam.csv"
df = pd.read_csv(csv)
print("CSV:", csv)
print("rows / matched tracklets:", len(df))
print("columns:")
print("  " + "\n  ".join(df.columns))
PY
```

------------------------------------------------------------------------

## 10) Creating relevant plots

### 10.1 CMS-like efficiency curve and fitted WP

After the HV pipeline has produced the JSON summaries, the CMS-like
curve can be created with:

``` bash
mkdir -p out_dynamic_regime_36BX/plots

python3 src/rpc_cms_style_eff_plot.py \
  --series "RxLW = 12 BX::out_dynamic_regime_36BX/summaries/summary_HV*beam*.json" \
  --textbox $'chamber 202 36BX\n1.4mm double gap iRPC\nFEB v2.3 Petiroc 2C\nthreshold ~ 40fC (dac10)' \
  --out out_dynamic_regime_36BX/plots/eff_cms_style_36BX_beam.png \
  --wp-offset-v 150 \
  --textbox-x 0.02 \
  --textbox-y 0.52
```

Use `--textbox-x` and `--textbox-y` to move the annotation block if it
overlaps the efficiency curve.

### 10.2 Map plots from ROOT histograms

The `rpc_plot_efficiency_maps.py` reads the standard ROOT map histograms
`h_total`, `h_tracklet`, and `h_eff`. It reproduces the three map-style
figures:

- Total hits
- Reconstructed tracklets
- Local Efficiency

Use it from the repository root:

``` bash
python3 src/rpc_plot_efficiency_maps.py \
  --input-dir out_dynamic_regime/maps \
  --plot-range 0 60 0 160 \
  --png
```

It writes outputs by default to:

``` text
out_dynamic_regime/map_plots/
```

For example:

``` text
out_dynamic_regime/map_plots/HV6_beam/
  Total_7kV_HV6_beam.pdf
  Tracklets_7kV_HV6_beam.pdf
  Local_Efficiency_7kV_HV6_beam.pdf
```

To keep only PDF output, omit `--png`:

``` bash
python3 src/rpc_plot_efficiency_maps.py \
  --input-dir out_dynamic_regime/maps \
  --plot-range 0 60 0 160
```

The plotting script also supports CMS/GIF++-style annotation. Example:

``` bash
python3 src/rpc_plot_efficiency_maps.py \
  --input-dir out_dynamic_regime/maps \
  --plot-range 0 60 0 160 \
  --cms-annotation "GIF++ Test Beam November 2025\n1.4 mm double gap iRPC\nFEB v2.3 Petiroc 2C\nthreshold ~ 40 fC\n\nGIF++ source off" \
  --png
```

### 10.3 Working-point validation plots

For approval packages, choose the HV point closest to the fitted working
point and generate the validation plots from the diagnostic ROOT file and
the per-tracklet CSV.

Example for HV6 beam:

``` bash
mkdir -p out_dynamic_regime_36BX/plots/wp_validation/HV6_beam

python3 src/rpc_wp_validation_plots.py \
  --root out_dynamic_regime_36BX/maps/tracklets_efficiency_HV6_beam.root \
  --csv out_dynamic_regime_36BX/maps/tracklets_full_HV6_beam.csv \
  --out-dir out_dynamic_regime_36BX/plots/wp_validation/HV6_beam \
  --label $'Chamber 202, 36BX\n1.4 mm double gap iRPC\nFEB v2.3 Petiroc 2C\nthreshold ~ 40 fC' \
  --png \
  --pdf
```

This produces:

``` text
03_delta_x_distribution.png/pdf
04_delta_y_distribution.png/pdf
05_theta_x_distribution.png/pdf
06_theta_y_distribution.png/pdf
07_theta_x_vs_theta_y_2d.png/pdf
08_tracklet_correlation_matrix.png/pdf
08_tracklet_correlation_matrix.csv
```

The first five plots are produced from the diagnostic ROOT histograms.
The correlation matrix is produced from the per-tracklet CSV supplied
with `--csv`.

These plots directly address the usual WP validation request:

``` text
Fit obtain WP, show beam profile for HV close to WP,
and show the correlation matrix and angles.
```

### 10.4 Creating the tracklet-level correlation matrix

The tracklet-level correlation matrix is created from the CSV file, not
from the final map histograms. The reason is that the matrix requires the
variables to be available row-by-row for each matched tracklet.

For the HV point closest to the fitted WP, use the CSV from the same HV
and regime as the ROOT diagnostic file:

``` text
ROOT: out_dynamic_regime_36BX/maps/tracklets_efficiency_HV6_beam.root
CSV : out_dynamic_regime_36BX/maps/tracklets_full_HV6_beam.csv
```

The simplest command is the full WP-validation command above, because it
creates both the angle plots and the correlation matrix:

``` bash
python3 src/rpc_wp_validation_plots.py \
  --root out_dynamic_regime_36BX/maps/tracklets_efficiency_HV6_beam.root \
  --csv out_dynamic_regime_36BX/maps/tracklets_full_HV6_beam.csv \
  --out-dir out_dynamic_regime_36BX/plots/wp_validation/HV6_beam \
  --label $'Chamber 202, 36BX\n1.4 mm double gap iRPC\nFEB v2.3 Petiroc 2C\nthreshold ~ 40 fC' \
  --png \
  --pdf
```

The matrix files are:

``` text
out_dynamic_regime_36BX/plots/wp_validation/HV6_beam/08_tracklet_correlation_matrix.png
out_dynamic_regime_36BX/plots/wp_validation/HV6_beam/08_tracklet_correlation_matrix.pdf
out_dynamic_regime_36BX/plots/wp_validation/HV6_beam/08_tracklet_correlation_matrix.csv
```

The plotted correlation coefficient is the Pearson coefficient. Typical
columns used are:

``` text
xA, yA, x_pred, y_pred, xT, yT,
dx_raw, dy_raw, dx, dy, dt,
theta_x_deg, theta_y_deg, theta_deg
```

Interpretation guide:

- strong `dx`--`xA` or `dy`--`yA` correlations can indicate residual
  geometry or alignment effects;
- strong `dt`--position correlations can indicate timing or readout
  effects across the chamber plane;
- `theta_x_deg` and `theta_y_deg` describe the angular spread of the
  accepted tracklets near the WP;
- a stable validation sample should not show unexpected large
  correlations between residuals, timing, and position variables.

To inspect the numerical matrix directly:

``` bash
column -s, -t < \
  out_dynamic_regime_36BX/plots/wp_validation/HV6_beam/08_tracklet_correlation_matrix.csv \
  | less -S
```

### 10.5 Suggested approval-plot order

A compact approval package can contain:

1. CMS-like efficiency vs HV curve with the fitted WP.
2. Beam profile or local efficiency map at the HV closest to WP.
3. `Δx` distribution at the HV closest to WP.
4. `Δy` distribution at the HV closest to WP.
5. `θx` distribution.
6. `θy` distribution.
7. `θx` vs `θy` 2D map.
8. Tracklet-level correlation matrix.

------------------------------------------------------------------------

## 11) Synthetic Dataset Generator

``` bash
python3 synthetic/make_synthetic_rpc_hvscan.py \
  --out-base synthetic_calib_3ch \
  --n-chambers 3 \
  --n-events 60000 \
  --z 0 50 100
```

------------------------------------------------------------------------

## 12) Troubleshooting

-   Ensure PyROOT is available.
-   Verify TTree name (`--tree`, default `events`).
-   For BX-specific studies, check that the selected tree exists, for
    example `events_36BX`.
-   Confirm correct `HV*/data.root` structure.
-   Check geometry (`--z`) and tolerances.
-   If beam and normal outputs look unexpectedly similar, verify whether
    the producer stage is truly using regime-aware selections or only
    regime-aware labels.
-   If dynamic cuts appear fixed, check whether `--auto-cuts`,
    `--auto-cuts-script`, or `--cuts-json` are being passed correctly by
    the pipeline.
-   If the WP validation script cannot find `h_dx`, `h_dy`,
    `h_theta_x`, or `h_theta_y`, rerun the analysis with the updated
    producer so that the diagnostic ROOT histograms are written.
-   If the correlation matrix cannot be produced, rerun the pipeline with
    `--write-csv`. The final map ROOT files alone do not contain enough
    event-by-event information for a full tracklet-level correlation
    matrix.
-   If a plot command fails because an output directory does not exist,
    create it with `mkdir -p` or use a script option that creates output
    directories automatically.

------------------------------------------------------------------------

## Contact

Klaudio Peqini --- University of Tirana, CMS-Albania UT-team\
Ilirjan Margjeka --- Luarasi University, CMS-Albania UT-team Associate
