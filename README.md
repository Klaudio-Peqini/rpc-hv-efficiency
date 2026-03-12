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
2.  **JSON summaries** with ROI and global efficiencies + binomial
    errors
3.  **Efficiency vs HV plots** (simple and CMS/RPC-group style)
4.  **3D tracklet visualizations** showing reconstructed trajectories
    through all chambers
5.  **Automatic cut estimation** for `Δx`, `Δy`, `Δt` and optional
    angular consistency, evaluated per working point when desired
6.  **Regime-aware analysis** for `normal` and `beam` running
    conditions

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
    │  └─ rpc_plot_tracklets_3d.py
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

### Install dependencies

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

------------------------------------------------------------------------

## 4) 3D Tracklet Visualization

The script:

    src/rpc_plot_tracklets_3d.py

provides a full 3D geometrical visualization of reconstructed tracklets across all chambers.

### 4.1 What the script does

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

### 4.2 Geometry Requirements

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

## 5) Running 3-Chamber Pipeline

### Fixed-cut mode

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

### Dynamic-cut mode with regime selection

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

### Using precomputed cut JSON

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

------------------------------------------------------------------------

## 6) Automatic Cut Estimation

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

## 7) Beam vs Normal Regimes

The repository now supports two analysis regimes:

-   **normal**: intended for natural background / broad detector acceptance studies
-   **beam**: intended for localized beam-region studies and ROI-centered analyses

In practice this means:

- output names include the regime label
- cut estimation can be performed separately for each regime
- summaries can be produced consistently without overwriting one another
- beam-region studies can be kept distinct from general detector-efficiency scans

------------------------------------------------------------------------

## 8) Synthetic Dataset Generator

``` bash
python3 synthetic/make_synthetic_rpc_hvscan.py \
  --out-base synthetic_calib_3ch \
  --n-chambers 3 \
  --n-events 60000 \
  --z 0 50 100
```

------------------------------------------------------------------------

## 9) Troubleshooting

-   Ensure PyROOT is available.
-   Verify TTree name (`--tree`, default `events`).
-   Confirm correct `HV*/data.root` structure.
-   Check geometry (`--z`) and tolerances.
-   If beam and normal outputs look unexpectedly similar, verify whether the producer stage is truly using regime-aware selections or only regime-aware labels.
-   If dynamic cuts appear fixed, check whether `--auto-cuts`, `--auto-cuts-script`, or `--cuts-json` are being passed correctly by the pipeline.

------------------------------------------------------------------------

## Contact

Klaudio Peqini --- University of Tirana, CMS-Albania UT-team\
Ilirjan Margjeka --- Luarasi University, CMS-Albania UT-team Associate
