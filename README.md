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

------------------------------------------------------------------------

# Download

``` bash
git clone https://github.com/Klaudio-Peqini/rpc-hv-efficiency
```

------------------------------------------------------------------------

## 1) Repository layout

    rpc-hv-efficiency/
    ├─ src/
    │  ├─ rpc_tracklets_efficiency.py
    │  ├─ rpc_hv_pipeline.py
    │  ├─ rpc_make_summary.py
    │  ├─ rpc_plot_eff_vs_hv.py
    │  ├─ rpc_plot_eff_vs_hv_cms.py
    │  └─ rpc_plot_tracklets_3d.py
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

------------------------------------------------------------------------

## 4) 3D Tracklet Visualization

The script:

    src/rpc_plot_tracklets_3d.py

provides a full 3D geometrical visualization of reconstructed tracklets across all chambers.

### 4.1 What the script does

For each event passing the matching criteria:

- Reconstructs a straight-line trajectory using the chamber hit coordinates and z-positions.
- Draws a thin black line representing the particle trajectory through space.
- Highlights the anchor chamber hit and the target chamber hit.
- Uses a consistent color for both endpoints of the same tracklet.
- Optionally displays intermediate chamber hits (for 3- or 4-chamber setups).
- Draws RPC chamber contours at their corresponding z-planes.

This allows:

- Direct visual validation of tracklet geometry
- Inspection of misalignment effects
- Verification of z-position configuration
- Presentation-ready detector visualization

### 4.2 Geometry Requirements

You must provide:

- ```--files```: one ROOT file per chamber (same HV point)
- ```--z```: one z-position (cm) per chamber
- ```--anchor-index```: index of the anchor chamber
- ```--target-index```: index of the target chamber

Also: 
- Thin black lines = reconstructed tracklets
- Colored anchor & target endpoints (same color per tracklet)
- Optional intermediate points
- RPC chamber contours drawn at each z-plane
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

------------------------------------------------------------------------

## 5) Running 3-Chamber Pipeline

``` bash
python3 src/rpc_hv_pipeline.py \
  --chambers ChA ChB ChC \
  --producer src/rpc_tracklets_efficiency.py \
  --producer-mode optionc \
  --summarizer src/rpc_make_summary.py \
  --hvscript src/rpc_plot_eff_vs_hv_cms.py \
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

------------------------------------------------------------------------

## 6) Synthetic Dataset Generator

``` bash
python3 synthetic/make_synthetic_rpc_hvscan.py \
  --out-base synthetic_calib_3ch \
  --n-chambers 3 \
  --n-events 60000 \
  --z 0 50 100
```

------------------------------------------------------------------------

## 7) Troubleshooting

-   Ensure PyROOT is available.
-   Verify TTree name (`--tree`, default `events`).
-   Confirm correct `HV*/data.root` structure.
-   Check geometry (`--z`) and tolerances.

------------------------------------------------------------------------

## Contact

Klaudio Peqini --- University of Tirana, CMS-Albania UT-team\
Ilirjan Margjeka --- University of Tirana, CMS-Albania UT-team
