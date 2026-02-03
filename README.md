# rpc-hv-efficiency
**Tracklet-based RPC efficiency vs High Voltage (HV) analysis**  
*(CMS/RPC-style workflow for laboratory or QC studies)*

---

## 0) What this repository does

This repository provides a **reproducible pipeline** to measure RPC efficiency as a function of HV using
**tracklet reconstruction**. It supports:

- **2-chamber** (A → B) efficiency (reference → target)
- **3-chamber** (A+B → C) efficiency (two references define a tracklet; third is target)
- **4-chamber** (A+B+C → D) efficiency (three references; fourth is target)

It produces, for each HV point:

1. **2D efficiency maps** (ROOT histograms):
   - `h_total` (denominator: eligible events)
   - `h_tracklet` (numerator: matched events)
   - `h_eff` (efficiency map)
2. **JSON summaries** with ROI and global efficiencies + binomial errors
3. **Efficiency vs HV plots** (simple and CMS/RPC-group style)


** Download**

To download the repository, type in the working directory

```bash
git clone https://github.com/Klaudio-Peqini/rpc-hv-efficiency
---

## 1) Repository layout

```
rpc-hv-efficiency/
├─ src/
│  ├─ rpc_tracklets_efficiency.py      # main producer (Option C, n-ch)
│  ├─ rpc_hv_pipeline.py               # HV-scan orchestrator (2/3/4 chambers)
│  ├─ rpc_make_summary.py              # ROOT -> JSON summaries (ROI + global)
│  ├─ rpc_plot_eff_vs_hv.py            # quick efficiency-vs-HV plot (debug/simple)
│  └─ rpc_plot_eff_vs_hv_cms.py        # CMS/RPC-group style plot (fit, plateau, WP)
├─ synthetic/
│  └─ make_synthetic_rpc_hvscan.py     # synthetic dataset generator (calibration/regression)
├─ README.md
├─ requirements.txt
├─ .gitignore
├─ CHANGELOG.md
└─ CITATION.cff
```

**What you run in practice:**
- Usually you run **`src/rpc_hv_pipeline.py`** (or the producer manually for one HV point).
- The pipeline automatically runs producer → summary → plot for all HV points (HV1…HV11).

---

## 2) Requirements and environment

### 2.1 Mandatory
- **Python ≥ 3.9**
- **ROOT with PyROOT enabled**  
  This must work:
```bash
python3 -c "import ROOT; print(ROOT.gROOT.GetVersion())"
```

### 2.2 Python packages (plotting + fitting)
Install:
```bash
python3 -m pip install -r requirements.txt
```

> Note: In most CERN/CMS setups, **PyROOT is not installed with pip**.  
> Use an environment where ROOT provides PyROOT (CMSSW, LCG, system ROOT, etc.).

---

## 3) Input data layout (required)

Each chamber must be arranged as:

```
ChamberXXX/
  HV1/data.root
  HV2/data.root
  ...
  HV11/data.root
```

### 3.1 HV mapping used by the pipeline
The pipeline assumes the following mapping (kV):

- HV1 → 6.0  
- HV2 → 6.5  
- HV3 → 6.7  
- HV4 → 6.8  
- HV5 → 6.9  
- HV6 → 7.0  
- HV7 → 7.1  
- HV8 → 7.2  
- HV9 → 7.3  
- HV10 → 7.4  
- HV11 → 7.5  

If your run uses different HV values, update the mapping inside the pipeline (the repo already supports this in code).

---

## 4) Physics definition (what “efficiency” means here)

We define:

**ε = N_matched / N_eligible**

- **Eligible events (denominator)** are defined *only* using the reference chamber(s).
- The **target/test chamber never contributes to the denominator**.
- **Matched events (numerator)** require a compatible hit in the target chamber, within:
  - spatial tolerances (`tol-x`, `tol-y`) and
  - time tolerance (`tol-t`)

This is the key condition that produces a physically meaningful sigmoid turn-on and stable plateau.

---

## 5) What the scripts produce

Assume an output directory `OUT_DIR` (e.g. `out_hv_3ch/`). You will get:

```
OUT_DIR/
  maps/
    tracklets_efficiency_HV1.root
    ...
    tracklets_efficiency_HV11.root

  summaries/
    summary_HV1.json
    ...
    summary_HV11.json

  plots/
    eff_scint_overlap_vs_hv.png
    eff_roi_cms_style.png          (if you run the CMS-style plotter)
    eff_global_cms_style.png       (optional)
```

### 5.1 ROOT map file content
Each `tracklets_efficiency_HVx.root` contains:
- `h_total`     : 2D denominator counts
- `h_tracklet`  : 2D numerator counts
- `h_eff`       : 2D efficiency map (numerator/denominator)

### 5.2 JSON summary content
Each `summary_HVx.json` contains:
- `hv_kv`
- ROI definition (typically scintillator overlap)
- ROI counts and ROI efficiency + binomial error
- Global counts and global efficiency + binomial error
- Timing metadata

---

## 6) Key parameters you must understand

### 6.1 Geometry: z positions
For n chambers you must provide one z-position per chamber (cm):

- **3 chambers:** `--z zA zB zC`
- **4 chambers:** `--z zA zB zC zD`

These are used to extrapolate tracklets between planes.

### 6.2 Chamber roles
- `--anchor-index`: which chamber is used as the anchor reference for defining the tracklet
- `--target-index`: which chamber is the *efficiency target* (the one tested)

In typical usage:
- For 3 chambers: `anchor-index=0`, `target-index=2`
- For 4 chambers: `anchor-index=0`, `target-index=3`

### 6.3 Matching tolerances
- `--tol-x` (cm): allowed |Δx| between predicted and observed hit
- `--tol-y` (cm): allowed |Δy|
- `--tol-t` (ns): allowed |Δt| from `Strip_Tdc_diff`

These strongly affect measured efficiency, especially near threshold HV.

### 6.4 Hit definition
- `--hit-mode centroid` (recommended): uses the centroid of the cluster points
- other modes may exist depending on your producer script (centroid is the stable choice)

### 6.5 Cluster selection
- `--require-ncluster-eq1`: enforce `nCluster == 1`
  - This matches typical RPC clean-hit selection logic.

---

## 7) Running the pipeline: 3-chamber analysis (BASH commands)

### 7.1 Example (A + B → C)
You have three chamber folders:

- `/path/to/ChA`
- `/path/to/ChB`
- `/path/to/ChC`

Each must contain `HV1..HV11/data.root`.

Choose z positions in cm, e.g.:
- A at 0 cm
- B at 50 cm
- C at 100 cm

Run:

```bash
python3 /path/to/src/rpc_hv_pipeline.py   --chambers /path/to/ChA /path/to/ChB /path/to/ChC   --producer /path/to/src/rpc_tracklets_efficiency.py   --producer-mode optionc   --summarizer /path/to/src/rpc_make_summary.py   --hvscript /path/to/src/rpc_plot_eff_vs_hv_cms.py   --z 0 50 100   --anchor-index 0   --target-index 2   --tol-x 7.0   --tol-y 20.0   --tol-t 3.0   --hit-mode centroid   --require-ncluster-eq1   --out-dir out_hv_3ch   --fit
```

### 7.2 Plot CMS/RPC-style ROI efficiency vs HV
```bash
python3 /path/to/src/rpc_plot_eff_vs_hv_cms.py   --use-roi   --series "ROI::out_hv_3ch/summaries/summary_HV*.json"   --out out_hv_3ch/plots/eff_roi_cms_style.pdf #.png
```

### 7.3 Plot CMS/RPC-style global efficiency vs HV (optional)
```bash
python3 /path/to/src/rpc_plot_eff_vs_hv_cms.py   --series "Global::out_hv_3ch/summaries/summary_HV*.json"   --out out_hv_3ch/plots/eff_global_cms_style.pdf #.png
```

---

## 8) Running the pipeline: 4-chamber analysis (BASH commands)

### 8.1 Example (A + B + C → D)
You have four chamber folders:

- `/path/to/ChA`
- `/path/to/ChB`
- `/path/to/ChC`
- `/path/to/ChD`

Choose z positions (cm), e.g.:
- 0, 40, 80, 120

Run:

```bash
python3 /path/to/src/rpc_hv_pipeline.py   --chambers /path/to/ChA /path/to/ChB /path/to/ChC /path/to/ChD   --producer /path/to/src/rpc_tracklets_efficiency.py   --producer-mode optionc   --summarizer /path/to/src/rpc_make_summary.py   --hvscript /path/to/src/rpc_plot_eff_vs_hv_cms.py   --z 0 40 80 120   --anchor-index 0   --target-index 3   --tol-x 7.0   --tol-y 20.0   --tol-t 3.0   --hit-mode centroid   --require-ncluster-eq1   --out-dir out_hv_4ch   --fit
```

### 8.2 Plot CMS/RPC-style ROI efficiency vs HV
```bash
python3 /path/to/src/rpc_plot_eff_vs_hv_cms.py   --use-roi   --series "ROI::out_hv_4ch/summaries/summary_HV*.json"   --out out_hv_4ch/plots/eff_roi_cms_style.pdf #.png
```

### 8.3 Plot CMS/RPC-style global efficiency vs HV (optional)
```bash
python3 /path/to/src/rpc_plot_eff_vs_hv_cms.py   --series "Global::out_hv_4ch/summaries/summary_HV*.json"   --out out_hv_4ch/plots/eff_global_cms_style.pdf #.png
```

---

## 9) Synthetic calibration dataset (recommended for regression)

The `/path/to/synthetic/make_synthetic_rpc_hvscan.py` generator creates ROOT files with:

- the same branch schema as real data
- controlled eligibility and matching
- enforced target efficiency points + binomial errors

### 9.1 Generate a 3-ch synthetic dataset
```bash
python3 synthetic/make_synthetic_rpc_hvscan.py   --out-base synthetic_calib_3ch   --n-chambers 3   --n-events 60000   --z 0 50 100
```

### 9.2 Analyze it with the real pipeline
```bash
python3 src/rpc_hv_pipeline.py   --chambers synthetic_calib_3ch/ChamberSim1 synthetic_calib_3ch/ChamberSim2 synthetic_calib_3ch/ChamberSim3   --producer src/rpc_tracklets_efficiency.py   --producer-mode optionc   --z 0 50 100   --anchor-index 0   --target-index 2   --tol-x 7.0 --tol-y 20.0 --tol-t 3.0   --hit-mode centroid   --require-ncluster-eq1   --out-dir out_syn_3ch   --fit
```

---

## 10) Troubleshooting checklist (most common issues)

### 10.1 PyROOT not available
Symptom: `ModuleNotFoundError: No module named ROOT`  
Fix: switch to a ROOT/CMSSW environment where PyROOT is provided.

### 10.2 Wrong TTree name
Symptom: “tree not found” or zero entries  
Fix: use the correct `--tree` option in the pipeline/producer (default `events`).

### 10.3 No outputs / empty folders
Common causes:
- input directory does not contain `HV*/data.root`
- wrong chamber path (points to parent directory)
- wrong permissions

### 10.4 Efficiency too high / too low
Check:
- `tol-x`, `tol-y`, `tol-t`
- z geometry (`--z`)
- whether `--require-ncluster-eq1` is appropriate for your dataset

---

## License / Contact
Update these fields as needed for your institution/group policy.
Klaudio Peqini: University of Tirana, CMS-Albania UT-team
Ilirjan Margjeka: University of Tirana, associated with, CMS-Albania UT-team
