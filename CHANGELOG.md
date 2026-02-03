# Changelog

All notable changes to this project will be documented in this file.

## [v1.0] - 2026-02-03
### Added
- Tracklet-based efficiency reconstruction (Option C) for 2/3/4 chambers.
- Automated HV scan pipeline over HV1..HV11 folder layout.
- ROOT→JSON summary extraction (ROI and global efficiencies).
- CMS/RPC-group style efficiency vs HV plotting with sigmoid fit and WP extraction.
- Synthetic calibration dataset generator with enforced efficiency points and binomial errors.

### Physics / Validation
- Denominator defined strictly from reference chambers (test chamber never enters denominator).
- Verified sigmoid turn-on and stable plateau behaviour.
- Synthetic regression dataset added to calibrate and validate future changes.

### Notes
- PyROOT/ROOT is required for production and histogram handling.
- Synthetic ROOT files are not tracked in git (only the generator script is tracked).

