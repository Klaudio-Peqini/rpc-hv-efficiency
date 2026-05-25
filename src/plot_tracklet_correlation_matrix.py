#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True)
    ap.add_argument("--out", default="tracklet_correlation_matrix.png")
    args = ap.parse_args()

    df = pd.read_csv(args.csv)

    cols = [
        "xA", "yA",
        "x_pred", "y_pred",
        "xT", "yT",
        "dx_raw", "dy_raw",
        "dx", "dy", "dt",
        "theta_x_deg", "theta_y_deg", "theta_deg",
    ]

    cols = [c for c in cols if c in df.columns]

    if len(cols) < 2:
        raise SystemExit(
            f"Not enough usable columns found. Available columns are:\n{list(df.columns)}"
        )

    corr = df[cols].corr(method="pearson")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(corr.values, vmin=-1, vmax=1, cmap="coolwarm")

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Pearson correlation coefficient")

    ax.set_xticks(range(len(cols)))
    ax.set_yticks(range(len(cols)))
    ax.set_xticklabels(cols, rotation=45, ha="right")
    ax.set_yticklabels(cols)

    for i in range(len(cols)):
        for j in range(len(cols)):
            ax.text(
                j, i,
                f"{corr.values[i, j]:.2f}",
                ha="center",
                va="center",
                fontsize=8,
            )

    ax.set_title("Tracklet-level correlation matrix near WP")
    fig.tight_layout()
    fig.savefig(out, dpi=250)
    corr.to_csv(out.with_suffix(".csv"))

    print(f"Wrote {out}")
    print(f"Wrote {out.with_suffix('.csv')}")


if __name__ == "__main__":
    main()

