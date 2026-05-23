#!/usr/bin/env python

import argparse
import gzip
import numpy as np
import cooler

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cool", required=True, help="e.g. sample.cool or sample.mcool::/resolutions/50000")
    ap.add_argument("--resolution", type=int, required=True)
    ap.add_argument("--lower", type=int, default=None)
    ap.add_argument("--upper", type=int, default=1000000)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--weight-name", default="weight")
    ap.add_argument("--bias-mode", choices=["inverse_weight", "coverage"], default="inverse_weight")
    args = ap.parse_args()

    lower = args.lower if args.lower is not None else args.resolution

    clr = cooler.Cooler(args.cool)
    bins = clr.bins()[:]
    pixels = clr.pixels()[:]

    chrom = bins["chrom"].astype(str).to_numpy()
    start = bins["start"].to_numpy().astype(int)
    end = bins["end"].to_numpy().astype(int)
    mid = ((start + end) // 2).astype(int)

    # Raw marginal coverage from cooler pixels.
    marginal = np.zeros(len(bins), dtype=float)
    b1 = pixels["bin1_id"].to_numpy()
    b2 = pixels["bin2_id"].to_numpy()
    cnt = pixels["count"].to_numpy().astype(float)

    np.add.at(marginal, b1, cnt)
    np.add.at(marginal, b2, cnt)
    same = b1 == b2
    if same.any():
        # Diagonal pixels were added twice above, but should contribute once to marginal.
        np.add.at(marginal, b1[same], -cnt[same])

    # 1. contactCounts: chr1 mid1 chr2 mid2 count
    contact_out = args.out_prefix + ".contactCounts.gz"
    with gzip.open(contact_out, "wt") as f:
        for i, j, c in zip(b1, b2, cnt):
            if chrom[i] != chrom[j]:
                continue
            d = abs(mid[j] - mid[i])
            if d < lower or d > args.upper:
                continue
            if c <= 0:
                continue
            if mid[i] <= mid[j]:
                f.write(f"{chrom[i]}\t{mid[i]}\t{chrom[j]}\t{mid[j]}\t{int(c)}\n")
            else:
                f.write(f"{chrom[j]}\t{mid[j]}\t{chrom[i]}\t{mid[i]}\t{int(c)}\n")

    # 2. fragmentsfile: chr start mid hitCount mappability
    # FitHiC mostly needs chromosome and midpoint; hitCount is useful metadata.
    frag_out = args.out_prefix + ".fragments.gz"
    with gzip.open(frag_out, "wt") as f:
        for c, s, m, cov in zip(chrom, start, mid, marginal):
            f.write(f"{c}\t{s}\t{m}\t{int(cov)}\t1\n")

    # 3. biasValues: chr mid bias
    bias_out = args.out_prefix + ".biasValues.gz"

    if args.bias_mode == "inverse_weight" and args.weight_name in bins.columns:
        w = bins[args.weight_name].astype(float).to_numpy()
        bias = np.full(len(bins), np.nan, dtype=float)
        ok = np.isfinite(w) & (w > 0)
        bias[ok] = 1.0 / w[ok]
    else:
        bias = marginal.astype(float)

    ok = np.isfinite(bias) & (bias > 0)

    # Remove extreme outliers before mean scaling.
    # This avoids one pathological bin dominating the whole bias vector.
    if ok.sum() > 0:
        vals = bias[ok]
        lo, hi = np.quantile(vals, [0.001, 0.999])
        ok = ok & (bias >= lo) & (bias <= hi)

    if ok.sum() == 0:
        raise RuntimeError("No valid bias values after filtering.")

    bias[ok] = bias[ok] / bias[ok].mean()

    with gzip.open(bias_out, "wt") as f:
        for c, m, b, is_ok in zip(chrom, mid, bias, ok):
            if is_ok:
                f.write(f"{c}\t{m}\t{b:.10g}\n")
            else:
                f.write(f"{c}\t{m}\t-1\n")

    print("Wrote:")
    print(contact_out)
    print(frag_out)
    print(bias_out)

if __name__ == "__main__":
    main()
