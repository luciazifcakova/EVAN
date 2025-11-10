#!/usr/bin/env python3
import argparse, os, csv, math
from collections import OrderedDict
from glob import glob

# BioPython
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

# Plotting
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from matplotlib import font_manager
import numpy as np

# --------------------------- helpers ---------------------------

def create_output_directory(directory, species, window):
    outdir = os.path.join(directory, f'{species}_output_{window}')
    os.makedirs(outdir, exist_ok=True)
    return outdir

def main_dir(fasta_file, species, window):
    outdir = create_output_directory(os.path.dirname(fasta_file), species, window)
    print("Output Directory:", outdir)
    print("Fasta File:", fasta_file)
    print("Species:", species)
    print("Window:", window)
    return outdir

def how_many(seq, alphabet):
    return sum(seq.count(a) for a in alphabet)

def valid_window(win, min_nonN_frac):
    L = len(win)
    if L == 0:
        return False
    nonN = sum(b in 'ACGTacgt' for b in win)
    return (nonN / L) >= min_nonN_frac

def save_grid(output_directory, filename, image_paths, columns=2):
    """Make a simple montage from existing PNGs."""
    if not image_paths:
        return
    rows = (len(image_paths) + columns - 1) // columns
    fig, axes = plt.subplots(rows, columns, figsize=(12*columns, 3*rows), squeeze=False)
    k = 0
    for r in range(rows):
        for c in range(columns):
            ax = axes[r][c]
            if k < len(image_paths):
                img = plt.imread(image_paths[k])
                ax.imshow(img)
                ax.axis('off')
                k += 1
            else:
                ax.axis('off')
    plt.tight_layout()
    plt.savefig(os.path.join(output_directory, filename), dpi=150, bbox_inches='tight', pad_inches=0)
    plt.close()

# --------------------------- args ---------------------------

parser = argparse.ArgumentParser(
    description="Plot/compute GC% in masked (lowercase) vs non-masked (uppercase) genome windows; skip high-N windows."
)
parser.add_argument("fasta_file", help="Path to the input FASTA file")
parser.add_argument("species", help="Species name")
parser.add_argument("window", type=int, help="Window size for GC% calculation")
parser.add_argument("--min-nonN", type=float, default=0.5,
                    help="Minimum fraction of non-N bases (A/C/G/T) required to keep a window [default: 0.5]")
args = parser.parse_args()

what_to_operate = [(args.fasta_file, args.species, [args.window])]

plt.rcParams["figure.figsize"] = [20.00, 20.00]

# --------------------------- core ---------------------------

for fasta_path, species, window_list in what_to_operate:
    for win_size in window_list:
        print(f"\n==> Processing: {fasta_path}  window={win_size}")

        output_directory = main_dir(fasta_path, species, win_size)

        per_window_gc = []           # floats for windows kept (all chromosomes)
        GC_chromosoms = OrderedDict()
        d_sec_values = dict()
        lengths = []
        chrom_index = 0

        for rec in SeqIO.parse(fasta_path, "fasta"):
            chrom_index += 1
            rec_id = rec.id
            seq = rec.seq
            Lrec = len(seq)
            lengths.append(Lrec)

            # chromosome-level
            soft_mask_chr = how_many(seq, 'acgt')
            ALL_all_chr   = how_many(seq, 'acgtACGT')
            GC_fraction_chr = gc_fraction(seq) if ALL_all_chr > 0 else 0.0
            softmask_frac_chr = (soft_mask_chr / ALL_all_chr) if ALL_all_chr > 0 else 0.0
            GC_chromosoms[f"{chrom_index}+{rec_id}"] = ((GC_fraction_chr, softmask_frac_chr), Lrec)

            for base in 'acgtACGTNn':
                d_sec_values[base] = d_sec_values.get(base, 0) + seq.count(base)

            # per-window for this chrom
            win_map = {}
            kept = 0
            for w in range(0, Lrec, win_size):
                win = seq[w:w+win_size]
                if not valid_window(win, args.min_nonN):
                    continue

                soft_mask = win.count('a') + win.count('t') + win.count('c') + win.count('g')
                ALL_all   = soft_mask + win.count('A') + win.count('T') + win.count('C') + win.count('G')
                if ALL_all == 0:
                    continue

                gcw = gc_fraction(win)
                per_window_gc.append(gcw)
                win_map[f"{w}|{rec_id}"] = (gcw, soft_mask / ALL_all)
                kept += 1

            print(f"{rec_id}: kept {kept} windows (win={win_size}, len={Lrec})")

            # profile plot (only if any windows)
            if win_map:
                names = list(win_map.keys())
                vals  = list(win_map.values())
                y_gc  = [100.0 * v[0] for v in vals]
                xpos  = [int(n.split('|')[0]) for n in names]
                colors = [v[1] for v in vals]   # 0..1 soft-mask fraction

                fig = plt.figure(figsize=(50, 20))
                ax = fig.add_subplot(111)
                sc = ax.scatter(xpos, y_gc, s=5, c=colors, cmap='RdYlGn', marker="o",
                                label='Color = soft-masked fraction')
                ax.grid(True)
                ax.set_xlim(0, Lrec)
                step = max(Lrec // 20, win_size)
                xticks = np.arange(0, Lrec + 1, step)
                ax.set_xticks(xticks)
                ax.set_xticklabels(['{:,.0f}'.format(x) for x in xticks])
                ax.set_facecolor("lightgrey")
                plt.title(f'GC% values of {species} with {win_size} - chromosome {rec_id}', fontsize=35)
                plt.ylabel('GC percent', fontsize=30)
                plt.xlabel(f'Chromosome {rec_id} windows (bp start)', fontsize=30)
                plt.xticks(fontsize=20, rotation=30)
                plt.yticks(fontsize=30)
                plt.legend()
                cbar = plt.colorbar(sc, label="Soft-masked fraction (0–1)", orientation="horizontal")
                cbar.ax.tick_params(labelsize=35)
                out_png = f'{output_directory}/{os.path.basename(species)}_{win_size}_profile_soft_unmask_{rec_id}.png'
                plt.savefig(out_png, dpi=150)
                plt.close()

            # small “nobar” plot
            if win_map:
                names = list(win_map.keys())
                vals  = list(win_map.values())
                y_gc  = [100.0 * v[0] for v in vals]
                xpos  = [int(n.split('|')[0]) for n in names]
                colors = [v[1] for v in vals]

                fig = plt.figure(figsize=(25, 10))
                ax = fig.add_subplot(111)
                sc = ax.scatter(xpos, y_gc, s=5, c=colors, cmap='RdYlGn', marker="o",
                                label='Color = soft-masked fraction')
                ax.set_xlim(0, Lrec)
                step = max(Lrec // 20, win_size)
                xticks = np.arange(0, Lrec + 1, step)
                ax.set_xticks(xticks)
                ax.set_xticklabels(['{:,.0f}'.format(x) for x in xticks])
                plt.title(f'GC% values of {species} - {win_size} chromosome {rec_id}', fontsize=30)
                plt.xticks(fontsize=20, rotation=30)
                plt.yticks(fontsize=20)
                plt.legend()
                out_png = f'{output_directory}/{os.path.basename(species)}_{win_size}_nobar_soft_unmask_{rec_id}.{Lrec}.png'
                plt.savefig(out_png, dpi=150)
                plt.close()

        # ---------------- write per-window CSV (filtered windows only) ----------------
        perwin_csv = os.path.join(output_directory, f'{os.path.basename(species)}_per_windows_{win_size}.csv')
        with open(perwin_csv, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['GC DNA'])
            for v in per_window_gc:
                w.writerow([v])

        # ---------------- per-chromosome summary CSV ----------------
        chrom_csv = os.path.join(output_directory, f'{os.path.basename(species)}_{win_size}.csv')
        GC_chromosoms = OrderedDict(sorted(GC_chromosoms.items(), key=lambda x: x[1][1], reverse=True))
        with open(chrom_csv, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['index','chromosome','mean_gc_fraction','softmask_fraction','length_bp'])
            for k, v in GC_chromosoms.items():
                w.writerow([k.split('+')[0], k.split('+')[1], v[0][0], v[0][1], v[1]])

        # ---------------- scatter of chromosome means ----------------
        if GC_chromosoms:
            names = [i.split('+')[0] for i in GC_chromosoms.keys()]
            vals  = list(GC_chromosoms.values())
            y_gc_mean = [100.0 * v[0][0] for v in vals]
            color_val = [v[0][1] for v in vals]

            fig = plt.figure(figsize=(15, 10))
            ax = fig.add_subplot(111)
            sc = ax.scatter(names, y_gc_mean, s=150, c=color_val, cmap='RdYlGn', marker="o")
            plt.grid(True)
            plt.title(f'GC% means of size-sorted chromosomes with {win_size} in {species}', fontsize=18)
            plt.ylabel('GC percent', fontsize=15)
            plt.xlabel('Chromosomes', fontsize=15)
            plt.xticks(rotation=30, fontsize=12)
            plt.yticks(fontsize=12)
            plt.colorbar(sc, label="Soft-masked fraction (0–1)", orientation="horizontal")
            out_png = f'{output_directory}/{os.path.basename(species)}_soft_unmask_per_chromosomes_{win_size}_{species}.png'
            plt.savefig(out_png, dpi=150)
            plt.close()

        # ---------------- base-counts & totals ----------------
        base_csv = os.path.join(output_directory, f'{os.path.basename(species)}_{win_size}_base_counts.csv')
        with open(base_csv, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['base','count'])
            for k in sorted(d_sec_values.keys()):
                w.writerow([k, d_sec_values[k]])

        denom = (d_sec_values.get('c', 0) + d_sec_values.get('g', 0) +
                 d_sec_values.get('C', 0) + d_sec_values.get('G', 0) +
                 d_sec_values.get('a', 0) + d_sec_values.get('t', 0) +
                 d_sec_values.get('A', 0) + d_sec_values.get('T', 0))
        cg_num = (d_sec_values.get('c', 0) + d_sec_values.get('g', 0) +
                  d_sec_values.get('C', 0) + d_sec_values.get('G', 0))
        gc_all = (cg_num / denom) if denom > 0 else 0.0
        with open(os.path.join(output_directory, f'{os.path.basename(species)}_all_{win_size}_GCproc.csv'), 'w') as f:
            f.write(f'%GC of {species} is {gc_all}\n')

        atgc_lower = d_sec_values.get('a', 0) + d_sec_values.get('t', 0) + d_sec_values.get('c', 0) + d_sec_values.get('g', 0)
        gc_lower   = d_sec_values.get('c', 0) + d_sec_values.get('g', 0)
        atgc_upper = d_sec_values.get('A', 0) + d_sec_values.get('T', 0) + d_sec_values.get('C', 0) + d_sec_values.get('G', 0)
        gc_upper   = d_sec_values.get('C', 0) + d_sec_values.get('G', 0)

        ratio_csv = os.path.join(output_directory, f'{os.path.basename(species)}_lower_upper_repeat_{win_size}.csv')
        with open(ratio_csv, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['metric','value'])
            w.writerow(['lower',   (gc_lower / atgc_lower) if atgc_lower > 0 else 0.0])
            w.writerow(['upper',   (gc_upper / atgc_upper) if atgc_upper > 0 else 0.0])
            w.writerow(['repeat%', (atgc_lower / (atgc_lower + atgc_upper)) if (atgc_lower + atgc_upper) > 0 else 0.0])

        # ---------------- all-chromosome montage (from nobar plots) ----------------
        nobars = sorted(glob(f"{output_directory}/{os.path.basename(species)}_{win_size}_nobar_soft_unmask_*.png"))
        save_grid(output_directory, f'{species}_ALL_CHROMS_{win_size}.png', nobars, columns=2)

print("\nDone.")
