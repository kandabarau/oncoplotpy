"""OncoPlotPy: Genomic waterfall plots in Matplotlib."""
import sys
import os
import argparse
import math
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt

from oncoplotpy.transformer import *
from oncoplotpy.plotter import *
from oncoplotpy.oncoplot import run_oncoplot

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='OncoPlotPy CLI Tool - Genomic Visualization')

    # --- Data Input (Required) ---
    parser.add_argument('--maf', required=True,
                        help='Path to input MAF file')

    # --- Sorting Order (Integer Columns) ---
    parser.add_argument('--genes', nargs='*', default=None,
                        help='List of genes in order (e.g., --genes TP53 EGFR)')
    parser.add_argument('--samples', nargs='*', default=None,
                        help='List of sample IDs in order')

    # --- Remove non-mutated ---
    parser.add_argument('--remove', action='store_true', dest='nonmut_rm', default=False,
                        help='Remove non-mutated samples')

    # --- Data & Annotation ---
    parser.add_argument('--ann', nargs='*', default=None,
                        help='Annotation columns to plot (e.g., --ann Stage Gender)')
    parser.add_argument('--ann_sort_by', action='store_true', default=False,
                        help='Does annotation precede in sample sort?')

    parser.add_argument('--print_ids', action='store_true', default=False,
                        help='Show sample IDs on the X-axis')
    parser.add_argument('--pct', action='store_true', dest='pct_ticks', default=False,
                        help='print percentage of mutated samples on gene histogram')
    parser.add_argument('--gene-ticks', action='store_true', default=False,
                        help='Show mutation count ticks on TMB axis')

    # --- Visual Layout ---
    parser.add_argument('--tile-width', type=float, default=1.0,
                        help='oncoplot tile width')
    parser.add_argument('--tile-height', type=float, default=0.5,
                        help='oncoplot tile height')
    parser.add_argument('--font-labels', type=float, default=20.0,
                        help='Font size for axis titles')
    parser.add_argument('--font-ticks', type=float, default=18.0,
                        help='Font size for tick labels')
    parser.add_argument('--font-legend', type=float, default=20.0,
                        help='Font size for legend text')
    parser.add_argument('--tmb-h', type=float, default=1.5,
                        help='Height of the TMB plot (inches)')
    parser.add_argument('--gene-w', type=float, default=1.0,
                        help='Width of the gene histogram (inches)')
    parser.add_argument('--leg-rows', type=int, default=4,
                        help='Number of rows per legend column')
    args = parser.parse_args()

    svg_path = str(Path(args.maf).with_suffix('').with_suffix('.svg'))

    run_oncoplot(
        maf_path=args.maf,
        svg_path=svg_path,
        remove_non_mutated=args.nonmut_rm,
        genes_order=args.genes,
        samples_order=args.samples,
        annotation_fields=args.ann,
        annotation_sort_by=args.ann_sort_by,
        tile_width=args.tile_width,
        tile_height=args.tile_height,
        show_sample_id=args.print_ids,
        labels_font_size=args.font_labels,
        tick_labels_font_size=args.font_ticks,
        legend_font_size=args.font_legend,
        sum_mut_genes_ticks=args.gene_ticks,
        pct_mut_samples_ticks=args.pct_ticks,
        tbm_height=args.tmb_h,
        gene_hist_width=args.gene_w,
        legend_nrows=args.leg_rows,
    )
