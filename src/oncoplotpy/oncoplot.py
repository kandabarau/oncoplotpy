"""OncoPlotPy: Genomic waterfall plots in Matplotlib."""
import math
import matplotlib.pyplot as plt

from oncoplotpy.transformer import OncoMatrix, SAMPLE_FIELD, VC_FIELD, VC_NONSYN, VC_MULTIHIT_LABEL
from oncoplotpy.plotter import (matrix_tmb, plot_tmb, plot_waterfall, matrix_genes, plot_genes,
                                legend_handles, plot_annot, MAFTOOLS_NONSYN_COLORS)

plt.rcParams['font.sans-serif'] = ['Arial', 'Liberation Sans', 'sans-serif']
plt.rcParams['font.family'] = 'sans-serif'


def run_oncoplot(
        maf_path: str,
        svg_path: str = None,
        genes_order: list[str] = None,
        samples_order: list[str] = None,
        annotation_fields: list[str] = None,
        annotation_sort_by: bool = False,
        remove_non_mutated: bool = False,
        tile_width: float = 1.0,
        tile_height: float = 0.5,
        sum_mut_genes_ticks: bool = False,
        pct_mut_samples_ticks: bool = True,
        show_sample_id: bool = False,
        labels_font_size: float = 20.0,
        tick_labels_font_size: float = 18.0,
        legend_font_size: float = 20.0,
        tbm_height: float = 1.5,
        gene_hist_width: float = 1.0,
        ann_height_set: float = 0.50,
        legend_height: float = 2.50,
        legend_nrows: int = 4
):
    """
    Instantiate OncoMatrix and assemble a waterfall heatmap with marginal mutation bar charts.

    This function processes a MAF file to create a publication-quality oncoplot,
    including Tumor Mutational Burden (TMB) on the top and gene mutation frequencies on the right.

    Parameters
    ----------
    maf_path : str
        Path to the input .maf file.
    svg_path : str, optional
        Path to save the .svg plot. If None, the plot is rendered to the
        active output device (e.g., Jupyter notebook).
    genes_order : list of str, optional
        Manual list of gene names to subset and order the vertical axis.
    samples_order : list of str, optional
        Manual list of sample IDs to subset and order the horizontal axis.
    annotation_fields : list of str, optional
        Column names in the MAF to use for clinical annotation tracks.
    annotation_sort_by : bool, default False
        Whether to prioritize sorting samples by annotation features before
        applying the mutation waterfall logic.
    remove_non_mutated : bool, default False
        If True, samples with zero mutations in the selected genes are excluded.
    tile_width : float, default 1.0
        Width of a single mutation cell (inches). Recommended < 0.05 for
        large cohorts (>500 samples).
    tile_height : float, default 0.5
        Height of a single mutation cell (inches).
    sum_mut_genes_ticks : bool, default False
        Show the raw count of mutated genes per sample on the TMB axis.
    pct_mut_samples_ticks : bool, default True
        Show the percentage of the cohort (initial .maf) mutated for each gene.
    show_sample_id : bool, default False
        Whether to print sample barcodes along the X-axis.
    labels_font_size : float, default 20.0
        Font size for primary axis labels.
    tick_labels_font_size : float, default 18.0
        Font size for gene names and scale ticks.
    legend_font_size : float, default 20.0
        Font size for the mutation and annotation legends.
    tbm_height : float, default 1.5
        Height of the top TMB bar plot (inches).
    gene_hist_width : float, default 1.0
        Width of the right-side gene histogram (inches).
    ann_height_set : float, default 0.5
        Height of each annotation track (inches).
    legend_height : float, default 2.5
        Total height reserved for the legend area.
    legend_nrows : int, default 4
        Maximum number of rows per legend category.

    Returns
    -------
    None
    Renders the plot to a file or the display device.
    """
    onco = OncoMatrix(
        maf_path=maf_path,
        annotation_fields=annotation_fields,
    )
    onco.prepare(remove_non_mutated=remove_non_mutated)
    if genes_order is not None:
        onco.sort_genes(genes_order=genes_order)
    else:
        onco.sort_genes()

    if samples_order is not None:
        onco.sort_samples(samples_order=samples_order)
    elif annotation_fields is not None and annotation_sort_by:
        onco.sort_samples(sortby=onco.annotation_fields)
    else:
        onco.sort_samples()

    tmb = matrix_tmb(onco.raw_df, onco.genes, onco.samples)
    gene_hist = matrix_genes(onco.data)

    n_samples: int = len(onco.samples)
    waterfall_width: float = n_samples * tile_width

    n_genes: int = len(onco.genes)
    waterfall_height: float = n_genes * tile_height


    # wider annotation axis in case sample labels
    ann_height: float = ann_height_set if not onco.annotation.empty else 0
    ann_height: float = 2 * ann_height if show_sample_id else ann_height


    if n_samples < 12:
        legend_nrows = 8
        legend_height = 5.00

    fig_width: float = waterfall_width + gene_hist_width
    width_space: float = 1 / fig_width
    fig_height: float = tbm_height + waterfall_height + ann_height + legend_height
    height_space: float = 1 / fig_height

    figsize: list[float] = [fig_width, fig_height]

    _fig, axes = plt.subplots(nrows=4, ncols=2, figsize=figsize,
                             gridspec_kw={'height_ratios': [tbm_height,
                                                            waterfall_height,
                                                            ann_height,
                                                            legend_height
                                                            ],
                                          'width_ratios': [waterfall_width,
                                                           gene_hist_width
                                                           ]
                                          },
                             )
    [[ax_tmb, ax2], [ax_wf, ax_genes], [ax_ann, ax6], [ax_lgd, ax8]] = axes

    # Create the TMB plot.
    plot_tmb(tmb, ax=ax_tmb)
    sum_mut_sample = tmb.T.sum(axis=0)
    ax_tmb.set_xticks([])
    if sum_mut_genes_ticks:
        height_space = height_space * 2
        sum_mut_sample = tmb.T.sum(axis=0)
        ax_tmb.set_xticks(range(len(sum_mut_sample)))
        ax_tmb.set_xticklabels([f'{v:.0f}' for v in sum_mut_sample], rotation=0)

    ax_tmb.spines['right'].set_visible(False)
    ax_tmb.spines['top'].set_visible(False)
    ax_tmb.spines['bottom'].set_visible(False)
    ax_tmb.set_ylabel('Mutation\nCount', fontsize=labels_font_size)
    ax_tmb.set_xlabel('')
    ax_tmb.set_yticks([0, sum_mut_sample.max()])
    ax_tmb.tick_params(which='major',
                       labelsize=tick_labels_font_size,
                       labelfontfamily='arial')
    ax_tmb.set_xlim(-0.5, len(onco.samples) - 0.5)

    ax2.remove()

    plot_waterfall(onco.data, ax=ax_wf, linewidths=0.5,
                   linecolor='#EAEAEA',
                   )
    ax_wf.set_xlabel('')
    ax_wf.tick_params(axis='y', which='major', length=0,
                      labelrotation=0, labelsize=tick_labels_font_size,
                      )
    for label in ax_wf.get_yticklabels():
        label.set_fontstyle('italic')

    plot_genes(gene_hist, ax=ax_genes, width=0.9)
    ax_genes.spines['right'].set_visible(False)
    ax_genes.spines['left'].set_visible(False)
    ax_genes.spines['top'].set_visible(False)
    ax_genes.set_xlabel('Patients', fontsize=labels_font_size)
    sum_mut_gene = gene_hist.sum(axis=1).iloc[::-1]
    ax_genes.set_xticks([0, sum_mut_gene.max()])
    ax_genes.set_yticks([])

    total_samples = onco.raw_df[SAMPLE_FIELD].nunique()
    if pct_mut_samples_ticks:
        width_space = width_space * 2
        ax_genes.set_yticks(range(len(sum_mut_gene)))
        ax_genes.set_yticklabels([f'{v:.0f}%' for v in sum_mut_gene / total_samples * 100])

    ax_genes.set_ylim(-0.5, len(onco.genes) - 0.5)
    ax_genes.tick_params(which='major',
                         labelsize=tick_labels_font_size,
                         labelfontfamily='arial',
                         )

    if not onco.annotation.empty:
        _, handles_ann = plot_annot(onco.annotation, group_col=onco.annotation.columns[0],
                                    xticklabels=show_sample_id,
                                    ax=ax_ann,
                                    )
        ax_ann.tick_params(which='major',
                           labelsize=tick_labels_font_size * 0.8,
                           labelfontfamily='arial',
                           )
        # keep legend labels italic
        # store annotation legend in case annotation exists and matters (more than one factor level)
        if len(handles_ann) > 1:
            leg_ann = ax_lgd.legend(handles=handles_ann, loc='upper left', title=onco.annotation.columns[0],
                                    ncol=math.ceil(len(handles_ann) / legend_nrows),
                                    prop={'style': 'italic', 'size': legend_font_size},
                                    fontsize=legend_font_size, title_fontsize=legend_font_size,
                                    frameon=False,
                                    )
            ax_ann.set_ylabel('')
            ax_lgd.add_artist(leg_ann)
        else:
            ax_ann.remove()
    else:
        ax_ann.remove()

    ax6.remove()

    handles_vc = legend_handles(VC_NONSYN + [VC_MULTIHIT_LABEL],
                                MAFTOOLS_NONSYN_COLORS + ['k'])
    handles_vc = [
        handle for handle in handles_vc if handle.get_label() in set(onco.raw_df[VC_FIELD]) | {VC_MULTIHIT_LABEL}
    ]
    leg_vc = ax_lgd.legend(handles=handles_vc,
                           loc='upper right',
                           title=VC_FIELD,
                           ncol=math.ceil(len(handles_vc) / legend_nrows),
                           labels=[
                               handle.get_label().replace("_",
                                                          " ").replace(" Mutation",
                                                                       "") for handle in handles_vc
                           ],
                           fontsize=legend_font_size, title_fontsize=legend_font_size, frameon=False,
                           )

    ax_lgd.add_artist(leg_vc)
    ax_lgd.axis('off')

    ax8.remove()

    plt.subplots_adjust(wspace=width_space, hspace=height_space)
    if svg_path:
        plt.savefig(svg_path, format="svg", bbox_inches='tight')
        plt.close()
    else:
        plt.show()
