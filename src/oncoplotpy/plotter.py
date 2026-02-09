"""OncoPlotPy: Genomic waterfall plots in Matplotlib."""
import numpy as np
import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.axes import Axes

from oncoplotpy.transformer import SAMPLE_FIELD, GENE_FIELD, VC_FIELD, VC_NONSYN

MAFTOOLS_COLORS = {
    "Nonstop_Mutation": "#A6CEE3FF",
    "Frame_Shift_Del": "#1F78B4FF",
    "IGR": "#B2DF8AFF",
    "Missense_Mutation": "#33A02CFF",
    "Silent": "#FB9A99FF",
    "Nonsense_Mutation": "#E31A1CFF",
    "RNA": "#FDBF6FFF",
    "Splice_Site": "#FF7F00FF",
    "Intron": "#CAB2D6FF",
    "Frame_Shift_Ins": "#6A3D9AFF",
    "In_Frame_Del": "#FFFF99FF",
    "ITD": "#9E0142FF",
    "In_Frame_Ins": "#D53E4FFF",
    "Translation_Start_Site": "#F46D43FF",
    "Multi_Hit": "#000000FF",
    "Amp": "#EE82EEFF",
    "Del": "#4169E1FF",
    "Complex_Event": "#7B7060FF",
    "pathway": "#535C68FF",
}

MAFTOOLS_NONSYN_COLORS = [MAFTOOLS_COLORS[name] for name in VC_NONSYN]
VC_NONMUTATED_LABEL = 'None'

def matrix_tmb(df: pd.DataFrame, genes: list[str], samples: list[str]) -> pd.DataFrame:
    """
    Compute the mutation burden matrix for the top barplot.

    Aggregates mutation counts per sample across the specified gene set
    before multi-hit logic is applied.

    Parameters
    ----------
    df : pd.DataFrame
        The raw mutation data (MAF-like format) before multi-hit aggregation.
    genes : list of str
        The list of genes to include in the count.
    samples : list of str
        The list of samples to include, defining the final row order.

    Returns
    -------
    pd.DataFrame
        A summary matrix where:
        - Rows: Sample IDs (sorted to match the waterfall horizontal order).
        - Columns: Mutation classifications (e.g., Missense, Nonsense).
        - Values: Integer counts of mutations per sample/type.
    """
    mask = (df[GENE_FIELD].isin(genes)) & (df[VC_FIELD].isin(VC_NONSYN))
    subset = df[mask].copy()

    subset[SAMPLE_FIELD] = pd.Categorical(subset[SAMPLE_FIELD], categories=samples)

    tmb = (
        subset.groupby([SAMPLE_FIELD, VC_FIELD], observed=False)
        .size()
        .unstack(fill_value=0)
    )

    tmb = tmb.reindex(index=samples, columns=VC_NONSYN, fill_value=0)

    return tmb.rename_axis(None, axis=1).astype(int)

def plot_tmb(df: pd.DataFrame, ax: Axes, width: float = 0.8, color: list = None, **kwargs) -> Axes:
    """
    Plot the mutation burden bar chart (TMB) above the waterfall heatmpap.

    Parameters
    ----------
    df : pd.DataFrame
        TMB matrix where columns are samples and rows represent mutation
        classifications (counts).
    ax : matplotlib.axes.Axes
        The axes object to draw the plot onto.
    width : float, default 0.8
        The width of individual bars.
    color : list of str, optional
        A list of colors corresponding to the mutation categories in the
        TMB matrix rows.
    **kwargs : dict
        Additional arguments passed to the pandas/matplotlib bar plot.

    Returns
    -------
    matplotlib.axes.Axes
        The axes object with the TMB plot rendered.
    """
    if color is None:
        color = MAFTOOLS_NONSYN_COLORS
    df.plot.bar(stacked=True, ax=ax, width=width, legend=False,
                color=color, **kwargs)
    ax.set_xlabel('Samples')
    ax.set_ylabel('Count')
    return ax

def plot_waterfall(df: pd.DataFrame, ax: Axes, cmap: list = None, **kwargs) -> Axes:
    """
    Plot the central mutation waterfall (heatmap).

    Parameters
    ----------
    df : pd.DataFrame
        The processed mutation matrix where rows are Genes and columns are Samples.
        Cell values contain categorical mutation types (e.g., 'Missense_Mutation').
    ax : matplotlib.axes.Axes
        The axes object to draw the plot onto.
    cmap : list of str, optional
        A list of colors defining the color map for the different
        mutation classifications.
    **kwargs : dict
        Additional arguments passed to `sns.heatmap` or the plotting backend
        (e.g., `linewidths`, `linecolor`).

    Returns
    -------
    matplotlib.axes.Axes
        The axes object with the waterfall heatmap rendered.
    """
    df = df.copy().fillna(VC_NONMUTATED_LABEL)
    l = reversed(VC_NONSYN + ['Multi_Hit', 'None'])
    d = {k: v for v, k in enumerate(l)}
    df = df.map(lambda x: d[x])
    if cmap is None:
        cmap = list(reversed(MAFTOOLS_NONSYN_COLORS + ['k', 'white']))
    sns.heatmap(df, ax=ax, xticklabels=False,
                cbar=False, vmin=0, vmax=len(cmap), cmap=cmap, **kwargs)
    ax.set_xlabel('Samples')
    ax.set_ylabel('')
    return ax

def matrix_genes(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute the gene mutation summary matrix for the right-side barplot.

    Calculates the distribution of mutation types (e.g., Missense, Nonsense,
    Multi-hit) for each gene across all included samples.

    Parameters
    ----------
    df : pd.DataFrame
        The processed mutation matrix (pivoted) where rows are Genes and
        columns are Samples.

    Returns
    -------
    pd.DataFrame
        A summary matrix where:
        - Rows: Genes (matches waterfall vertical axis).
        - Columns: Mutation Variant Classifications (plus 'Multi_Hit').
        - Values: Frequency counts of samples containing each mutation type.
    """
    df = df.apply(lambda r: r.value_counts(), axis=1)
    for varcls in VC_NONSYN + ['Multi_Hit']:
        if varcls not in df.columns:
            df[varcls] = np.nan
    df = df[VC_NONSYN + ['Multi_Hit']]
    df = df.fillna(0)
    return df

def plot_genes(df: pd.DataFrame, ax: Axes, color: list = None, **kwargs) -> Axes:
    """
    Plot the gene mutation summary bar chart to the right of the waterfall heatmap.

    Uses a horizontal stacked bar chart to display the frequency of mutation
    types for each gene. The order is reversed internally to align with the
    top-to-bottom layout of the heatmap.

    Parameters
    ----------
    df : pd.DataFrame
        Summary matrix where:
        - Rows: Genes (must match waterfall vertical order).
        - Columns: Mutation classifications (e.g., Missense, Multi_Hit).
        - Values: Counts or proportions of mutated samples.
    ax : matplotlib.axes.Axes
        The axes object to draw the plot onto.
    color : list of str, optional
        Colors for the mutation categories. Defaults to MAFTOOLS_NONSYN_COLORS
        with black ('k') appended for Multi_Hit.
    **kwargs : dict
        Additional arguments passed to `df.plot()`.

    Returns
    -------
    matplotlib.axes.Axes
        The axes object with the horizontal bar chart rendered.
    """
    if color is None:
        color = MAFTOOLS_NONSYN_COLORS
    color = color + ['k']
    df = df.iloc[::-1]
    kind = 'barh'
    xlabel, ylabel = 'Count', ''
    df.plot(
        kind=kind, ax=ax, stacked=True, legend=False,
        color=color, **kwargs
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    return ax

def legend_handles(labels: list, colors='tab10') -> list:
    """
    Create a list of matplotlib patch handles for the plot legend.

    Generates colored rectangles (patches) that associate mutation types
    or annotation categories with their respective colors.

    Parameters
    ----------
    labels : list of str
        The text labels to display in the legend (e.g., mutation classifications).
    colors : str or list, default 'tab10'
        If a string, the name of a matplotlib colormap to use.
        If a list, the specific colors (hex, RGB, or names) to map to the labels.

    Returns
    -------
    list of matplotlib.patches.Patch
        A list of handle objects ready to be passed to `ax.legend()`.

    Raises
    ------
    TypeError
        If the colors argument is neither a string (colormap name) nor a list.
    """
    if isinstance(colors, str):
        colors = plt.get_cmap(colors).colors
    elif isinstance(colors, list):
        pass
    else:
        raise TypeError(f'Incorrect type of colors: {type(colors)}')
    handles = []
    for i, label in enumerate(labels):
        handles.append(mpatches.Patch(color=colors[i], label=label))
    return handles

def plot_annot(df, group_col: str, ax: Axes, colors: str ='tab10',
               xticklabels: bool = True) -> tuple[Axes, list]:
    """
        Plot a categorical annotation track below the waterfall heatmap.

        Converts categorical labels into a color-coded bar, ensuring the
        horizontal alignment matches the samples in the waterfall.

        Parameters
        ----------
        df : pd.DataFrame
            Annotation data where the index corresponds to Sample IDs.
        group_col : str
            The specific column name in `df` to be visualized.
        ax : matplotlib.axes.Axes
            The axes object to draw the annotation track onto.
        colors : str, default 'tab10'
            The name of the matplotlib colormap to use for categories.
        xticklabels : bool, default True
            Whether to display sample IDs below the annotation bar.

        Returns
        -------
        ax : matplotlib.axes.Axes
            The axes object with the annotation heatmap.
        handles : list of matplotlib.patches.Patch
            A list of legend handles for the specific categories in `group_col`.
        """
    s = df[group_col]
    s = s.fillna('N/A')
    group_order = s.unique()

    d = {k: v for v, k in enumerate(group_order)}
    df = s.to_frame().map(lambda x: d[x])

    cmap = plt.get_cmap(colors)
    color_list = [cmap(i / len(group_order)) for i in range(len(group_order))]

    sns.heatmap(
        df.T, ax=ax, cmap=color_list, xticklabels=xticklabels, cbar=False,
        linewidths=0.5
    )
    # twice thinner heatmap in case xticklabels
    ba = 0.5 / len(s)
    if xticklabels:
        ax.set_box_aspect(ba)
        ax.set_anchor('N')
    else:
        ax.set_box_aspect(None)
        ax.set_anchor('C')
    ax.set_xlabel('')
    ax.set_ylabel(group_col)
    ax.set_yticks([])

    handles = legend_handles(group_order, color_list)
    return ax, handles
