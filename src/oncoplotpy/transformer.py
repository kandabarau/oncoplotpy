"""OncoPlotPy: Genomic waterfall plots in Matplotlib."""
from typing import Callable, Any
import pandas as pd
import numpy as np


# Constants
SAMPLE_FIELD = 'Tumor_Sample_Barcode'
GENE_FIELD = 'Hugo_Symbol'
VC_FIELD = 'Variant_Classification'
MAF_REQUIRED_FIELDS = [SAMPLE_FIELD, GENE_FIELD, VC_FIELD]
VC_NONSYN = [
    "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
    "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site",
    "Translation_Start_Site", "Missense_Mutation"
]
VC_MULTIHIT_LABEL = 'Multi_Hit'


def custom_agg(x, nonsyn: list[str] = None, multihit: str = None) -> str:
    """
    A function to aggregate mutations while pivoting maf into gene/subject matrix
    :param x:
    :param nonsyn: a list of non-synonymous mutation types in focus
    :param multihit: a label for multiple mutations in one gene
    :return: a string to populate gene/subject matrix
    """
    if nonsyn is None:
        nonsyn = VC_NONSYN
    if multihit is None:
        multihit = VC_MULTIHIT_LABEL
    filtered = [item for item in x if item in nonsyn]
    if len(filtered) == 1:
        return filtered[0]

    if len(filtered) > 1:
        return multihit

    return np.nan


def read_maf(path: str, fields: list[str] = None,
             annotation_fields: list[str] = None) -> pd.DataFrame:
    """
    A function to read the system file into pandas df
    Controls the required columns to be present and loads additional non-required data
    :param path: a path to .maf file
    :param fields: a minimal list of fields from .maf needed to pivot into gene/subject matrix
    :param annotation_fields: a list of subject annotation fields in .maf, in order they will be used in sorts
    :return: a dataframe read from raw file
    """
    # Check for all columns in .maf
    available_columns = pd.read_table(path, nrows=0).columns.tolist()
    if hasattr(path, 'seek'):
        path.seek(0)

    # either provided via parameter of default required fields needed for gene/subject matrix
    load_fields = list(fields) if fields else list(MAF_REQUIRED_FIELDS)
    # subject annotation fields in .maf if any
    extra_annotations = annotation_fields if annotation_fields else []

    # target columns, given in parameters
    target_cols = list(dict.fromkeys(load_fields + extra_annotations))
    # which of target columns are present in .maf
    valid_cols = [c for c in target_cols if c in available_columns]

    # can't proceed if not all required columns present in .maf
    missing_required = [f for f in MAF_REQUIRED_FIELDS if f not in available_columns]
    if missing_required:
        raise ValueError(
            f"Required {missing_required} columns missing from file."
            f"Check your file headers."
        )

    dtype_dict = {SAMPLE_FIELD: str, GENE_FIELD: str, VC_FIELD: str}
    df = pd.read_table(path, usecols=valid_cols, dtype=dtype_dict)
    if hasattr(path, 'seek'):
        path.seek(0)

    # provided in parameters but not present in .maf - ignored in matrix sort
    for ann in extra_annotations:
        if ann and ann not in df:
            print(f"no '{ann}' read from maf, will not be used to annotate subjects")

    return df


def units_sorter(df: pd.DataFrame, column_to_order: str, orderby: str) -> Any | None:
    """
    A function that orders a target features based on external index
    :param df: input data frame
    :param column_to_order: a string column with features to order
    :param orderby: an integer column to order the string column with features
    :return: unique features ordered
    """
    if not orderby or orderby not in df.columns:
        return None
    valid_df = df.dropna(subset=[column_to_order, orderby])
    sorted_units = (
        valid_df[[column_to_order, orderby]]
        .drop_duplicates()
        .sort_values(orderby)[column_to_order]
        .tolist()
    )
    return sorted_units


class OncoMatrix:
    """
        Pivots raw MAF data into a matrix format suitable for oncoplots.

        This class handles data aggregation, filtering and sorting.
        """
    def __init__(self, maf_path: str, fields: list[str] = None,
                 annotation_fields: list[str] = None
                 ):
        self.raw_df = read_maf(path=maf_path, fields=fields,
                               annotation_fields=annotation_fields)
        self.matrix: pd.DataFrame = pd.DataFrame()

        self.annotation_fields = [
            f for f in (annotation_fields or [])
            if f in self.raw_df.columns
        ]
        self.annotation = (
            self.raw_df[[SAMPLE_FIELD] + self.annotation_fields]
            .groupby(SAMPLE_FIELD)
            .first()  # This ensures 1 row per sample ID
        )
        self.samples = self.raw_df[SAMPLE_FIELD].dropna().unique().tolist()
        self.genes = self.raw_df[GENE_FIELD].dropna().unique().tolist()

    def prepare(self, genes: str = GENE_FIELD, samples: str = SAMPLE_FIELD,
                classification: str = VC_FIELD, agg_fn: Callable = custom_agg, remove_non_mutated: bool = False):
        """
        Pivot the raw MAF data into a Gene x Sample mutation matrix.

        This method populates the internal matrix by aggregating variant
        classifications and cleaning null indices.

        Parameters
        ----------
        genes : str
            Column name for gene symbols.
        samples : str
            Column name for sample/patient identifiers.
        classification : str
            Column name for variant types (e.g., Missense).
        agg_fn : Callable
            Function to handle multiple mutations in the same gene/sample cell.
        remove_non_mutated : bool, default False
            If True, samples with zero mutations across the gene set are dropped.

        Returns
        -------
        self : OncoMatrix
            The current instance with .matrix, .genes, and .samples attributes populated.
        """

        if not remove_non_mutated:
            print(f"keeping all {self.raw_df[SAMPLE_FIELD].count()} samples")

        self.matrix = self.raw_df.pivot_table(
            index=genes,
            columns=samples,
            values=classification,
            aggfunc=agg_fn,
            dropna=remove_non_mutated
        )
        if None in self.matrix.index or pd.isna(self.matrix.index).any():
            self.matrix = self.matrix.loc[self.matrix.index.notna()]
        if None in self.matrix.columns or pd.isna(self.matrix.columns).any():
            self.matrix = self.matrix.loc[:, self.matrix.columns.notna()]
        self.genes, self.samples = self.matrix.index.tolist(), self.matrix.columns.tolist()
        return self

    def sort_genes(self, genes_order: list[str] = None):
        """
        ort rows in the mutation matrix by mutation frequency or a provided list.

        If no order is specified, genes are ranked by the total number of
        mutated samples in descending order.

        Parameters
        ----------
        genes_order : list of str, optional
            A custom list of genes to force a specific vertical order.

        Returns
        -------
        self : OncoMatrix
            The instance with updated .matrix and .genes attributes.
        """
        if genes_order is not None:
            self.matrix = self.matrix.reindex(index=genes_order)
        else:
            sorted_indices = (
                self.matrix.count(axis=1)
                .sort_values(ascending=False)
                .index
            )
            self.matrix = self.matrix.reindex(index=sorted_indices)

        self.genes = self.matrix.index.tolist()
        return self

    def sort_samples(self, samples_order: list[str] = None, sortby: list[str]|str = None):
        """
        Sort columns in the mutation matrix to create a waterfall effect.

        Organizes samples based on a hierarchical priority: custom list first,
        followed by clinical annotations (if requested), and finally the
        standard mutation "staircase" logic.

        Parameters
        ----------
        samples_order : list of str, optional
            A custom list of sample IDs to force a specific horizontal order.
        sortby : list[str] or str, optional
            Annotation column(s) to prioritize for sorting. Samples are grouped
            by category size (largest groups first) within these fields.

        Returns
        -------
        self : OncoMatrix
            The instance with updated .matrix, .samples, and .annotation order.
        """
        if samples_order is not None:
            self.matrix = self.matrix.reindex(columns=samples_order)
        else:
            if self.matrix.empty:
                return self

            # 1. Calculate the Mutation Waterfall (Staircase)
            staircase_columns = (
                self.matrix.notna()
                .sort_values(by=list(self.matrix.index), axis=1, ascending=False)
                .columns.tolist()
            )
            staircase_rank = {s: i for i, s in enumerate(staircase_columns)}

            # 2. Build the Sorting Dataframe
            current_samples = self.matrix.columns.tolist()
            sort_df = self.annotation.loc[current_samples].copy()

            sort_keys = []
            ascending_flags = []

            # 3. Handle Annotation Sorting
            if sortby and isinstance(sortby, list):
                for field in sortby:
                    if field not in self.annotation.columns:
                        continue
                    freq_col = f"{field}_group_size"
                    sort_df[freq_col] = sort_df.groupby(field)[field].transform('count')
                    sort_keys.extend([freq_col, field])
                    ascending_flags.extend([False, True])
            elif sortby and isinstance(sortby, str) and sortby in self.annotation.columns:
                freq_col = f"{sortby}_group_size"
                sort_df[freq_col] = sort_df.groupby(sortby)[sortby].transform('count')
                sort_keys.extend([freq_col, sortby])
                ascending_flags.extend([False, True])

            # 4. Always add the staircase rank as the final tie-breaker
            sort_df['staircase_rank'] = sort_df.index.map(staircase_rank).astype(int)
            sort_keys.append('staircase_rank')
            ascending_flags.append(True)

            final_order = sort_df.sort_values(
                by=sort_keys,
                ascending=ascending_flags
            ).index.tolist()

            self.matrix = self.matrix[final_order].rename_axis(None, axis=1)

        # Sync metadata
        self.samples = self.matrix.columns.tolist()
        if not self.annotation.empty:
            self.annotation = self.annotation.reindex(self.samples)

        return self

    @property
    def data(self):
        """Returns the processed mutation matrix."""
        return self.matrix
