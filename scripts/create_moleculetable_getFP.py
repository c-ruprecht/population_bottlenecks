#!/usr/bin/env python3
"""
Creates molecule tables for getFP analysis from UMI clustering output.

This script aggregates individual sample molecule tables into wide-format tables
grouped by experimental conditions (cage, dose, etc.) and filters for barcodes
present in gavage samples.
"""

import pandas as pd
import os
import argparse
import sys


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Create molecule tables for getFP analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Simple project (P4C2T4T5)
  python create_moleculetable_getFP.py \\
    --project P4C2T4T5 \\
    --ill-sheets ILL130_P4C2T4.xlsx,ILL135_sample_sheet.xlsx \\
    --metadata P4C2T4_cr_mouse_sample.tsv,P4C2T5_cr_mouse_sample.tsv \\
    --basepath /path/to/data \\
    --output-dir /path/to/output \\
    --strains ST1,ST2,ST4,ST5,ST6

  # Dose-response project (P4C3T1)
  python create_moleculetable_getFP.py \\
    --project P4C3T1 \\
    --ill-sheets ILL134_sample_sheet.xlsx \\
    --metadata P4C3T1_cr_mouse_sample.tsv \\
    --basepath /path/to/data \\
    --output-dir /path/to/output \\
    --strains ST1,ST2 \\
    --grouping-column gavage_dose_200ul
        """
    )

    # Required arguments
    parser.add_argument(
        "--project",
        required=True,
        help="Project name (e.g., P4C1T8, P4C2T4T5, P4C3T1, P4C3T2, P4C4T1)"
    )
    parser.add_argument(
        "--ill-sheets",
        required=True,
        help="Comma-separated paths to ILL Excel sheet(s)"
    )
    parser.add_argument(
        "--metadata",
        required=True,
        help="Comma-separated paths to metadata TSV file(s)"
    )
    parser.add_argument(
        "--basepath",
        required=True,
        help="Base directory containing {ILL}/readtables/ folders"
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Output directory for molecule tables"
    )
    parser.add_argument(
        "--strains",
        required=True,
        help="Comma-separated list of strains (e.g., ST1,ST2,ST4)"
    )

    # Optional arguments
    parser.add_argument(
        "--gavage-sample-type",
        default=None,
        help="Sample type identifier for gavage samples (default: auto-detect)"
    )
    parser.add_argument(
        "--ntc-pattern",
        default="WATER|NTC",
        help="Regex pattern to identify NTC/water samples (default: 'WATER|NTC')"
    )
    parser.add_argument(
        "--grouping-column",
        default=None,
        help="Column to group samples by (e.g., 'cage_id', 'gavage_dose_200ul')"
    )
    parser.add_argument(
        "--special-mode",
        default=None,
        choices=['c1gavage', 'invitro'],
        help="Special processing mode for specific project variants"
    )
    parser.add_argument(
        "--gavage-project-path",
        default=None,
        help="Path to import gavage data from (for in vitro experiments)"
    )

    return parser.parse_args()


def gavage_only_readtables(indir, outdir):
    """
    Filter molecule tables to keep only barcodes present in gavage samples.

    Processes all CSV files in indir, keeping only rows where at least one
    gavage column has non-zero values. Applies downsampling if totals exceed
    R integer limits.

    Args:
        indir (str): Directory containing input CSV files
        outdir (str): Directory to write filtered output files
    """
    # Create output directory
    os.makedirs(outdir, exist_ok=True)

    # R integer maximum
    MAX_INT = 2147483647
    # Use a safe threshold (2 billion) to stay well below the limit
    SAFE_THRESHOLD = 2e9

    for file in os.listdir(indir):
        if file.endswith('.csv') and not file.startswith('._') and 'onlygavage' not in file:
            print(f"\nProcessing: {file}")
            df = pd.read_csv(os.path.join(indir, file))

            # Get gavage columns (case-insensitive)
            li_gavage = [col for col in df.columns if 'gavage' in col.lower()]

            if not li_gavage:
                print(f"  WARNING: No gavage columns found in {file}, skipping...")
                continue

            # Drop rows where all gavage columns are 0
            df = df[df[li_gavage].sum(axis=1) > 0].copy()

            # Check total reads across all numeric columns
            numeric_cols = df.select_dtypes(include=['int64', 'float64']).columns

            # Calculate max total for any single sample column
            max_sample_total = df[numeric_cols].sum(axis=0).max()

            print(f"  Max sample total: {max_sample_total:,.0f}")

            if max_sample_total > SAFE_THRESHOLD:
                # Calculate scaling factor to bring below threshold
                scale_factor = max_sample_total / SAFE_THRESHOLD
                print(f"  WARNING: Exceeds safe threshold!")
                print(f"  Scaling down by factor: {scale_factor:.2f}")

                # Downsample all numeric columns proportionally
                df[numeric_cols] = (df[numeric_cols] / scale_factor).round().astype(int)

                # Verify new max
                new_max = df[numeric_cols].sum(axis=0).max()
                print(f"  New max sample total: {new_max:,.0f}")

            # Convert numeric columns to integers to save space
            df[numeric_cols] = df[numeric_cols].astype(int)

            # Save output
            output_file = os.path.join(
                outdir,
                file.replace('_moleculestable.csv', '_onlygavage_moleculestable.csv')
            )
            df.to_csv(output_file, index=False)
            print(f"  Saved: {output_file}")


def load_ill_metadata(ill_sheet_paths, metadata_paths, project):
    """
    Load and merge ILL sequencing sheets with sample metadata.

    Args:
        ill_sheet_paths (list): List of paths to ILL Excel files
        metadata_paths (list): List of paths to metadata TSV files
        project (str): Project name for filtering

    Returns:
        tuple: (merged_df, ntc_df) - merged sample info and NTC samples
    """
    all_ill_dfs = []
    all_meta_dfs = []
    all_ntc_dfs = []

    # Load ILL sheets
    for ill_path in ill_sheet_paths:
        print(f"Loading ILL sheet: {ill_path}")
        df_ill = pd.read_excel(
            ill_path,
            sheet_name='Individual Library Information',
            header=1
        )
        df_ill = df_ill.dropna(subset=['Sample Name/Pool Name*'])

        # Parse sample information
        df_ill['sample_id'] = df_ill['Sample Name/Pool Name*'].apply(
            lambda x: str(x).split('_')[0]
        )
        df_ill['project'] = df_ill['Library Name*'].apply(
            lambda x: str(x).split('_')[1]
        )

        # Handle ILL138/ILL139 naming correction
        df_ill['ill'] = df_ill['Library Name*'].apply(
            lambda x: 'ILL138' if str(x).split('_')[0] == 'ILL139'
            else str(x).split('_')[0]
        )

        df_ill['file_name'] = df_ill['Sample Name/Pool Name*'].apply(
            lambda x: str(x).replace('_', '-')
        )

        all_ill_dfs.append(df_ill)

        # Extract NTC samples
        df_ntc = df_ill[df_ill['Sample Name/Pool Name*'].str.contains("WATER|NTC")]
        all_ntc_dfs.append(df_ntc)

    # Combine all ILL sheets
    df_ill_combined = pd.concat(all_ill_dfs, axis=0, ignore_index=True)
    df_ntc_combined = pd.concat(all_ntc_dfs, axis=0, ignore_index=True)

    # Load metadata sheets
    for meta_path in metadata_paths:
        print(f"Loading metadata: {meta_path}")
        df_meta = pd.read_csv(meta_path, sep='\t')
        # Filter for DNA samples only
        df_meta = df_meta[df_meta['sample_type'].str.contains('DNA', na=False)].copy()
        all_meta_dfs.append(df_meta)

    # Combine all metadata
    df_meta_combined = pd.concat(all_meta_dfs, axis=0, ignore_index=True)

    # Merge ILL and metadata
    df_merge = pd.merge(
        df_ill_combined[['sample_id', 'ill', 'file_name', 'project']],
        df_meta_combined,
        how='right',
        on=['sample_id', 'project']
    )

    # Drop rows without file_name (not sequenced)
    df_merge = df_merge.dropna(subset=['file_name'])

    print(f"\nMerged data: {len(df_merge)} samples")

    return df_merge, df_ntc_combined


def load_molecule_table(basepath, ill, file_name):
    """
    Load a single molecule table CSV file.

    Args:
        basepath (str): Base directory path
        ill (str): ILL run name
        file_name (str): Sample file name

    Returns:
        pd.DataFrame: Molecule table with columns [strain, umi_seq, molecules]
    """
    filepath = os.path.join(
        basepath, ill, 'readtables', file_name,
        f"{file_name}_total_moleculestable.csv"
    )
    return pd.read_csv(filepath)


def process_ntc_samples(df_ntc, basepath):
    """
    Load and pivot NTC sample data.

    Args:
        df_ntc (pd.DataFrame): DataFrame with NTC sample info
        basepath (str): Base directory path

    Returns:
        pd.DataFrame: Pivoted NTC data indexed by [strain, umi_seq]
    """
    df_ntc_reads = pd.DataFrame()

    for idx, row in df_ntc.iterrows():
        try:
            df = load_molecule_table(basepath, row['ill'], row['file_name'])
            df['sample'] = f"{row['Sample Name/Pool Name*']}_{row['sample_id']}"
            df_ntc_reads = pd.concat([df_ntc_reads, df])
        except FileNotFoundError:
            print(f"  WARNING: NTC file not found: {row['file_name']}")
            continue

    if df_ntc_reads.empty:
        print("  WARNING: No NTC data found, returning empty dataframe")
        return pd.DataFrame()

    # Pivot to wide format
    df_ntc_pivot = df_ntc_reads.pivot(
        index=['strain', 'umi_seq'],
        columns='sample',
        values='molecules'
    ).fillna(0)

    # Drop problematic columns if they exist
    cols_to_drop = [col for col in df_ntc_pivot.columns if 'WATER_D2_WATER' in col]
    if cols_to_drop:
        df_ntc_pivot.drop(columns=cols_to_drop, inplace=True)

    return df_ntc_pivot


def process_sample_group(df_samples, basepath, gavage_sample_type, project_name):
    """
    Load and pivot sample data for a group (gavage or experimental).

    Args:
        df_samples (pd.DataFrame): DataFrame with sample info
        basepath (str): Base directory path
        gavage_sample_type (str): Sample type identifier for gavage
        project_name (str): Project name for column naming

    Returns:
        pd.DataFrame: Pivoted data indexed by [strain, umi_seq]
    """
    df_reads = pd.DataFrame()

    for idx, row in df_samples.iterrows():
        try:
            df = load_molecule_table(basepath, row['ill'], row['file_name'])

            # Create sample column name
            if gavage_sample_type and row['sample_type'] == gavage_sample_type:
                df['sample'] = f"{row['sample_id']}_{row['sample_type']}"
            elif 'project' in row and pd.notna(row['project']):
                df['sample'] = f"{row['project']}_{row['file_name']}_{row['sample_type']}"
            else:
                df['sample'] = f"{row['sample_id']}_{row['sample_type']}"

            df_reads = pd.concat([df_reads, df])
        except FileNotFoundError:
            print(f"  WARNING: File not found: {row['file_name']}")
            continue

    if df_reads.empty:
        return pd.DataFrame()

    # Pivot to wide format
    df_pivot = df_reads.pivot(
        index=['strain', 'umi_seq'],
        columns='sample',
        values='molecules'
    ).fillna(0)

    return df_pivot


def create_and_export_tables(df_gavage_pivot, df_samples_pivot, df_ntc_pivot,
                             strains, output_dir, group_prefix=""):
    """
    Combine pivoted dataframes and export by strain.

    Args:
        df_gavage_pivot (pd.DataFrame): Pivoted gavage data
        df_samples_pivot (pd.DataFrame): Pivoted sample data
        df_ntc_pivot (pd.DataFrame): Pivoted NTC data
        strains (list): List of strain names
        output_dir (str): Output directory
        group_prefix (str): Prefix for output filenames (e.g., "c1_", "dose1_")
    """
    # Concatenate all dataframes
    dfs_to_concat = []

    if not df_gavage_pivot.empty:
        dfs_to_concat.append(df_gavage_pivot)
    if not df_samples_pivot.empty:
        dfs_to_concat.append(df_samples_pivot)
    if not df_ntc_pivot.empty:
        dfs_to_concat.append(df_ntc_pivot)

    if not dfs_to_concat:
        print("  ERROR: No data to export")
        return

    df_combined = pd.concat(dfs_to_concat, axis=1).fillna(0)

    # Clean column names
    df_combined.columns = df_combined.columns.str.replace(' ', '-')
    df_combined.reset_index(inplace=True)

    # Convert float columns to int
    float_cols = df_combined.select_dtypes(include=['float']).columns
    df_combined = df_combined.astype({col: 'int' for col in float_cols})

    # Export by strain
    os.makedirs(output_dir, exist_ok=True)

    for strain in strains:
        df_strain = df_combined[df_combined['strain'] == strain].drop('strain', axis=1)
        output_file = os.path.join(
            output_dir,
            f"{group_prefix}{strain}_moleculestable.csv"
        )
        df_strain.to_csv(output_file, index=False)
        print(f"  Exported: {output_file}")


def process_simple_project(df_merge, df_ntc, args):
    """
    Process projects without grouping (P4C2T4T5).

    Args:
        df_merge (pd.DataFrame): Merged sample metadata
        df_ntc (pd.DataFrame): NTC sample info
        args: Command-line arguments
    """
    print("\n=== Processing simple project (no grouping) ===")

    # Process NTC samples
    print("\nProcessing NTC samples...")
    df_ntc_pivot = process_ntc_samples(df_ntc, args.basepath)

    # Process gavage samples
    print("\nProcessing gavage samples...")
    gavage_type = args.gavage_sample_type or 'Gavage DNA'
    df_gavage = df_merge[df_merge['sample_type'] == gavage_type].copy()
    df_gavage_pivot = process_sample_group(
        df_gavage, args.basepath, gavage_type, args.project
    )

    # Process experimental samples
    print("\nProcessing experimental samples...")
    df_samples = df_merge[df_merge['sample_type'] != gavage_type].copy()
    df_samples_pivot = process_sample_group(
        df_samples, args.basepath, gavage_type, args.project
    )

    # Export tables
    print("\nExporting tables...")
    strains = args.strains.split(',')
    create_and_export_tables(
        df_gavage_pivot, df_samples_pivot, df_ntc_pivot,
        strains, args.output_dir
    )


def process_grouped_project(df_merge, df_ntc, args):
    """
    Process projects with grouping by cage_id or gavage_dose.

    Args:
        df_merge (pd.DataFrame): Merged sample metadata
        df_ntc (pd.DataFrame): NTC sample info
        args: Command-line arguments
    """
    print(f"\n=== Processing grouped project (grouping by {args.grouping_column}) ===")

    # Process NTC samples once
    print("\nProcessing NTC samples...")
    df_ntc_pivot = process_ntc_samples(df_ntc, args.basepath)

    # Get grouping values
    groups = df_merge[args.grouping_column].unique()
    print(f"\nFound {len(groups)} groups: {groups}")

    gavage_type = args.gavage_sample_type or 'Gavage DNA'
    strains = args.strains.split(',')

    # Process each group
    for group in groups:
        print(f"\n--- Processing group: {group} ---")
        df_group = df_merge[df_merge[args.grouping_column] == group].copy()

        # Process gavage for this group
        print(f"  Processing gavage samples...")
        df_gavage = df_group[df_group['sample_type'] == gavage_type].copy()
        df_gavage_pivot = process_sample_group(
            df_gavage, args.basepath, gavage_type, args.project
        )

        # Process experimental samples for this group
        print(f"  Processing experimental samples...")
        df_samples = df_group[df_group['sample_type'] != gavage_type].copy()
        df_samples_pivot = process_sample_group(
            df_samples, args.basepath, gavage_type, args.project
        )

        # Export tables with group prefix
        print(f"  Exporting tables...")
        create_and_export_tables(
            df_gavage_pivot, df_samples_pivot, df_ntc_pivot,
            strains, args.output_dir, group_prefix=f"{group}_"
        )


def process_mouse_grouped_project(df_merge, df_ntc, args):
    """
    Process P4C1T8 with mouse-level gavage grouping.

    Args:
        df_merge (pd.DataFrame): Merged sample metadata
        df_ntc (pd.DataFrame): NTC sample info
        args: Command-line arguments
    """
    print("\n=== Processing mouse-grouped project (P4C1T8) ===")

    # Process NTC samples once
    print("\nProcessing NTC samples...")
    df_ntc_pivot = process_ntc_samples(df_ntc, args.basepath)

    # Get cages
    cages = df_merge['cage_id'].unique()
    print(f"\nFound {len(cages)} cages: {cages}")

    gavage_type = args.gavage_sample_type or 'gavage DNA'

    # Define strain lists for different cages
    strains_mA = ['ST1', 'ST2']
    strains_mB = ['ST5', 'ST6']
    strains_c5 = ['ST1', 'ST2', 'ST4']

    # Process each cage
    for cage in cages:
        print(f"\n--- Processing cage: {cage} ---")
        df_cage = df_merge[df_merge['cage_id'] == cage].copy()

        if cage != 'c5':
            # Process mouse A and B separately
            df_mA_gavage = df_cage[
                (df_cage['mouse_id'] == 'mA') &
                (df_cage['sample_type'] == gavage_type)
            ]
            df_mB_gavage = df_cage[
                (df_cage['mouse_id'] == 'mB') &
                (df_cage['sample_type'] == gavage_type)
            ]

            print(f"  Processing mouse A gavage...")
            df_gavage_mA = process_sample_group(
                df_mA_gavage, args.basepath, gavage_type, args.project
            )

            print(f"  Processing mouse B gavage...")
            df_gavage_mB = process_sample_group(
                df_mB_gavage, args.basepath, gavage_type, args.project
            )

            # Process experimental samples (shared between mice)
            print(f"  Processing experimental samples...")
            df_samples = df_cage[df_cage['sample_type'] != gavage_type].copy()
            df_samples_pivot = process_sample_group(
                df_samples, args.basepath, gavage_type, args.project
            )

            # Export for mouse A strains
            print(f"  Exporting mouse A tables...")
            create_and_export_tables(
                df_gavage_mA, df_samples_pivot, df_ntc_pivot,
                strains_mA, args.output_dir, group_prefix=f"{cage}_"
            )

            # Export for mouse B strains
            print(f"  Exporting mouse B tables...")
            create_and_export_tables(
                df_gavage_mB, df_samples_pivot, df_ntc_pivot,
                strains_mB, args.output_dir, group_prefix=f"{cage}_"
            )
        else:
            # Cage c5 - single gavage
            print(f"  Processing cage c5 gavage...")
            df_gavage = df_cage[df_cage['sample_type'] == gavage_type].copy()
            df_gavage_pivot = process_sample_group(
                df_gavage, args.basepath, gavage_type, args.project
            )

            print(f"  Processing experimental samples...")
            df_samples = df_cage[df_cage['sample_type'] != gavage_type].copy()
            df_samples_pivot = process_sample_group(
                df_samples, args.basepath, gavage_type, args.project
            )

            print(f"  Exporting tables...")
            create_and_export_tables(
                df_gavage_pivot, df_samples_pivot, df_ntc_pivot,
                strains_c5, args.output_dir, group_prefix=f"{cage}_"
            )


def process_c1gavage_mode(df_merge, df_ntc, args):
    """
    Process P4C1T8 with only c1 cage as gavage source.

    Args:
        df_merge (pd.DataFrame): Merged sample metadata
        df_ntc (pd.DataFrame): NTC sample info
        args: Command-line arguments
    """
    print("\n=== Processing c1gavage special mode ===")

    # Process NTC samples once
    print("\nProcessing NTC samples...")
    df_ntc_pivot = process_ntc_samples(df_ntc, args.basepath)

    gavage_type = args.gavage_sample_type or 'gavage DNA'

    # Get c1 gavage only
    print("\nProcessing c1 gavage samples...")
    df_mA_gavage = df_merge[
        (df_merge['mouse_id'] == 'mA') &
        (df_merge['sample_type'] == gavage_type) &
        (df_merge['cage_id'] == 'c1')
    ]
    df_mB_gavage = df_merge[
        (df_merge['mouse_id'] == 'mB') &
        (df_merge['sample_type'] == gavage_type) &
        (df_merge['cage_id'] == 'c1')
    ]

    df_gavage_mA = process_sample_group(
        df_mA_gavage, args.basepath, gavage_type, args.project
    )
    df_gavage_mB = process_sample_group(
        df_mB_gavage, args.basepath, gavage_type, args.project
    )

    # Define strain lists
    strains_mA = ['ST1', 'ST2']
    strains_mB = ['ST5', 'ST6']

    # Process each cage with c1 gavage
    cages = [c for c in df_merge['cage_id'].unique() if c != 'c5']

    for cage in cages:
        print(f"\n--- Processing cage: {cage} (with c1 gavage) ---")
        df_cage = df_merge[df_merge['cage_id'] == cage].copy()

        # Process experimental samples only (no gavage)
        df_samples = df_cage[df_cage['sample_type'] != gavage_type].copy()
        df_samples_pivot = process_sample_group(
            df_samples, args.basepath, gavage_type, args.project
        )

        # Export with c1 gavage
        print(f"  Exporting mouse A tables...")
        create_and_export_tables(
            df_gavage_mA, df_samples_pivot, df_ntc_pivot,
            strains_mA, args.output_dir, group_prefix=f"{cage}_"
        )

        print(f"  Exporting mouse B tables...")
        create_and_export_tables(
            df_gavage_mB, df_samples_pivot, df_ntc_pivot,
            strains_mB, args.output_dir, group_prefix=f"{cage}_"
        )


def process_invitro_mode(df_merge, df_ntc, args):
    """
    Process P4C4T1 in vitro experiment (no gavage, imported later).

    Args:
        df_merge (pd.DataFrame): Merged sample metadata
        df_ntc (pd.DataFrame): NTC sample info
        args: Command-line arguments
    """
    print("\n=== Processing in vitro special mode ===")

    # Process samples without gavage
    print("\nProcessing in vitro samples...")
    df_samples_pivot = process_sample_group(
        df_merge, args.basepath, None, args.project
    )

    # Export by strain (no gavage, no NTC)
    strains = args.strains.split(',')
    empty_df = pd.DataFrame()

    print("\nExporting tables...")
    create_and_export_tables(
        empty_df, df_samples_pivot, empty_df,
        strains, args.output_dir
    )

    # Import gavage from another project if specified
    if args.gavage_project_path:
        print(f"\n=== Importing gavage from {args.gavage_project_path} ===")

        for strain in strains:
            # Read invitro table
            invitro_file = os.path.join(
                args.output_dir,
                f"{strain}_moleculestable.csv"
            )
            df_invitro = pd.read_csv(invitro_file)
            df_invitro.set_index('umi_seq', inplace=True)

            # Read gavage table
            gavage_file = os.path.join(
                args.gavage_project_path,
                f"{strain}_onlygavage_moleculestable.csv"
            )

            if not os.path.exists(gavage_file):
                print(f"  WARNING: Gavage file not found: {gavage_file}")
                continue

            df_gavage = pd.read_csv(gavage_file)
            df_gavage.set_index('umi_seq', inplace=True)

            # Get gavage columns
            gavage_cols = [col for col in df_gavage.columns if 'Gavage' in col]

            # Combine
            df_combined = pd.concat([df_gavage[gavage_cols], df_invitro], axis=1)
            df_combined.reset_index(inplace=True)
            df_combined.fillna(0, inplace=True)

            # Filter for rows with gavage present
            df_combined = df_combined[df_combined[gavage_cols].sum(axis=1) > 0]

            # Save
            output_file = os.path.join(
                args.output_dir,
                f"{strain}_onlygavage_moleculestable.csv"
            )
            df_combined.to_csv(output_file, index=False)
            print(f"  Exported: {output_file}")


def main():
    """Main execution function."""
    args = parse_arguments()

    print("=" * 80)
    print(f"Create Molecule Tables for getFP Analysis")
    print(f"Project: {args.project}")
    print("=" * 80)

    # Parse comma-separated inputs
    ill_sheets = [s.strip() for s in args.ill_sheets.split(',')]
    metadata_files = [s.strip() for s in args.metadata.split(',')]

    # Load and merge data
    print("\n=== Loading data ===")
    df_merge, df_ntc = load_ill_metadata(ill_sheets, metadata_files, args.project)

    # Route to appropriate processing function
    if args.special_mode == 'c1gavage':
        process_c1gavage_mode(df_merge, df_ntc, args)
    elif args.special_mode == 'invitro':
        process_invitro_mode(df_merge, df_ntc, args)
    elif args.project == 'P4C1T8':
        process_mouse_grouped_project(df_merge, df_ntc, args)
    elif args.grouping_column:
        process_grouped_project(df_merge, df_ntc, args)
    else:
        process_simple_project(df_merge, df_ntc, args)

    # Apply gavage filtering
    if args.special_mode != 'invitro' or not args.gavage_project_path:
        print("\n=== Applying gavage filtering ===")
        gavage_only_readtables(args.output_dir, args.output_dir)

    print("\n" + "=" * 80)
    print("Processing complete!")
    print(f"Output directory: {args.output_dir}")
    print("=" * 80)


if __name__ == "__main__":
    main()
