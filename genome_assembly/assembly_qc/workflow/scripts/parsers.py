import re
import numpy as np
import pandas as pd
from scipy.signal import find_peaks, peak_widths
import portion as P
import argparse
import plotly.graph_objects as go

############################
###                      ###
### FUNCTIONS FOR KMERS  ###
###                      ###
############################


def calculate_overlap_ratios(df, interval1, interval2):
    """
    Calculates overlap ratios between 'read-only' kmers and assembly kmers.

    Parameters:
        df (pd.DataFrame): Kmer histogram DataFrame with 'Assembly', 'kmer_multiplicity', and 'Count' columns.
        interval1 (list): Intervals representing 'read-only' kmer peaks (start, end).
        interval2 (portion.Interval): Interval representing a specific peak in the assembly.

    Returns:
        float: Average overlap ratio for the given intervals.
    """
    overlap_ratios = []

    for i, read_interv in enumerate(interval1):
        # Calculate overlap between intervals
        overlap = interval2 & P.open(read_interv[0], read_interv[1])

        if not overlap.empty:
            # Sum counts for the overlapping kmers in the assembly
            assembly = df[
                (df["Assembly"] != "read-only")
                & (df["kmer_multiplicity"].isin(list(P.iterate(overlap, step=1))))
            ]
            sum_assembly = assembly["Count"].sum()

            # Sum counts for overlapping kmers in the read-only data
            read = df[
                (df["Assembly"] == "read-only")
                & (df["kmer_multiplicity"].isin(list(P.iterate(overlap, step=1))))
            ]
            sum_read = read["Count"].sum()

            # Calculate ratio of read to assembly counts
            overlap_ratio = sum_read / sum_assembly if sum_assembly != 0 else 0
            overlap_ratios.append(overlap_ratio)
        else:
            overlap_ratios.append(0)  # No overlap for this interval

    # Return average overlap ratio
    return sum(overlap_ratios) / len(overlap_ratios) if overlap_ratios else 0


def calculate_peak_ratio(first_peak_count, second_peak_count):
    """
    Computes the ratio of the second peak height to the first peak height.

    Parameters:
        first_peak_count (float): Count at the first peak.
        second_peak_count (float): Count at the second peak.

    Returns:
        float: Ratio of the second peak to the first.
    """
    return second_peak_count / first_peak_count


def analyze_peaks(df, estimated_coverage):
    """
    Identifies key peaks in a kmer histogram and validates their properties.

    Parameters:
        df (pd.DataFrame): Histogram with 'Assembly', 'kmer_multiplicity', and 'Count' columns.
        estimated_coverage (float): Coverage estimate to locate homozygous peaks.

    Returns:
        tuple: Analysis results including overlap ratio, presence of peaks, and peak ratio.
    """
    # Extract kmers for analysis
    non_read_only = df[df["Assembly"] != "read-only"]
    read_only = df[df["Assembly"] == "read-only"]

    kmer_multiplicity = non_read_only["kmer_multiplicity"].values[:100]
    counts = non_read_only["Count"].values[:100]
    read_only_counts = read_only["Count"].astype(int).values[:100]

    # Detect peaks in the assembly and read-only data
    peaks, _ = find_peaks(counts)
    read_only_peaks, _ = find_peaks(read_only_counts)
    read_only_widths = peak_widths(read_only_counts, read_only_peaks)

    # Construct intervals for read-only peaks
    start_intervals = np.vectorize(round)(read_only_widths[2])
    end_intervals = np.vectorize(round)(read_only_widths[3])
    read_only_intervals = [
        (start, end)
        for start, end in zip(start_intervals, end_intervals)
        if start != end
    ]

    # Identify homozygous peak based on coverage z-score
    mu, std_dev = (
        estimated_coverage,
        np.std(np.random.normal(estimated_coverage, 1, size=1000)),
    )
    second_peak_idx = next(
        (
            i
            for i, mult in enumerate(kmer_multiplicity[peaks])
            if abs((mult - mu) / std_dev) <= 2
        ),
        None,
    )

    if second_peak_idx is None:
        raise ValueError("Could not find a homozygous peak near the coverage estimate.")

    # Peak details
    second_peak_count = counts[peaks[second_peak_idx]]
    second_peak_interval = P.open(
        *np.vectorize(round)(peak_widths(counts, [peaks[second_peak_idx]])[2:4])
    )

    # Heterozygous peak analysis
    if second_peak_idx > 0:
        first_peak_idx = second_peak_idx - 1
        first_peak_count = counts[peaks[first_peak_idx]]
        first_peak_interval = P.open(
            *np.vectorize(round)(peak_widths(counts, [peaks[first_peak_idx]])[2:4])
        )
        homozygous_peak, heterozygous_peak = 1, 1
        peak_ratio = calculate_peak_ratio(first_peak_count, second_peak_count)
    else:
        first_peak_interval = None
        homozygous_peak, heterozygous_peak = 1, 0
        peak_ratio = 1

    # Overlap ratios
    overlap_ratio_second_peak = calculate_overlap_ratios(
        df, read_only_intervals, second_peak_interval
    )
    overlap_ratio_first_peak = (
        calculate_overlap_ratios(df, read_only_intervals, first_peak_interval)
        if first_peak_interval
        else 0
    )
    overlap_ratio = (overlap_ratio_second_peak + overlap_ratio_first_peak) / 2

    return overlap_ratio, homozygous_peak, heterozygous_peak, peak_ratio


def kmer_extractor(file_path, estimated_coverage):
    """
    Main function to extract kmer profile metrics from a histogram file.

    Parameters:
        file_path (str): Path to the kmer histogram file.
        estimated_coverage (float): Estimated genome coverage.

    Returns:
        dict: Metrics including overlap ratio, peak presence, and peak ratios.
    """
    kmer = pd.read_csv(file_path, sep="\t")

    try:
        # Analyze peaks in the kmer data
        overlap_ratio, homozygous_peak, heterozygous_peak, peak_ratio = analyze_peaks(
            kmer, estimated_coverage
        )
    except ValueError as e:
        # Handle missing peaks
        return {"Error": str(e)}

    return {
        "Missing read kmer from assembly": overlap_ratio,
        "Homozygous peak found": homozygous_peak,
        "Heterozygous peak found": heterozygous_peak,
        "Peak ratio": peak_ratio,
    }


##################################
###                            ###
### PARSERS OF DIFFERENT TOOLS ###
###                            ###
##################################


# Extract BUSCO scores
def parse_busco_report(file_path):
    busco_metrics = {}

    # Define regex patterns to match each metric
    patterns = {
        "Complete BUSCOs": r"(\d+)\s+Complete BUSCOs \(C\)",
        "Complete single-copy BUSCOs": r"(\d+)\s+Complete and single-copy BUSCOs \(S\)",
        "Complete duplicated BUSCOs": r"(\d+)\s+Complete and duplicated BUSCOs \(D\)",
        "Fragmented BUSCOs": r"(\d+)\s+Fragmented BUSCOs \(F\)",
        "Missing BUSCOs": r"(\d+)\s+Missing BUSCOs \(M\)",
    }

    # Process each file found
    with open(file_path) as f:
        for line in f:
            # Check each pattern in the line
            for metric, pattern in patterns.items():
                match = re.search(pattern, line)
                if match:
                    # Store the metric as an integer in the dictionary
                    busco_metrics[metric] = int(match.group(1))

    return busco_metrics


# Parse Inspector output summary statistics
def parse_summary_statistics(file_path, estimated_genome_size):
    stats = {
        "Total length": None,
        "Number of contigs": None,
        "Longest contig": None,
        "N50": None,
        "Mapping rate %": None,
        "Depth": None,
        "Mapping rate in large contigs %": None,
        "Expansion": None,
        "Collapse": None,
        "Haplotype switch": None,
        "Inversion": None,
        "Small-scale assembly error per Mbp": None,
    }

    patterns = {
        "Total length": r"Total length\s+(\d+)",
        "Number of contigs": r"Number of contigs\s+(\d+)",
        "Longest contig": r"Longest contig\s+(\d+)",
        "N50": r"N50\s+(\d+)",
        "Mapping rate %": r"Mapping rate /%\s+([\d.]+)",
        "Depth": r"Depth\s+([\d.]+)",
        "Mapping rate in large contigs %": r"Mapping rate in large contigs /%\s+([\d.]+)",
        "Expansion": r"Expansion\s+(\d+)",
        "Collapse": r"Collapse\s+(\d+)",
        "Haplotype switch": r"Haplotype switch\s+(\d+)",
        "Inversion": r"Inversion\s+(\d+)",
        "Small-scale assembly error per Mbp": r"Small-scale assembly error /per Mbp\s+([\d.]+)",
    }

    with open(file_path) as f:
        for line in f:
            for metric, pattern in patterns.items():
                match = re.search(pattern, line)
                if match:
                    value = (
                        float(match.group(1))
                        if "." in match.group(1)
                        else int(match.group(1))
                    )
                    stats[metric] = value

    # Calculate how close the total length is to the estimated genome size
    if stats["Total length"] is not None:
        stats["Difference from estimated genome size"] = abs(
            stats["Total length"] - estimated_genome_size
        )

    return stats


# Extract QV values from Merqury run
def parse_merqury_qv(file_path):
    qv_value = None

    with open(file_path) as f:
        for line in f:
            columns = line.split()
            if len(columns) >= 4:
                try:
                    qv_value = float(columns[3])
                    break
                except ValueError:
                    pass  # In case the 4th column is not a number

    return qv_value


# We will introduce a function to find the priority index.
# Then we provide the attributes data to this function.
# Source: https://www.analyticsvidhya.com/blog/2023/05/multi-criteria-decision-making-using-ahp-in-python/
def ahp_attributes(ahp_df):
    # Creating an array of sum of values in each column
    sum_array = np.array(ahp_df.sum(numeric_only=True))
    # Creating a normalized pairwise comparison matrix.
    # By dividing each column cell value with the sum of the respective column.
    cell_by_sum = ahp_df.div(sum_array, axis=1)
    # Creating Priority index by taking avg of each row
    priority_df = pd.DataFrame(
        cell_by_sum.mean(axis=1), index=ahp_df.index, columns=["priority index"]
    )
    priority_df = priority_df.transpose()
    return priority_df


def consistency_ratio(priority_index, ahp_df):
    random_matrix = {
        1: 0,
        2: 0,
        3: 0.58,
        4: 0.9,
        5: 1.12,
        6: 1.24,
        7: 1.32,
        8: 1.14,
        9: 1.45,
        10: 1.49,
        11: 1.51,
        12: 1.48,
        13: 1.56,
        14: 1.57,
        15: 1.59,
        16: 1.605,
        17: 1.61,
        18: 1.615,
        19: 1.62,
        20: 1.625,
    }
    # Check for consistency
    consistency_df = ahp_df.multiply(
        np.array(priority_index.loc["priority index"]), axis=1
    )
    consistency_df["sum_of_col"] = consistency_df.sum(axis=1)
    # To find lambda max
    lambda_max_df = consistency_df["sum_of_col"].div(
        np.array(priority_index.transpose()["priority index"]), axis=0
    )
    lambda_max = lambda_max_df.mean()
    # To find the consistency index
    consistency_index = round(
        (lambda_max - len(ahp_df.index)) / (len(ahp_df.index) - 1), 3
    )
    print(f"The Consistency Index is: {consistency_index}")
    # To find the consistency ratio
    consistency_ratio = round(consistency_index / random_matrix[len(ahp_df.index)], 3)
    print(f"The Consistency Ratio is: {consistency_ratio}")
    if consistency_ratio < 0.1:
        print("The model is consistent")
    else:
        print("The model is not consistent")


###################################
###                             ###
### FUNCTION FOR SANKEY DIAGRAM ###
###                             ###
###################################

def parse_args():
    parser = argparse.ArgumentParser(description="Collect assembly statistics and generate BUSCO Sankey diagram")
    parser.add_argument("--assemblies", required=True, help="Comma-separated list of assembly names")
    return parser.parse_args()


import pandas as pd
import plotly.graph_objects as go
import plotly.colors as pc

def create_busco_sankey(full_table_paths, output_file="assembly_evaluation/busco_sankey.html"):
    """
    Generates an interactive Sankey diagram showing shared BUSCO gene statuses across assemblies.

    Args:
        full_table_paths (dict): Dictionary of assembly names to paths to the full_table.tsv file.
        output_file (str): Path to save the interactive HTML output.
    """

    assembly_data = {}
    all_busco_ids = set()

    # Parse each full_table.tsv
    for assembly, path in full_table_paths.items():
        try:
            df = pd.read_csv(path, comment="#", sep="\t", usecols=[0, 1], names=['Busco id', 'Status'])
            if 'Busco id' not in df.columns or 'Status' not in df.columns:
                print(f"Warning: Required columns missing in {path}. Skipping.")
                continue
            assembly_data[assembly] = df
            all_busco_ids.update(df['Busco id'].tolist())
        except FileNotFoundError:
            print(f"Error: File not found at {path}. Skipping.")
            continue
        except pd.errors.ParserError as e:
            print(f"Error: Could not parse {path}: {e}. Skipping.")
            continue
        except Exception as e:
            print(f"Error: Could not process {path}: {e}. Skipping.")
            continue

    all_busco_ids = sorted(list(all_busco_ids))
    assemblies = sorted(list(full_table_paths.keys()))

    # Create a dictionary to store shared BUSCO statuses between assemblies
    shared_statuses = {(a1, a2): {'Complete': 0, 'Duplicated': 0, 'Fragmented': 0, 'Missing': 0}
                       for a1 in assemblies for a2 in assemblies if a1 < a2}

    # Compare statuses for each BUSCO ID across assemblies
    for busco_id in all_busco_ids:
        statuses_per_assembly = {}
        for assembly in assemblies:
            if assembly in assembly_data and busco_id in assembly_data[assembly]['Busco id'].values:
                status = assembly_data[assembly].loc[assembly_data[assembly]['Busco id'] == busco_id, 'Status'].iloc[0]
                statuses_per_assembly[assembly] = status

        # Check pairwise combinations of assemblies for shared statuses
        for a1, a2 in shared_statuses.keys():
            if a1 in statuses_per_assembly and a2 in statuses_per_assembly:
                if statuses_per_assembly[a1] == statuses_per_assembly[a2]:  # Same status for this BUSCO ID
                    shared_statuses[(a1, a2)][statuses_per_assembly[a1]] += 1

    # Create source, target, value lists for the Sankey diagram
    source_indices = []
    target_indices = []
    values = []
    colors = []

    # Generate colors for BUSCO categories
    categories = ['Complete', 'Duplicated', 'Fragmented', 'Missing']
    color_scale = pc.sequential.Viridis
    category_colors = {cat: color_scale[i * (len(color_scale) // len(categories))] for i, cat in enumerate(categories)}

    # Populate Sankey data based on shared statuses
    for (a1, a2), status_counts in shared_statuses.items():
        for status, count in status_counts.items():
            if count > 0:  # Only include links with shared statuses
                source_indices.append(assemblies.index(a1))
                target_indices.append(assemblies.index(a2))
                values.append(count)
                colors.append(category_colors[status])

    print(shared_statuses.items())

    # Create the Sankey diagram
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=assemblies,
            color="lightblue"
        ),
        link=dict(
            source=source_indices,
            target=target_indices,
            value=values,
            color=colors
        ))])

    # Add color legend
    for category, color in category_colors.items():
        fig.add_trace(go.Scatter(
            x=[None],
            y=[None],
            mode='markers',
            marker=dict(size=10, color=color),
            showlegend=True,
            name=category
        ))

    fig.update_layout(
        title_text="Shared BUSCO Gene Statuses Across Assemblies",
        font_size=10,
        legend_title_text='BUSCO Categories',
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        )
    )

    fig.write_html(output_file)
    print(f"Sankey diagram saved to {output_file}")