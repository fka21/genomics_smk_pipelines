import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt 
from pathlib import Path
from glob import glob
from scipy.stats import zscore
from scipy.spatial.distance import pdist
import os
import argparse
import numpy as np
from parsers import *  # noqa: E402, F403

def parse_args():
    parser = argparse.ArgumentParser(description='Assembly data collection script')
    parser.add_argument('--assemblies', required=True, help='List of assemblies')
    parser.add_argument('--genome_size', type=int, required=True, help='Genome size in base pairs')
    parser.add_argument('--coverage', type=int, required=True, help='Estimated coverage')  # Note the parameter name
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    assemblies = args.assemblies.split(' ')
    genome_size = args.genome_size
    estimated_coverage = args.coverage

    # Build the expected file paths
    file_paths = {s : os.path.join("busco_report",s,"run_metazoa_odb10","full_table.tsv") for s in assemblies}

    # -------------------------------------------------------
    # BEGIN Assembly Statistics Collection and Plot Generation
    # -------------------------------------------------------

    # Initialize list to hold results
    results = []

    for sample in assemblies:
        # Define paths for each sample
        busco_path = glob(os.path.join(f"busco_report/{sample}/short_summary.specific.metazoa_odb10.{sample}.txt"))
        kmer_profile_path = glob(os.path.join(f"kmer_report/{sample}/{sample}.spectra-asm.hist"))
        kmer_coverage_path = glob(os.path.join(f"kmer_report/{sample}/{sample}.qv"))
        summary_stats_path = glob(os.path.join(f"inspector_report/{sample}/summary_statistics"))

        # Check if required files are found for each sample
        if busco_path and kmer_coverage_path and summary_stats_path:
            busco_score = parse_busco_report(busco_path[0])
            kmer_coverage = parse_merqury_qv(kmer_coverage_path[0])
            summary_statistics = parse_summary_statistics(summary_stats_path[0], estimated_genome_size=genome_size)
            #kmer_metrics = kmer_extractor(kmer_profile_path[0], estimated_coverage=estimated_coverage)

            # Aggregate results into a dictionary
            result = {
                "Sample": sample,
                "BUSCO Score": busco_score,
                #"Kmer Coverage": kmer_coverage,
                "Summary Statistics": summary_statistics,
            }

            # Add detailed metrics to the result dictionary
            result.update(busco_score)
            result["Kmer Coverage"] = kmer_coverage
            result.update(summary_statistics)
            #result.update(kmer_metrics)

            # Append to results list
            results.append(result)
        else:
            print(f"Missing required files for sample: {sample}")

    # Convert list of results to DataFrame
    results_df = pd.DataFrame(results)

    # Dropping unused columns
    if not results_df.empty and all(col in results_df.columns for col in ["BUSCO Score", "Summary Statistics", "Total length"]):
        results_df = results_df.drop(["BUSCO Score", "Summary Statistics", "Total length"], axis=1, errors='ignore')

    # Convert list of results to DataFrame
    results_df = pd.DataFrame(results_df)
    
    # Define a custom color palette
    cmap = sns.diverging_palette(220, 20, as_cmap=True)

    # Drop columns with zero variance
    non_zero_variance_cols = results_df.set_index("Sample").var() != 0
    filtered_df = results_df.set_index("Sample").loc[:, non_zero_variance_cols]

    # Apply Z-score normalization
    z_transformed = filtered_df.apply(zscore, axis=0).T

    # Plot data
    clustermap = sns.clustermap(
    z_transformed,
    cmap=cmap,
    figsize=(12, 8),  # Set figure size
    dendrogram_ratio=(0.1, 0.1),  # Adjust dendrogram size ratio
    cbar_pos=(-0.05, 0.6, 0.03, 0.3),  # Legend position
    center=0,
    )

    # Reduce the size of the color bar's tick labels
    clustermap.ax_cbar.set_title("Z-Score", fontsize=8)
    clustermap.ax_cbar.tick_params(labelsize=6)

    # Adjust layout for better aesthetics
    plt.setp(clustermap.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    plt.setp(clustermap.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)

    # Save the heatmap as a high-resolution image
    output_path = "assembly_evaluation/clustermap_heatmap.pdf"
    clustermap.savefig(output_path, dpi=300, format="pdf")

    # Save the DataFrame
    if not results_df.empty:
        results_df.to_csv("assembly_evaluation/assembly_metrics_summary.csv", index=False)
    else:
        print("Warning: No data to save to CSV.")

    # -------------------------------------------------------
    # END Assembly Statistics Collection and Plot Generation
    # -------------------------------------------------------

    # Generate BUSCO Sankey diagram
    create_busco_sankey(file_paths, "assembly_evaluation/busco_sankey.html")

    print("Assembly statistics collection, interactive heatmap, radar plot, BUSCO Sankey diagram, and CSV file have been generated.")