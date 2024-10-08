# metadomain_peak_caller/plotting.py
import argparse
import logging
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cooler
import glob

def parse_args():
    parser = argparse.ArgumentParser(description='Plot Top Hits of Metadomain Peak Caller')
    parser.add_argument('--output_dir', required=True, help='Directory containing output files and log file')
    parser.add_argument('--top_n', type=int, default=5, help='Number of top metadomains to plot')
    parser.add_argument('--intermediate_dir', default=None, help='Directory containing intermediate files')
    parser.add_argument('--save_plots', action='store_true', help='Save plots to files instead of showing')
    args = parser.parse_args()
    return args

def main_plot_cli():
    args = parse_args()
    output_dir = args.output_dir
    top_n = args.top_n
    intermediate_dir = args.intermediate_dir or os.path.join(output_dir, 'intermediate_files/')
    save_plots = args.save_plots

    log_file = os.path.join(output_dir, 'InterDomain.log')
    results_file_pattern = os.path.join(output_dir, 'intra_metadomains.*.tsv')
    results_files = glob.glob(results_file_pattern)
    if len(results_files) == 0:
        print(f"No results files found in {output_dir}")
        return
    results_file = results_files[0]

    # Read the log file to get the necessary parameters
    cool_file = None
    resolution = None
    filter_n = None
    filter_width = None
    label = None
    with open(log_file, 'r') as f:
        for line in f:
            if 'Cool file:' in line:
                cool_file = line.strip().split('Cool file:')[1].strip()
            elif 'Resolution:' in line:
                resolution = int(line.strip().split('Resolution:')[1].strip())
            elif 'filter_n:' in line:
                filter_n = int(line.strip().split('filter_n:')[1].strip())
            elif 'filter_width:' in line:
                filter_width = int(line.strip().split('filter_width:')[1].strip())
            elif 'Label:' in line:
                label = line.strip().split('Label:')[1].strip()

    if cool_file is None or resolution is None or filter_n is None or filter_width is None:
        print("Could not find necessary parameters in log file.")
        return

    # Load the cool file
    cool = cooler.Cooler(cool_file)

    # Read the results file
    results = pd.read_csv(results_file, sep='\t')
    if 'log_p_value' not in results.columns:
        print("The results file does not contain 'log_p_value'. Cannot determine top hits.")
        return

    # Sort by log_p_value (higher log_p_value means more significant)
    results = results.sort_values(by='log_p_value', ascending=False)
    top_results = results.head(top_n)

    # Plot the top N metadomains
    for idx, row in top_results.iterrows():
        chrom = row['chrom']
        start = int(row['start'])
        end = int(row['end'])
        log_p_value = row['log_p_value']

        # Fetch the matrix around the metadomain
        # Adjust the window size as needed
        window_size = filter_n * 2  # Extend beyond the outside filter
        bin_start = max(0, min(start, end) - window_size)
        bin_end = max(start, end) + window_size

        region = (chrom, bin_start * resolution, bin_end * resolution)
        mat = cool.matrix(balance=True).fetch(region, region)

        plt.figure(figsize=(6,6))
        plt.matshow(np.log1p(mat), fignum=False, cmap='coolwarm')
        plt.title(f'Chromosome {chrom}, Start {start}, End {end}, Log_p_value {log_p_value}')
        ax = plt.gca()

        # Calculate positions for the filters
        center_x = (start - bin_start)  # Adjusted for the fetched region
        center_y = (end - bin_start)

        # Inside filter rectangle
        inside_size = filter_width
        rect_inside = plt.Rectangle(
            (center_y - inside_size//2, center_x - inside_size//2),
            inside_size,
            inside_size,
            edgecolor='green',
            facecolor='none',
            linewidth=2,
            label='Inside Filter'
        )
        ax.add_patch(rect_inside)

        # Outside filter rectangle
        outside_size = filter_n
        rect_outside = plt.Rectangle(
            (center_y - outside_size//2, center_x - outside_size//2),
            outside_size,
            outside_size,
            edgecolor='red',
            facecolor='none',
            linewidth=2,
            label='Outside Filter'
        )
        ax.add_patch(rect_outside)

        plt.legend(handles=[rect_inside, rect_outside], loc='upper right')

        if save_plots:
            plot_file = os.path.join(output_dir, f'metadomain_{chrom}_{start}_{end}.png')
            plt.savefig(plot_file)
            print(f"Plot saved to {plot_file}")
        else:
            plt.show()
        plt.close()

        # If intermediate files are saved, plot them
        if os.path.exists(intermediate_dir):
            # Load intermediate matrices if they exist
            prefix = f"{label}_chrL={chrom}_chrR={chrom}_"
            try:
                peak_smooth_X1 = np.load(os.path.join(intermediate_dir, f"{prefix}peak_smooth_X1.npy"))
                peak_smooth_Y1 = np.load(os.path.join(intermediate_dir, f"{prefix}peak_smooth_Y1.npy"))
                full_logp_mat = np.load(os.path.join(intermediate_dir, f"{prefix}full_logp_mat.npy"))
                collapsed_logp_mat = np.load(os.path.join(intermediate_dir, f"{prefix}collapsed_logp_mat.npy"))
                bal_intra = np.load(os.path.join(intermediate_dir, f"{prefix}bal_intra.npy"))
                raw_intra = np.load(os.path.join(intermediate_dir, f"{prefix}raw_intra.npy"))
                oe = np.load(os.path.join(intermediate_dir, f"{prefix}oe.npy"))

                # Plot each intermediate matrix
                intermediate_matrices = {
                    'peak_smooth_X1': peak_smooth_X1,
                    'peak_smooth_Y1': peak_smooth_Y1,
                    'full_logp_mat': full_logp_mat,
                    'collapsed_logp_mat': collapsed_logp_mat,
                    'bal_intra': bal_intra,
                    'raw_intra': raw_intra,
                    'oe': oe
                }

                fig, axes = plt.subplots(1, len(intermediate_matrices), figsize=(20, 5))
                for ax, (name, matrix) in zip(axes, intermediate_matrices.items()):
                    sub_matrix = matrix[bin_start:bin_end, bin_start:bin_end]
                    im = ax.matshow(np.log1p(sub_matrix), cmap='coolwarm')
                    ax.set_title(name)
                    fig.colorbar(im, ax=ax)

                if save_plots:
                    plot_file = os.path.join(output_dir, f'metadomain_{chrom}_{start}_{end}_intermediates.png')
                    plt.savefig(plot_file)
                    print(f"Intermediate plots saved to {plot_file}")
                else:
                    plt.show()
                plt.close()

            except FileNotFoundError as e:
                print(f"Intermediate file not found: {e}")
                continue

if __name__ == '__main__':
    main_plot_cli()
