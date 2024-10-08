# InterDomain/plotting.py

import argparse
import logging
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cooler
import glob

def parse_args():
    parser = argparse.ArgumentParser(description='Plot Top Hits of InterDomain Peak Caller')
    parser.add_argument('--output_dir', default='bedfile_output/', help='Directory containing output files and log file')
    parser.add_argument('--top_n', type=int, default=5, help='Number of top metadomains to plot')
    parser.add_argument('--intermediate_dir', default=None, help='Directory containing intermediate files')
    parser.add_argument('--save_plots', default=True, action='store_true', help='Save plots to files instead of showing')
    parser.add_argument('--type', choices=['intra', 'inter'], required=True, help='Type of metadomains to plot: "intra" or "inter"')
    args = parser.parse_args()
    return args

def main_plot_cli():
    args = parse_args()
    output_dir = args.output_dir
    top_n = args.top_n
    plot_type = args.type  # 'intra' or 'inter'
    intermediate_dir = args.intermediate_dir
    save_plots = args.save_plots

    # Set default intermediate directory if not provided
    if intermediate_dir is None:
        intermediate_dir = os.path.join(output_dir, 'intermediate_files/')

    plots_dir = os.path.join(output_dir, 'plots')
    os.makedirs(plots_dir, exist_ok=True)  # Create the directory if it doesn't exist

    # Adjusted to handle specified type of metadomains
    results_files = []
    results_dir = os.path.join(output_dir, plot_type)
    if os.path.exists(results_dir):
        pattern = f'{plot_type}_metadomains.*.tsv'
        results_files.extend(glob.glob(os.path.join(results_dir, pattern)))

    if len(results_files) == 0:
        print(f"No results files found for the specified type '{plot_type}' in {output_dir}")
        return

    # Concatenate all results files
    results = pd.concat([pd.read_csv(f, sep='\t') for f in results_files], ignore_index=True)

    # Read the log file to get the necessary parameters
    log_file = os.path.join(output_dir, plot_type, 'InterDomain.log')
    if not os.path.exists(log_file):
        print(f"No log file found for the specified type '{plot_type}' in {output_dir}")
        return

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

    if cool_file is None or resolution is None or filter_n is None or filter_width is None or label is None:
        print("Could not find necessary parameters in log file.")
        return

    # Load the cool file
    cool = cooler.Cooler(cool_file)

    if 'log_p_value' not in results.columns:
        print("The results file does not contain 'log_p_value'. Cannot determine top hits.")
        return

    # Sort by log_p_value (higher log_p_value means more significant)
    results = results.sort_values(by='log_p_value', ascending=False)
    top_results = results.head(top_n)

    # Plot the top N metadomains
    for idx, row in top_results.iterrows():
        # Determine if it's intra- or inter-chromosomal based on plot_type
        if plot_type == 'inter':
            # Inter-chromosomal metadomain
            chrom1 = row['chrom1']
            start1 = int(row['start1'])
            end1 = int(row['end1'])
            chrom2 = row['chrom2']
            start2 = int(row['start2'])
            end2 = int(row['end2'])
        else:
            # Intra-chromosomal metadomain
            chrom1 = chrom2 = row['chrom']
            start1 = int(row['start'])
            end1 = int(row['start'])
            start2 = int(row['end'])
            end2 = int(row['end'])

        log_p_value = row['log_p_value']

        # Fetch the interaction matrix between region_L and region_R
        window_size = filter_n * 2  # Extend beyond the outside filter

        # Define the left and right regions
        bin_start_L = max(0, start1 - window_size)
        bin_end_L = end1 + window_size
        region_L = (chrom1, bin_start_L * resolution, bin_end_L * resolution)

        bin_start_R = max(0, start2 - window_size)
        bin_end_R = end2 + window_size
        region_R = (chrom2, bin_start_R * resolution, bin_end_R * resolution)

        # Fetch the interaction matrix between region_L and region_R
        mat = cool.matrix(balance=True).fetch(region_L, region_R)

        # Determine if regions cross over into each other (diagonal is included)
        diagonal_included = (chrom1 == chrom2) and (bin_end_L > bin_start_R) and (bin_end_R > bin_start_L)

        # If diagonal is included, take logarithm of raw Hi-C data before plotting
        if diagonal_included:
            # Fetch raw Hi-C data
            mat_to_plot = np.log(mat)
        else:
            # Use balanced data
            mat_to_plot = np.log1p(mat)

        plt.figure(figsize=(6,6))
        plt.matshow(mat_to_plot, fignum=False, cmap='gist_heat_r')
        plt.title(f'Interaction between {chrom1}:{start1}-{end1} and {chrom2}:{start2}-{end2}\nLog_p_value {log_p_value}')
        ax = plt.gca()

        # Calculate positions for the filters
        center_x = (start1 - bin_start_L)
        center_y = (start2 - bin_start_R)

        # Inside filter rectangle
        inside_size = filter_width
        rect_inside = plt.Rectangle(
            (center_y - inside_size//2 - .5, center_x - inside_size//2 - .5),
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
            (center_y - outside_size//2 - .5, center_x - outside_size//2 - .5),
            outside_size,
            outside_size,
            edgecolor='black',
            facecolor='none',
            linewidth=2,
            label='Outside Filter'
        )
        ax.add_patch(rect_outside)

        plt.legend(handles=[rect_inside, rect_outside], loc='upper right')

        if save_plots:
            plot_filename = f'metadomain_{chrom1}_{start1}_{chrom2}_{start2}.png'
            plot_file = os.path.join(plots_dir, plot_filename)
            plt.savefig(plot_file)
            print(f"Plot saved to {plot_file}")
        else:
            plt.show()
        plt.close()

        # If intermediate files are present, plot them
        if intermediate_dir:
            # Construct the prefix for the intermediate files
            prefix = f'{label}_chrL={chrom1}_chrR={chrom2}_'
            intermediate_dir_full = os.path.join(output_dir, plot_type, 'intermediate_files/')
            if not os.path.exists(intermediate_dir_full):
                print(f"Intermediate directory {intermediate_dir_full} does not exist.")
                continue

            # Paths to intermediate files
            peak_smooth_X1_file = os.path.join(intermediate_dir_full, prefix + 'peak_smooth_X1.npy')
            peak_smooth_Y1_file = os.path.join(intermediate_dir_full, prefix + 'peak_smooth_Y1.npy')
            full_logp_mat_file = os.path.join(intermediate_dir_full, prefix + 'full_logp_mat.npy')
            collapsed_logp_mat_file = os.path.join(intermediate_dir_full, prefix + 'collapsed_logp_mat.npy')
            bal_mat_file = os.path.join(intermediate_dir_full, prefix + ('bal_intra.npy' if plot_type == 'intra' else 'bal_inter.npy'))
            raw_mat_file = os.path.join(intermediate_dir_full, prefix + ('raw_intra.npy' if plot_type == 'intra' else 'raw_inter.npy'))
            oe_file = os.path.join(intermediate_dir_full, prefix + 'oe.npy')

            # Check if all intermediate files exist
            intermediate_files = [
                peak_smooth_X1_file,
                peak_smooth_Y1_file,
                full_logp_mat_file,
                collapsed_logp_mat_file,
                bal_mat_file,
                raw_mat_file,
                oe_file,
            ]

            if all(os.path.exists(f) for f in intermediate_files):
                try:
                    # Load intermediate matrices
                    full_logp_mat = np.load(full_logp_mat_file)
                    collapsed_logp_mat = np.load(collapsed_logp_mat_file)
                    bal_mat = np.load(bal_mat_file)
                    oe = np.load(oe_file)

                    # Extract the relevant submatrices
                    slice_L = slice(bin_start_L, bin_end_L)
                    slice_R = slice(bin_start_R, bin_end_R)

                    full_logp_mat_sub = full_logp_mat[slice_L, slice_R]
                    collapsed_logp_mat_sub = collapsed_logp_mat[slice_L, slice_R]
                    bal_mat_sub = bal_mat[slice_L, slice_R]
                    oe_sub = oe[slice_L, slice_R]

                    # Prepare the intermediate matrices for plotting
                    intermediate_matrices = {
                        'bal_mat': bal_mat_sub,
                        'oe': oe_sub,
                        'full_logp_mat': full_logp_mat_sub,
                        'collapsed_logp_mat': collapsed_logp_mat_sub,
                    }
                    kwargs = {
                        'bal_mat' : {
                            'cmap': 'gist_heat_r',
                            },
                        'oe' : {
                            'cmap': 'coolwarm',
                            },
                        'full_logp_mat' : {
                            'cmap': 'coolwarm',
                            },  
                        'collapsed_logp_mat' : {
                            'cmap': 'coolwarm',
                            },
                        }
                    titles = {
                        'bal_mat': 'Balanced Matrix',
                        'oe': 'Observed/Expected Matrix',
                        'full_logp_mat': 'Full -Log P-value Matrix',
                        'collapsed_logp_mat': 'Collapsed -Log P-value Matrix',
                    }

                    # Apply logarithmic transformation if diagonal is included
                    for name in ['bal_mat']:
                        mat = intermediate_matrices[name]
                        if diagonal_included:
                            mat[mat <= 0] = np.nan
                            intermediate_matrices[name] = np.log(mat)
                        else:
                            intermediate_matrices[name] = np.log1p(mat)

                    fig, axes = plt.subplots(2, 2, figsize=(10, 10))
                    axes = np.ravel(axes)
                    for ax, (name, matrix) in zip(axes, intermediate_matrices.items()):
                        im = ax.matshow(matrix, 
                                        **kwargs[name])
                        ax.set_title(titles[name])
                        fig.colorbar(im, ax=ax, shrink=.4)
                        inside_size = filter_width
                        rect_inside = plt.Rectangle(
                            (center_y - .5 - inside_size//2, center_x - .5 - inside_size//2),
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
                            (center_y - .5 - outside_size//2, center_x - .5 - outside_size//2),
                            outside_size,
                            outside_size,
                            edgecolor='black',
                            facecolor='none',
                            linewidth=2,
                            label='Outside Filter'
                        )
                        ax.add_patch(rect_outside)

                    plt.legend(handles=[rect_inside, rect_outside], loc='upper right')

                    if save_plots:
                        plot_filename = f'metadomain_{chrom1}_{start1}_{chrom2}_{start2}_intermediates.png'
                        plot_file = os.path.join(plots_dir, plot_filename)
                        plt.savefig(plot_file)
                        print(f"Intermediate plots saved to {plot_file}")
                    else:
                        plt.show()
                    plt.close()

                except Exception as e:
                    print(f"Error plotting intermediate files for {chrom1}-{chrom2}: {e}")
                    continue
            else:
                print(f"Intermediate files not found for {chrom1}-{chrom2}. Skipping intermediate plots.")

if __name__ == '__main__':
    main_plot_cli()
