# InterDomain/processing_inter.py

import logging
import os
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
from .utils import prepare_inputs, compute_prominent_peaks, get_filter_pvalue, collapse_filter_and_peak_pvalues

# import logging
# logging.basicConfig(level=logging.INFO)
# logger = logging.getLogger(__name__)

def main_inter(
    cool,
    n_workers=1,
    label='test',
    save_intermediates=True,
    intermediate_dir='intermediate_files/',
    output_dir='bedfile_output/',
    minweight=0.0025,
    pc=1e-6,
    useSigma=True,
    sigma=2.0,
    prominence=4,
    filter_n=15,
    filter_width=3,
    inter_pseudocount=0.5,
    pco=5,
    pmin=1e-300,
    frac_min_valid=0.0,
    logger=None,
    minChromSize=2_000_000,
    chroms_to_ignore=[],
):
    if logger is None:
        logger = logging.getLogger(__name__)

    chroms = cool.chromsizes
    chroms = chroms[chroms > minChromSize]
    chroms = chroms.index
    chroms = [chrom for chrom in chroms if chrom not in chroms_to_ignore]
    resolution = cool.binsize

    logger.info(f"Running Interchromosomal")
    logger.info(f"Resolution: {resolution}")
    logger.info(f"Chromosomes to use: {chroms}")

    all_inputs = make_inputs_inter(
        cool,
        label,
        chroms,
        save_all=save_intermediates,
        save_dir=intermediate_dir,
        resolution=resolution,
        minweight=minweight,
        pc=pc,
        useSigma=useSigma,
        sigma=sigma,
        prominence=prominence,
        filter_n=filter_n,
        filter_width=filter_width,
        inter_pseudocount=inter_pseudocount,
        pco=pco,
        pmin=pmin,
        frac_min_valid=frac_min_valid,
    )
    logger.info(f"Starting inter-chromosomal processing with {n_workers} workers")
    with ProcessPoolExecutor(max_workers=n_workers) as e:
        results = list(e.map(run_inter, all_inputs))
    results = pd.concat(results, axis=0)
    num_detected = len(results)
    logger.info(f"Number of inter-chromosomal metadomains detected: {num_detected}")
    output_file = os.path.join(output_dir, f'inter_metadomains.{label}.res_{resolution}.tsv')
    results.to_csv(output_file, sep='\t', index=False)
    logger.info(f"Results saved to {output_file}")

def make_inputs_inter(cool, label, chroms, **kwargs):
    inps = []
    for i, chrom1 in enumerate(chroms):
        for chrom2 in chroms[i+1:]:  # Ensure each pair is unique and chrom1 != chrom2
            d = {
                'cool': cool,
                'chrom1': chrom1,
                'chrom2': chrom2,
                'label': label,
            }
            d.update(kwargs)
            inps.append(d)
    return inps

def run_inter(kwargs):
    cool = kwargs.get('cool')
    chrom1 = kwargs.get('chrom1')
    chrom2 = kwargs.get('chrom2')
    label = kwargs.get('label')
    save_all = kwargs.get('save_all', True)
    save_dir = kwargs.get('save_dir', 'intermediate_files/')
    minweight = kwargs.get('minweight', 0.0025)
    pc = kwargs.get('pc', 1e-6)
    useSigma = kwargs.get('useSigma', True)
    sigma = kwargs.get('sigma', 2.0)
    prominence = kwargs.get('prominence', 4)
    filter_n = kwargs.get('filter_n', 15)
    filter_width = kwargs.get('filter_width', 3)
    inter_pseudocount = kwargs.get('inter_pseudocount', 0.5)
    pco = kwargs.get('pco', 5)
    pmin = kwargs.get('pmin', 1e-300)
    frac_min_valid = kwargs.get('frac_min_valid', 0.0)
    print("Running inter-chromosomal processing for", chrom1, "-", chrom2, flush=True)
    print("Using parameters: useSigma={}, sigma={}, prominence={}".format(useSigma, sigma, prominence), flush=True)
    print("Using parameters: filter_n={}, filter_width={}, inter_pseudocount={}, pco={}, pmin={}, frac_min_valid={}".format(filter_n, filter_width, inter_pseudocount, pco, pmin, frac_min_valid), flush=True)
    try:
        # Prepare inputs
        bal_inter, raw_inter, oe = prepare_inputs(
            cool, chrom1, chrom2, type='inter', minweight=minweight, pc=pc
        )
        oe[oe > 200] = np.nan
        bal_inter[np.isnan(oe)] = np.nan
        raw_inter[np.isnan(bal_inter)] = 0

        # Compute prominent peaks
        peak_smooth_X1, peak_smooth_Y1, z = compute_prominent_peaks(
            oe, useSigma=useSigma, sigma=sigma, prominence=prominence
        )

        # Get filter p-values
        full_logp_mat = get_filter_pvalue(
            raw_inter,
            bal_inter,
            type='inter',
            filter_n=filter_n,
            filter_width=filter_width,
            inter_pseudocount=inter_pseudocount,
            pmin=pmin,
            frac_min_valid=frac_min_valid,
        )

        # Collapse filter and peak p-values
        collapsed_logp_mat = collapse_filter_and_peak_pvalues(
            full_logp_mat, peak_smooth_X1, peak_smooth_Y1, oe, pco=pco
        )

        # Apply OE cutoff to collapsed logp matrix
        collapsed_logp_mat[oe < 6] = 0

        # Save intermediate files
        if save_all:
            prefix = os.path.join(save_dir, f'{label}_chrL={chrom1}_chrR={chrom2}_')
            np.save(prefix + 'peak_smooth_X1.npy', peak_smooth_X1)
            np.save(prefix + 'peak_smooth_Y1.npy', peak_smooth_Y1)
            np.save(prefix + 'full_logp_mat.npy', full_logp_mat)
            np.save(prefix + 'collapsed_logp_mat.npy', collapsed_logp_mat)
            np.save(prefix + 'bal_inter.npy', bal_inter)
            np.save(prefix + 'raw_inter.npy', raw_inter)
            np.save(prefix + 'oe.npy', oe)

        # Extract significant interactions
        X, Y = np.where(collapsed_logp_mat > 0)
        log_p_values = collapsed_logp_mat[X, Y]

        # Convert bin indices to genomic coordinates
        bins_chrom1 = cool.bins().fetch(chrom1).reset_index(drop=True)
        bins_chrom2 = cool.bins().fetch(chrom2).reset_index(drop=True)

        results = pd.DataFrame()
        results['start1'] = bins_chrom1.iloc[X]['start'].values
        results['end1'] = bins_chrom1.iloc[X]['end'].values
        results['chrom1'] = chrom1
        results['chrom2'] = chrom2
        results['start2'] = bins_chrom2.iloc[Y]['start'].values
        results['end2'] = bins_chrom2.iloc[Y]['end'].values
        results['log_p_value'] = log_p_values
        results = results[['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'log_p_value']]
        return results

    except Exception as e:
        logging.error(f"Error processing {chrom1}-{chrom2}: {e}")
        return pd.DataFrame()  # Return empty DataFrame on error

