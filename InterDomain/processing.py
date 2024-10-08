import logging
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
from .utils import prepare_inputs, compute_prominent_peaks, get_filter_pvalue, collapse_filter_and_peak_pvalues

def main(
    cool,
    n_workers=1,
    cutoff=2000000,
    label='test',
    save_intermediates=True,
    intermediate_dir='metadomain_prominent_peak_output/',
    output_dir='bedfile_output/',
    minweight=0.0025,
    pc=1e-6,
    useSigma=False,
    sigma=0.75,
    prominence=4,
    filter_n=35,
    filter_width=3,
    inter_pseudocount=0.5,
    pco=20,
    pmin=1e-300,
    frac_min_valid=0.0,
    logger=None,
):
    if logger is None:
        logger = logging.getLogger(__name__)

    chroms = cool.chromsizes
    resolution = cool.binsize
    chroms = chroms[chroms > cutoff]
    logger.info(f"Resolution: {resolution}")
    logger.info(f"Chromosomes to use: {list(chroms.index)}")

    all_inputs = make_inputs(
        cool,
        label,
        chroms,
        cutoff=cutoff,
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
    logger.info(f"Starting processing with {n_workers} workers")
    with ProcessPoolExecutor(max_workers=n_workers) as e:
        results = list(e.map(run, all_inputs))
    results = pd.concat(results, axis=0)
    num_detected = len(results)
    logger.info(f"Number of metadomains detected: {num_detected}")
    output_file = output_dir + f'intra_metadomains.{label}.res_{resolution}.tsv'
    results.to_csv(output_file, sep='\t', index=False)
    logger.info(f"Results saved to {output_file}")

def make_inputs(cool, label, CHROMS_TO_USE, **kwargs):
    inps = []
    for chrom in CHROMS_TO_USE.index:
        d = {
            'cool': cool,
            'chrom': chrom,
            'label': label,
        }
        d.update(kwargs)
        inps.append(d)
    return inps

def run(kwargs):
    cool = kwargs.get('cool')
    chrom = kwargs.get('chrom')
    label = kwargs.get('label')
    save_all = kwargs.get('save_all', True)
    save_dir = kwargs.get('save_dir', 'metadomain_prominent_peak_output/')
    minweight = kwargs.get('minweight', 0.0025)
    pc = kwargs.get('pc', 1e-6)
    useSigma = kwargs.get('useSigma', False)
    sigma = kwargs.get('sigma', 0.75)
    prominence = kwargs.get('prominence', 4)
    filter_n = kwargs.get('filter_n', 35)
    filter_width = kwargs.get('filter_width', 3)
    inter_pseudocount = kwargs.get('inter_pseudocount', 0.5)
    pco = kwargs.get('pco', 20)
    pmin = kwargs.get('pmin', 1e-300)
    frac_min_valid = kwargs.get('frac_min_valid', 0.0)

    bal_intra, raw_intra, oe = prepare_inputs(
        cool, chrom, chrom, minweight=minweight, pc=pc
    )
    peak_smooth_X1, peak_smooth_Y1, z = compute_prominent_peaks(
        oe, useSigma=useSigma, sigma=sigma, prominence=prominence
    )
    full_logp_mat = get_filter_pvalue(
        raw_intra,
        bal_intra,
        filter_n=filter_n,
        filter_width=filter_width,
        inter_pseudocount=inter_pseudocount,
        pmin=pmin,
        frac_min_valid=frac_min_valid,
    )
    collapsed_logp_mat = collapse_filter_and_peak_pvalues(
        full_logp_mat, peak_smooth_X1, peak_smooth_Y1, oe, pco=pco
    )

    if save_all:
        np.save(save_dir + f'{label}_chrL={chrom}_chrR={chrom}_peak_smooth_X1.npy', peak_smooth_X1)
        np.save(save_dir + f'{label}_chrL={chrom}_chrR={chrom}_peak_smooth_Y1.npy', peak_smooth_Y1)
        np.save(save_dir + f'{label}_chrL={chrom}_chrR={chrom}_full_logp_mat.npy', full_logp_mat)
        np.save(save_dir + f'{label}_chrL={chrom}_chrR={chrom}_collapsed_logp_mat.npy', collapsed_logp_mat)
        np.save(save_dir + f'{label}_chrL={chrom}_chrR={chrom}_bal_intra.npy', bal_intra)
        np.save(save_dir + f'{label}_chrL={chrom}_chrR={chrom}_raw_intra.npy', raw_intra)
        np.save(save_dir + f'{label}_chrL={chrom}_chrR={chrom}_oe.npy', oe)

    X, Y = np.where(collapsed_logp_mat > 0)
    log_p_values = collapsed_logp_mat[X, Y]

    results = pd.DataFrame()
    results['chrom'] = [chrom] * len(X)
    results['start'] = X
    results['end'] = Y
    results['log_p_value'] = log_p_values
    return results
