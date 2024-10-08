import numpy as np
import scipy
import time
import logging
import scipy.ndimage
import scipy.stats
import scipy.signal

logger = logging.getLogger(__name__)

def prepare_inputs(cool, chrom1, chrom2, type='intra', minweight=0.0025, pc=1e-6):
    bal_intra = cool.matrix().fetch(chrom1, chrom2)
    raw_intra = cool.matrix(balance=False).fetch(chrom1, chrom2)

    w1 = cool.bins().fetch(chrom1)['weight']
    w2 = cool.bins().fetch(chrom2)['weight']
    nanind1 = np.zeros_like(w1)
    nanind2 = np.zeros_like(w2)

    nanind1[w1 <= minweight] = np.nan
    nanind2[w2 <= minweight] = np.nan
    nanind1[w1 > minweight] = 0
    nanind2[w2 > minweight] = 0

    extra_nan_filter = np.outer(nanind1, nanind2)
    bal_intra[np.isnan(extra_nan_filter)] = np.nan

    oe, exp = make_obs_exp_nolog(bal_intra, mat_type=type, pc=pc)
    return bal_intra, raw_intra, oe

def make_obs_exp_nolog(balanced_mat, mat_type='intra', pc=1e-6):
    exp = make_expected(balanced_mat, mat_type=mat_type)
    f = (balanced_mat + pc) / (exp + pc)
    return f, exp

def is_symmetric(mat, tol=1e-8):
    if mat.shape[0] != mat.shape[1]:
        return False
    return np.nansum(np.abs(mat - mat.T)) < tol

def make_expected(balanced_mat, mat_type='intra'):
    if (mat_type == 'inter') or (not is_symmetric(balanced_mat)):
        return np.nanmean(balanced_mat)

    exp = np.zeros(balanced_mat.shape)
    iu_x, iu_y = np.diag_indices(len(balanced_mat))
    last_m = 0
    for off_diag_k in range(len(balanced_mat)):
        m = np.nanmean(np.diag(balanced_mat, k=off_diag_k))
        if off_diag_k > 0 and np.isnan(m):
            exp[(iu_x, iu_y)] = last_m
            exp[(iu_y, iu_x)] = last_m
        else:
            exp[(iu_x, iu_y)] = m
            exp[(iu_y, iu_x)] = m
        iu_x, iu_y = iu_x[:-1], iu_y[1:]
        last_m = m
    return exp

def compute_prominent_peaks(oe, useSigma=False, sigma=0.75, prominence=4):
    z, peak_X1, peak_Y1 = compute_peaks(oe, useSigma=useSigma, sigma=sigma, prominence=prominence)
    ker = np.ones((3, 3))
    peak_smooth_X1 = scipy.ndimage.correlate(peak_X1, ker) > 0
    peak_smooth_Y1 = scipy.ndimage.correlate(peak_Y1, ker) > 0
    return peak_smooth_X1, peak_smooth_Y1, z

def get_filter_pvalue(raw_intra, bal_intra, type='intra', filter_n=35, filter_width=3, inter_pseudocount=0.5, pmin=1e-300, frac_min_valid=0.0):
    raw_intra = raw_intra.copy().astype(float)
    if type == 'inter':
        _, full_logp_mat, _, _, _, _, _ = call_peak_with_poisson_and_peak_cutoff(
            raw_intra + inter_pseudocount,
            filter_n=filter_n,
            filter_width=filter_width,
            pmin=pmin,
            frac_min_valid=frac_min_valid,
        )
    elif type == 'intra':
        _, full_logp_mat, _, _, _, _, _ = call_peak_with_poisson_and_peak_cutoff(
            raw_intra,
            filter_n=filter_n,
            filter_width=filter_width,
            pmin=pmin,
            frac_min_valid=frac_min_valid,
        )
    elif type == 'mega':
        _, full_logp_mat, _, _, _, _, _ = call_peak_with_poisson_and_peak_cutoff(
            raw_intra,
            filter_n=filter_n,
            filter_width=filter_width,
            pmin=pmin,
            frac_min_valid=frac_min_valid,
        )

    full_logp_mat[np.isnan(bal_intra)] = 0
    return full_logp_mat

def collapse_filter_and_peak_pvalues(full_logp_mat, peak_smooth_X1, peak_smooth_Y1, oe, pco=20):
    megaloops_of_interest = (full_logp_mat > pco) & peak_smooth_X1 & peak_smooth_Y1
    label_mat, _ = scipy.ndimage.label(megaloops_of_interest)
    obj_locations = scipy.ndimage.find_objects(label_mat)
    collapsed_logp_mat = np.zeros(megaloops_of_interest.shape) * 2

    for c, i in enumerate(obj_locations):
        l, r = i
        submat = oe[l, r].copy() * (label_mat[l, r] == (c + 1))
        S, E = np.unravel_index(np.nanargmax(np.ravel(submat)), submat.shape)
        s1, s2 = l.start, r.start
        collapsed_logp_mat[s1 + S, s2 + E] = full_logp_mat[s1 + S, s2 + E].copy()
    return collapsed_logp_mat

def call_peak_with_poisson_and_peak_cutoff(
    image,
    sigma=1,
    pmin=1e-300,
    logp_co=10,
    filter_n=15,
    filter_width=1,
    frac_min_valid=0.0,
    verbose=False,
):
    if verbose:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)

    C, O, L, R, U, D = set_filters_custom(filter_n, filter_width)
    nanfilter_image = image.copy()
    nanfilter_image[np.isnan(nanfilter_image)] = 0
    filt_image = scipy.ndimage.gaussian_filter(nanfilter_image, sigma=sigma)
    n_C = C.sum()
    n_O = O.sum()

    logger.info("Correlating")
    frac_valid_C = 1 - scipy.ndimage.correlate(np.isnan(image).astype(float), C) / n_C
    frac_valid_O = 1 - scipy.ndimage.correlate(np.isnan(image).astype(float), O) / n_O

    counts_C = (scipy.ndimage.correlate(nanfilter_image, C) / frac_valid_C).astype(int)
    counts_O = (scipy.ndimage.correlate(nanfilter_image, O) / frac_valid_O).astype(int)
    logger.info("Done correlating")

    logger.info("Poisson")
    pval_mat = scipy.stats.poisson(counts_O * n_C / n_O).sf(counts_C) + pmin
    logp_mat = -np.log10(pval_mat)

    frac_valid_bool = (frac_valid_C < frac_min_valid) | (frac_valid_O < frac_min_valid)
    logp_mat[frac_valid_bool] = 0

    is_sig = logp_mat > logp_co
    logger.info("Done poisson")
    row, col = np.where(is_sig)

    ind_mat = is_sig

    no_edge_ind_mat = np.zeros_like(ind_mat)
    sl = slice(10, -10)
    no_edge_ind_mat[sl, sl] = ind_mat[sl, sl]
    peak_row, peak_col = np.where(no_edge_ind_mat)

    resultsdict = {}
    resultsdict['counts_inner'] = counts_C
    resultsdict['counts_outer'] = counts_O
    return filt_image, logp_mat, row, col, peak_row, peak_col, resultsdict

def compute_peaks(z, useSigma=True, sigma=0.75, prominence=5):
    z = z.copy()
    z[np.isnan(z)] = np.nanmedian(z)
    if useSigma:
        z = scipy.ndimage.gaussian_filter(z, sigma=sigma)
    peak_X = np.zeros_like(z)
    peak_Y = np.zeros_like(z)
    for c, row in enumerate(z):
        peaks = scipy.signal.find_peaks(row, prominence=prominence)[0]
        peak_X[c, peaks] = 1

    for c, row in enumerate(z.T):
        peaks = scipy.signal.find_peaks(row, prominence=prominence)[0]
        peak_Y[peaks, c] = 1
    return z, peak_X, peak_Y

def set_filters_custom(n, w):
    inner_filter = np.zeros((n, n))
    outer_filter = np.ones((n, n))

    mid = inner_filter.shape[1] // 2
    inner_filter[mid - w : mid + w + 1, mid - w : mid + w + 1] = 1
    outer_filter -= inner_filter

    left_filter = np.zeros((n, n))
    left_filter[mid - w : mid + w + 1, : n // 2] = 1
    left_filter[inner_filter == 1] = 0

    right_filter = np.zeros((n, n))
    right_filter[mid - w : mid + w + 1, n // 2 :] = 1
    right_filter[inner_filter == 1] = 0

    up_filter = left_filter.T
    down_filter = right_filter.T

    return inner_filter, outer_filter, left_filter, right_filter, up_filter, down_filter

