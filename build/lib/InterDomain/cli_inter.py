# metadomain_peak_caller/cli.py
import argparse
import logging
import os
from .processing_inter import main_inter
import cooler

def parse_args():
    parser = argparse.ArgumentParser(description='InterDomain inter-chromosomal Metadomain Caller')
    parser.add_argument('cool_file', help='Path to the .cool file')
    parser.add_argument('--n_workers', type=int, default=1, help='Number of workers for parallel processing')
    parser.add_argument('--label', default='test', help='Label for output files')
    parser.add_argument('--save_intermediates', action='store_true', help='Save intermediate files used for calculation')
    parser.add_argument('--output_dir', default='bedfile_output/', help='Directory to save output files')

    parser.add_argument('--minweight', type=float, default=0.0025, help='Minimum weight for bin filtering in prepare_inputs')
    parser.add_argument('--pc', type=float, default=1e-6, help='Pseudocount for Obs/Exp calculation (usually not an important parameter needed)')
    parser.add_argument('--pmin', type=float, default=1e-300, help='Minimum p-value in call_peak_with_poisson_and_peak_cutoff')
    parser.add_argument('--frac_min_valid', type=float, default=0.0, help='Minimum fraction of valid data required in call_peak_with_poisson_and_peak_cutoff')
    parser.add_argument('--minChromSize', type=float, default= 2_000_000, help='Minimum chromosome size')
    parser.add_argument('--chroms_to_ignore', type=str, default='', help='Comma-separated list of chromosomes to ignore (e.g., 10,12,M,Y or chr10,chr12, etc.)')

    parser.add_argument('--useSigma', action='store_true', help='Use Gaussian smoothing of Hi-C matrix to find peaks (in case of sparsity)')
    parser.add_argument('--sigma', type=float, default=0.75, help='Sigma value Gaussian smoothing')
    parser.add_argument('--prominence', type=float, default=4, help='Prominence value to identify peaks (how enriched above baseline should a metadomain)')
    parser.add_argument('--filter_n', type=int, default=35, help='The width for the outside filter used as the local control')
    parser.add_argument('--filter_width', type=int, default=3, help='The width for the inside filter used for calling metadomains')
    parser.add_argument('--inter_pseudocount', type=float, default=0.5, help='Pseudocount for interchromosomal metadomains')
    parser.add_argument('--pco', type=float, default=20, help='P-value cutoff for calling metadomains')
    args = parser.parse_args()
    return args

def main_cli_inter():
    args = parse_args()
    # Ensure output directories exist
    args.output_dir = os.path.join(args.output_dir, 'inter/')
    os.makedirs(args.output_dir, exist_ok=True)
    intermediate_dir = os.path.join(args.output_dir, 'intermediate_files/')
    os.makedirs(intermediate_dir, exist_ok=True)

    # Set up logging
    log_file = os.path.join(args.output_dir, 'InterDomain.log')
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    # Create handlers
    c_handler = logging.StreamHandler()
    f_handler = logging.FileHandler(log_file)
    c_handler.setLevel(logging.INFO)
    f_handler.setLevel(logging.INFO)

    # Create formatters and add them to handlers
    c_format = logging.Formatter('%(message)s')
    f_format = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    c_handler.setFormatter(c_format)
    f_handler.setFormatter(f_format)

    # Add handlers to the logger
    logger.addHandler(c_handler)
    logger.addHandler(f_handler)

    # Log basic information

    logger.info(f"Cool file: {args.cool_file}")
    logger.info(f"Number of workers: {args.n_workers}")
    logger.info(f"Label: {args.label}")
    logger.info(f"Save intermediates: {args.save_intermediates}")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info(f"filter_n: {args.filter_n}")
    logger.info(f"filter_width: {args.filter_width}")
    # Log other parameters if necessary

    # Load the cool file
    cool = cooler.Cooler(args.cool_file)

    # Call the main function with arguments
    main_inter(
        cool=cool,
        n_workers=args.n_workers,
        label=args.label,
        save_intermediates=args.save_intermediates,
        intermediate_dir=intermediate_dir,
        output_dir=args.output_dir,
        minweight=args.minweight,
        pc=args.pc,
        useSigma=args.useSigma,
        sigma=args.sigma,
        prominence=args.prominence,
        filter_n=args.filter_n,
        filter_width=args.filter_width,
        inter_pseudocount=args.inter_pseudocount,
        pco=args.pco,
        pmin=args.pmin,
        frac_min_valid=args.frac_min_valid,
        minChromSize=args.minChromSize,
        chroms_to_ignore=args.chroms_to_ignore,
        logger=logger,  # Pass logger to main
    )

if __name__ == '__main__':
    main_cli_inter()
