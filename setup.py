# setup.py
from setuptools import setup, find_packages

setup(
    name='InterDomain',
    version='0.2.0',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'pandas',
        'cooler',
        'scipy',
        'matplotlib',
    ],
    entry_points={
        'console_scripts': [
            'call_metadomains_intra=InterDomain.cli_intra:main_cli_intra',
            'call_metadomains_inter=InterDomain.cli_inter:main_cli_inter',
            'interdomain_plot=InterDomain.plotting:main_plot_cli',
        ],
    },
)


# call_metadomains_intra ~/pritlab/jupys/tregs/Treg_all.mcool::/resolutions/50000 --n_workers 16 --save_intermediates --filter_width 1 --filter_n 15 --cutoff 0
# call_metadomains_inter ~/pritlab/jupys/tregs/Treg_all.mcool::/resolutions/50000 --n_workers 16 --save_intermediates --filter_width 3 --filter_n 35