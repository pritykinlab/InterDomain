from setuptools import setup, find_packages
setup(
    name='metadomain_peak_caller',
    version='0.1.0',
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
            'metadomain_peak_caller=metadomain_peak_caller.cli:main_cli',
            'plot_top_hits=metadomain_peak_caller.plotting:main_plot_cli',
        ],
    },
)
