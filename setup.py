from setuptools import setup, find_packages
setup(
    name='InterDomain',
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
            'InterDomain=InterDomain.cli:main_cli',
            'plot_top_hits=InterDomain.plotting:main_plot_cli',
        ],
    },
)
