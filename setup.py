import io
from os.path import dirname, join
from setuptools import setup

# read the contents of your README file
from os import path

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()


def get_version(relpath):
    """Read version info from a file without importing it"""
    for line in io.open(join(dirname(__file__), relpath), encoding="cp437"):
        if "__version__" in line:
            if '"' in line:
                # __version__ = "0.9"
                return line.split('"')[1]
            elif "'" in line:
                return line.split("'")[1]


setup(
    name='cfphage_pipe',
    version=get_version("cfphage_pipe/__init__.py"),
    url='https://github.com/julianzaugg/CFphage_pipe',
    license='GPL-3.0',
    author='Julian Zaugg',
    author_email='j.zaugg@uq.edu.au',
    description='cfphage_pipe - pipeline to process Cystic Fibrosis isolates for phage discovery and evaluation',
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=['cfphage_pipe'],
    package_data={'': [
        "cfphage_pipe/*",
    ]},
    data_files=[(".", ["README.md", "LICENSE"])],
    include_package_data=True,
    install_requires=[
        "snakemake",
        "numpy",
        "ruamel.yaml>=0.15.99",
        "pandas",
        "biopython",
    ],
    # install via conda: click, pandas, pyyaml, snakemake
    entry_points={
        'console_scripts': [
            'cfphage_pipe = cfphage_pipe.cfphage_pipe:main'
        ]
    },
    classifiers=["Topic :: Scientific/Engineering :: Bio-Informatics"],
)
