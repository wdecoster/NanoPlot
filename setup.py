# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))
exec(open('nanoplot/version.py').read())

setup(
    name='NanoPlot',
    version=__version__,
    description='Plotting suite for Oxford Nanopore sequencing data and alignments',
    long_description=open(path.join(here, "README.md"), encoding='utf-8').read(),
    long_description_content_type="text/markdown",
    url='https://github.com/wdecoster/NanoPlot',
    author='Wouter De Coster',
    author_email='decosterwouter@gmail.com',
    license='GPLv3',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    keywords='nanopore sequencing plotting quality control',
    packages=find_packages() + ['extra', 'scripts'],
    python_requires='>=3',
    install_requires=['biopython',
                      'pysam>0.10.0.0',
                      'pandas>=1.1.0',
                      'numpy>=1.16.5',
                      'scipy',
                      'python-dateutil',
                      'nanoget>=1.14.0',
                      'nanomath>=1.0.0',
                      'plotly>=4.1.0',
                      'pyarrow',
                      'kaleido'
                      ],
    package_data={'NanoPlot': []},
    package_dir={'nanoplot': 'nanoplot'},
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'NanoPlot=nanoplot.NanoPlot:main',
        ],
    },
    data_files=[("", ["LICENSE"])]
)
