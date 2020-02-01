from setuptools import setup
from kobe.__version__ import version

setup(
    name='kobe',
    version=version,
    python_requires='>=2.7',
    packages=[
        'kobe',
        'kobe.pipeline',
        'kobe.view',
        'kobe.circulate',
        'kobe.interface',
    ],
    install_requires=[
        'meander',
        'astropy',
        'numpy',
        'healpy',
        'matplotlib',
        'astroquery',
    ],
    classifiers=[                     
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Astronomy',       
    ], 
    zip_safe = False
)
