from setuptools import setup
from kobe.__version__ import version, description

setup(
    name='astrokobe',
    version=version,
    description=description,
    python_requires='>=2.7',
    packages=['kobe'],
    install_requires=[
        'astropy',
        'numpy',
        'healpy',
        'matplotlib',
    ],
    classifiers=[                     
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Astronomy',       
    ], 
    zip_safe = False
)
