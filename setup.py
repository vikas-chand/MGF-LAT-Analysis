from setuptools import setup, find_packages

setup(
    name='fermi_lat_transients',
    version='0.1.0',
    description='Fermi-LAT analysis pipeline for transient sources (MGFs, FXTs, GRBs)',
    author='Vikas Chand',
    url='https://github.com/vikas-chand/MGF-LAT-Analysis',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'pandas',
        'astropy',
        'matplotlib',
        'pyyaml',
        'scipy',
    ],
    entry_points={
        'console_scripts': [
            'lat-transients-run = fermi_lat_transients.client:main',
        ],
    },
    package_data={
        'fermi_lat_transients': ['inputs.yaml'],
    },
    python_requires='>=3.8',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Programming Language :: Python :: 3',
    ],
)
