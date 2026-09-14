from setuptools import setup

setup(
    name="mkits",
    description="multi-DFT codes assistant program.",
    long_description="""mkits is a Python library that provides pre-processing and 
post-processing utilities for first-principles calculation software (such as VASP, 
Quantum ESPRESSO, and ORCA). It is designed to be called directly (for instance, 
within a Jupyter environment) rather than executed as a standalone command-line 
interface (CLI) tool. Its features encompass structure processing (supporting 
multiple read/write formats, spglib-based symmetry analysis, and k-point mesh/path 
generation), input file generation, automated multi-step workflows (including 
staged relaxation, opt→scf→band→dos calculation chains, convergence scans, Born 
effective charge calculations, and deformation potential, elastic modulus, and 
effective mass scans), and post-processing analysis (extraction of DOS, PROCAR, 
and band data, total energy analysis, and effective mass fitting). Furthermore, 
the library contains a carrier mobility estimation tool based on deformation 
potential theory, which has been cross-validated with mobility results calculated 
using the Boltzmann transport equation via epw.x.""",
    version="1.0.0",
    author="Leon Ma",
    author_email="blustery.med@hotmail.com",
    license="GPLv3+",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
    keywords="DFT code assistant",
    packages=["mkits"],
    entry_points={"console_scripts": ["mkits = mkits.main:mkits_main"]},
    
    install_requires=[
        "spglib", "numpy", "shapely", "scipy", "matplotlib", "pandas"
    ],
    python_requires=">=3.7",
    url="https://github.com/leon-venir/mkits",
    download_url=(
        "https://github.com/leon-venir/mkits"
    )
)
