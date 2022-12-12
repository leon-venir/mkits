from setuptools import setup

setup(
    name="mkits",
    description="multi-DFT codes assistant program.",
    long_description=
    """mkits is a python written tool containing many helpful initial- or post-processing commands for some popular first-principles calculation codes.
    """,
    version="0.83",
    author="Leon Ma",
    author_email="blustery.med@hotmail.com",
    license="GPLv3+",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    keywords="DFT code assistant",
    packages=["mkits"],
    entry_points={"console_scripts": ["mkits = mkits.main:mkits_main"]},
    
    install_requires=[
        "spglib", "numpy", "ase", "pandas"
    ],
    python_requires=">=3.5",
    url="https://github.com/leon-venir/mkits",
    download_url=(
        "https://github.com/leon-venir/mkits"
    )
)
