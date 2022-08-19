from setuptools import setup

setup(
    name="mkits",
    description="multi-DFT codes assistant program.",
    long_description=
    """mkits 
    """,
    version="0.1",
    author="Leon Ma",
    author_email="blustery.med@hotmail.com",
    license="GPLv3+",
    classifiers=[
        "Development Status :: 1 - Beta",
        "Programming Language :: Python :: 3.6",
        "License :: OSI Approved :: GNU General Public License v3 or"
        " later (GPLv3+)",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    keywords="DFT code assistant",
    packages=["mkits"],
    entry_points={"console_scripts": ["mkits = mkits.main:mkits_main"]},
    
    install_requires=[
        "spglib", "numpy", "ase"
    ],
    python_requires=">=3.5",
    url="https://xxxx",
    download_url=(
        "https://xxxxx"
    )
)