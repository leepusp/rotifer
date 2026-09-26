from setuptools import setup, find_packages

setup(
    name="rotifer",
    version="0.1.0",  # Update as appropriate
    description=(
        "High-level libraries and command line tools for comparative genomics "
        "and computational analysis of biological sequences"
    ),
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="Robson F. de Souza and contributors",
    url="https://github.com/leepusp/rotifer",
    project_urls={
        "Documentation": "https://leepusp.github.io/rotifer/",
        "Source": "https://github.com/leepusp/rotifer",
        "Issues": "https://github.com/leepusp/rotifer/issues",
    },
    packages=find_packages("lib"),  # Finds packages in the "lib" folder
    package_dir={"": "lib"},  # Root directory for the packages
    python_requires=">=3.6",  # Specify minimum Python version
    install_requires=[
        "numpy",   # Add your dependencies here
        "pandas",
    ],
    license="BSD-3-Clause",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)

