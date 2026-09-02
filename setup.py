from setuptools import setup, find_packages

setup(
    name="rotifer",
    version="0.1.0",  # Update as appropriate
    description="A Python library for rotifer bioinformatics",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="Your Name",  # Replace with your name
    author_email="your.email@example.com",  # Replace with your email
    url="https://github.com/leepbioinfo/rotifer",
    packages=find_packages("lib"),  # Finds packages in the "lib" folder
    package_dir={"": "lib"},  # Root directory for the packages
    python_requires=">=3.6",  # Specify minimum Python version
    install_requires=[
        "numpy",   # Add your dependencies here
        "pandas",
    ],
    extras_require={
        # pyarrow speeds up the scans in rotifer.db.uniprot.io, which
        # fall back to the standard library when it is absent
        "speedups": ["pyarrow"],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",  # Update license if different
        "Operating System :: OS Independent",
    ],
)

