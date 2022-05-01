import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="custardpy",
    version="0.1.0",
    license="GPL3.0",
    install_requires=[
        "numpy<=1.17",
        "pandas>=0.22.0",
        "scipy>=1.3",
        "matplotlib>=3.2.2",
        "seaborn>=0.11.1"
        "h1d>=0.2.0"
    ],
    author="Ryuichiro Nakato",
    author_email="rnakato@iqb.u-tokyo.ac.jp",
    description="Hi-C analysis tools by Python3",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/rnakato/custardpy",
    keywords="Hi-C analysis, 3D genome, NGS",
  ####  scripts=['scripts/eeisp', 'eeisp/eeisp_add_genename_from_geneid'],
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)
