import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="custardpy",
    version="0.1.10",
    license="GPL3.0",
    install_requires=[
        "numpy>=1.18",
        "pandas>=1.3.0",
        "scipy>=1.3",
        "scikit-learn>=1.0.0",
        "matplotlib>=3.2.2",
        "seaborn>=0.11.1",
        "h1d>=0.2.0",
        "hic-straw>=1.3.0",
    ],
    author="Ryuichiro Nakato",
    author_email="rnakato@iqb.u-tokyo.ac.jp",
    description="Hi-C analysis tools by Python3",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/rnakato/custardpy",
    keywords="Hi-C analysis, 3D genome, NGS",
    scripts=['custardpy/checkHiCfile.py',
             'custardpy/drawSquareMulti.py',
             'custardpy/drawSquarePair.py',
             'custardpy/drawSquareRatioMulti.py',
             'custardpy/drawSquareRatioPair.py',
             'custardpy/drawTriangleMulti.py',
             'custardpy/drawTrianglePair.py',
             'custardpy/drawTriangleRatioMulti.py',
             'custardpy/HMMDRF.py',
             'custardpy/plotCompartmentGenome.py',
             'custardpy/plotHiCfeature.py'],
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)
