import setuptools
import pathlib

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

HERE = pathlib.Path(__file__).parent
install_requires = (HERE / "devtools/requirements.txt").read_text().splitlines()

setuptools.setup(
    name="pycbh-colliner",
    version="0.0.0",
    author="Eric M. Collins",
    author_email="colliner@indiana.edu",
    description="pyCBH",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/colliner/pyCBH",
    packages=['pycbh'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: MIT License",
    ],
    python_requires='>=3.6',
    install_requires=install_requires,
    dependency_links=['http://github.com/colliner/xyz2mol#egg=xyz2mol'])
