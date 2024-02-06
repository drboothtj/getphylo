from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    description = fh.read()

setup(
    name="getphylo",
    version="0.2.0",
    author="Thomas J. Booth",
    author_email="thoboo@biosustain.dtu.dk",
    packages=find_packages(),
    description="a python package for automated generation of heuristic phylogenetic trees from genbank files",
    long_description=description,
    long_description_content_type="text/markdown",
    url="https://github.com/DrBoothTJ/getphylo",
    license='GNU General Public License v3.0',
    python_requires='>=3.7',
    install_requires=['Bio'],
    entry_points={'console_scripts': ["getphylo=getphylo.__main__:entrypoint"]}
)
