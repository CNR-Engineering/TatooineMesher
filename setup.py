#!/usr/bin/env python

from glob import glob
from setuptools import find_packages, setup

from tatooinemesher import VERSION


EXCLUDE_FROM_PACKAGES = ["cli", "media"]


with open("requirements.txt") as f:
    requirements = f.read().splitlines()

cli_files = []
for file in glob("cli/*.py"):
    if not file.endswith("__init__.py"):
        cli_files.append(file)

setup(
    name="TatooineMesher",
    version=VERSION,
    author="Luc Duron",
    author_email="l.duron@cnr.tm.fr",
    packages=find_packages(exclude=EXCLUDE_FROM_PACKAGES),
    scripts=cli_files,
    install_requires=requirements,
    description="Channel mesher and interpolator from 1D cross-sections and constraint lines",
    url="https://github.com/CNR-Engineering/TatooineMesher",
)
