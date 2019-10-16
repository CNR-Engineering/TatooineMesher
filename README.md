TatooineMesher
==============

<img style="float: right" src="https://github.com/CNR-Engineering/TatooineMesher/raw/master/media/logo_tatooinemesher_256px.png" width="168px" />

> Channel mesher and interpolator from 1D cross-sections and constraint lines

## Description

Set of 4 command line scripts:
* `densify_profiles.py` : interpolate intermediate cross-sections
* `mesher_and_interpolator.py`: mesh and interpolate from a set of cross-sections
* `mesh_crue10_run.py` : visualize Crue10 model geometry and results
* `mesh_mascaret_run.py` : visualize Mascaret model geometry and results

The commande line scripts of TatooineMesher are located in the `cli` folder.
 
See [wiki pages](https://github.com/CNR-Engineering/TatooineMesher/wiki) to learn how to use them.

## Pre-requisites

**Python 3** with **packages** listed in [requirements.txt](requirements.txt).

[PyTelTools](https://github.com/CNR-Engineering/PyTelTools) and [Crue10_tools](https://github.com/CNR-Engineering/Crue10_tools) are dependencies that can be installed in another directory
and linked with the `PYTHONPATH` environment variable.
These packages are used to read/write some specific file formats of Telemac, Mascaret and Crue10 tools.
