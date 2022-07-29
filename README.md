TatooineMesher
==============

![Python package](https://github.com/CNR-Engineering/TatooineMesher/workflows/Python%20package/badge.svg)

Tested versions: 3.6, 3.7, 3.8, 3.9 and 3.10.

<img style="float: right" src="https://github.com/CNR-Engineering/TatooineMesher/raw/master/media/logo_tatooinemesher_256px.png" width="168px" />

> Channel mesher and interpolator from 1D cross-sections and constraint lines

## Description

Set of 4 command line scripts:
* `densify_cross_sections.py` : interpolate initial and intermediate cross-sections
* `mesh_and_interpolate.py`: mesh and/or interpolate from a set of cross-sections
* `mesh_mascaret_run.py` : visualize Mascaret model geometry and results
* `mesh_crue10_run.py` : visualize Crue10 model geometry and results

The command line scripts of TatooineMesher are located in the `cli` folder.
 
See [wiki pages](https://github.com/CNR-Engineering/TatooineMesher/wiki) to learn how to use them.

## Pre-requisites and installation

**Python 3** with **packages** listed in [requirements.txt](requirements.txt). To install them use:
```bash
pip install -r requirements.txt
```

[PyTelTools](https://github.com/CNR-Engineering/PyTelTools) and [Crue10_tools](https://github.com/CNR-Engineering/Crue10_tools) are dependencies that can be installed in another directory
and linked with the `PYTHONPATH` environment variable.
These packages are used to read/write some specific file formats of Telemac, Mascaret and Crue10 tools.

## Reference
_TatooineMesher: Anisotropic interpolation from 1D cross-sections and 2D channel mesher_, L. Duron, F.-X. Cierco and K. Saad, **XXVIth TELEMAC-MASCARET User Conference 2019**, Toulouse (France). See [paper](https://github.com/CNR-Engineering/TatooineMesher/raw/master/media/publi/article/TUC2019_paper_Duron-et-al_TatooineMesher.pdf) and [presentation](https://github.com/CNR-Engineering/TatooineMesher/raw/master/media/publi/presentation/TUC2019_slides_Duron-et-al_TatooineMesher.pdf).
