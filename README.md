# future-wind

[![DOI](https://zenodo.org/badge/497835472.svg)](https://zenodo.org/badge/latestdoi/497835472)

Processing six hourly CMIP6 model output for a limited region and interpolating winds to constant heights above ground.
The search of the necessary fields (u, v, T, q, and ps) is done using [pyesgf](https://esgf-pyclient.readthedocs.io/en/latest/index.html) and ESGF OPeNDAP servers.

We provide a set of common subroutines and some examples for a few types of CMIP6 model structures. More example scripts will be provided in the future.

| model | calendar | horizontal grid | vertical grid | script name | extra comments |
|-------|----------|-----------------|---------------|-------------|----------------|
| ACCESS-CM2 | proleptic_gregorian | arakawa-C | height-level| ACCESS.py|  |
| CanESM5| noLeap| gaussian | sigma-pressure | CanESM5.py | |
| HadGEM3-GC31-LL | 360_day| arakawa-C | height-level | HadGEM.py| very slow server |

This code accompanies the article Current and future wind energy resources in the North Sea according to CMIP6 by Andrea N. Hahmann, Oscar García-Santiago and Alfredo Peña, published in Wind Energy Science in 2022.
