# future-wind
Processing of CMIP6 data for wind interpolation using pyesgf and ESGF OPeNDAP servers.

We provide a set of common subroutines and some examples for a few types of CMIP6 model structures. 

| model | calendar | horizontal grid | vertical grid | script name | extra comments |
|-------|----------|-----------------|---------------|-------------|----------------|
| ACCESS-CM2 | proleptic_gregorian | arakawa-C | height-level| ACCESS.py|  |
| CanESM5	| noLeap	| gaussian | sigma-pressure | CanESM5.py | |
| HadGEM3-GC31-LL | 360_day	| arakawa-C | height-level | HadGEM.py| very slow server |

More example scripts will be provided in the future.

The code has been tested using:


