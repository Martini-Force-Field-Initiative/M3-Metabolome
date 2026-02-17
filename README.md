# The Martini 3 Metabolome

Parameters, structure files, and parameterisation reference files for the Martini 3 Metabolome.

The currently available models are listed in the [database](database.md).

If you are interested in using the parameters, the structure and parameter files are contained
in the [gros](gros) and [itps](itps) subdirectories respectively.

If you are interested in more of the details of the parameterisation process, mapping files and 
figures comparing simulated and reference distributions are contained in the [mappings](mappings)
and [distributions](distributions) subdirectories.

The [misc](misc) folder contains any other files that may be relevant or of use.  For the moment, 
it contains the [database.md](database.md) file in csv format, for ease of use.

## Citation
Citation forthcoming

## Repo organisation
```commandline
.
├── distributions
│   ├── nucleotide_plots
│   │   └── *.png
│   └── *.png
├── gros
│   └── *.gro
├── itps
│   ├── bonded
│   │   └── martini_v3.0.0_nucleotide_ffbonded_v1.itp
│   └── molecules
│       └── *.itp
├── mappings
│   ├── map_files
│   │       └── *.map
│   └── mappings.pdf
└── misc
```

## Contributions

If there's a metabolite not covered by the parameters in the current set, we welcome all contributions subject 
to suitable validation. Contributions can be made through [issues](https://github.com/Martini-Force-Field-Initiative/M3-Metabolome/issues)
or [pull requests](https://github.com/Martini-Force-Field-Initiative/M3-Metabolome/pulls). Please note that the
decision of whether contributions can give authorship on the resulting academic paper is left to our sole discretion.

## Full dataset

The full dataset (including reference atomistic trajectories) for the metabolome is available on Zenodo [LINK]
