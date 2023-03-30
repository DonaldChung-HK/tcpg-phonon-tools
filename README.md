# tcpg-phonon-tools
This is a python package that contains easy python tools I wrote for:
- easier setup of INPUT for `CASTEP`, `CRYSTAL` and `DFTB+`
- one-line execution of optimisation and phonon
- wrapper for `DFTB+` to use ASE `SocketIOCalculator` for optimisation and parallel job `phonopy` calculation
- wrapper for data processing and visualisation for `ABINS` and Earth Mover Distance `scipy`
## install
```
pip install .
```
### notes
#### Requires `Python >= 3.10` but MANTID is locked to `Python 3.8`.
This tools requires MANTID but there are no `> 3.8` code in the part for mantid and analysis of data. The main `> 3.10` part is for importing non python file as this package employs quite a few text template files for setting up input for abinitio softwares
#### Euphonic mantid issue
Check that the `numpy` version is matching the one euphonic uses, if not the C wrapper will fail to run and fallback to pure Python calculation 
## Usage
Here are the command line tools 
### `dftbp-socket-setup`
Set up the `dftb_in.hsd` for optimisation
### `ase-dftbp-optimise`
Optimise using DFTB+ with ASE optimisers
### `castep-phonopy-setup`
Set up phonon calculation with `phonopy` using `CASTEP` with `run.slurm` for SCARF
### `dftbp-phonopy-setup`
Set up phonon calculation with `phonopy` using `CASTEP` with `run.slurm` for SCARF or `run.py` for single machine execution
### `crystal-opt-setup`
Setup `CRYSTAL` calculation from `.cif` files as only `cif2cell` can get the symmetry of a system with `run.slurm` for SCARF
### `crystal-phonopy-setup`
Set up phonon calculation with `phonopy` using `CRYSTAL` with `run.slurm` for SCARF
### `abins-mantid-runner`
Calculate INS spectra using ABINS and get data into `.csv` for further analysis
### `emd-chart` and `emd-chart-multi`
Plot the spectra with the calculated spectra normalised with reference spectra with earth mover distance. `emd-chart-multi` is for overplotting multiple calculated spectra.