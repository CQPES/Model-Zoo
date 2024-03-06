# Model-Zoo

Dataset and PES in Published Papers.

Maintainer: mizu-bai

## Prerequests

If you want to use the Gaussian interface, make sure you have `Gaussian-PES` installed.

Run

```shell
$ pip install git+git@github.com:CQPES/Gaussian-PES.git
```

or

```shell
$ git clone git@github.com:CQPES/Gaussian-PES.git
$ cd Gaussian-PES
$ pip install .
```

## PES

### `CH4-MRCI-2015`

species

- CH<sub>4</sub>

Majumder, M.; Hegger, S. E.; Dawes, R.; Manzhos, S.; Wang, X.-G.; Tucker, C., Jr; Li, J.; Guo, H. Explicitly Correlated MRCI-F12 Potential Energy Surfaces for Methane Fit with Several Permutation Invariant Schemes and Full-Dimensional Vibrational Calculations. Mol. Phys. 2015, 113 (13–14), 1823–1833. https://doi.org/10.1080/00268976.2015.1015642.

### `H4O-CC-MRCI-2022`

species

- H<sub>2</sub>
- H<sub>2</sub>O

reactions

- H<sub>2</sub> + H′<sub>2</sub> → HH′ + HOH′

Li, J.; Liu, Y.; Guo, H.; Li, J. An Accurate Full-Dimensional H4O Potential Energy Surface and Dynamics of an Exchange Reaction. Phys. Chem. Chem. Phys. 2022, 24 (44), 27548–27557. https://doi.org/10.1039/d2cp04521d.

### `H5CO-CC-2020`

species

- H
- H<sub>2</sub>
- H<sub>2</sub>O
- CH<sub>3</sub>
- CH<sub>3</sub>O
- CH<sub>2</sub>OH
- CH<sub>3</sub>OH

reactions

- H + CH<sub>3</sub>OH → H<sub>2</sub> + CH<sub>3</sub>O
- H + CH<sub>3</sub>OH → H<sub>2</sub> + CH<sub>2</sub>OH
- H + CH<sub>3</sub>OH → H<sub>2</sub>O + CH<sub>3</sub>

Lu, D.; Behler, J.; Li, J. Accurate Global Potential Energy Surfaces for the H + CH3OH Reaction by Neural Network Fitting with Permutation Invariance. J. Phys. Chem. A 2020, 124 (28), 5737–5745. https://doi.org/10.1021/acs.jpca.0c04182.
