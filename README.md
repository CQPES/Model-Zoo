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

### [`CH4-MRCI-2015`](https://github.com/CQPES/Model-Zoo/tree/main/CH4-MRCI-2015)

#### Species

- CH<sub>4</sub>

#### Reference

(1) Majumder, M.; Hegger, S. E.; Dawes, R.; Manzhos, S.; Wang, X.-G.; Tucker, C., Jr; Li, J.; Guo, H. Explicitly Correlated MRCI-F12 Potential Energy Surfaces for Methane Fit with Several Permutation Invariant Schemes and Full-Dimensional Vibrational Calculations. _Mol. Phys._ **2015**, _113_ (13–14), 1823–1833. https://doi.org/10.1080/00268976.2015.1015642.

### [`H4CFO-CC-2019`](https://github.com/CQPES/Model-Zoo/tree/main/H4CFO-CC-2019)

#### Species

- H
- HF
- CH<sub>3</sub>O
- CH<sub>2</sub>OH
- CH<sub>3</sub>OH

#### Channels

- F + CH<sub>3</sub>OH → HF + CH<sub>3</sub>O
- F + CH<sub>3</sub>OH → HF + CH<sub>2</sub>OH

#### Reference

(1) Lu, D.-D.; Xie, C.-J.; Li, J.; Guo, H. Rate Coefficients and Branching Ratio for Multi-Channel Hydrogen Abstractions from CH<sub>3</sub>OH by F. Chin. _J. Chem. Phys._ **2019**, _32_ (1), 84–88. https://doi.org/10.1063/1674-0068/cjcp1811256.

(2) Lu, D.; Li, J.; Guo, H. Stereodynamical Control of Product Branching in Multi-Channel Barrierless Hydrogen Abstraction of CH<sub>3</sub>OH by F. _Chem. Sci._ **2019**, _10_ (34), 7994–8001. https://doi.org/10.1039/c9sc02445j.

(3) Weichman, M. L.; DeVine, J. A.; Babin, M. C.; Li, J.; Guo, L.; Ma, J.; Guo, H.; Neumark, D. M. Feshbach Resonances in the Exit Channel of the F + CH<sub>3</sub>OH → HF + CH<sub>3</sub>O Reaction Observed Using Transition-State Spectroscopy. _Nat. Chem._ **2017**, _9_ (10), 950–955. https://doi.org/10.1038/nchem.2804.

### [`H4CClO-CC-2020`](https://github.com/CQPES/Model-Zoo/tree/main/H4O-CC-MRCI-2022)

#### Species

- H
- HCl
- CH<sub>3</sub>O
- CH<sub>2</sub>OH
- CH<sub>3</sub>OH

#### Channels

- Cl + CH<sub>3</sub>OH → HCl + CH<sub>3</sub>O
- Cl + CH<sub>3</sub>OH → HCl + CH<sub>2</sub>OH

#### Reference

(1) Lu, D.; Li, J.; Guo, H. Comprehensive Investigations of the Cl + CH<sub>3</sub>OH → HCl + CH<sub>3</sub>O/CH<sub>2</sub>OH Reaction: Validation of Experiment and Dynamic Insights. _CCS Chem._ **2020**, _2_ (5), 882–894. https://doi.org/10.31635/ccschem.020.202000195.

### [`H5CO-CC-2020`](https://github.com/CQPES/Model-Zoo/tree/main/H5CO-CC-2020)

#### Species

- H
- H<sub>2</sub>
- H<sub>2</sub>O
- CH<sub>3</sub>
- CH<sub>3</sub>O
- CH<sub>2</sub>OH
- CH<sub>3</sub>OH

#### Channels

- H + CH<sub>3</sub>OH → H<sub>2</sub> + CH<sub>3</sub>O
- H + CH<sub>3</sub>OH → H<sub>2</sub> + CH<sub>2</sub>OH
- H + CH<sub>3</sub>OH → H<sub>2</sub>O + CH<sub>3</sub>

#### Reference

(1) Lu, D.; Behler, J.; Li, J. Accurate Global Potential Energy Surfaces for the H + CH<sub>3</sub>OH Reaction by Neural Network Fitting with Permutation Invariance. _J. Phys. Chem. A_ **2020**, _124_ (28), 5737–5745. https://doi.org/10.1021/acs.jpca.0c04182.

### [`NH4-anion-CC-2022`](https://github.com/CQPES/Model-Zoo/tree/main/NH4-anion-CC-2022)

#### Species

- H<sup>-</sub>
- H<sub>2</sub>
- NH<sub>2</sub><sup>-</sup>
- NH<sub>3</sub>

#### Channels

- H<sub>2</sub> + NH<sub>2</sub><sup>-</sup> → H<sup>-</sup> + NH<sub>3</sub>
- H<sup>-</sup> + NH<sub>3</sub> → H<sup>-</sup> + NH<sub>3</sub>
- H<sup>-</sup> + NH<sub>3</sub> → NH<sub>4</sub><sup>-</sup>

#### Reference

(1) Song, K.; Song, H.; Li, J. Validating Experiments for the Reaction H<sub>2</sub> + NH<sub>2</sub><sup>−</sup> by Dynamical Calculations on an Accurate Full-Dimensional Potential Energy Surface. _Phys. Chem. Chem. Phys._ **2022**, _24_ (17), 10160–10167. https://doi.org/10.1039/d2cp00870j.

### [`H4SiCl-CC-2022`](https://github.com/CQPES/Model-Zoo/tree/main/H4SiCl-CC-2022)

#### Species

- H
- Cl
- H<sub>2</sub>
- HCl
- SiH<sub>3</sub>
- SiH<sub>2</sub>Cl
- SiH<sub>4</sub>
- SiH<sub>3</sub>Cl

#### Channels

- Cl + SiH<sub>4</sub> → HCl + SiH<sub>3</sub>
- H + SiH<sub>3</sub>Cl → H<sub>2</sub> + SiH<sub>2</sub>Cl
- HCl + SiH<sub>3</sub> → H + SiH<sub>3</sub>Cl

#### Reference

(1) Xu, X.; Li, J. Deciphering Dynamics of the Cl + SiH<sub>4</sub> → H + SiH<sub>3</sub>Cl Reaction on a Machine Learning Made Globally Accurate Full-Dimensional Potential Energy Surface. _J. Phys. Chem. A_ **2022**, _126_ (37), 6456–6466. https://doi.org/10.1021/acs.jpca.2c05417.

### [`H4O-CC-MRCI-2022`](https://github.com/CQPES/Model-Zoo/tree/main/H4O-CC-MRCI-2022)

#### Species

- H<sub>2</sub>
- H<sub>2</sub>O

#### Channels

- H<sub>2</sub> + H′<sub>2</sub>O → HH′ + HOH′

#### Reference

(1) Li, J.; Liu, Y.; Guo, H.; Li, J. An Accurate Full-Dimensional H<sub>4</sub>O Potential Energy Surface and Dynamics of an Exchange Reaction. _Phys. Chem. Chem. Phys._ **2022**, _24_ (44), 27548–27557. https://doi.org/10.1039/d2cp04521d.

### [`H4S2-CC-2024`](https://github.com/CQPES/Model-Zoo/tree/main/H4S2-CC-2024)

#### Species

- H<sub>2</sub>S

#### Channels

- H<sub>2</sub>S + H′<sub>2</sub>S′ → HH′S + HH′S′

#### Reference

(1) Deng, F.; Liu, Y.; Li, J.; Li, J. An Accurate Full-Dimensional Potential Energy Surface of the Hydrogen Sulfide Dimer System and Its Kinetics for the Hydrogen Exchange Channel. _Mol. Phys._ **2024**. https://doi.org/10.1080/00268976.2024.2438329.
