import numpy as np
from gau_pes import BasePES
from h5co_pes import H5COPES

# constants
_NUM_ATOMS = 3
_EQ_CH3 = np.array([
    [-0.0000000000, -0.0000004806, -1.0777339555],  # H
    [0.9333449540, -0.0000004806, 0.5388669718],  # H
    [-0.9333449779, -0.0000004806, 0.5388669408],  # H
    [0.0000000020, 0.0000000000, 0.0000000036],  # C
])  # Angstrom
_DISPLACE = 20.0  # Angstrom


class H2OPES(BasePES):
    def __init__(self) -> None:
        self.h5co_pes = H5COPES()

    def calc_energy(
        self,
        coords: np.array,
    ) -> float:
        """Calculate relative potenial energy of

        1. H2O in H + CH3OH -> H2O + CH3 channel

        Order: H H O
        """
        self._check_coords(_NUM_ATOMS, coords)

        # concat a fixed CH3 molecule with input H2O
        # now that order of the atoms is H H | H H H C | O
        _eq_ch3 = _EQ_CH3 + coords.mean(axis=0) + _DISPLACE
        new_coords = np.concatenate((
            coords[[0, 1], :],  # H H
            _eq_ch3,  # H H H C
            coords[[2, ], :],  # O
        ))

        return self.h5co_pes.calc_energy(new_coords)


if __name__ == "__main__":
    from gau_pes import GauDriver

    driver = GauDriver()
    pes = H2OPES()

    driver.write(
        energy=pes.calc_energy(driver.coords),
        gradients=(pes.calc_gradients(driver.coords)
                   if driver.derivs in [1, 2] else None),
        force_constants=(pes.calc_force_constants(
            driver.coords) if driver.derivs == 2 else None),
    )
