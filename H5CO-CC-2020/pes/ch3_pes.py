import numpy as np
from gau_pes import BasePES
from h5co_pes import H5COPES

# constants
_NUM_ATOMS = 4
_EQ_H2O = np.array([
    [0.0000000000, -0.7579202604,  0.5215934515],
    [0.0000000000, -0.0000000001, -0.0657193274],
    [0.0000000000,  0.7579202612,  0.5215934508],
])  # Angstrom
_DISPLACE = 20.0  # Angstrom


class CH3PES(BasePES):
    def __init__(self) -> None:
        self.h5co_pes = H5COPES()

    def calc_energy(
        self,
        coords: np.array,
    ) -> float:
        """Calculate relative potenial energy of

        1. CH3 in H + CH3OH -> H2O + CH3 channel

        Order: H H H C
        """
        self._check_coords(_NUM_ATOMS, coords)

        # concat a fixed H2O molecule with input CH3
        # now that order of the atoms is H H | H H H C | O
        _eq_h2o = _EQ_H2O + coords.mean(axis=0) + _DISPLACE
        new_coords = np.concatenate((
            _eq_h2o[[0, 1], :],  # H H
            coords,  # H H H C
            _eq_h2o[[2, ], :],  # O
        ))

        return self.h5co_pes.calc_energy(new_coords)


if __name__ == "__main__":
    from gau_pes import GauDriver

    driver = GauDriver()
    pes = CH3PES()

    driver.write(
        energy=pes.calc_energy(driver.coords),
        gradients=(pes.calc_gradients(driver.coords)
                   if driver.derivs in [1, 2] else None),
        force_constants=(pes.calc_force_constants(
            driver.coords) if driver.derivs == 2 else None),
    )
