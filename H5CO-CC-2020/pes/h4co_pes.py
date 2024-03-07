import numpy as np
from gau_pes import BasePES
from h5co_pes import H5COPES

# constants
_NUM_ATOMS = 6
_EQ_H = np.array([
    [0.0000000000, 0.0000000000, -0.3708705760],  # H
])  # Angstrom
_DISPLACE = 20.0  # Angstrom


class H4COPES(BasePES):
    def __init__(self) -> None:
        self.h5co_pes = H5COPES()

    def calc_energy(
        self,
        coords: np.array,
    ) -> float:
        """Calculate relative potenial energy of

        1. CH3OH in H + CH3OH -> H2 + CH3O channel
        2. CH3OH in H + CH3OH -> H2 + CH2OH channel

        Order: H H H H C O
        """
        self._check_coords(_NUM_ATOMS, coords)

        # concat input coords with a fixed H atom
        # now that order of the atoms is H | H H H H C O
        _eq_h = _EQ_H + np.ones_like(_EQ_H) * _DISPLACE
        new_coords = np.concatenate((
            _eq_h,
            coords,
        ))

        return self.h5co_pes.calc_energy(new_coords)


if __name__ == "__main__":
    from gau_pes import GauDriver

    driver = GauDriver()
    pes = H4COPES()

    driver.write(
        energy=pes.calc_energy(driver.coords),
        gradients=(pes.calc_gradients(driver.coords)
                   if driver.derivs in [1, 2] else None),
        force_constants=(pes.calc_force_constants(
            driver.coords) if driver.derivs == 2 else None),
    )
