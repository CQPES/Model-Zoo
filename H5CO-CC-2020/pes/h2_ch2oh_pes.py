import numpy as np
from gau_pes import BasePES
from h5co_pes import H5COPES

# constants
_NUM_ATOMS = 2
_EQ_CH2OH = np.array([
    [-0.8160287642, -0.0652099757, -1.0051374955],
    [0.9245894030,  0.1315186489,  1.2231422666],
    [-0.9568081607,  0.2045999444,  1.2003029743],
    [-0.0162644273, -0.0396818636,  0.7269265325],
    [0.0656483770,  0.0127229291, -0.6350665457],
])  # Angstrom
_DISPLACE = 20.0  # Angstrom


class H2PES(BasePES):
    def __init__(self) -> None:
        self.h5co_pes = H5COPES()

    def calc_energy(
        self,
        coords: np.array,
    ) -> float:
        """Calculate relative potenial energy of H2 in H2 + CH2OH channel"""

        self._check_coords(_NUM_ATOMS, coords)

        # concat a fixed CH2OH molecule with input H2
        # now that order of the atoms is H H | H H H C O
        _eq_ch2oh = _EQ_CH2OH + np.ones_like(_EQ_CH2OH) * _DISPLACE
        new_coords = np.concatenate((
            coords,
            _eq_ch2oh,
        ))

        return self.h5co_pes.calc_energy(new_coords)


if __name__ == "__main__":
    from gau_pes import GauDriver

    driver = GauDriver()
    pes = H2PES()

    driver.write(
        energy=pes.calc_energy(driver.coords),
        gradients=(pes.calc_gradients(driver.coords)
                   if driver.derivs in [1, 2] else None),
        force_constants=(pes.calc_force_constants(
            driver.coords) if driver.derivs == 2 else None),
    )
