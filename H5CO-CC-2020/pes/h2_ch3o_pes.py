import numpy as np
from gau_pes import BasePES
from h5co_pes import H5COPES

# constants
_NUM_ATOMS = 2
_EQ_CH3O = np.array([
    [0.0001793492,  1.0575127554,  0.9587567996],
    [-0.9061984573, -0.4549143972,  1.0914602192],
    [0.9060632119, -0.4551919446,  1.0914477449],
    [-0.0000026895, -0.0075814712,  0.6704905509],
    [-0.0000007594, -0.0035948704, -0.7012682718],
])  # Angstrom
_DISPLACE = 20.0  # Angstrom


class H2PES(BasePES):
    def __init__(self) -> None:
        self.h5co_pes = H5COPES()

    def calc_energy(
        self,
        coords: np.array,
    ) -> float:
        """Calculate relative potenial energy of H2 in H2 + CH3O channel"""

        self._check_coords(_NUM_ATOMS, coords)

        # concat a fixed CH3O molecule with input H2
        # now that order of the atoms is H H | H H H C O
        _eq_ch3o = _EQ_CH3O + np.ones_like(_EQ_CH3O) * _DISPLACE
        new_coords = np.concatenate((
            coords,
            _eq_ch3o,
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
