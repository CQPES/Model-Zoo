import numpy as np
from gau_pes import BasePES
from h4cfo_pes import H4CFOPES

# constatns
_NUM_ATOMS = 2
_EQ_CH3O = np.array([
    [0.0001793492,  1.0575127554, 0.9587567996],  # H
    [-0.9061984573, -0.4549143972, 1.0914602192],  # H
    [0.9060632119, -0.4551919446, 1.0914477449],  # H
    [-0.0000026895, -0.0075814712, 0.6704905509],  # C
    [-0.0000007594, -0.0035948704, -0.7012682718],  # O
])
_DISPLACE = 20.0  # Angstrom


class HFPES(BasePES):
    def __init__(self) -> None:
        self.h4cfo_pes = H4CFOPES()

    def calc_energy(
        self,
        coords: np.array,
    ) -> float:
        """Calculate relative potenial energy of

        1. HF in F + CH3OH -> HF + CH3O channel
        2. HF in F + CH3OH -> HF + CH2OH channel

        Order of atoms: H F
        """

        self._check_coords(_NUM_ATOMS, coords)

        # concat the input coords with a fixed F atom
        # now that order of the atoms is H H H | H | C | F | O
        _eq_ch3o = _EQ_CH3O + coords.mean(axis=0) + _DISPLACE
        new_coords = np.concatenate((
            _eq_ch3o[[0, 1, 2], :],  # H H H
            coords[[0], :],  # H in HF
            _eq_ch3o[[3], :],  # C
            coords[[1], :],  # F in HF
            _eq_ch3o[[4], :],  # O
        ))

        return self.h4cfo_pes.calc_energy(new_coords)


if __name__ == "__main__":
    from gau_pes import GauDriver

    driver = GauDriver()
    pes = HFPES()

    driver.write(
        energy=pes.calc_energy(driver.coords),
        gradients=(pes.calc_gradients(driver.coords)
                   if driver.derivs in [1, 2] else None),
        force_constants=(pes.calc_force_constants(
            driver.coords) if driver.derivs == 2 else None),
    )
