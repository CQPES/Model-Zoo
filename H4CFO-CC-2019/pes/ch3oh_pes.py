import numpy as np
from gau_pes import BasePES
from h4cfo_pes import H4CFOPES

# constants
_NUM_ATOMS = 6
_EQ_F = np.array([
    [0.0000000000, 0.0000000000, 0.0000000000],  # F
])  # Angstrom
_DISPLACE = 20.0  # Angstrom


class CH3OHPES(BasePES):
    def __init__(self) -> None:
        self.h4cfo_pes = H4CFOPES()

    def calc_energy(
        self,
        coords: np.array,
    ) -> float:
        """Calculate relative potenial energy of

        1. CH3OH in F + CH3OH -> HF + CH3O channel
        2. CH3OH in F + CH3OH -> HF + CH2OH channel

        Order of atoms: H H H H C O
        """

        self._check_coords(_NUM_ATOMS, coords)

        # concat the input coords with a fixed F atom
        # now that order of the atoms is H H H H C | F | O
        _eq_f = _EQ_F + coords.mean(axis=0) + _DISPLACE
        new_coords = np.concatenate((
            coords[[0, 1, 2, 3, 4], :],  # H H H H C
            _eq_f,
            coords[[5], :],
        ))

        return self.h4cfo_pes.calc_energy(new_coords)


if __name__ == "__main__":
    from gau_pes import GauDriver

    driver = GauDriver()
    pes = CH3OHPES()

    driver.write(
        energy=pes.calc_energy(driver.coords),
        gradients=(pes.calc_gradients(driver.coords)
                   if driver.derivs in [1, 2] else None),
        force_constants=(pes.calc_force_constants(
            driver.coords) if driver.derivs == 2 else None),
    )
