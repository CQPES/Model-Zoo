import numpy as np
from gau_pes import BasePES
from h4sicl_pes import H4SiClPES

# constants
_NUM_ATOMS = 4
_EQ_H2 = np.array([
    [0.0000000000, 0.0000000000, -0.3708705760],
    [0.0000000000, 0.0000000000,  0.3708705760],
])  # Angstrom
_DISPLACE = 20.0  # Angstrom


class SiH2ClPES(BasePES):
    def __init__(self) -> None:
        self.h4sicl_pes = H4SiClPES()

    def calc_energy(
        self,
        coords: np.array,
    ) -> float:
        """Calculate relative potential energy of

        1. SiH2Cl in H + SiH3Cl -> H2 + SiH2Cl channel

        Order: H H Si Cl
        """

        self._check_coords(_NUM_ATOMS, coords)

        # concat a fixed H2 molecule with input SiH2Cl molecule
        # now the order of atoms is H H | H H Si Cl
        new_coords = np.concatenate((
            _EQ_H2 + coords.mean(axis=0) + _DISPLACE,  # H H
            coords,  # H H Si Cl
        ))

        return self.h4sicl_pes.calc_energy(new_coords)


if __name__ == "__main__":
    from gau_pes import GauDriver

    driver = GauDriver()
    pes = SiH2ClPES()

    driver.write(
        energy=pes.calc_energy(driver.coords),
        gradients=(pes.calc_gradients(driver.coords)
                   if driver.derivs in [1, 2] else None),
        force_constants=(pes.calc_force_constants(
            driver.coords) if driver.derivs == 2 else None),
    )
