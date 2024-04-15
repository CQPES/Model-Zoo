import numpy as np
from gau_pes import BasePES
from h4sicl_pes import H4SiClPES

# constants
_NUM_ATOMS = 5
_EQ_Cl = np.array([
    [0.00000000, 0.00000000, 0.00000000],
])  # Angstrom
_DISPLACE = 20.0  # Angstrom


class SiH4PES(BasePES):
    def __init__(self) -> None:
        self.h4sicl_pes = H4SiClPES()

    def calc_energy(
        self,
        coords: np.array,
    ) -> float:
        """Calculate relative potential energy of

        1. SiH4 in Cl + SiH4 -> H + SiH3Cl channel

        Order: H H H H Si
        """

        self._check_coords(_NUM_ATOMS, coords)

        # concat a fixed Cl atom with input SiH4 molecule
        # now the order of atoms is H H H H Si | Cl
        new_coords = np.concatenate((
            coords,  # H H H H Si
            _EQ_Cl + coords.mean(axis=0) + _DISPLACE,  # Cl
        ))

        return self.h4sicl_pes.calc_energy(new_coords)


if __name__ == "__main__":
    from gau_pes import GauDriver

    driver = GauDriver()
    pes = SiH4PES()

    driver.write(
        energy=pes.calc_energy(driver.coords),
        gradients=(pes.calc_gradients(driver.coords)
                   if driver.derivs in [1, 2] else None),
        force_constants=(pes.calc_force_constants(
            driver.coords) if driver.derivs == 2 else None),
    )
