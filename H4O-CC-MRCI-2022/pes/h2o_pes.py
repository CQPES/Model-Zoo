import numpy as np
from gau_pes import BasePES
from h4o_pes import H4OPES

# constants
_NUM_ATOMS = 3
_EQ_H2 = np.array([
    [0.00000000, 0.00000000, 0.00000000],
    [0.00000000, 0.00000000, 0.74172494],
])  # Angstrom
_DISPLACE = 20.0  # Angstrom


class H2OPES(BasePES):
    def __init__(self) -> None:
        self.h4o_pes = H4OPES()

    def calc_energy(
        self,
        coords: np.array,
    ) -> float:
        """Calculate relative potential energy of H2O in H2 + H2O system"""

        self._check_coords(_NUM_ATOMS, coords)

        # concat a fixed H2 molecule with input H2O
        # now the order of atoms is H H | H H O
        # and the H2 molecule is displaced by _DISPLACE in each axis
        new_coords = np.concatenate((
            _EQ_H2 + coords.mean(axis=0) + _DISPLACE,  # H H
            coords,  # H H O
        ))

        return self.h4o_pes.calc_energy(new_coords)


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
