import numpy as np
from gau_pes import BasePES
from h4o_pes import H4OPES

# constants
_NUM_ATOMS = 2
_EQ_H2O = np.array([
    [-0.00000002, 0.00000000, -0.00000012],
    [0.92850569, 0.00000000,  1.19804603],
    [-0.00000008, 0.00000000,  0.95882669],
])  # Angstrom
_DISPLACE = 20.0  # Angstrom


class H2PES(BasePES):
    def __init__(self) -> None:
        self.h4o_pes = H4OPES()

    def calc_energy(
        self,
        coords: np.array,
    ) -> float:
        """Calculate relative potential energy of H2O in H2 + H2O system"""

        self._check_coords(_NUM_ATOMS, coords)

        # concat a fixed H2O molecule with input H2
        # now the order of atoms is H H | H H O
        new_coords = np.concatenate((
            coords,  # H H
            _EQ_H2O + coords.mean(axis=0) + _DISPLACE,  # H H O
        ))

        return self.h4o_pes.calc_energy(new_coords)


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
