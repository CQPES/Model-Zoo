import numpy as np
from gau_pes import BasePES
from h4s2_pes import H4S2PES

# constants
_NUM_ATOMS = 3
_EQ_H2S = np.array([
    [0.00000000, -0.96543800, 0.87283600],  # H
    [0.00000000,  0.96543800, 0.87283600],  # H
    [0.00000000,  0.00000000,-0.05483400],  # S
])  # Angstrom
_DISPLACE = 20.0  # Angstrom


class H2SPES(BasePES):
    def __init__(self) -> None:
        self.h4s2_pes = H4S2PES()

    def calc_energy(
        self,
        coords: np.ndarray,
    ):
        """Calculate relative potential energy of H2S in H2S + H2S system"""

        self._check_coords(_NUM_ATOMS, coords)

        # concat a fixed H2S molecule with input H2S
        # now the order of atoms is H H | H H S | S
        new_coords = np.concatenate((
            coords[0: 2], # H H
            _EQ_H2S + coords.mean(axis=0) + _DISPLACE,  # H H S
            coords[2:],  # S
        ))

        return self.h4s2_pes.calc_energy(new_coords)


if __name__ == "__main__":
    from gau_pes import GauDriver

    driver = GauDriver()
    pes = H2SPES()

    driver.write(
        energy=pes.calc_energy(driver.coords),
        gradients=(pes.calc_gradients(driver.coords)
                   if driver.derivs in [1, 2] else None),
        force_constants=(pes.calc_force_constants(
            driver.coords) if driver.derivs == 2 else None),
    )
