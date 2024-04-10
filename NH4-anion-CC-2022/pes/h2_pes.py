import numpy as np
from gau_pes import BasePES
from h4n_anion_pes import H4NAnionPES

# constants
_NUM_ATOMS = 2
_EQ_NH2_ANION = np.array([
    [0.0000000000,  0.0000000031, -0.0817721211],  # H
    [0.0000000000,  0.8004020112,  0.5681675516],  # H
    [0.0000000000, -0.8004020549,  0.5681675162],  # N
])  # Angstrom
_DISPLACE = 20.0  # Angstrom


class H2PES(BasePES):
    def __init__(self) -> None:
        self.h4n_anion_pes = H4NAnionPES()

    def calc_energy(
        self,
        coords: np.array,
    ) -> float:
        """Calculate relative potenial energy of

        1. H2 in H2 + NH2^- -> H^- + NH3 channel

        Order: H H N
        """

        self._check_coords(_NUM_ATOMS, coords)

        # concat a fixed NH2 anion with input H2
        # now the order of atoms is H H | H H N
        new_coords = np.concatenate((
            coords,  # H H
            _EQ_NH2_ANION + coords.mean(axis=0) + _DISPLACE,  # H H N
        ))

        return self.h4n_anion_pes.calc_energy(new_coords)


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
