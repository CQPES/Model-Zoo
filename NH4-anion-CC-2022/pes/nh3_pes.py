import numpy as np
from gau_pes import BasePES
from h4n_anion_pes import H4NAnionPES

# constants
_NUM_ATOMS = 4
_EQ_H = np.array([
    [0.0000000000, 0.0000000000, 0.0000000000],  # H
])  # Angstrom
_DISPLACE = 20.0  # Angstrom


class NH3PES(BasePES):
    def __init__(self) -> None:
        self.h4n_anion_pes = H4NAnionPES()

    def calc_energy(
        self,
        coords: np.array,
    ) -> float:
        """Calculate relative potenial energy of

        1. NH3 in H2 + NH2^- -> H^- + NH3 channel

        Order: H H H N
        """

        self._check_coords(_NUM_ATOMS, coords)

        # concat a fixed H atom with input NH3
        # now the order of atoms is H | H H H N
        # and the H atom is displaced by _DISPLACE in each axis
        new_coords = np.concatenate((
            _EQ_H + coords.mean(axis=0) + _DISPLACE,  # H
            coords,  # H H H N
        ))

        return self.h4n_anion_pes.calc_energy(new_coords)


if __name__ == "__main__":
    from gau_pes import GauDriver

    driver = GauDriver()
    pes = NH3PES()

    driver.write(
        energy=pes.calc_energy(driver.coords),
        gradients=(pes.calc_gradients(driver.coords)
                   if driver.derivs in [1, 2] else None),
        force_constants=(pes.calc_force_constants(
            driver.coords) if driver.derivs == 2 else None),
    )
