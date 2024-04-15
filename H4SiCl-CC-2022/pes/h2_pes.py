import numpy as np
from gau_pes import BasePES
from h4sicl_pes import H4SiClPES

# constants
_NUM_ATOMS = 2
_EQ_SiH2Cl = np.array([
    [-1.2227527348,  0.6176043148,  1.5906107660],  # H
    [ 1.2227315992,  0.6176150217,  1.5906490954],  # H
    [ 0.0000001226, -0.0549877359,  1.0964472683],  # Si
    [ 0.0000005038,  0.0084430395, -0.9590384119],  # Cl
])  # Angstrom
_DISPLACE = 20.0  # Angstrom


class H2PES(BasePES):
    def __init__(self) -> None:
        self.h4sicl_pes = H4SiClPES()

    def calc_energy(
        self,
        coords: np.array,
    ) -> float:
        """Calculate relative potential energy of

        1. H2 in H + SiH3Cl -> H2 + SiH2Cl channel

        Order: H H
        """

        self._check_coords(_NUM_ATOMS, coords)

        # concat a fixed SiH2Cl molecule with input H2 molecule
        # now the order of atoms is H H | H H Si Cl
        new_coords = np.concatenate((
            coords,  # H H
            _EQ_SiH2Cl + coords.mean(axis=0) + _DISPLACE,  # H H Si Cl
        ))

        return self.h4sicl_pes.calc_energy(new_coords)


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
