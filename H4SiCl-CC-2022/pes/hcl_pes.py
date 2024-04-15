import numpy as np
from gau_pes import BasePES
from h4sicl_pes import H4SiClPES

# constants
_NUM_ATOMS = 2
_EQ_SiH3 = np.array([
    [ 0.076396,  1.407745,  0.370732],  # H
    [ 1.180967, -0.769932,  0.370791],  # H
    [-1.257347, -0.637626,  0.370792],  # H
    [-0.000001, -0.000013, -0.079451],  # Si
])  # Angstrom
_DISPLACE = 20.0  # Angstrom


class HClPES(BasePES):
    def __init__(self) -> None:
        self.h4sicl_pes = H4SiClPES()

    def calc_energy(
        self,
        coords: np.array,
    ) -> float:
        """Calculate relative potential energy of

        1. HCl in H + SiH3Cl -> HCl + SiH3 channel

        Order: H Cl
        """

        self._check_coords(_NUM_ATOMS, coords)

        # concat a fixed SiH3 molecule with input HCl molecule
        # now the order of atoms is H | H H H Si | Cl
        new_coords = np.concatenate((
            coords[0: 1],  # H
            _EQ_SiH3 + coords.mean(axis=0) + _DISPLACE,  # H H H Si
            coords[1:],  # Cl
        ))

        return self.h4sicl_pes.calc_energy(new_coords)


if __name__ == "__main__":
    from gau_pes import GauDriver

    driver = GauDriver()
    pes = HClPES()

    driver.write(
        energy=pes.calc_energy(driver.coords),
        gradients=(pes.calc_gradients(driver.coords)
                   if driver.derivs in [1, 2] else None),
        force_constants=(pes.calc_force_constants(
            driver.coords) if driver.derivs == 2 else None),
    )
