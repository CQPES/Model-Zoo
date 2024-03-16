import libh4cclo
import numpy as np
from gau_pes import BasePES

_NUM_ATOMS = 7


class H4CClOPES(BasePES):
    def __init__(self) -> None:
        libh4cclo.init()

    def calc_energy(
        self,
        coords: np.array,
    ) -> float:
        """Calculate potential energy of Cl + CH3OH system

        Order of atoms: H H H H C Cl O
        """
        self._check_coords(_NUM_ATOMS, coords)

        return libh4cclo.calc_energy(coords)


if __name__ == "__main__":
    from gau_pes import GauDriver

    driver = GauDriver()
    pes = H4CClOPES()

    driver.write(
        energy=pes.calc_energy(driver.coords),
        gradients=(pes.calc_gradients(driver.coords)
                   if driver.derivs in [1, 2] else None),
        force_constants=(pes.calc_force_constants(
            driver.coords) if driver.derivs == 2 else None),
    )
