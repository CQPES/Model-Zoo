import libh4o
import numpy as np
from gau_pes import BasePES

_NUM_ATOMS = 5


class H4OPES(BasePES):
    def __init__(self) -> None:
        libh4o.init()

    def calc_energy(
        self,
        coords: np.array,
    ) -> float:
        """Calculate potential energy of H2 + H2O system"""
        self._check_coords(_NUM_ATOMS, coords)

        return libh4o.calc_energy(coords)


if __name__ == "__main__":
    from gau_pes import GauDriver

    driver = GauDriver()
    pes = H4OPES()

    driver.write(
        energy=pes.calc_energy(driver.coords),
        gradients=(pes.calc_gradients(driver.coords)
                   if driver.derivs in [1, 2] else None),
        force_constants=(pes.calc_force_constants(
            driver.coords) if driver.derivs == 2 else None),
    )
