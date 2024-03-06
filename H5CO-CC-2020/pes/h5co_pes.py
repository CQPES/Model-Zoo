import libh5co
import numpy as np
from gau_pes import BasePES

_NUM_ATOMS = 7


class H5COPES(BasePES):
    def __init__(self) -> None:
        libh5co.init()

    def calc_energy(
        self,
        coords: np.array,
    ) -> float:
        """Calculate potential energy of H + CH3OH system"""
        self._check_coords(_NUM_ATOMS, coords)

        return libh5co.calc_energy(coords)


if __name__ == "__main__":
    from gau_pes import GauDriver

    driver = GauDriver()
    pes = H5COPES()

    driver.write(
        energy=pes.calc_energy(driver.coords),
        gradients=(pes.calc_gradients(driver.coords)
                   if driver.derivs in [1, 2] else None),
        force_constants=(pes.calc_force_constants(
            driver.coords) if driver.derivs == 2 else None),
    )
