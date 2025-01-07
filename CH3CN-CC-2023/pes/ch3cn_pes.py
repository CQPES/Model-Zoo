import libch3cn
import numpy as np
from gau_pes import BasePES

_NUM_ATOMS = 6


class CH3CNPES(BasePES):
    def __init__(self) -> None:
        libch3cn.init()

    def calc_energy(
            self,
            coords: np.ndarray,
        ) -> float:
        """Calculate potential energy of CH3CN <-> CH3NC system"""
        self._check_coords(_NUM_ATOMS, coords)

        return libch3cn.calc_energy(coords)


if __name__ == "__main__":
    from gau_pes import GauDriver

    driver = GauDriver()
    pes = CH3CNPES()

    driver.write(
        energy=pes.calc_energy(driver.coords),
        gradients=(pes.calc_gradients(driver.coords)
                   if driver.derivs in [1, 2] else None),
        force_constants=(pes.calc_force_constants(
            driver.coords) if driver.derivs == 2 else None),
    )
