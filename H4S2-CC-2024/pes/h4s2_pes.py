import libh4s2
import numpy as np
from gau_pes import BasePES

_NUM_ATOMS = 6


class H4S2PES(BasePES):
    def __init__(self):
        libh4s2.init()

    def calc_energy(
        self,
        coords: np.ndarray,
    ) -> float:
        """Calculate potential energy of H2S + H2S system"""
        self._check_coords(_NUM_ATOMS, coords)

        return libh4s2.calc_energy(coords)

if __name__ == "__main__":
    from gau_pes import GauDriver

    driver = GauDriver()
    pes = H4S2PES()

    driver.write(
        energy=pes.calc_energy(driver.coords),
        gradients=(pes.calc_gradients(driver.coords)
                   if driver.derivs in [1, 2] else None),
        force_constants=(pes.calc_force_constants(
            driver.coords) if driver.derivs == 2 else None),
    )
