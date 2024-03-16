import numpy as np
from gau_pes import BasePES

_NUM_ATOMS = 1


class ClPES(BasePES):
    def calc_energy(
        self,
        coords: np.array,
    ) -> float:
        """Calculate relative potential energy of F atom"""
        self._check_coords(_NUM_ATOMS, coords)

        return 0.0


if __name__ == "__main__":
    from gau_pes import GauDriver

    driver = GauDriver()
    pes = ClPES()

    driver.write(
        energy=pes.calc_energy(driver.coords),
        gradients=(pes.calc_gradients(driver.coords)
                   if driver.derivs in [1, 2] else None),
        force_constants=(pes.calc_force_constants(
            driver.coords) if driver.derivs == 2 else None),
    )
