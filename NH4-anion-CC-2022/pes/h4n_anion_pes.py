import libh4n_anion
import numpy as np
from gau_pes import BasePES

_NUM_ATOMS = 5


class H4NAnionPES(BasePES):
    def __init__(self) -> None:
        libh4n_anion.init()

    def calc_energy(
        self,
        coords: np.array,
    ) -> float:
        """Calculate potential energy of H2 + NH2 anion system"""
        self._check_coords(_NUM_ATOMS, coords)

        return libh4n_anion.calc_energy(coords)


if __name__ == "__main__":
    from gau_pes import GauDriver

    driver = GauDriver()
    pes = H4NAnionPES()

    driver.write(
        energy=pes.calc_energy(driver.coords),
        gradients=(pes.calc_gradients(driver.coords)
                   if driver.derivs in [1, 2] else None),
        force_constants=(pes.calc_force_constants(
            driver.coords) if driver.derivs == 2 else None),
    )
