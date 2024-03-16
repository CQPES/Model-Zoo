import numpy as np
from gau_pes import BasePES
from h4cclo_pes import H4CClOPES

# constants
_NUM_ATOMS = 5
_EQ_HCl = np.array([
    [0.0000000000, 0.0000000000, -1.2417634173],
    [0.0000000000, 0.0000000000, 0.0353037266],
])  # Angstrom
_DISPLACE = 20.0  # Angstrom


class H3COPES(BasePES):
    def __init__(self) -> None:
        self.h4cclo_pes = H4CClOPES()

    def calc_energy(
        self,
        coords: np.array,
    ) -> float:
        """Calculate relative potenial energy of

        1. CH3O in Cl + CH3OH -> HCl + CH3O channel
        2. CH2OH in Cl + CH3OH -> HCl + CH2OH channel

        Order: H H H C O
        """

        self._check_coords(_NUM_ATOMS, coords)

        # concat input coords with a fixed HF molecule
        # now that order of the atoms is H H H | H | C | Cl | O
        _eq_hcl = _EQ_HCl + np.ones_like(_EQ_HCl) * _DISPLACE
        new_coords = np.concatenate((
            coords[[0, 1, 2], :],  # H H H
            _eq_hcl[[0], :],  # H in HCl
            coords[[3], :],  # C
            _eq_hcl[[1], :],  # Cl in HCl
            coords[[4], :],  # O
        ))

        return self.h4cclo_pes.calc_energy(new_coords)


if __name__ == "__main__":
    from gau_pes import GauDriver

    driver = GauDriver()
    pes = H3COPES()

    driver.write(
        energy=pes.calc_energy(driver.coords),
        gradients=(pes.calc_gradients(driver.coords)
                   if driver.derivs in [1, 2] else None),
        force_constants=(pes.calc_force_constants(
            driver.coords) if driver.derivs == 2 else None),
    )
