subroutine init()
    use nnparam
    implicit none

    call pes_init()
end subroutine init

subroutine calc_energy(coords, natom, vpes)
    implicit none

    integer, intent(in) :: natom                 !> number of atoms
    real(kind=8), intent(in) :: coords(natom, 3) !> Angstrom
    real(kind=8), intent(out) :: vpes            !> Hartree

    real(kind=8), parameter :: Eh2eV = 27.211399d0

    real(kind=8) :: vpesa, vpesb, vpesc

    call ch4pipNN(transpose(coords), vpes, vpesa, vpesb, vpesc)
    vpes = vpes / Eh2eV
end subroutine calc_energy
