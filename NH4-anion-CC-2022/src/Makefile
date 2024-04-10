FC = gfortran
F2PY = f2py

all: LIBH4NANION

h4n_anion_pipnn.o: h4n_anion_pipnn.f
	$(FC) -c h4n_anion_pipnn.f

LIBH4NANION: h4n_anion_pipnn.o libh4n_anion.f90
	$(F2PY) libh4n_anion.f90 -m libh4n_anion -h libh4n_anion.pyf
	$(F2PY) -c libh4n_anion.f90 h4n_anion_pipnn.f libh4n_anion.pyf

.PHONY: clean
clean:
	rm -rf *.o *.mod *.pyf *.so
