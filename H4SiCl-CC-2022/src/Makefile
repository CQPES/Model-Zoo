FC = gfortran
F2PY = f2py

all: LIBH4SICL

h4sicl_pipnn.o: bemsa411_5.o h4sicl_pipnn.f
	$(FC) -c h4sicl_pipnn.f

bemsa411_5.o: bemsa411_5.f
	$(FC) -c bemsa411_5.f

LIBH4SICL: bemsa411_5.o h4sicl_pipnn.o libh4sicl.f90
	$(F2PY) libh4sicl.f90 -m libh4sicl -h libh4sicl.pyf
	$(F2PY) -c libh4sicl.f90 bemsa411_5.f h4sicl_pipnn.f libh4sicl.pyf

.PHONY: clean
clean:
	rm -rf *.o *.mod *.pyf *.so
