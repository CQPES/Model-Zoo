FC = gfortran
F2PY = f2py

all: LIBH4S2

h4s2_pipnn.o: bemsa42_6.o h4s2_pipnn.f
	$(FC) -c h4s2_pipnn.f

bemsa42_6.o: bemsa42_6.f
	$(FC) -c bemsa42_6.f

LIBH4S2: bemsa42_6.o h4s2_pipnn.o libh4s2.f90
	$(F2PY) libh4s2.f90 -m libh4s2 -h libh4s2.pyf
	$(F2PY) -c libh4s2.f90 bemsa42_6.f h4s2_pipnn.f libh4s2.pyf

.PHONY: clean
clean:
	rm -rf *.o *.mod *.pyf *.so
