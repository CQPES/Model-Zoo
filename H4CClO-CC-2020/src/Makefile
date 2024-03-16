FC = gfortran
F2PY = f2py

all: LIBH4CCLO

h4cclo_pipnn.o: h4cclo_pipnn.f
	$(FC) -c h4cclo_pipnn.f

LIBH4CCLO: h4cclo_pipnn.o libh4cclo.f90
	$(F2PY) libh4cclo.f90 -m libh4cclo -h libh4cclo.pyf
	$(F2PY) -c libh4cclo.f90 h4cclo_pipnn.f libh4cclo.pyf

.PHONY: clean
clean:
	rm -rf *.o *.mod *.pyf *.so
