FC = gfortran
F2PY = f2py

all: LIBH4CFO

h4cfo_pipnn.o: h4cfo_pipnn.f
	$(FC) -c h4cfo_pipnn.f

LIBH4CFO: h4cfo_pipnn.o libh4cfo.f90
	$(F2PY) libh4cfo.f90 -m libh4cfo -h libh4cfo.pyf
	$(F2PY) -c libh4cfo.f90 h4cfo_pipnn.f libh4cfo.pyf

.PHONY: clean
clean:
	rm -rf *.o *.mod *.pyf *.so
