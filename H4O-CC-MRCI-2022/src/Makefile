FC = gfortran
F2PY = f2py

all: LIBH4O

h4o_pipnn.o: h4o_pipnn.f
	$(FC) -c h4o_pipnn.f

LIBH4O: h4o_pipnn.o libh4o.f90
	$(F2PY) libh4o.f90 -m libh4o -h libh4o.pyf
	$(F2PY) -c libh4o.f90 h4o_pipnn.f libh4o.pyf

.PHONY: clean
clean:
	rm -rf *.o *.mod *.pyf *.so
