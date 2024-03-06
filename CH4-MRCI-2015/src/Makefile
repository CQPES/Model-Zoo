FC = gfortran
F2PY = f2py

all: LIBCH4

ch4_pipnn.o: ch4_pipnn.f
	$(FC) -c ch4_pipnn.f

LIBCH4: ch4_pipnn.o libch4.f90
	$(F2PY) libch4.f90 -m libch4 -h libch4.pyf
	$(F2PY) -c libch4.f90 ch4_pipnn.f libch4.pyf

.PHONY: clean
clean:
	rm -rf *.o *.mod *.pyf *.so
