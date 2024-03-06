FC = gfortran
F2PY = f2py

all: LIBH5CO

h5co_pipnn.o: bemsa511_5.o h5co_pipnn.f
	$(FC) -c h5co_pipnn.f

bemsa511_5.o: bemsa511_5.f90
	$(FC) -c bemsa511_5.f90

LIBH5CO: bemsa511_5.o h5co_pipnn.o libh5co.f90
	$(F2PY) libh5co.f90 -m libh5co -h libh5co.pyf
	$(F2PY) -c libh5co.f90 bemsa511_5.f90 h5co_pipnn.f libh5co.pyf

.PHONY: clean
clean:
	rm -rf *.o *.mod *.pyf *.so
