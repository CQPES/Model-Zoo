FC = gfortran
F2PY = f2py

all: LIBCH3CN

ch3cn_pipnn.o: ch3cn_pipnn.f90
	$(FC) -c ch3cn_pipnn.f90

LIBCH3CN: ch3cn_pipnn.o libch3cn.f90
	$(F2PY) libch3cn.f90 -m libch3cn -h libch3cn.pyf
	$(F2PY) -c libch3cn.f90 ch3cn_pipnn.f90 libch3cn.pyf

.PHONY: clean
clean:
	rm -rf *.o *.mod *.pyf *.so
