FC = gfortran

all: sgns means

sgns:
	$(FC) -o sgns sgns.f90 -I../build/include ../build/libpcapack_static.a -lgomp

means:
	$(FC) -o means means.f90 -I../build/include ../build/libpcapack_static.a -lgomp


clean:
	rm -rf ./sgns
	rm -rf ./means

