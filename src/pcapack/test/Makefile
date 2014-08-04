FC = gfortran

all: sgns

sgns:
	$(FC) -o sgns sgns.f90 -I../build/include ../build/libpcapack_static.a -lgomp

clean:
	rm -rf ./sgns

