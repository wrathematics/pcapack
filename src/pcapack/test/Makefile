FC = gfortran

all: sgns means float_equals distances

sgns:
	$(FC) -o sgns sgns.f90 -I../build/include ../build/libpcapack_static.a -lgomp

means:
	$(FC) -o means means.f90 -I../build/include ../build/libpcapack_static.a -lgomp

float_equals:
	$(FC) -o float_equals float_equals.f90 -I../build/include ../build/libpcapack_static.a -lgomp -llapack

distances:
	$(FC) -o distances distances.f90 -I../build/include ../build/libpcapack_static.a -lgomp -lblas -llapack

covariance:
	$(FC) -o covariance covariance.f90 -I../build/include ../build/libpcapack_static.a -lgomp -lblas -llapack

pca:
	$(FC) -o pca pca.f90 -I../build/include ../build/libpcapack_static.a -lgomp -lblas -llapack


clean:
	rm -rf ./sgns
	rm -rf ./means
	rm -rf ./float_equals
	rm -rf ./distances
	rm -rf ./covariance
	rm -rf ./pca

