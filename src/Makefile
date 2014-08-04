all: mk_pcapack mk_pkg 


mk_pcapack:
	cd pcapack && make && cd ..


mk_pkg:
	mkdir -p build;
	cd build; \
	cmake ..; \
	make install


mk_pkg: mk_pcapack


clean:
	clean:
	rm -rf ./pcapack/build

