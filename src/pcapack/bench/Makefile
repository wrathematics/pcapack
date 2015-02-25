CC=gcc -std=c99
CFLAGS=-L../build -lpcapack_static -lm -lblas -llapack -lgomp

all: cov prcomp xpose

cov:
	$(CC) -o cov cov.c $(CFLAGS)

prcomp:
	$(CC) -o prcomp prcomp.c $(CFLAGS)

xpose:
	$(CC) -o xpose xpose.c $(CFLAGS)

.PHONY: cov prcomp xpose

clean:
	rm -rf cov prcomp xpose *~
