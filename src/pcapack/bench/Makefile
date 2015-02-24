CC=gcc -std=c99
CFLAGS=-L../build -lpcapack_static -lm -lblas -llapack -lgomp

all: cov xpose

cov:
	$(CC) -o cov cov.c $(CFLAGS)

xpose:
	$(CC) -o xpose xpose.c $(CFLAGS)

.PHONY: cov xpose

clean:
	rm -rf cov xpose *~
