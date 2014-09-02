CC=		gcc
CFLAGS=	-std=c99 -Wall

all: demo

demo: ec_c_p_mul.o
	$(CC) $(CFLAGS) ec_c_p_mul.c -o demo

.PHONY: clean
clean:
	rm -rf *o demo