# Use GNU Make to process this file
CC = gcc
CFLAGS=-O2 -Wall -fPIC -std=c99

files = adaptloop.o

all: libamtoolbox.a

libamtoolbox.a: $(files) 
	ar rvu libamtoolbox.a $(files_unix)
	ranlib libamtoolbox.a

%.o: %.c Makefile
	$(CC) $(CFLAGS) -c $<

clean:
	rm *.o *.a 
