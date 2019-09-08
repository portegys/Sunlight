#
# makefile for sunlight program
#

all: sunlight

sunlight: sunlight.o
	gcc  -o sunlight -L/usr/openwin/lib \
            sunlight.o -lXol -lXt -lX11 -lsunmath -lm

sunlight.o: sunlight.C
	gcc -I/usr/openwin/include -c sunlight.C
