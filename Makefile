#
# Makefile for direct.
#
CFLAGS	=   -O2
LIBS	=   -lm

default:	direct

clean:
	rm -f *.o

direct: main.o kd.o grav.o ewald.o
	$(CC) $(CFLAGS) -o direct main.o kd.o grav.o ewald.o $(LIBS)

main.o: main.c kd.h
	$(CC) $(CFLAGS) -c main.c

kd.o: kd.c kd.h grav.h ewald.h tipsydefs.h
	$(CC) $(CFLAGS) -c kd.c

grav.o: grav.c grav.h kd.h spline.h 
	$(CC) $(CFLAGS) -c grav.c

ewald.o: ewald.c ewald.h kd.h spline.h 
	$(CC) $(CFLAGS) -c ewald.c

