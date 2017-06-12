all:
	make linux

linux:
	mpif77  -I/usr/include -O3 -w -o pd2 pd2.f
	rm -f *.o

clean:
	rm -f hostfile *.data *.o *.dat pd2 fort.* *.err *.out *.xyz *.png
