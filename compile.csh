gcc -I/usr/include/python2.7/ -c model.c -o model.o
gcc -I/usr/include/python2.7/ -c fitsio.c -o fitsio.o
gcc -g -I./ -I/usr/include/python2.7/  -L/DATA/disks/rjharris/fftw3/lib64/ -L/home/rjharris/lib64/ fitsio.o  driver_model.c model.o  -o driver -lm -lfftw3  -ldl -lpython2.7
