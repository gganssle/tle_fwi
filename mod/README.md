#Modeling tools

This directory contains a simple seismic modeling algorithm. Really
this is just a little program which convolves a wavelet, in this
case a Ricker wavelet, with a series of normally incident reflection
coefficients calculated from a 3D velocity model.

The output is a 3D "poststack" dataset. Though that's really a
misnomer, because there was no prestack domain to start with.
