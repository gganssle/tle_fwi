# This program generates a 3D synthetic seismic section by convolution with
# a the Ricker wavelet defined in Seismic.jl's ShotProfileWEM. To do this
# it reads a velocity model in, calculates reflection coeficients, and uses
# a 1D convolution per trace in the frequency domain.
#
# by GRAM | 8 Apr 2016

