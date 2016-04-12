# This is the simplest example of full waveform inversion.
#
# by GRAM | 12 Apr 2016

using Seismic

# import vel model
vel,vel_h = SeisRead("../vel/vel_1")
nz = size(vel)[1]
nx = size(vel)[2]
ny = size(vel)[3]

# initialize


# write out the wavefield
h = Array(Header,nx*ny)
ex = Seismic.Extent(convert(Int32,nz), convert(Int32,nx), convert(Int32,ny), 
	1, 1, 0, 0, 0, 0, 0, convert(Float32,dz), convert(Float32,dx), 
	convert(Float32,dy), 1, 1, "Depth", "mx", "my", "", "", "", "", "", 
	"", "", "")

#SeisWrite("image",vel,vel_h,ex)







