# This program generates a 3D synthetic seismic section by convolution with
# a the Ricker wavelet defined in Seismic.jl's ShotProfileWEM. To do this
# it reads a velocity model in, calculates reflection coeficients, and uses
# a 1D convolution per trace in the frequency domain.
#
# function call vars:
# peakF is the Ricker peak frequency, and samp is the sampling interval
#
# by GRAM | 8 Apr 2016

using Seismic

function modeler(vel="", out="", peakF=100, samp=0.001, verbose=0)

	# import vel model
	vel,vel_h = SeisRead(vel)
	nz = size(vel)[1]
	nx = size(vel)[2]
	ny = size(vel)[3]

	# initialize
	dz = 10
	dx = 20
	dy = 20

	# build a wavelet for modeling
	temp = ricker(peakF, samp)
	rick = zeros(Float32, size(temp)[1])
	for i = 1 : size(temp)[1]
		rick[i] = convert(Float32, temp[i])
	end
	
	if verbose > 0
		print("\nThis Ricker wavelet is ",size(rick)[1]," samples long.\n\n")
	end

	#write out wavelet for use with finite difference modeling
	ex = Seismic.Extent(convert(Int32,size(rick)[1]), convert(Int32,1), convert(Int32,1), 
		1, 1, 0, 0, 0, 0, 0, convert(Float32,dz), convert(Float32,1), 
		convert(Float32,1), 1, 1, "Depth", "mx", "0", "", "", "", "", "", 
		"", "", "")
	h_w = Array(Header,1);
	h_w[1] = Seismic.InitSeisHeader();
	h_w[1].tracenum = 1;
	h_w[1].n1 = size(rick)[1];
	h_w[1].d1 = samp;
	SeisWrite("../dat/wav", rick[:], h_w, ex);

	# calculate (normal incidence) reflection coeficients
	r = zeros(Float32,nz,nx,ny)
	d1 = d2 = 1 # acoustic assumption

	for k = 1 : ny
		for j = 1 : nx
			for i = 1 : nz - 1
				r[i,j,k] = ((d2 * vel[i+1,j,k] - d1 * vel[i,j,k]) /
					   (d2 * vel[i+1,j,k] + d1 * vel[i,j,k])) ^ 2
			end
		end
	end

#=
	# write out reflection coefficients
	h = Array(Header,nx*ny)
	ex = Seismic.Extent(convert(Int32,nz), convert(Int32,nx), convert(Int32,ny), 
		1, 1, 0, 0, 0, 0, 0, convert(Float32,dz), convert(Float32,dx), 
		convert(Float32,dy), 1, 1, "Depth", "mx", "my", "", "", "", "", "", 
		"", "", "")

	for i = 1 : ny
		for j = 1 : nx
			h[j + nx * (i - 1)] = Seismic.InitSeisHeader()
			h[j + nx * (i - 1)].tracenum = convert(Int32, j + nx * (i - 1))
			h[j + nx * (i - 1)].o1 = convert(Float32, 0)		
			h[j + nx * (i - 1)].n1 = convert(Int32, nz)
			h[j + nx * (i - 1)].d1 = convert(Float32, dz)
			h[j + nx * (i - 1)].mx = convert(Float32, j)
			h[j + nx * (i - 1)].my = convert(Float32, i)
		end
	end

	SeisWrite("../dat/refl",r,vel_h,ex)
=#

	# convolve RC with Ricker wavelet to generate 0 offset section
	nz_new = size(rick)[1] + nz - 1
	m = zeros(Float32, nz_new, nx, ny)
	for j = 1 : ny
		for i = 1 : nx
			m[:,i,j] = conv(r[:,i,j], rick)
		end
	end

	# window out middle samples of convolution for physical reality
	extra = size(m)[1] - nz
	if isodd(extra) == true
		for j = 1 : ny
			for i = 1 : nx
				vel[:,i,j] = m[Int(ceil(extra / 2)) : nz - 1 + Int(ceil(extra / 2)), i, j]
			end
		end
	else
		for j = 1 : ny
			for i = 1 : nx
				vel[:,i,j] = m[Int(extra / 2) : nz - 1 + Int(extra / 2), i, j]
			end
		end
	end

	# write out the wavefield
	h = Array(Header,nx*ny)
	ex = Seismic.Extent(convert(Int32,nz), convert(Int32,nx), convert(Int32,ny), 
		1, 1, 0, 0, 0, 0, 0, convert(Float32,dz), convert(Float32,dx), 
		convert(Float32,dy), 1, 1, "Depth", "mx", "my", "", "", "", "", "", 
		"", "", "")

	SeisWrite(out,vel,vel_h,ex)

end # end function





