# This is the simplest example of full waveform inversion. It is
# built in the poststack domain using a simple convolution, which 
# highlights the theory of FWI, but is useless commercially.
#
# Because this is a tutorial I've written the code for clarity, not
# speed. So if you try to use some huge velocity model, it'll take
# a ton of memory and crash. Probably.
#
# by GRAM | 23 Apr 2016

using Seismic

# read initial vel model, seismic image, and Ricker wavelet
vel,vel_h = SeisRead("../vel/vel_1")
sei, sei_h = SeisRead("../vel/vel_1")
rick, rick_h = SeisRead("../mod/wav")

# initialize
nz = size(vel)[1]
nx = size(vel)[2]
ny = size(vel)[3]
dz = vel_h[1].d1
dx = vel_h[1].d2
dy = vel_h[1].d3
inc = .10		# velocity update percentage
d1 = d2 = 1		# acoustic assumption
v_p = zeros(Float32,nz,nx,ny)
v_n = zeros(Float32,nz,nx,ny)
r_p = zeros(Float32,nz,nx,ny)
r_n = zeros(Float32,nz,nx,ny)
r_o = zeros(Float32,nz,nx,ny)
s_p = zeros(Float32,nz,nx,ny)
s_n = zeros(Float32,nz,nx,ny)
s_o = zeros(Float32,nz,nx,ny)

#### main FWI loop ########################
for i = 1 : ny, j = 1 : nx, k = 1 : nz
	
	print("\n Now updating trace x = ", j, ", y = ", i, " \n")

	# build perturbed velocity models
	v_p = vel[:,:,:]
	v_p[k,i,j] += v_p[k,i,j] * (1 + inc)
	
	v_n = vel[:,:,:]
	v_n[k,i,j] += v_n[k,i,j] * (1 - inc)

	# calculate reflection coeficients
	for kk = 1 : ny, jj = 1 : nx, ii = 1 : nz - 1

		r_p[ii,jj,kk] = ((d2 * v_p[ii+1,jj,kk] - d1 * v_p[ii,jj,kk]) /
				   (d2 * v_p[ii+1,jj,kk] + d1 * v_p[ii,jj,kk])) ^ 2
		
		r_n[ii,jj,kk] = ((d2 * v_n[ii+1,jj,kk] - d1 * v_n[ii,jj,kk]) /
				   (d2 * v_n[ii+1,jj,kk] + d1 * v_n[ii,jj,kk])) ^ 2

		r_o[ii,jj,kk] = ((d2 * vel[ii+1,jj,kk] - d1 * vel[ii,jj,kk]) /
				   (d2 * vel[ii+1,jj,kk] + d1 * vel[ii,jj,kk])) ^ 2

	end

	# model perturbed seismic data
		# increase trace length for convolved samples
	nz_new = size(rick)[1] + nz - 1

	m_p = zeros(Float32, nz_new, nx, ny)
	m_n = zeros(Float32, nz_new, nx, ny)
	m_o = zeros(Float32, nz_new, nx, ny)

		# model seismic with a convolution
	for jj = 1 : ny, ii = 1 : nx
		m_p[:,ii,jj] = conv(r_p[:,ii,jj], rick)
		m_n[:,ii,jj] = conv(r_n[:,ii,jj], rick)
		m_o[:,ii,jj] = conv(r_o[:,ii,jj], rick)
	end

	# clip off unused (non-physical) ends of convolution
	extra = size(m_p)[1] - nz

	if isodd(extra) == true
		for jj = 1 : ny, ii = 1 : nx
			s_p[:,ii,jj] = m_p[ceil(extra / 2) : nz - 1 + ceil(extra / 2), ii, jj]
			s_n[:,ii,jj] = m_n[ceil(extra / 2) : nz - 1 + ceil(extra / 2), ii, jj]
			s_o[:,ii,jj] = m_o[ceil(extra / 2) : nz - 1 + ceil(extra / 2), ii, jj]
		end
	else
		for jj = 1 : ny, ii = 1 : nx
			s_p[:,ii,jj] = m_p[extra / 2 : nz - 1 + extra / 2, ii, jj]
			s_n[:,ii,jj] = m_n[extra / 2 : nz - 1 + extra / 2, ii, jj]
			s_o[:,ii,jj] = m_o[extra / 2 : nz - 1 + extra / 2, ii, jj]
		end
	end

	# compare images
		# correlate signals
		
		# difference and weight

		# sum the differences
	comp = [,,]
	
	# move in the direction of improvement
	dir = find(comp .== max(comp[1],comp[2],comp[3]))
	
	if dir == 1 
		vel[k,j,i] = v_p[k,j,i]
	elseif dir == 2
		vel[k,j,i] = v_n[k,j,i]
	end
	
end

# write out the wavefield
ex = Seismic.Extent(convert(Int32,nz), convert(Int32,nx), convert(Int32,ny), 
	1, 1, 0, 0, 0, 0, 0, convert(Float32,dz), convert(Float32,dx), 
	convert(Float32,dy), 1, 1, "Depth", "mx", "my", "", "", "", "", "", 
	"", "", "")

SeisWrite("updated_vel",vel,vel_h,ex)







