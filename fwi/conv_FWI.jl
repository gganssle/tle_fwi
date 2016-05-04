# This is the simplest example of full waveform inversion. It is
# built in the 3D poststack domain using a simple convolution, which 
# highlights the theory of FWI, but is inaccurate on real data.
#
# by GRAM | 23 Apr 2016

print("\n3D Poststack full waveform inversion via convolution\n\n
In this program FWI updates a velocity model one sample at a time,
one trace at a time. The velocity sample is perturbed positively,
negatively, and not at all. Then the algorithm models new seismic
from all three of the perturbed velocity models, and compares them to
the input (real) seismic data. Whichever perturbation generates the
least amount of error is chosen to be correct, and input into the 
updated velocity model.
This algorithm accepts six inputs in the following order:\n
	1. Initial velocity model
	2. Poststack 3D seismic volume in Seis format
	3. One dimensional wavelet in Seis format
	4. Output file name (updated vel model in Seis format)
	5. Velocity update increment percentage in decimal
	6. Maximum number of velocity update iterations\n
Here is an example to get you started on the syntax:\n\n
fwi(\"../dat/vel_incorrect\", \"../dat/image\", \"../dat/wav\", \"../dat/updated_vel\", .1, 5)\n\n
Please note this algorithm is written for clarity not speed, so
use forgivingly small velocity models.\n\n")

using Seismic

function fwi(vel="",sei="",wav="",out="",inc=.1,itr_max=5)

	# read initial vel model, seismic image, and wavelet
	vel,vel_h = SeisRead(vel)
	sei, sei_h = SeisRead(sei)
	wav, wav_h = SeisRead(wav)

	# initialize
	nz = size(vel)[1]
	nx = size(vel)[2]
	ny = size(vel)[3]
	dz = vel_h[1].d1
	dx = 20
	dy = 20
#	itr_max = 5		# max number of vel iterations
#	inc = .1		# velocity update percentage
	d1 = d2 = 1		# acoustic assumption
	r_p = zeros(Float32,nz)
	r_n = zeros(Float32,nz)
	r_o = zeros(Float32,nz)

	#### main FWI loop ########################
	for i = 1 : ny, j = 1 : nx, k = 1 : nz
	
		print("\nNow updating trace x = ", j, ", y = ", i, " \n")

		for itr = 1 : itr_max
	
			# build perturbed velocity models
			v_p = vel[:,j,i]
			v_p[k] += v_p[k] * (1 + inc)
		
			v_n = vel[:,j,i]
			v_n[k] += v_n[k] * (1 - inc)
		
			v_o = vel[:,j,i]

			# calculate reflection coeficients
			for kk = 1 : nz - 1

				r_p[kk] = ((d2 * v_p[kk + 1] - d1 * v_p[kk]) /
						   (d2 * v_p[kk + 1] + d1 * v_p[kk])) ^ 2
			
				r_n[kk] = ((d2 * v_n[kk + 1] - d1 * v_n[kk]) /
						   (d2 * v_n[kk + 1] + d1 * v_n[kk])) ^ 2

				r_o[kk] = ((d2 * v_o[kk + 1] - d1 * v_o[kk]) /
						   (d2 * v_o[kk + 1] + d1 * v_o[kk])) ^ 2

			end

			# model perturbed seismic data
				# increase trace length for convolved samples
			nz_new = size(wav)[1] + nz - 1

			m_p = zeros(Float32, nz_new, nx, ny)
			m_n = zeros(Float32, nz_new, nx, ny)
			m_o = zeros(Float32, nz_new, nx, ny)

				# model seismic with a convolution
			m_p = conv(r_p, wav)
			m_n = conv(r_n, wav)
			m_o = conv(r_o, wav)

				# clip off unused (non-physical) ends of convolution
			extra = size(m_p)[1] - nz

			if isodd(extra) == true
				s_p = m_p[ceil(extra / 2) : nz - 1 + ceil(extra / 2)]
				s_n = m_n[ceil(extra / 2) : nz - 1 + ceil(extra / 2)]
				s_o = m_o[ceil(extra / 2) : nz - 1 + ceil(extra / 2)]
			else
				s_p = m_p[div(extra, 2) : nz - 1 + div(extra, 2)]
				s_n = m_n[div(extra, 2) : nz - 1 + div(extra, 2)]
				s_o = m_o[div(extra, 2) : nz - 1 + div(extra, 2)]
			end

			# compare images
				# correlate signals
			cor_p = ifft(conj(fft(s_p)) .* fft(sei[:,j,i]))
			cor_n = ifft(conj(fft(s_n)) .* fft(sei[:,j,i]))
			cor_o = ifft(conj(fft(s_o)) .* fft(sei[:,j,i]))
			acor = ifft(conj(fft(sei[:,j,i])) .* fft(sei[:,j,i]))

				# convolution builds in neglible reals 
				# instead of zeros, so elim trailing imaginaries
			cor_p = real(cor_p)
			cor_n = real(cor_n)
			cor_o = real(cor_o)
			acor = real(acor)

				# normalize amplitudes
			cor_p = cor_p ./ maximum(cor_p)
			cor_n = cor_n ./ maximum(cor_n)
			cor_o = cor_o ./ maximum(cor_o)
			acor = acor ./ maximum(acor)

				# difference and weight
			dif_p = abs(acor - cor_p)
			dif_n = abs(acor - cor_n)
			dif_o = abs(acor - cor_o)

				# sum the difference vectors
			comp = [sum(dif_p),sum(dif_n),sum(dif_o)]
		
			print("\n",comp,"\n")

			# move in the direction of improvement
			dir = 0
			dir = find(comp .== minimum(comp))[1]
		
			print("\n",dir,"\n")

			if dir == 1
				vel[k,j,i] = v_p[k]
			elseif dir == 2
				vel[k,j,i] = v_n[k]
			elseif dir == 3
				break
			else
				print("\n xcor diff symmetric @", k,"\n")
			end
	
		end
	
	end

	# write out the wavefield
	ex = Seismic.Extent(convert(Int32,nz), convert(Int32,nx), convert(Int32,ny), 
		1, 1, 0, 0, 0, 0, 0, convert(Float32,dz), convert(Float32,dx), 
		convert(Float32,dy), 1, 1, "Depth", "mx", "my", "", "", "", "", "", 
		"", "", "")

	SeisWrite("../dat/updated_vel",vel,vel_h,ex)

end





