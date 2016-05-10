# This is the simplest example of full waveform inversion. It is
# built in the 3D poststack domain using a simple convolution, which 
# highlights the theory of FWI, but is inaccurate on real data.
#
# by GRAM | 23 Apr 2016

print("\n3D Poststack full waveform inversion algorithm\n\n
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
	6. Maximum number of velocity update iterations
	7. Verbosity (see note)\n
If verbose is 0 operation is silent. If verbose is 1 updates will
print to stdout. If verbose is 2 debugging info will print to stdout.
Here is an example to get you started on the syntax:\n\n
fwi(\"../dat/vel_incorrect\", \"../dat/image_correct\", \"../dat/wav\", \"../dat/updated_vel\", .1, 5, 0)\n\n
Please note this algorithm is written for clarity not speed, so
use forgivingly small velocity models.\n\n")

using Seismic
include("compare_wavelet.jl")	# wavefield comparison algorithm
include("model_conv.jl")	# seismic modeling algorithm

function fwi(vel="",sei="",wav="",out="",inc=.1,itr_max=5,verbose=0)

	# read initial vel model, seismic image, and wavelet
	vel, vel_h, ex = SeisRead(vel)
	sei, sei_h = SeisRead(sei)
	wav, wav_h = SeisRead(wav)

	# initialize
	nz = size(vel)[1]
	nx = size(vel)[2]
	ny = size(vel)[3]
	dz = vel_h[1].d1
	dx = 20
	dy = 20
	thrs = 0.01		# change in vel threshold
	wvln = size(wav)[1]	# wavelength of the wavelet

	#### start FWI ##########################
	for i = 1 : ny, j = 1 : nx

		# tell the user where we are in the volume
		if verbose > 0
			print("\nNow updating trace x = ", j, ", y = ", i, " \n")
		end

		# calc velocity intervals
		layer = [0]
		for kk = 2 : nz
			current = vel[kk - 1, j, i]
			if vel[kk, j, i] > current + thrs
				append!(layer, [kk - 1])
			end
		end

		if verbose == 2 # debugging information
			print("\nlayer = ", layer,"\n")
		end

		# work over all vel intervals
		for k = 1 : size(layer)[1] - 1
	
			# iteratively update the velocity model
			for itr = 1 : itr_max

				# build perturbed velocity models
				v_p = vel[:,j,i]
				v_p[layer[k] + 1 : layer[k + 1]] = v_p[layer[k] + 1] * (1 + inc)

				v_n = vel[:,j,i]
				v_n[layer[k] + 1 : layer[k + 1]] = v_n[layer[k] + 1] * (1 - inc)

				v_o = vel[:,j,i]

				# model perturbed seismic data
				s_p, s_n, s_o = model(v_p, v_n, v_o, nx, ny, nz, wav)

				# compare images
				comp = compare(s_p, s_n, s_o, sei, i, j, k, wvln, layer, verbose)

				# move in the direction of improvement
				dir = 0
				dir = find(comp .== minimum(comp))[1]

				if verbose == 2 # debugging information
					print("\n",comp,"\n")
					print("\n",dir,"\n")
				end

				if dir == 1
					vel[:,j,i] = v_p
				elseif dir == 2
					vel[:,j,i] = v_n
				elseif dir == 3
					break
				else
					print("\n diff symmetric @", k,"\n")
				end

			end # end of vel iterations

		end # end of vel intervals

	end #### end of FWI ###########################

	# write out the FWI updated velocity model
	SeisWrite(out,vel,vel_h,ex)

end # end of function





