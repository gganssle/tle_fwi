# This is an image comparison via time domain
# subtraction in a sliding wavelet sized window. 
#
# author: GRAM 5/1/16

function compare(s_p, s_n, s_o, sei, i, j, k, wvln, layer, verbose)

	# split the interface sample
	ind_1 = Int(layer[k + 1] + ceil(wvln / 2))
	ind_2 = Int(layer[k + 1] - floor(wvln / 2))

	if verbose == 2
		print("\nwindow goes $ind_2 to $ind_1\n")
	end

	# find error
	dif_p = abs(sei[ind_2 : ind_1, i, j] - s_p[ind_2 : ind_1])
	dif_n = abs(sei[ind_2 : ind_1, i, j] - s_n[ind_2 : ind_1])
	dif_o = abs(sei[ind_2 : ind_1, i, j] - s_o[ind_2 : ind_1])

	# sum the difference vectors
	return [sum(dif_p),sum(dif_n),sum(dif_o)]

end
