# This is an image comparison via frequency domain
# correlation. The image correlations are differenced
# and weighted against the original signal's 
# autocorrelation.
#
# author: GRAM 5/1/16

function compare(s_p, s_n, s_o, sei, i, j)

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
	return [sum(dif_p),sum(dif_n),sum(dif_o)]

end
