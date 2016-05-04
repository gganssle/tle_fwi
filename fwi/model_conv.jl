# This is a simple seismic modeling algorithm in the
# 3D poststack domain. It just convolves a wavelet
# with normal incident reflection coefficients
# generated from a velocity model.
# 
# author: GRAM 5/1/16

function model(v_p, v_n, v_o, nx, ny, nz, wav)

	# initialize
	r_p = zeros(Float32,nz)
	r_n = zeros(Float32,nz)
	r_o = zeros(Float32,nz)
	d1 = d2 = 1		# acoustic assumption

	# calculate reflection coeficients
	for kk = 1 : nz - 1

		r_p[kk] = ((d2 * v_p[kk + 1] - d1 * v_p[kk]) /
				   (d2 * v_p[kk + 1] + d1 * v_p[kk])) ^ 2

		r_n[kk] = ((d2 * v_n[kk + 1] - d1 * v_n[kk]) /
				   (d2 * v_n[kk + 1] + d1 * v_n[kk])) ^ 2

		r_o[kk] = ((d2 * v_o[kk + 1] - d1 * v_o[kk]) /
				   (d2 * v_o[kk + 1] + d1 * v_o[kk])) ^ 2

	end

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

	return s_p, s_n, s_o

end
