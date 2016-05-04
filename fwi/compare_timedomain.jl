# This is an image comparison via simple subtraction.
# All we do here is take the difference of each 
# modeled trace with the real data trace. Whichever
# has the lowest error (and I'm using the term loosely,
# will be chosen by the FWI algorithm.
#
# author: GRAM 5/1/16

function compare(s_p, s_n, s_o, sei, i, j)

	# difference and weight
	dif_p = abs(sei[:,i,j] - s_p)
	dif_n = abs(sei[:,i,j] - s_n)
	dif_o = abs(sei[:,i,j] - s_o)

#=
	print("\n",dif_p,"\n")
	print("\n",dif_n,"\n")
	print("\n",dif_o,"\n")
=#

	# sum the difference vectors
	return [sum(dif_p),sum(dif_n),sum(dif_o)]

end
