# This program creates a simple 3D vel model for FD seismic modeling
#
# by GRAM | 8 Apr 2016

# def size
nx = 3
ny = 2
nz = 10

vel = zeros(Float32,nz,nx,ny)

#build model
for i = 1:nz, j = 1:nx, k = 1:ny
	if i < nz / 2
		vel[i,j,k] = 1
	else
		vel[i,j,k] = 2
	end
end

#= display
print(size(vel), "\n")
print(vel,"\n")
=#

#write out text file for QC
out = open("vel_1.txt","w")
for k = 1:ny
	for i = 1:nz
		for j = 1:nx
			showcompact(out, vel[i,j,k])
			write(out," ")
		end
		write(out, "\r\n")
	end
	write(out, "\r\n")
end
close(out)

#write out model (data) in binary for Seismic.jl
out = open("vel_1.bin","w")
for k = 1:ny, i = 1:nz, j = 1:nx
			write(out, vel[i,j,k])
end
close(out)

#write out model (headers) in binary for Seismic.jl

#write out model (text info) in binary for Seismic.jl

