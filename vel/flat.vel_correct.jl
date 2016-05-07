# This program creates a simple 3D vel model for FD seismic modeling
#
# by GRAM | 8 Apr 2016

using Seismic

# def size
nx = 50
ny = 50
nz = 200
dz = 10 # depth increment
dx = 20
dy = 20

vel = zeros(Float32,nz,nx,ny)

# build model
for i = 1:nz, j = 1:nx, k = 1:ny
	if i < nz / 3
		vel[i,j,k] = 2000
	elseif ((i < 2 * nz / 3) & (i > nz / 3))
		vel[i,j,k] = 3000
	else
		vel[i,j,k] = 4000
	end
end

#= display
print(size(vel), "\n")
print(vel,"\n")
=#

#= write out text file for QC
out = open("../dat/vel_1.txt","w")
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
=#

# write out vel model for Seismic.jl
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

SeisWrite("../dat/vel_correct",vel,h,ex)







