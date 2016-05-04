#FWI tools

This directory contains (or did/will contain) some alternative tools 
for the FWI algorithm. You can swap out seismic modeling algorithms
or wavefield comparison algorithms to suit your needs, and there
are a couple non-fancy ones in here for you to try out if you'd like.

To change algorithms just change the the include names in the FWI
program. It should look something like this:

include("xxx.jl")