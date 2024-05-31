module PROB

global pbm
include("s2xlib.jl")

function setup(n)
#   println( "In setup" )
   pbm = n
end

function eval()
   println( pbm )
end

end
