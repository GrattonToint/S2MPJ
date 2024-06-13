#  Tests Julia problems
#
#  Programming: Ph. Toint and S. Gratton
#  This version : 10 VI 2024
#
########################################################################

include("s2mpjlib.jl")

start = 1
stop  = 10000

second_evaluation = true
eval_matvec       = true
fullcheck         = false

usejulia  = true

#########################################################################

using Printf
using LinearAlgebra

filename = "fullproblist"

problem = String[]
invalue = String[]

open(filename, "r") do file

    # Read all lines

    lines = readlines(file)
    for line in lines
        parts = split(line, r"\s+")
        if length(parts) >= 2
            push!(problem, parts[1])
            push!(invalue, parts[2])
        end
    end
end

#problem = [ "ALLINITC" ]
#invalue = [ "4"]

for iprob = min(start,length(problem)):min(stop,length(problem))

    prob  = problem[iprob]
    if startswith(prob, "#" )
       continue
    end

    inval = invalue[iprob]
    println( "=== Checking problem ", iprob, ": ", prob )
    
    #  Julia
    

    jlfile = "./julia_problems/"*prob*".jl"
    run( `touch $jlfile` )
    
    include(jlfile)
    theprob = eval( Meta.parse( prob ) )
    theargs = "\"setup\","*@sprintf("%s",inval)

    @time begin
    pb,pbm  = eval( Meta.parse( prob*"("*theargs*")" ) )
    print("   Julia setup       :")
    end

    if fullcheck
       dump(pb)
       dump(pbm)
    end

    yjl  = ones( pb.m, 1 )
    x0jl = pb.x0

    @time begin
    if fullcheck
       fxjl  = theprob( "fx", pbm, x0jl )
       println( "fx Julia = ", fxjl )#D
       cxjl  = theprob( "cx", pbm, x0jl )
       Lxyjl = theprob( "Lxy", pbm, x0jl, yjl )
    else
       Lxyjl,Lgxyjl,LHxyjl = theprob( "LgHxy", pbm, x0jl, yjl )
    end
    print("   Julia evaluation  :" )
    end

    if second_evaluation
       @time begin
       if fullcheck
          fxjl  = theprob( "fx", pbm, x0jl )
          cxjl  = theprob( "cx", pbm, x0jl )
          Lxyjl,Lgxyjl,LHxyjl = theprob( "LgHxy", pbm, x0jl, yjl )
       else
          Lxyjl,Lgxyjl,LHxyjl = theprob( "LgHxy", pbm, x0jl, yjl )
       end
       print("   Julia evaluation  :" )
       end
    end

    if eval_matvec
       vjl =  ones( pb.n, 1 )
       @time begin
       Hxyvjl = theprob("LHxyv", pbm, x0jl, yjl, vjl )
       print("   Julia prod H*v    :" )
       end
    end

    if fullcheck
       println( "Lxy Julia = ", Lxyjl )
       println( "fx Julia = ", fxjl )
       println( "y'*c = ", only( yjl'*cxjl ) )
       Lbjl = fxjl + only(yjl'*cxjl)
       println( "Lbjl = ", Lbjl )
#       println( "g Lag Julia = ", Lgxyjl )
       println( " Verif matvec Julia = ", norm( LHxyjl * vjl - Hxyvjl ) )
    end
    
end

