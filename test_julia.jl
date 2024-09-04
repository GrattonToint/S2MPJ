#  Tests Julia problems
#
#  Programming: Ph. Toint and S. Gratton
#  This version: 2 IX 2024
#
########################################################################

include("s2mpjlib.jl")

start = 1071
stop  = 10000


second_evaluation = true
eval_matvec       = true
fullcheck         = false
printtimes        = true  # print timings in test_julia.times
printres          = true  # print values of x0, L0, g0, and Hv in test_julia.data

#########################################################################

using Printf
using LinearAlgebra
using TimerOutputs

problem = String[]
invalue = String[]

filename = "fullproblist"
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
timer    = TimerOutput()

for iprob = min(start,length(problem)):min(stop,length(problem))

    prob  = problem[iprob]
    if startswith(prob, "#" )
       continue
    end
    
    reset_timer!(timer::TimerOutput)
    
    inval = invalue[iprob]
    println( "=== Checking problem ", iprob, ": ", prob )
    if printtimes
        fidtimes = open("test_julia.times","a")
        @printf( fidtimes, "=== Checking problem %d: %s\n", iprob, prob )
    end
    jlfile = "./julia_problems/"*prob*".jl"
    run( `touch $jlfile` )
    
    include(jlfile)
    theprob = eval( Meta.parse( prob ) )
    theargs = "\"setup\","*@sprintf("%s",inval)

    @timeit timer "setup"  pb,pbm  = eval( Meta.parse( prob*"("*theargs*")" ) )
    thetime =  10^(-9) * TimerOutputs.time(timer["setup"])
    @printf(stdout, "   Julia setup       :  %f8 seconds\n", thetime)
    if printtimes
        @printf( fidtimes, "   Julia setup       :  %f8 seconds\n", thetime )
    end

    if printres
        fidres = open("test_julia.data","a")
        @printf( fidres, "=== Checking problem  %d :  %s ( n = %d  m = %d )\n", iprob, prob, pb.n, pb.m )
        j = 0
        while j < pb.n
             j += 1
             @printf( fidres, "x0 = %+18.15e",  pb.x0[j]  )
             for k in [ 1 2 3 4 ]
                 j += 1
                 if j > pb.n
                    break
                 end
                 @printf( fidres, " %+18.15e",  pb.x0[j] )
             end
             println( fidres, "" )
        end
    end

    if fullcheck
       dump(pb)
       dump(pbm)
    end

    yjl  = ones( pb.m, 1 )
    x0jl = pb.x0
    ifixed = findall(x->x==0,pb.xupper-pb.xlower)
    if length(ifixed) > 0
       x0jl[ifixed] = pb.xlower[ifixed]
    end

    @timeit timer "eval1" begin
    if fullcheck
       fxjl  = theprob( "fx", pbm, x0jl )
       println( "fx Julia = ", fxjl )#D
       cxjl  = theprob( "cx", pbm, x0jl )
       Lxyjl = theprob( "Lxy", pbm, x0jl, yjl )
    else
       Lxyjl,Lgxyjl,LHxyjl = theprob( "LgHxy", pbm, x0jl, yjl )
    end
    end
    thetime =  10^(-9) * TimerOutputs.time(timer["eval1"])
    @printf( stdout, "   Julia evaluation  :  %f8 seconds\n", thetime )
    if printtimes
        @printf( fidtimes, "   Julia evaluation  :  %f8 seconds\n", thetime )
    end
    if printres
        @printf( fidres, "L0 = %+18.15e\n", Lxyjl  )
        j = 0
        while j < pb.n
             j += 1
             @printf( fidres, "g0 = %+18.15e", Lgxyjl[j]  )
             for k in [ 1 2 3 4 ]
                 j += 1
                 if j > pb.n
                    break
                 end
                 @printf( fidres, " %+18.15e",  Lgxyjl[j] )
             end
             println( fidres, "" )
        end
    end
    

    if second_evaluation
       @timeit timer "eval2"  begin
       if fullcheck
          fxjl  = theprob( "fx", pbm, x0jl )
          cxjl  = theprob( "cx", pbm, x0jl )
          Lxyjl,Lgxyjl,LHxyjl = theprob( "LgHxy", pbm, x0jl, yjl )
       else
          Lxyjl,Lgxyjl,LHxyjl = theprob( "LgHxy", pbm, x0jl, yjl )
       end
       end
       thetime =  10^(-9) * TimerOutputs.time(timer["eval2"])
       @printf( stdout, "   Julia evaluation  :  %f8 seconds\n", thetime )
       if printtimes
           @printf( fidtimes, "   Julia evaluation  :  %f8 seconds\n", thetime )
       end
    end

    if eval_matvec
       vjl =  ones( pb.n, 1 )
       if  length(ifixed) > 0
          vjl[ifixed] = zeros( length(ifixed), 1 )
       end
       @timeit timer "prod"  Hxyvjl = theprob("LHxyv", pbm, x0jl, yjl, vjl )
       thetime =  10^(-9) * TimerOutputs.time(timer["prod"])
       @printf( stdout, "   Julia prod H*v    :  %f8 seconds\n", thetime )
       if printtimes
           @printf( fidtimes, "   Julia prod H*v    :  %f8 seconds\n", thetime )
           close( fidtimes )
       end
       if printres
           j = 0
           while j < pb.n
               j += 1
               @printf( fidres, "Hv = %+18.15e", Hxyvjl[j]  )
               for k in [ 1 2 3 4 ]
                   j += 1
                   if j > pb.n
                      break
                   end
                   @printf( fidres, " %+18.15e",  Hxyvjl[j] )
               end
               println( fidres, "" )
           end
           close( fidres )
       end
    end

    if fullcheck
       println( "Lxy Julia = ", Lxyjl )
       println( "fx Julia = ", fxjl )
       println( "y'*c = ", only( yjl'*cxjl ) )
       Lbjl = fxjl + only(yjl'*cxjl)
       println( "Lbjl = ", Lbjl )
       println( " Verif matvec Julia = ", norm( LHxyjl * vjl - Hxyvjl ) )
    end
    
end

