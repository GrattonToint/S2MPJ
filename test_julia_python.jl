#  Tests Julia and Python problems
#
#  Programming: Ph. Toint and S. Gratton
#  This version : 17 V 2024
#
########################################################################

include("s2xlib.jl")

start = 1
stop  = 10000

second_evaluation = true
eval_matvec       = true
fullcheck         = false

usejulia  = true
usepython = true

#########################################################################

if ( usepython )
   ENV["PYTHON"] = "/usr/bin/python3"
   using PyCall
   pushfirst!(PyVector(pyimport("sys")."path"), "/home/philippe/s2x")
end
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
#problem = [ "QUARTC" ]
#problem = [ "ROTDISC" ]
#problem = [ "SIMBQP" ]
#problem = [ "TRIDIA" ]
#problem = [ "ZAMB2m10" ]
#invalue = [ "4"]
#invalue = [ "10"]

for iprob = min(start,length(problem)):min(stop,length(problem))

    prob  = problem[iprob]
    if startswith(prob, "#" )
       continue
    end

    inval = invalue[iprob]
    println( "=== Checking problem ", iprob, ": ", prob )
    
    #  Julia
    
if( usejulia )

    jlfile = prob*".jl"
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
    
    if( usepython )
#        println( " Julia ok!")
    end
end 

    #  Python

if ( usepython )

    pyfile = prob*".py"
    run( `touch $pyfile` )
    
    inval   = eval( Meta.parse(inval) ) 
    Mod     = pyimport( prob )
    pyfunc  = Mod[prob]

    @time begin
    pyprob  = pyfunc(inval...)
    print("   Python setup      :" )
    end
  
    ypy  = ones(pyprob.pb.m,1)
    x0py = pyprob.pb.x0
    
    @time begin
    if fullcheck
       fxpy    = theprob( "fx", pbm, x0py )
       println( "fx python = ", fxpy )
       cxpy    = theprob( "cx", pbm, x0py )
       Lxypy   = pyprob.Lxy( x0py, ypy )
    else
       Lxypy,Lgxypy,LHxypy = pyprob.LgHxy( x0py, ypy )
    end
    print("   Python evaluation :" )
    end

    if second_evaluation
       @time begin
       if fullcheck
          fxpy = theprob( "fx", pbm, x0py )
          cxpy = theprob( "cx", pbm, x0py )
          Lxypy,Lgxypy,LHxypy  = pyprob.LgHxy( x0py, ypy )
       else
          Lxypy,Lgxypy,LHxypy  = pyprob.LgHxy( x0py, ypy )
       end
       print("   Python evaluation :" )
       end
    end

    if eval_matvec
       vpy =  ones( pyprob.pb.n, 1 )
       @time begin
       Hxyvpy = pyprob.LHxyv( x0py, ypy, vpy )
       print("   Python prod H*v   :" )
       end
    end

    if fullcheck
       fxypy,fgxypy,fHxypy  = pyprob.fgHx( x0py )#D
       println( "true Hv = ", fHxypy*vpy )#D
       cxypy,cJpy, cHpy  = pyprob.cJHx( x0py )#D
       yHiv = zeros( pyprob.pb.n,1)
       for ig in 1:pyprob.pb.m
          yHiv += ypy[ig] * cHpy[ig]* vpy
       end
       println(" true yHiv = ", yHiv )
       println("true LHv = ", fHxypy*vpy+yHiv)
    
       println( "Lxy Python = ", Lxypy )
       println( "fx python = ", fxpy )
       println( "y'*c = ", only( ypy'*cxpy ) )
       Lbpy = fxpy + only(ypy'*cxpy)
       println( "Lbpy = ", Lbpy )
#       println( "g Lag Python = ", Lgxypy )
       println( "Verif matvec Python = ", norm( LHxypy * vpy - Hxyvpy ) )
    end

    if ( usejulia )
#        println( " Python ok!" )
    end
end

    #  Print error

if ( usejulia && usepython )
   if fullcheck
       println( " Error on the start point: absolute = ", norm(x0jl-x0py), ", relative = ", norm(x0jl-x0py)/(1e-15+norm(x0py) ) )
       println( " Error on the objective  : absolute = ", norm(fxjl-fxpy), ", relative = ", norm(fxjl-fxpy)/(1e-15+norm(fxpy) ) )
       println( " Error on the constraint : absolute = ", norm(cxjl-cxpy), ", relative = ", norm(cxjl-cxpy)/(1e-15+norm(cxpy) ) )
       println( " Error on the multiplier : absolute = ", norm(yjl -ypy ), ", relative = ", norm(yjl-ypy)/(1e-15+norm(ypy) ) )
       println( " Error on the Lagrangian : absolute = ", norm(Lxyjl-Lxypy), ", relative = ", norm(Lxyjl-Lxypy)/(1e-15+norm(Lxypy) ) )
       println( " Error on the Lag grad   : absolute = ", norm(Lgxyjl-Lgxypy), ", relative = ", norm(Lgxyjl-Lgxypy)/(1e-15+norm(Lgxypy) ) )
       if eval_matvec
          println( " Error on the product H*v: absolute = ", norm(Hxyvjl-Hxyvpy), ", relative = ", norm(Hxyvjl-Hxyvpy)/(1e-15+norm(Hxyvpy) ) )
       end
    else
       println( " Error on the Lag gradient: absolute = ", norm(Lgxyjl-Lgxypy), ", relative = ", norm(Lgxyjl-Lgxypy)/(1e-15+norm(Lgxypy) ) )
       if eval_matvec
          println( " Error on the product H*v : absolute = ", norm(Hxyvjl-Hxyvpy), ", relative = ", norm(Hxyvjl-Hxyvpy)/(1e-15+norm(Hxyvpy) ) )
       end
    end
end

end
