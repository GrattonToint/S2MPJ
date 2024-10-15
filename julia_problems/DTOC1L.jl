function DTOC1L(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DTOC1L
#    *********
# 
#    This is a discrete time optimal control (DTOC) problem.  
#    The system has N time periods, NX control variables and NY state variables.
#    The parameter mu in the original problem formulation is set to zero, 
#    yielding linear transition functions, hence the L in the problem's name.
# 
#    The problem is convex.
# 
#    Sources: problem 1 (with mu = 0) in
#    T.F. Coleman and A. Liao,
#    "An Efficient Trust Region Method for Unconstrained Discret-Time Optimal
#    Control Problems",
#    Tech. Report, ctc93tr144,  Advanced Computing Research Institute, 
#    Cornell University, 1992.
# 
#    L.Z. Liao and C.A. Shoemaker,
#    "Advantages of differential dynamic programming over Newton's method for
#    discrete-time optimal control problems",
#    Tech. Report ctc92tr97, Advanced Computing Research Institute, 
#    Cornell University, 1992.
# 
#    SIF input: Ph. Toint, August 1993
# 
#    classification = "C-OLR2-AN-V-V"
# 
#    Problem variants: they are identified by the values of
#    the parameter vector ( N, NX, NY )
# 
#    The problem has (N-1)*NX+N*NY  variables (of which NY are fixed),
#    and (N-1)*NY constraints
# 
#       Alternative values for the SIF file parameters:
# IE N                   10             $-PARAMETER # periods  } original value
# IE NX                  2              $-PARAMETER # controls } n=   58, m=  36
# IE NY                  4              $-PARAMETER # states   }
# 
# IE N                   50             $-PARAMETER # periods  }
# IE NX                  2              $-PARAMETER # controls } n=  298, m= 196
# IE NY                  4              $-PARAMETER # states   }
# 
# IE N                   100            $-PARAMETER # periods  }
# IE NX                  2              $-PARAMETER # controls } n=  598, m= 396
# IE NY                  4              $-PARAMETER # states   }
# 
# IE N                   500            $-PARAMETER # periods  }
# IE NX                  2              $-PARAMETER # controls } n= 2998, m=1996
# IE NY                  4              $-PARAMETER # states   }
# 
# IE N                   1000           $-PARAMETER # periods  }
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "DTOC1L"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        if nargin<1
            v_["N"] = Int64(10);  #  SIF file default value
        else
            v_["N"] = Int64(args[1]);
        end
# IE NX                  2              $-PARAMETER # controls } n= 5998, m=3996
        if nargin<2
            v_["NX"] = Int64(2);  #  SIF file default value
        else
            v_["NX"] = Int64(args[2]);
        end
# IE NY                  4              $-PARAMETER # states   }
        if nargin<3
            v_["NY"] = Int64(4);  #  SIF file default value
        else
            v_["NY"] = Int64(args[3]);
        end
# IE N                   10             $-PARAMETER # periods  }
# IE NX                  5              $-PARAMETER # controls } n=  145, m=  90
# IE NY                  10             $-PARAMETER # states   }
# IE N                   50             $-PARAMETER # periods  }
# IE NX                  5              $-PARAMETER # controls } n=  745, m= 490
# IE NY                  10             $-PARAMETER # states   }
# IE N                   100            $-PARAMETER # periods  }
# IE NX                  5              $-PARAMETER # controls } n= 1495, m= 990
# IE NY                  10             $-PARAMETER # states   }
# IE N                   500            $-PARAMETER # periods  }
# IE NX                  5              $-PARAMETER # controls } n= 7495, m=4990
# IE NY                  10             $-PARAMETER # states   }
# IE N                   1000           $-PARAMETER # periods  }
# IE NX                  5              $-PARAMETER # controls } n=14995, m=9990
# IE NY                  10             $-PARAMETER # states   }
        v_["N-1"] = -1+v_["N"]
        v_["1"] = 1
        v_["2"] = 2
        v_["NY-1"] = -1+v_["NY"]
        v_["NX+NY"] = v_["NX"]+v_["NY"]
        v_["RXY"] = Float64(v_["NX+NY"])
        v_["1/RXY"] = 1.0/v_["RXY"]
        for J = Int64(v_["1"]):Int64(v_["NX"])
            for I = Int64(v_["1"]):Int64(v_["NY"])
                v_["I-J"] = I-J
                v_["RI-J"] = Float64(v_["I-J"])
                v_["B"*string(I)*","*string(J)] = v_["RI-J"]*v_["1/RXY"]
            end
        end
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for T = Int64(v_["1"]):Int64(v_["N-1"])
            for I = Int64(v_["1"]):Int64(v_["NX"])
                iv,ix_,_ = s2mpj_ii("X"*string(T)*","*string(I),ix_)
                arrset(pb.xnames,iv,"X"*string(T)*","*string(I))
            end
        end
        for T = Int64(v_["1"]):Int64(v_["N"])
            for I = Int64(v_["1"]):Int64(v_["NY"])
                iv,ix_,_ = s2mpj_ii("Y"*string(T)*","*string(I),ix_)
                arrset(pb.xnames,iv,"Y"*string(T)*","*string(I))
            end
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for T = Int64(v_["1"]):Int64(v_["N-1"])
            for I = Int64(v_["1"]):Int64(v_["NX"])
                ig,ig_,_ = s2mpj_ii("OX"*string(T)*","*string(I),ig_)
                arrset(gtype,ig,"<>")
                iv = ix_["X"*string(T)*","*string(I)]
                pbm.A[ig,iv] += Float64(1.0)
            end
        end
        for T = Int64(v_["1"]):Int64(v_["N"])
            for I = Int64(v_["1"]):Int64(v_["NY"])
                ig,ig_,_ = s2mpj_ii("OY"*string(T)*","*string(I),ig_)
                arrset(gtype,ig,"<>")
                iv = ix_["Y"*string(T)*","*string(I)]
                pbm.A[ig,iv] += Float64(1.0)
            end
        end
        for T = Int64(v_["1"]):Int64(v_["N-1"])
            v_["T+1"] = 1+T
            ig,ig_,_ = s2mpj_ii("TT"*string(T)*","*string(Int64(v_["1"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"TT"*string(T)*","*string(Int64(v_["1"])))
            iv = ix_["Y"*string(Int64(v_["T+1"]))*","*string(Int64(v_["1"]))]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("TT"*string(T)*","*string(Int64(v_["1"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"TT"*string(T)*","*string(Int64(v_["1"])))
            iv = ix_["Y"*string(T)*","*string(Int64(v_["1"]))]
            pbm.A[ig,iv] += Float64(0.5)
            iv = ix_["Y"*string(T)*","*string(Int64(v_["2"]))]
            pbm.A[ig,iv] += Float64(0.25)
            for I = Int64(v_["1"]):Int64(v_["NX"])
                ig,ig_,_ = s2mpj_ii("TT"*string(T)*","*string(Int64(v_["1"])),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"TT"*string(T)*","*string(Int64(v_["1"])))
                iv = ix_["X"*string(T)*","*string(I)]
                pbm.A[ig,iv] += Float64(v_["B"*string(Int64(v_["1"]))*","*string(I)])
            end
            for J = Int64(v_["2"]):Int64(v_["NY-1"])
                v_["J-1"] = -1+J
                v_["J+1"] = 1+J
                ig,ig_,_ = s2mpj_ii("TT"*string(T)*","*string(J),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"TT"*string(T)*","*string(J))
                iv = ix_["Y"*string(Int64(v_["T+1"]))*","*string(J)]
                pbm.A[ig,iv] += Float64(-1.0)
                iv = ix_["Y"*string(T)*","*string(J)]
                pbm.A[ig,iv] += Float64(0.5)
                iv = ix_["Y"*string(T)*","*string(Int64(v_["J-1"]))]
                pbm.A[ig,iv] += Float64(-0.25)
                iv = ix_["Y"*string(T)*","*string(Int64(v_["J+1"]))]
                pbm.A[ig,iv] += Float64(0.25)
                for I = Int64(v_["1"]):Int64(v_["NX"])
                    ig,ig_,_ = s2mpj_ii("TT"*string(T)*","*string(J),ig_)
                    arrset(gtype,ig,"==")
                    arrset(pb.cnames,ig,"TT"*string(T)*","*string(J))
                    iv = ix_["X"*string(T)*","*string(I)]
                    pbm.A[ig,iv] += Float64(v_["B"*string(J)*","*string(I)])
                end
            end
            ig,ig_,_ = s2mpj_ii("TT"*string(T)*","*string(Int64(v_["NY"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"TT"*string(T)*","*string(Int64(v_["NY"])))
            iv = ix_["Y"*string(Int64(v_["T+1"]))*","*string(Int64(v_["NY"]))]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("TT"*string(T)*","*string(Int64(v_["NY"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"TT"*string(T)*","*string(Int64(v_["NY"])))
            iv = ix_["Y"*string(T)*","*string(Int64(v_["NY"]))]
            pbm.A[ig,iv] += Float64(0.5)
            iv = ix_["Y"*string(T)*","*string(Int64(v_["NY-1"]))]
            pbm.A[ig,iv] += Float64(-0.25)
            for I = Int64(v_["1"]):Int64(v_["NX"])
                ig,ig_,_ = s2mpj_ii("TT"*string(T)*","*string(Int64(v_["NY"])),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"TT"*string(T)*","*string(Int64(v_["NY"])))
                iv = ix_["X"*string(T)*","*string(I)]
                pbm.A[ig,iv] += Float64(v_["B"*string(Int64(v_["NY"]))*","*string(I)])
            end
        end
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        legrps = findall(x->x=="<=",gtype)
        eqgrps = findall(x->x=="==",gtype)
        gegrps = findall(x->x==">=",gtype)
        pb.nle = length(legrps)
        pb.neq = length(eqgrps)
        pb.nge = length(gegrps)
        pb.m   = pb.nle+pb.neq+pb.nge
        pbm.congrps = [[legrps;eqgrps];gegrps]
        pb.nob = ngrp-pb.m
        pbm.objgrps = findall(x->x=="<>",gtype)
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        for T = Int64(v_["1"]):Int64(v_["N-1"])
            for I = Int64(v_["1"]):Int64(v_["NX"])
                pbm.gconst[ig_["OX"*string(T)*","*string(I)]] = Float64(-0.5)
            end
        end
        for T = Int64(v_["1"]):Int64(v_["N"])
            for I = Int64(v_["1"]):Int64(v_["NY"])
                pbm.gconst[ig_["OY"*string(T)*","*string(I)]] = Float64(-0.25)
            end
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        for I = Int64(v_["1"]):Int64(v_["NY"])
            pb.xlower[ix_["Y"*string(Int64(v_["1"]))*","*string(I)]] = 0.0
            pb.xupper[ix_["Y"*string(Int64(v_["1"]))*","*string(I)]] = 0.0
        end
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gL4",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for T = Int64(v_["1"]):Int64(v_["N-1"])
            for I = Int64(v_["1"]):Int64(v_["NX"])
                ig = ig_["OX"*string(T)*","*string(I)]
                arrset(pbm.grftype,ig,"gL4")
            end
        end
        for T = Int64(v_["1"]):Int64(v_["N"])
            for I = Int64(v_["1"]):Int64(v_["NY"])
                ig = ig_["OY"*string(T)*","*string(I)]
                arrset(pbm.grftype,ig,"gL4")
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
# LO SOLUTION(  10,2, 4) 0.0735931360
# LO SOLUTION(  50,2, 4) 0.2299411960
# LO SOLUTION( 100,2, 4) 0.4253681120
# LO SOLUTION( 500,2, 4) 1.9887794988
# LO SOLUTION(1000,2, 4) 3.9430507151
# LO SOLUTION(  10,5,10) 1.1498579294
# LO SOLUTION(  50,5,10) 6.1678479713
# LO SOLUTION( 100,5,10) 12.439954329
# LO SOLUTION( 500,5,10) 62.616843379
# LO SOLUTION(1000,5,10) 125.33793359
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons   = collect(1:length(pbm.congrps))
        pb.pbclass = "C-OLR2-AN-V-V"
        pb.x0          = zeros(Float64,pb.n)
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2]
        pb.conderlvl  = pbm.conderlvl;
        return pb, pbm

# ********************
#  SET UP THE GROUPS *
#  ROUTINE           *
# ********************

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    elseif action == "gL4"

        GVAR_   = args[1]
        igr_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_= GVAR_^4
        if nargout>1
            g_ = 4.0*GVAR_^3
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_ = 12.0*GVAR_^2
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    #%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    elseif action in  ["fx","fgx","fgHx","cx","cJx","cJHx","cIx","cIJx","cIJHx","cIJxv","fHxv",
                       "cJxv","cJtxv","cIJtxv","Lxy","Lgxy","LgHxy","LIxy","LIgxy","LIgHxy",
                       "LHxyv","LIHxyv"]

        pbm = args[1]
        if pbm.name == name
            pbm.has_globs = [0,0]
            return s2mpj_eval(action,args...)
        else
            println("ERROR: please run "*name*" with action = setup")
            return ntuple(i->undef,args[end])
        end

    else
        println("ERROR: action "*action*" unavailable for problem "*name*".jl")
        return ntuple(i->undef,args[end])
    end

end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

