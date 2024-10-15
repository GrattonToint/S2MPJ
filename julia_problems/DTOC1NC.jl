function DTOC1NC(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DTOC1NC
#    *********
# 
#    This is a discrete time optimal control (DTOC) problem.  
#    The system has N time periods, NX control variables and NY state variables.
#    The nonlinearity parameter mu is set to 0.5.
# 
#    The problem is not convex.
# 
#    Sources: problem 1 in
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
#    classification = "C-OQR2-AN-V-V"
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

    name = "DTOC1NC"

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
        if nargin<4
            v_["MU"] = Float64(0.5);  #  SIF file default value
        else
            v_["MU"] = Float64(args[4]);
        end
        v_["N-1"] = -1+v_["N"]
        v_["0"] = 0
        v_["1"] = 1
        v_["2"] = 2
        v_["NY-1"] = -1+v_["NY"]
        v_["NX+NY"] = v_["NX"]+v_["NY"]
        v_["RXY"] = Float64(v_["NX+NY"])
        v_["1/RXY"] = 1.0/v_["RXY"]
        v_["MU/RXY"] = v_["MU"]*v_["1/RXY"]
        v_["NYNX"] = v_["NX"]*v_["NY"]
        v_["NYNX-1"] = -1+v_["NYNX"]
        for J = Int64(v_["1"]):Int64(v_["NX"])
            for I = Int64(v_["1"]):Int64(v_["NY"])
                v_["I-J"] = I-J
                v_["RI-J"] = Float64(v_["I-J"])
                v_["B"*string(I)*","*string(J)] = v_["RI-J"]*v_["1/RXY"]
                v_["I+J"] = I+J
                v_["RI+J"] = Float64(v_["I+J"])
                v_["C"*string(I)*","*string(J)] = v_["RI+J"]*v_["MU/RXY"]
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
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "ePR", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"MUC")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for T = Int64(v_["1"]):Int64(v_["N-1"])
            for K = Int64(v_["0"]):Int64(v_["NYNX-1"])
                v_["I"] = trunc(Int,(K/v_["NX"]))
                v_["INX"] = v_["I"]*v_["NX"]
                v_["J"] = K-v_["INX"]
                v_["I"] = 1+v_["I"]
                v_["J"] = 1+v_["J"]
                ename = "E"*string(T)*","*string(K)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"ePR")
                arrset(ielftype,ie,iet_["ePR"])
                vname = "Y"*string(T)*","*string(Int64(v_["I"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(T)*","*string(Int64(v_["J"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="X",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                posep = findfirst(x->x=="MUC",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["C"*string(Int64(v_["I"]))*","*string(Int64(v_["J"]))]))
            end
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
        for T = Int64(v_["1"]):Int64(v_["N-1"])
            for J = Int64(v_["1"]):Int64(v_["NY"])
                for K = Int64(v_["0"]):Int64(v_["NYNX-1"])
                    ig = ig_["TT"*string(T)*","*string(J)]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["E"*string(T)*","*string(K)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,1.)
                end
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
# LO S(  10,2, 4)        0.31231015014
# LO S(  50,2, 4)        1.72189403445
# LO S( 100,2, 4)        3.48385849052
# LO S( 500,2, 4)        17.5796246912
# LO S(1000,2, 4)        35.1993544122
# LO S(  10,5,10)        2.15407839920
# LO S(  50,5,10)        12.2942814975
# LO S( 100,5,10)        24.9697864449
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
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-OQR2-AN-V-V"
        pb.x0          = zeros(Float64,pb.n)
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2]
        pb.conderlvl  = pbm.conderlvl;
        return pb, pbm

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "ePR"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = pbm.elpar[iel_][1]*EV_[1]*EV_[2]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = pbm.elpar[iel_][1]*EV_[2]
            g_[2] = pbm.elpar[iel_][1]*EV_[1]
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = pbm.elpar[iel_][1]
                H_[2,1] = H_[1,2]
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

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

