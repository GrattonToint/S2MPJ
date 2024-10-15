function DTOC2(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DTOC2
#    *********
# 
#    This is a discrete time optimal control (DTOC) problem.  
#    The system has N time periods, 2 control variables and 4 state variables.
# 
#    The problem is not convex.
# 
#    Sources: problem 2 in
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
#    classification = "C-OOR2-AN-V-V"
# 
#    Problem variants: they are identified by the value of the parameter N.
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

    name = "DTOC2"

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
        v_["3"] = 3
        v_["4"] = 4
        v_["NY-1"] = -1+v_["NY"]
        v_["2NY"] = v_["NY"]+v_["NY"]
        v_["R2NY"] = Float64(v_["2NY"])
        v_["1/2NY"] = 1.0/v_["R2NY"]
        for J = Int64(v_["1"]):Int64(v_["NX"])
            for I = Int64(v_["1"]):Int64(v_["NY"])
                v_["I+J"] = I+J
                v_["RI+J"] = Float64(v_["I+J"])
                v_["C"*string(I)*","*string(J)] = v_["RI+J"]*v_["1/2NY"]
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
        for T = Int64(v_["1"]):Int64(v_["N"])
            ig,ig_,_ = s2mpj_ii("O"*string(T),ig_)
            arrset(gtype,ig,"<>")
        end
        for T = Int64(v_["1"]):Int64(v_["N-1"])
            v_["T+1"] = 1+T
            for J = Int64(v_["1"]):Int64(v_["NY"])
                ig,ig_,_ = s2mpj_ii("TT"*string(T)*","*string(J),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"TT"*string(T)*","*string(J))
                iv = ix_["Y"*string(Int64(v_["T+1"]))*","*string(J)]
                pbm.A[ig,iv] += Float64(-1.0)
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
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        for I = Int64(v_["1"]):Int64(v_["NY"])
            v_["RI"] = Float64(I)
            v_["TMP"] = v_["RI"]*v_["1/2NY"]
            pb.xlower[ix_["Y"*string(Int64(v_["1"]))*","*string(I)]] = v_["TMP"]
            pb.xupper[ix_["Y"*string(Int64(v_["1"]))*","*string(I)]] = v_["TMP"]
        end
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        for I = Int64(v_["1"]):Int64(v_["NY"])
            v_["RI"] = Float64(I)
            v_["TMP"] = v_["RI"]*v_["1/2NY"]
            pb.x0[ix_["Y"*string(Int64(v_["1"]))*","*string(I)]] = Float64(v_["TMP"])
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eOEL", iet_)
        loaset(elftv,it,1,"YY1")
        loaset(elftv,it,2,"YY2")
        loaset(elftv,it,3,"YY3")
        loaset(elftv,it,4,"YY4")
        loaset(elftv,it,5,"XX1")
        loaset(elftv,it,6,"XX2")
        it,iet_,_ = s2mpj_ii( "eSQ", iet_)
        loaset(elftv,it,1,"YY")
        it,iet_,_ = s2mpj_ii( "eSINE", iet_)
        loaset(elftv,it,1,"ZZ")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for T = Int64(v_["1"]):Int64(v_["N-1"])
            ename = "EO"*string(T)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eOEL")
            arrset(ielftype,ie,iet_["eOEL"])
            vname = "Y"*string(T)*","*string(Int64(v_["1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="YY1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "Y"*string(T)*","*string(Int64(v_["2"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="YY2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "Y"*string(T)*","*string(Int64(v_["3"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="YY3",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "Y"*string(T)*","*string(Int64(v_["4"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="YY4",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(T)*","*string(Int64(v_["1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="XX1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(T)*","*string(Int64(v_["2"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="XX2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            for J = Int64(v_["1"]):Int64(v_["NY"])
                ename = "SY"*string(T)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eSINE")
                arrset(ielftype,ie,iet_["eSINE"])
                vname = "Y"*string(T)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="ZZ",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
            for I = Int64(v_["1"]):Int64(v_["NX"])
                ename = "SX"*string(T)*","*string(I)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eSINE")
                arrset(ielftype,ie,iet_["eSINE"])
                vname = "X"*string(T)*","*string(I)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="ZZ",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
        end
        for J = Int64(v_["1"]):Int64(v_["NY"])
            ename = "YNSQ"*string(J)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSQ")
            arrset(ielftype,ie,iet_["eSQ"])
            vname = "Y"*string(Int64(v_["N"]))*","*string(J)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="YY",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for T = Int64(v_["1"]):Int64(v_["N-1"])
            ig = ig_["O"*string(T)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["EO"*string(T)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        for J = Int64(v_["1"]):Int64(v_["NY"])
            ig = ig_["O"*string(Int64(v_["N"]))]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["YNSQ"*string(J)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        for T = Int64(v_["1"]):Int64(v_["N-1"])
            for J = Int64(v_["1"]):Int64(v_["NY"])
                ig = ig_["TT"*string(T)*","*string(J)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["SY"*string(T)*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,1.)
                for I = Int64(v_["1"]):Int64(v_["NX"])
                    ig = ig_["TT"*string(T)*","*string(J)]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["SX"*string(T)*","*string(I)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["C"*string(J)*","*string(I)]))
                end
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
# LO SOLUTION(  10)      0.485983918948
# LO SOLUTION(  20)      0.486212213803
# LO SOLUTION(  30)      0.486383392574
# LO SOLUTION(  40)      0.486572686778
# LO SOLUTION(  50)      0.486884900389
# LO SOLUTION( 100)      0.487532342563
# LO SOLUTION( 500)      0.490996540460
# LO SOLUTION(1000)      0.490200910983
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
        pb.pbclass = "C-OOR2-AN-V-V"
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

    elseif action == "eSQ"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]*EV_[1]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[1]+EV_[1]
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 2.0
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eSINE"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        SZ = sin(EV_[1])
        f_   = SZ
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = cos(EV_[1])
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = -SZ
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eOEL"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        XN2 = EV_[5]*EV_[5]+EV_[6]*EV_[6]
        YN2 = EV_[1]*EV_[1]+EV_[2]*EV_[2]+EV_[3]*EV_[3]+EV_[4]*EV_[4]
        SZ = sin(0.5*XN2)
        CZ = cos(0.5*XN2)
        SZ2 = SZ*SZ+1.0
        SC = SZ*CZ
        CCSS = CZ*CZ-SZ*SZ
        f_   = YN2*SZ2
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[5] = 2.0*YN2*SC*EV_[5]
            g_[6] = 2.0*YN2*SC*EV_[6]
            g_[1] = 2.0*EV_[1]*SZ2
            g_[2] = 2.0*EV_[2]*SZ2
            g_[3] = 2.0*EV_[3]*SZ2
            g_[4] = 2.0*EV_[4]*SZ2
            if nargout>2
                H_ = zeros(Float64,6,6)
                H_[5,5] = 2.0*YN2*(SC+EV_[5]*EV_[5]*CCSS)
                H_[5,6] = 2.0*YN2*EV_[5]*EV_[6]*CCSS
                H_[6,5] = H_[5,6]
                H_[5,1] = 4.0*EV_[1]*SC*EV_[5]
                H_[1,5] = H_[5,1]
                H_[5,2] = 4.0*EV_[2]*SC*EV_[5]
                H_[2,5] = H_[5,2]
                H_[5,3] = 4.0*EV_[3]*SC*EV_[5]
                H_[3,5] = H_[5,3]
                H_[5,4] = 4.0*EV_[4]*SC*EV_[5]
                H_[4,5] = H_[5,4]
                H_[6,6] = 2.0*YN2*(SC+EV_[6]*EV_[6]*CCSS)
                H_[6,1] = 4.0*EV_[1]*SC*EV_[6]
                H_[1,6] = H_[6,1]
                H_[6,2] = 4.0*EV_[2]*SC*EV_[6]
                H_[2,6] = H_[6,2]
                H_[6,3] = 4.0*EV_[3]*SC*EV_[6]
                H_[3,6] = H_[6,3]
                H_[6,4] = 4.0*EV_[4]*SC*EV_[6]
                H_[4,6] = H_[6,4]
                H_[1,1] = 2.0*SZ2
                H_[2,2] = 2.0*SZ2
                H_[3,3] = 2.0*SZ2
                H_[4,4] = 2.0*SZ2
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

