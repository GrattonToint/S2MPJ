function DTOC4(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DTOC4
#    *********
# 
#    This is a discrete time optimal control (DTOC) problem.  
#    The system has N time periods, 1 control variable and 2 state variables.
# 
#    The problem is not convex.
# 
#    Sources: problem 4 in
#    T.F. Coleman and A. Liao,
#    "An Efficient Trust Region Method for Unconstrained Discret-Time Optimal
#    Control Problems",
#    Tech. Report, ctc93tr144,  Advanced Computing Research Institute, 
#    Cornell University, 1992.
# 
#    G. Di Pillo, L. Grippo and F. Lampariello,
#    "A class of structures quasi-Newton algorithms for optimal control
#    problems",
#    in H.E. Rauch, ed., IFAC Applications of nonlinear programming to
#    optimization and control, pp. 101-107, IFAC, Pergamon Press, 1983.
# 
#    SIF input: Ph. Toint, August 1993
# 
#    classification = "C-QOR2-AN-V-V"
# 
#    Problem variants: they are identified by the value of the parameter N.
# 
#    The problem has 3N-1  variables (of which 2 are fixed),
#    and 2(N-1) constraints
# 
#       Alternative values for the SIF file parameters:
# IE N                   10             $-PARAMETER  n=   29,m= 18 original value
# IE N                   50             $-PARAMETER  n=  149,m= 98
# IE N                   100            $-PARAMETER  n=  299,m=198
# IE N                   500            $-PARAMETER  n= 1499,m=998
# IE N                   1000           $-PARAMETER  n= 2999,m=1998
# IE N                   1500           $-PARAMETER  n= 4499,m=2998
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "DTOC4"

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
# IE N                   5000           $-PARAMETER  n=14999,m=9998
        v_["N-1"] = -1+v_["N"]
        v_["1"] = 1
        v_["2"] = 2
        v_["RN"] = Float64(v_["N"])
        v_["H"] = 1.0/v_["RN"]
        v_["5H"] = 5.0*v_["H"]
        v_["1/5H"] = 1.0/v_["5H"]
        v_["1+5H"] = 1.0+v_["5H"]
        v_["-5H"] = -1.0*v_["5H"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for T = Int64(v_["1"]):Int64(v_["N-1"])
            iv,ix_,_ = s2mpj_ii("X"*string(T),ix_)
            arrset(pb.xnames,iv,"X"*string(T))
        end
        for T = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("Y"*string(T)*","*string(Int64(v_["1"])),ix_)
            arrset(pb.xnames,iv,"Y"*string(T)*","*string(Int64(v_["1"])))
            iv,ix_,_ = s2mpj_ii("Y"*string(T)*","*string(Int64(v_["2"])),ix_)
            arrset(pb.xnames,iv,"Y"*string(T)*","*string(Int64(v_["2"])))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        arrset(pbm.gscale,ig,Float64(v_["1/5H"]))
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
            pbm.A[ig,iv] += Float64(v_["1+5H"])
            ig,ig_,_ = s2mpj_ii("TT"*string(T)*","*string(Int64(v_["1"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"TT"*string(T)*","*string(Int64(v_["1"])))
            iv = ix_["Y"*string(T)*","*string(Int64(v_["2"]))]
            pbm.A[ig,iv] += Float64(v_["-5H"])
            ig,ig_,_ = s2mpj_ii("TT"*string(T)*","*string(Int64(v_["1"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"TT"*string(T)*","*string(Int64(v_["1"])))
            iv = ix_["X"*string(T)]
            pbm.A[ig,iv] += Float64(v_["5H"])
            ig,ig_,_ = s2mpj_ii("TT"*string(T)*","*string(Int64(v_["2"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"TT"*string(T)*","*string(Int64(v_["2"])))
            iv = ix_["Y"*string(Int64(v_["T+1"]))*","*string(Int64(v_["2"]))]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["Y"*string(T)*","*string(Int64(v_["2"]))]
            pbm.A[ig,iv] += Float64(1.0)
            ig,ig_,_ = s2mpj_ii("TT"*string(T)*","*string(Int64(v_["2"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"TT"*string(T)*","*string(Int64(v_["2"])))
            iv = ix_["Y"*string(T)*","*string(Int64(v_["1"]))]
            pbm.A[ig,iv] += Float64(v_["5H"])
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
        pb.xlower[ix_["Y"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))]] = 0.0
        pb.xupper[ix_["Y"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))]] = 0.0
        pb.xlower[ix_["Y"*string(Int64(v_["1"]))*","*string(Int64(v_["2"]))]] = 1.0
        pb.xupper[ix_["Y"*string(Int64(v_["1"]))*","*string(Int64(v_["2"]))]] = 1.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        pb.x0[ix_["Y"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))]]  = (
              Float64(0.0))
        pb.x0[ix_["Y"*string(Int64(v_["1"]))*","*string(Int64(v_["2"]))]]  = (
              Float64(1.0))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eSQ", iet_)
        loaset(elftv,it,1,"Z")
        it,iet_,_ = s2mpj_ii( "eAAB", iet_)
        loaset(elftv,it,1,"A")
        loaset(elftv,it,2,"B")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "Y1SQ"*string(Int64(v_["1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        ename = "Y1SQ"*string(Int64(v_["1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "Y"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "Y2SQ"*string(Int64(v_["1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        ename = "Y2SQ"*string(Int64(v_["1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "Y"*string(Int64(v_["1"]))*","*string(Int64(v_["2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "XSQ"*string(Int64(v_["1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        ename = "XSQ"*string(Int64(v_["1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        for T = Int64(v_["2"]):Int64(v_["N-1"])
            ename = "Y1SQ"*string(T)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSQ")
            arrset(ielftype,ie,iet_["eSQ"])
            vname = "Y"*string(T)*","*string(Int64(v_["1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "Y2SQ"*string(T)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSQ")
            arrset(ielftype,ie,iet_["eSQ"])
            vname = "Y"*string(T)*","*string(Int64(v_["2"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "XSQ"*string(T)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSQ")
            arrset(ielftype,ie,iet_["eSQ"])
            vname = "X"*string(T)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        ename = "Y1SQ"*string(Int64(v_["N"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        ename = "Y1SQ"*string(Int64(v_["N"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "Y"*string(Int64(v_["N"]))*","*string(Int64(v_["1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "Y2SQ"*string(Int64(v_["N"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        ename = "Y2SQ"*string(Int64(v_["N"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "Y"*string(Int64(v_["N"]))*","*string(Int64(v_["2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        for T = Int64(v_["1"]):Int64(v_["N-1"])
            ename = "E"*string(T)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eAAB")
            arrset(ielftype,ie,iet_["eAAB"])
            vname = "Y"*string(T)*","*string(Int64(v_["2"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="A",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "Y"*string(T)*","*string(Int64(v_["1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="B",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["OBJ"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["Y1SQ"*string(Int64(v_["1"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(0.5))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["Y2SQ"*string(Int64(v_["1"]))])
        loaset(pbm.grelw,ig,posel,Float64(0.5))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["XSQ"*string(Int64(v_["1"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(1.0))
        for T = Int64(v_["2"]):Int64(v_["N-1"])
            ig = ig_["OBJ"]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["Y1SQ"*string(T)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["Y2SQ"*string(T)])
            loaset(pbm.grelw,ig,posel, 1.)
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["XSQ"*string(T)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        ig = ig_["OBJ"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["Y1SQ"*string(Int64(v_["N"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(0.5))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["Y2SQ"*string(Int64(v_["N"]))])
        loaset(pbm.grelw,ig,posel,Float64(0.5))
        for T = Int64(v_["1"]):Int64(v_["N-1"])
            ig = ig_["TT"*string(T)*","*string(Int64(v_["1"]))]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["E"*string(T)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["-5H"]))
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
# LO SOLUTION(  10)      3.75078392210
# LO SOLUTION(  50)      3.02963141755
# LO SOLUTION( 100)      2.94726711402
# LO SOLUTION( 500)      2.87827434035
# LO SOLUTION(1000)      2.87483889886
# LO SOLUTION(5000)      2.86386891514
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
        pb.pbclass = "C-QOR2-AN-V-V"
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

    elseif action == "eAAB"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]*EV_[1]*EV_[2]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 2.0*EV_[1]*EV_[2]
            g_[2] = EV_[1]*EV_[1]
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = EV_[2]+EV_[2]
                H_[1,2] = EV_[1]+EV_[1]
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

