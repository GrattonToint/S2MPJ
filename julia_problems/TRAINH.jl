function TRAINH(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : TRAINH
#    *********
# 
#    The problem is to minimize the energy spent to move a train 
#    from the beginning of a track to its end in a given time.  The train
#    is slowed down by some drag (assumed to be quadratic in the the velocity).
#    The track follows the slope of a hill.  The track geometry is given
#    by the equation
# 
#           1                    1  ns-1                          x - z_i
#    g(x) = - ( a_1 + a_{ns} ) + -- SUM ( s_{i+1} - s_i ) arctan( ------- )
#           2                    pi  1                              eps
# 
#    where the z_i are the breakpoints between sections of the track, and where 
#    the s_i are the "slopes" on these sections (eps is a regularization
#    parameter). Here we have a track of the overall shape
# 
#                      ______
#                     /      \      z0 = 0, z1 = 2, z2 = 4, z3 = 6
#                    /        \     s1 = 2, s2 = 0, s3 = -2
#                   /          \    eps = 0.05
# 
#    The control variables are the acceleration force (UA) and the braking
#    force (UB) applied on the train.
# 
#    Source: adapted from
#    J. Kautsky and N. K. Nichols,
#    "OTEP-2: Optimal Train Energy Programme, mark 2",
#    Numerical Analysis Report NA/4/83,
#    Department of Mathematics, University of Reading, 1983.
# 
#    SIF input: N. Nichols and Ph. Toint, April 1993
# 
#    classification = "C-QOR2-MN-V-V"
# 
#    Number of discretized points in the interval
# 
#       Alternative values for the SIF file parameters:
# IE N                   11             $-PARAMETER n=48, m=22
# IE N                   51             $-PARAMETER n=208, m=102
# IE N                   101            $-PARAMETER n=408, m=202  original value
# IE N                   201            $-PARAMETER n=808, m=402
# IE N                   501            $-PARAMETER n=2008, m=1002 
# IE N                   1001           $-PARAMETER n=4008, m=2002
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "TRAINH"

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
            v_["N"] = Int64(11);  #  SIF file default value
        else
            v_["N"] = Int64(args[1]);
        end
# IE N                   5001           $-PARAMETER n=20008, m=10002
        if nargin<2
            v_["TIME"] = Float64(4.8);  #  SIF file default value
        else
            v_["TIME"] = Float64(args[2]);
        end
        if nargin<3
            v_["LENGTH"] = Float64(6.0);  #  SIF file default value
        else
            v_["LENGTH"] = Float64(args[3]);
        end
        if nargin<4
            v_["NS"] = Int64(3);  #  SIF file default value
        else
            v_["NS"] = Int64(args[4]);
        end
        if nargin<5
            v_["Z1"] = Float64(2.0);  #  SIF file default value
        else
            v_["Z1"] = Float64(args[5]);
        end
        if nargin<6
            v_["Z2"] = Float64(4.0);  #  SIF file default value
        else
            v_["Z2"] = Float64(args[6]);
        end
        if nargin<7
            v_["S1"] = Float64(2.0);  #  SIF file default value
        else
            v_["S1"] = Float64(args[7]);
        end
        if nargin<8
            v_["S2"] = Float64(0.0);  #  SIF file default value
        else
            v_["S2"] = Float64(args[8]);
        end
        if nargin<9
            v_["S3"] = Float64(-2.0);  #  SIF file default value
        else
            v_["S3"] = Float64(args[9]);
        end
        v_["N-1"] = -1+v_["N"]
        v_["RN"] = Float64(v_["N"])
        v_["H"] = v_["TIME"]/v_["RN"]
        v_["H/2"] = 0.5*v_["H"]
        v_["-H"] = -1.0*v_["H"]
        v_["-H/2"] = -1.0*v_["H/2"]
        v_["UAMAX"] = 10.0
        v_["UBMIN"] = -2.0
        v_["VMAX"] = 10.0
        v_["A"] = 0.3
        v_["B"] = 0.14
        v_["C"] = 0.16
        v_["EPS"] = 0.05
        v_["0"] = 0
        v_["1"] = 1
        v_["PI"] = 3.1415926535
        v_["NS-1"] = -1+v_["NS"]
        v_["BH/2"] = v_["B"]*v_["H/2"]
        v_["1+BH/2"] = 1.0+v_["BH/2"]
        v_["BH/2-1"] = -1.0+v_["BH/2"]
        v_["-AH"] = v_["A"]*v_["-H"]
        v_["LENGTH/N"] = v_["LENGTH"]/v_["RN"]
        v_["CH/2"] = v_["C"]*v_["H/2"]
        v_["SUMS"] = v_["S"*string(Int64(v_["1"]))]+v_["S"*string(Int64(v_["NS"]))]
        v_["-AVS"] = -0.5*v_["SUMS"]
        v_["-AVSH"] = v_["-AVS"]*v_["H"]
        v_["CNST"] = v_["-AH"]+v_["-AVSH"]
        v_["H/2PI"] = v_["H/2"]/v_["PI"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["0"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        for I = Int64(v_["0"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("V"*string(I),ix_)
            arrset(pb.xnames,iv,"V"*string(I))
        end
        for I = Int64(v_["0"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("UA"*string(I),ix_)
            arrset(pb.xnames,iv,"UA"*string(I))
        end
        for I = Int64(v_["0"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("UB"*string(I),ix_)
            arrset(pb.xnames,iv,"UB"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("ENERGY",ig_)
        arrset(gtype,ig,"<>")
        for I = Int64(v_["0"]):Int64(v_["N-1"])
            v_["I+1"] = 1+I
            ig,ig_,_ = s2mpj_ii("XEQ"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"XEQ"*string(I))
            iv = ix_["X"*string(Int64(v_["I+1"]))]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["V"*string(Int64(v_["I+1"]))]
            pbm.A[ig,iv] += Float64(v_["-H/2"])
            iv = ix_["V"*string(I)]
            pbm.A[ig,iv] += Float64(v_["-H/2"])
            ig,ig_,_ = s2mpj_ii("VEQ"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"VEQ"*string(I))
            iv = ix_["V"*string(Int64(v_["I+1"]))]
            pbm.A[ig,iv] += Float64(v_["1+BH/2"])
            iv = ix_["V"*string(I)]
            pbm.A[ig,iv] += Float64(v_["BH/2-1"])
            iv = ix_["UA"*string(Int64(v_["I+1"]))]
            pbm.A[ig,iv] += Float64(v_["-H/2"])
            iv = ix_["UA"*string(I)]
            pbm.A[ig,iv] += Float64(v_["-H/2"])
            iv = ix_["UB"*string(Int64(v_["I+1"]))]
            pbm.A[ig,iv] += Float64(v_["-H/2"])
            iv = ix_["UB"*string(I)]
            pbm.A[ig,iv] += Float64(v_["-H/2"])
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
        for I = Int64(v_["0"]):Int64(v_["N-1"])
            pbm.gconst[ig_["VEQ"*string(I)]] = Float64(v_["CNST"])
        end
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower[ix_["X"*string(Int64(v_["0"]))]] = 0.0
        pb.xupper[ix_["X"*string(Int64(v_["0"]))]] = 0.0
        pb.xlower[ix_["V"*string(Int64(v_["0"]))]] = 0.0
        pb.xupper[ix_["V"*string(Int64(v_["0"]))]] = 0.0
        pb.xlower[ix_["UA"*string(Int64(v_["0"]))]] = v_["UAMAX"]
        pb.xupper[ix_["UA"*string(Int64(v_["0"]))]] = v_["UAMAX"]
        pb.xlower[ix_["UB"*string(Int64(v_["0"]))]] = 0.0
        pb.xupper[ix_["UB"*string(Int64(v_["0"]))]] = 0.0
        for I = Int64(v_["1"]):Int64(v_["N-1"])
            pb.xlower[ix_["X"*string(I)]] = -Inf
            pb.xupper[ix_["X"*string(I)]] = +Inf
            pb.xlower[ix_["V"*string(I)]] = -Inf
            pb.xupper[ix_["V"*string(I)]] = +Inf
            pb.xlower[ix_["UA"*string(I)]] = 0.0
            pb.xupper[ix_["UA"*string(I)]] = v_["UAMAX"]
            pb.xlower[ix_["UB"*string(I)]] = v_["UBMIN"]
            pb.xupper[ix_["UB"*string(I)]] = 0.0
        end
        pb.xlower[ix_["X"*string(Int64(v_["N"]))]] = v_["LENGTH"]
        pb.xupper[ix_["X"*string(Int64(v_["N"]))]] = v_["LENGTH"]
        pb.xlower[ix_["V"*string(Int64(v_["N"]))]] = 0.0
        pb.xupper[ix_["V"*string(Int64(v_["N"]))]] = 0.0
        pb.xlower[ix_["UA"*string(Int64(v_["N"]))]] = 0.0
        pb.xupper[ix_["UA"*string(Int64(v_["N"]))]] = 0.0
        pb.xlower[ix_["UB"*string(Int64(v_["N"]))]] = v_["UBMIN"]
        pb.xupper[ix_["UB"*string(Int64(v_["N"]))]] = v_["UBMIN"]
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        pb.x0[ix_["X"*string(Int64(v_["0"]))]] = Float64(0.0)
        pb.x0[ix_["V"*string(Int64(v_["0"]))]] = Float64(0.0)
        pb.x0[ix_["UA"*string(Int64(v_["0"]))]] = Float64(v_["UAMAX"])
        pb.x0[ix_["UB"*string(Int64(v_["0"]))]] = Float64(0.0)
        for I = Int64(v_["1"]):Int64(v_["N-1"])
            v_["RI"] = Float64(I)
            v_["PI"] = v_["LENGTH/N"]*v_["RI"]
            pb.x0[ix_["X"*string(I)]] = Float64(v_["PI"])
            pb.x0[ix_["V"*string(I)]] = Float64(v_["LENGTH/N"])
            pb.x0[ix_["UA"*string(I)]] = Float64(0.0)
            pb.x0[ix_["UB"*string(I)]] = Float64(0.0)
        end
        pb.x0[ix_["X"*string(Int64(v_["N"]))]] = Float64(v_["LENGTH"])
        pb.x0[ix_["V"*string(Int64(v_["N"]))]] = Float64(0.0)
        pb.x0[ix_["UA"*string(Int64(v_["N"]))]] = Float64(0.0)
        pb.x0[ix_["UB"*string(Int64(v_["N"]))]] = Float64(v_["UBMIN"])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "ePROD", iet_)
        loaset(elftv,it,1,"UU")
        loaset(elftv,it,2,"VV")
        it,iet_,_ = s2mpj_ii( "eSQ", iet_)
        loaset(elftv,it,1,"VVV")
        it,iet_,_ = s2mpj_ii( "eATAN", iet_)
        loaset(elftv,it,1,"XX")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"ZZ")
        loaset(elftp,it,2,"E")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["0"]):Int64(v_["N"])
            v_["I+1"] = 1+I
            ename = "VISQ"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSQ")
            arrset(ielftype,ie,iet_["eSQ"])
            vname = "V"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="VVV",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            for J = Int64(v_["1"]):Int64(v_["NS-1"])
                ename = "A"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eATAN")
                arrset(ielftype,ie,iet_["eATAN"])
                vname = "X"*string(I)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="XX",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                posep = findfirst(x->x=="ZZ",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["Z"*string(J)]))
                posep = findfirst(x->x=="E",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["EPS"]))
            end
        end
        for I = Int64(v_["1"]):Int64(v_["N-1"])
            ename = "UV"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"ePROD")
            arrset(ielftype,ie,iet_["ePROD"])
            vname = "UA"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="UU",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "V"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="VV",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["0"]):Int64(v_["N-1"])
            v_["I+1"] = 1+I
            ig = ig_["VEQ"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["VISQ"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["CH/2"]))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["VISQ"*string(Int64(v_["I+1"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["CH/2"]))
            for J = Int64(v_["1"]):Int64(v_["NS-1"])
                v_["J+1"] = 1+J
                v_["DS"] = v_["S"*string(Int64(v_["J+1"]))]-v_["S"*string(J)]
                v_["WJ"] = v_["DS"]*v_["H/2PI"]
                ig = ig_["VEQ"*string(I)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A"*string(I)*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["WJ"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A"*string(Int64(v_["I+1"]))*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["WJ"]))
            end
        end
        for I = Int64(v_["1"]):Int64(v_["N-1"])
            ig = ig_["ENERGY"]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["UV"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["H"]))
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
# LO SOLUTION(11)        12.423025536
# LO SOLUTION(51)        12.295777964
# LO SOLUTION(101)       12.306399739
# LO SOLUTION(201)       12.309848614
# LO SOLUTION(1001)      12.307801327
# LO SOLUTION(5001)      12.221148056
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
        pb.pbclass = "C-QOR2-MN-V-V"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2]
        pb.conderlvl  = pbm.conderlvl;
        return pb, pbm

# ********************
#  SET UP THE GROUPS *
#  ROUTINE           *
# ********************

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

    elseif action == "ePROD"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]*EV_[2]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[2]
            g_[2] = EV_[1]
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = 1.0
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

    elseif action == "eATAN"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        DX = EV_[1]-pbm.elpar[iel_][1]
        E2 = pbm.elpar[iel_][2]*pbm.elpar[iel_][2]
        f_   = atan(DX/pbm.elpar[iel_][2])
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = pbm.elpar[iel_][2]/(E2+DX*DX)
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = -2.0*DX*pbm.elpar[iel_][2]/(E2+DX*DX)^2
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

