function WATSONNE(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : WATSONNE
#    *********
# 
#    Watson problem in 12 variables. This is a nonlinear equation version
#    of problem WATSON.
# 
#    This function  is a nonlinear least squares with 31 groups.  Each
#    group has 1 nonlinear and 1 linear elements.
# 
#    Source:  problem 20 in
#    J.J. More', B.S. Garbow and K.E. Hillstrom,
#    "Testing Unconstrained Optimization Software",
#    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
# 
#    See also Buckley#128 (p. 100).
# 
#    SIF input: Ph. Toint, Dec 1989.
#    Modification as a set of nonlinear equations: Nick Gould, Oct 2015.
# 
#    classification = "C-NOR2-AN-V-31"
# 
#    The number of variables can be varied, but should be smaller than
#    31
# 
#    Number of variables
# 
#       Alternative values for the SIF file parameters:
# IE N                   12             $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "WATSONNE"

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
            v_["N"] = Int64(12);  #  SIF file default value
        else
            v_["N"] = Int64(args[1]);
        end
# IE N                   31             $-PARAMETER
        v_["M"] = 31
        v_["1"] = 1
        v_["2"] = 2
        v_["29"] = 29
        v_["30"] = 30
        v_["1/29"] = 0.024482759
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["1"]):Int64(v_["29"])
            v_["RI"] = Float64(I)
            v_["TI"] = v_["RI"]*v_["1/29"]
            v_["LNTI"] = log(v_["TI"])
            for J = Int64(v_["2"]):Int64(v_["N"])
                v_["RJ"] = Float64(J)
                v_["RJ-1"] = -1.0+v_["RJ"]
                v_["RJ-2"] = -2.0+v_["RJ"]
                v_["AE"] = v_["RJ-2"]*v_["LNTI"]
                v_["C0"] = exp(v_["AE"])
                v_["C"] = v_["C0"]*v_["RJ-1"]
                ig,ig_,_ = s2mpj_ii("G"*string(I),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"G"*string(I))
                iv = ix_["X"*string(J)]
                pbm.A[ig,iv] += Float64(v_["C"])
            end
        end
        ig,ig_,_ = s2mpj_ii("G"*string(Int64(v_["30"])),ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"G"*string(Int64(v_["30"])))
        iv = ix_["X"*string(Int64(v_["1"]))]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("G"*string(Int64(v_["M"])),ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"G"*string(Int64(v_["M"])))
        iv = ix_["X"*string(Int64(v_["2"]))]
        pbm.A[ig,iv] += Float64(1.0)
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
        #%%%%%%%%%%%%%%%%%%  CONSTANTS %%%%%%%%%%%%%%%%%%%
        pbm.gconst = fill(1.0,ngrp)
        pbm.gconst[ig_["G"*string(Int64(v_["30"]))]] = Float64(0.0)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eMWSQ", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        loaset(elftv,it,3,"V3")
        loaset(elftv,it,4,"V4")
        loaset(elftv,it,5,"V5")
        loaset(elftv,it,6,"V6")
        loaset(elftv,it,7,"V7")
        loaset(elftv,it,8,"V8")
        loaset(elftv,it,9,"V9")
        loaset(elftv,it,10,"V10")
        loaset(elftv,it,11,"V11")
        loaset(elftv,it,12,"V12")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"T1")
        loaset(elftp,it,2,"T2")
        loaset(elftp,it,3,"T3")
        loaset(elftp,it,4,"T4")
        loaset(elftp,it,5,"T5")
        loaset(elftp,it,6,"T6")
        loaset(elftp,it,7,"T7")
        loaset(elftp,it,8,"T8")
        loaset(elftp,it,9,"T9")
        loaset(elftp,it,10,"T10")
        loaset(elftp,it,11,"T11")
        loaset(elftp,it,12,"T12")
        it,iet_,_ = s2mpj_ii( "eMSQ", iet_)
        loaset(elftv,it,1,"V1")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["29"])
            v_["RI"] = Float64(I)
            v_["TI"] = v_["RI"]*v_["1/29"]
            v_["LNTI"] = log(v_["TI"])
            ename = "E"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eMWSQ")
            arrset(ielftype,ie,iet_["eMWSQ"])
            vname = "X1"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X2"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X3"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X4"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V4",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X5"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V5",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X6"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V6",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X7"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V7",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X8"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V8",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X9"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V9",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X10"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V10",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X11"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V11",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X12"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V12",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            for J = Int64(v_["1"]):Int64(v_["N"])
                v_["J-1"] = -1+J
                v_["RJ-1"] = Float64(v_["J-1"])
                v_["CE0"] = v_["RJ-1"]*v_["LNTI"]
                v_["CE"*string(J)] = exp(v_["CE0"])
            end
            ename = "E"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            posep = findfirst(x->x=="T1",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["CE1"]))
            posep = findfirst(x->x=="T2",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["CE2"]))
            posep = findfirst(x->x=="T3",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["CE3"]))
            posep = findfirst(x->x=="T4",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["CE4"]))
            posep = findfirst(x->x=="T5",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["CE5"]))
            posep = findfirst(x->x=="T6",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["CE6"]))
            posep = findfirst(x->x=="T7",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["CE7"]))
            posep = findfirst(x->x=="T8",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["CE8"]))
            posep = findfirst(x->x=="T9",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["CE9"]))
            posep = findfirst(x->x=="T10",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["CE10"]))
            posep = findfirst(x->x=="T11",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["CE11"]))
            posep = findfirst(x->x=="T12",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["CE12"]))
        end
        ename = "E"*string(Int64(v_["M"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eMSQ")
        arrset(ielftype,ie,iet_["eMSQ"])
        ename = "E"*string(Int64(v_["M"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["29"])
            ig = ig_["G"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["E"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        ig = ig_["G"*string(Int64(v_["M"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["M"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
        pb.objlower = 0.0
#    Solution
# LO SOLTN(12)           2.27559922D-9
# LO SOLTN(31)           1.53795068D-9
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
        pb.pbclass = "C-NOR2-AN-V-31"
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

    elseif action == "eMSQ"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = -EV_[1]*EV_[1]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = -EV_[1]-EV_[1]
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = -2.0
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eMWSQ"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U  = (
              pbm.elpar[iel_][1]*EV_[1]+pbm.elpar[iel_][2]*EV_[2]+pbm.elpar[iel_][3]*EV_[3]+pbm.elpar[iel_][4]*EV_[4]+pbm.elpar[iel_][5]*EV_[5]+pbm.elpar[iel_][6]*EV_[6]+pbm.elpar[iel_][7]*EV_[7]+pbm.elpar[iel_][8]*EV_[8]+pbm.elpar[iel_][9]*EV_[9]+pbm.elpar[iel_][10]*EV_[10]+pbm.elpar[iel_][11]*EV_[11]+pbm.elpar[iel_][12]*EV_[12])
        TWOT1 = pbm.elpar[iel_][1]+pbm.elpar[iel_][1]
        TWOT2 = pbm.elpar[iel_][2]+pbm.elpar[iel_][2]
        TWOT3 = pbm.elpar[iel_][3]+pbm.elpar[iel_][3]
        TWOT4 = pbm.elpar[iel_][4]+pbm.elpar[iel_][4]
        TWOT5 = pbm.elpar[iel_][5]+pbm.elpar[iel_][5]
        TWOT6 = pbm.elpar[iel_][6]+pbm.elpar[iel_][6]
        TWOT7 = pbm.elpar[iel_][7]+pbm.elpar[iel_][7]
        TWOT8 = pbm.elpar[iel_][8]+pbm.elpar[iel_][8]
        TWOT9 = pbm.elpar[iel_][9]+pbm.elpar[iel_][9]
        TWOT10 = pbm.elpar[iel_][10]+pbm.elpar[iel_][10]
        TWOT11 = pbm.elpar[iel_][11]+pbm.elpar[iel_][11]
        TWOT12 = pbm.elpar[iel_][12]+pbm.elpar[iel_][12]
        f_   = -U*U
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = -TWOT1*U
            g_[2] = -TWOT2*U
            g_[3] = -TWOT3*U
            g_[4] = -TWOT4*U
            g_[5] = -TWOT5*U
            g_[6] = -TWOT6*U
            g_[7] = -TWOT7*U
            g_[8] = -TWOT8*U
            g_[9] = -TWOT9*U
            g_[10] = -TWOT10*U
            g_[11] = -TWOT11*U
            g_[12] = -TWOT12*U
            if nargout>2
                H_ = zeros(Float64,12,12)
                H_[1,1] = -TWOT1*pbm.elpar[iel_][1]
                H_[1,2] = -TWOT1*pbm.elpar[iel_][2]
                H_[2,1] = H_[1,2]
                H_[1,3] = -TWOT1*pbm.elpar[iel_][3]
                H_[3,1] = H_[1,3]
                H_[1,4] = -TWOT1*pbm.elpar[iel_][4]
                H_[4,1] = H_[1,4]
                H_[1,5] = -TWOT1*pbm.elpar[iel_][5]
                H_[5,1] = H_[1,5]
                H_[1,6] = -TWOT1*pbm.elpar[iel_][6]
                H_[6,1] = H_[1,6]
                H_[1,7] = -TWOT1*pbm.elpar[iel_][7]
                H_[7,1] = H_[1,7]
                H_[1,8] = -TWOT1*pbm.elpar[iel_][8]
                H_[8,1] = H_[1,8]
                H_[1,9] = -TWOT1*pbm.elpar[iel_][9]
                H_[9,1] = H_[1,9]
                H_[1,10] = -TWOT1*pbm.elpar[iel_][10]
                H_[10,1] = H_[1,10]
                H_[1,11] = -TWOT1*pbm.elpar[iel_][11]
                H_[11,1] = H_[1,11]
                H_[1,12] = -TWOT1*pbm.elpar[iel_][12]
                H_[12,1] = H_[1,12]
                H_[2,2] = -TWOT2*pbm.elpar[iel_][2]
                H_[2,3] = -TWOT2*pbm.elpar[iel_][3]
                H_[3,2] = H_[2,3]
                H_[2,4] = -TWOT2*pbm.elpar[iel_][4]
                H_[4,2] = H_[2,4]
                H_[2,5] = -TWOT2*pbm.elpar[iel_][5]
                H_[5,2] = H_[2,5]
                H_[2,6] = -TWOT2*pbm.elpar[iel_][6]
                H_[6,2] = H_[2,6]
                H_[2,7] = -TWOT2*pbm.elpar[iel_][7]
                H_[7,2] = H_[2,7]
                H_[2,8] = -TWOT2*pbm.elpar[iel_][8]
                H_[8,2] = H_[2,8]
                H_[2,9] = -TWOT2*pbm.elpar[iel_][8]
                H_[9,2] = H_[2,9]
                H_[2,10] = -TWOT2*pbm.elpar[iel_][10]
                H_[10,2] = H_[2,10]
                H_[2,11] = -TWOT2*pbm.elpar[iel_][11]
                H_[11,2] = H_[2,11]
                H_[2,12] = -TWOT2*pbm.elpar[iel_][12]
                H_[12,2] = H_[2,12]
                H_[3,3] = -TWOT3*pbm.elpar[iel_][3]
                H_[3,4] = -TWOT3*pbm.elpar[iel_][4]
                H_[4,3] = H_[3,4]
                H_[3,5] = -TWOT3*pbm.elpar[iel_][5]
                H_[5,3] = H_[3,5]
                H_[3,6] = -TWOT3*pbm.elpar[iel_][6]
                H_[6,3] = H_[3,6]
                H_[3,7] = -TWOT3*pbm.elpar[iel_][7]
                H_[7,3] = H_[3,7]
                H_[3,8] = -TWOT3*pbm.elpar[iel_][8]
                H_[8,3] = H_[3,8]
                H_[3,9] = -TWOT3*pbm.elpar[iel_][8]
                H_[9,3] = H_[3,9]
                H_[3,10] = -TWOT3*pbm.elpar[iel_][10]
                H_[10,3] = H_[3,10]
                H_[3,11] = -TWOT3*pbm.elpar[iel_][11]
                H_[11,3] = H_[3,11]
                H_[3,12] = -TWOT3*pbm.elpar[iel_][12]
                H_[12,3] = H_[3,12]
                H_[4,4] = -TWOT4*pbm.elpar[iel_][4]
                H_[4,5] = -TWOT4*pbm.elpar[iel_][5]
                H_[5,4] = H_[4,5]
                H_[4,6] = -TWOT4*pbm.elpar[iel_][6]
                H_[6,4] = H_[4,6]
                H_[4,7] = -TWOT4*pbm.elpar[iel_][7]
                H_[7,4] = H_[4,7]
                H_[4,8] = -TWOT4*pbm.elpar[iel_][8]
                H_[8,4] = H_[4,8]
                H_[4,9] = -TWOT4*pbm.elpar[iel_][8]
                H_[9,4] = H_[4,9]
                H_[4,10] = -TWOT4*pbm.elpar[iel_][10]
                H_[10,4] = H_[4,10]
                H_[4,11] = -TWOT4*pbm.elpar[iel_][11]
                H_[11,4] = H_[4,11]
                H_[4,12] = -TWOT4*pbm.elpar[iel_][12]
                H_[12,4] = H_[4,12]
                H_[5,5] = -TWOT5*pbm.elpar[iel_][5]
                H_[5,6] = -TWOT5*pbm.elpar[iel_][6]
                H_[6,5] = H_[5,6]
                H_[5,7] = -TWOT5*pbm.elpar[iel_][7]
                H_[7,5] = H_[5,7]
                H_[5,8] = -TWOT5*pbm.elpar[iel_][8]
                H_[8,5] = H_[5,8]
                H_[5,9] = -TWOT5*pbm.elpar[iel_][8]
                H_[9,5] = H_[5,9]
                H_[5,10] = -TWOT5*pbm.elpar[iel_][10]
                H_[10,5] = H_[5,10]
                H_[5,11] = -TWOT5*pbm.elpar[iel_][11]
                H_[11,5] = H_[5,11]
                H_[5,12] = -TWOT5*pbm.elpar[iel_][12]
                H_[12,5] = H_[5,12]
                H_[6,6] = -TWOT6*pbm.elpar[iel_][6]
                H_[6,7] = -TWOT6*pbm.elpar[iel_][7]
                H_[7,6] = H_[6,7]
                H_[6,8] = -TWOT6*pbm.elpar[iel_][8]
                H_[8,6] = H_[6,8]
                H_[6,9] = -TWOT6*pbm.elpar[iel_][8]
                H_[9,6] = H_[6,9]
                H_[6,10] = -TWOT6*pbm.elpar[iel_][10]
                H_[10,6] = H_[6,10]
                H_[6,11] = -TWOT6*pbm.elpar[iel_][11]
                H_[11,6] = H_[6,11]
                H_[6,12] = -TWOT6*pbm.elpar[iel_][12]
                H_[12,6] = H_[6,12]
                H_[7,7] = -TWOT7*pbm.elpar[iel_][7]
                H_[7,8] = -TWOT7*pbm.elpar[iel_][8]
                H_[8,7] = H_[7,8]
                H_[7,9] = -TWOT7*pbm.elpar[iel_][8]
                H_[9,7] = H_[7,9]
                H_[7,10] = -TWOT7*pbm.elpar[iel_][10]
                H_[10,7] = H_[7,10]
                H_[7,11] = -TWOT7*pbm.elpar[iel_][11]
                H_[11,7] = H_[7,11]
                H_[7,12] = -TWOT7*pbm.elpar[iel_][12]
                H_[12,7] = H_[7,12]
                H_[8,8] = -TWOT8*pbm.elpar[iel_][8]
                H_[8,9] = -TWOT8*pbm.elpar[iel_][8]
                H_[9,8] = H_[8,9]
                H_[8,10] = -TWOT8*pbm.elpar[iel_][10]
                H_[10,8] = H_[8,10]
                H_[8,11] = -TWOT8*pbm.elpar[iel_][11]
                H_[11,8] = H_[8,11]
                H_[8,12] = -TWOT8*pbm.elpar[iel_][12]
                H_[12,8] = H_[8,12]
                H_[9,9] = -TWOT9*pbm.elpar[iel_][9]
                H_[9,10] = -TWOT9*pbm.elpar[iel_][10]
                H_[10,9] = H_[9,10]
                H_[9,11] = -TWOT9*pbm.elpar[iel_][11]
                H_[11,9] = H_[9,11]
                H_[9,12] = -TWOT9*pbm.elpar[iel_][12]
                H_[12,9] = H_[9,12]
                H_[10,10] = -TWOT10*pbm.elpar[iel_][10]
                H_[10,11] = -TWOT10*pbm.elpar[iel_][11]
                H_[11,10] = H_[10,11]
                H_[10,12] = -TWOT10*pbm.elpar[iel_][12]
                H_[12,10] = H_[10,12]
                H_[11,11] = -TWOT11*pbm.elpar[iel_][11]
                H_[11,12] = -TWOT11*pbm.elpar[iel_][12]
                H_[12,11] = H_[11,12]
                H_[12,12] = -TWOT12*pbm.elpar[iel_][12]
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

