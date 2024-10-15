function CRESC4(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : CRESC4
#    *********
# 
#    This problem consists in finding the crescent of smallest area containing
#    a set of points given in the plane.   This problem arises as a subproblem
#    in pattern recognition and has been suggested by J.P. Rasson.  It
#    originates in the detection of "salt domes" (with the potential presence of
#    oil pockets!) from geological data.
# 
#    The present problem is a simplified version where the crescent is entirely
#    determined by the only four data points.
# 
#    The problem is not convex.
# 
#    A crescent is defined as follows.  Assume one has two circles of center
#    C1 and C2 and of radii r1 and r2 respectively. Assume furthermore that
#    r1 >= r2 and that C2 is within C1.  Assume finally that the distance from
#    C1 to C2 is >= r1 - r2.  Then the crescent is the part of the plane
#    contained in circle 2 but not in circle 1.
# 
#    In order to preserve feasibility at all stages (ensuring that the
#    crescent exists and that its area can be computed), the following
#    parametrization is used:
# 
#    ( C2x, C2y ) = ( C1x, C1y ) + a * d * ( cos(t), sin(t) )
# 
#    r1 = a * d + r 
# 
#    r2 = ( a + 1 ) * d + r
# 
#    with the bounds
# 
#    a >= 1, 0 <= t <= 2 * pi, r2 >= 0 , 0 <= d <= 1.
# 
#    SIF input: Ph. Toint, June 1993.
# 
#    classification = "C-OOR2-MY-6-8"
# 
#    number of points to be included in the crescent.
#    the number of constraints is 2*NP
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "CRESC4"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["NP"] = 4
        v_["X1"] = 1.0
        v_["Y1"] = 0.0
        v_["X2"] = 0.0
        v_["Y2"] = 1.0
        v_["X3"] = 0.0
        v_["Y3"] = -1.0
        v_["X4"] = 0.5
        v_["Y4"] = 0.0
        v_["1"] = 1
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        iv,ix_,_ = s2mpj_ii("V1",ix_)
        arrset(pb.xnames,iv,"V1")
        iv,ix_,_ = s2mpj_ii("W1",ix_)
        arrset(pb.xnames,iv,"W1")
        iv,ix_,_ = s2mpj_ii("D",ix_)
        arrset(pb.xnames,iv,"D")
        iv,ix_,_ = s2mpj_ii("A",ix_)
        arrset(pb.xnames,iv,"A")
        iv,ix_,_ = s2mpj_ii("T",ix_)
        arrset(pb.xnames,iv,"T")
        iv,ix_,_ = s2mpj_ii("R",ix_)
        arrset(pb.xnames,iv,"R")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        for I = Int64(v_["1"]):Int64(v_["NP"])
            ig,ig_,_ = s2mpj_ii("IS2"*string(I),ig_)
            arrset(gtype,ig,"<=")
            arrset(pb.cnames,ig,"IS2"*string(I))
            ig,ig_,_ = s2mpj_ii("OS1"*string(I),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"OS1"*string(I))
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
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower[ix_["V1"]] = -Inf
        pb.xupper[ix_["V1"]] = +Inf
        pb.xlower[ix_["W1"]] = -Inf
        pb.xupper[ix_["W1"]] = +Inf
        pb.xlower[ix_["R"]] = 0.39
        pb.xlower[ix_["A"]] = 1.0
        pb.xlower[ix_["T"]] = 0.0
        pb.xupper[ix_["T"]] = 6.2831852
        pb.xlower[ix_["D"]] = 1.0e-8
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        pb.x0[ix_["V1"]] = Float64(-40.0)
        pb.x0[ix_["W1"]] = Float64(5.0)
        pb.x0[ix_["R"]] = Float64(0.75)
        pb.x0[ix_["A"]] = Float64(2.0)
        pb.x0[ix_["T"]] = Float64(1.5)
        pb.x0[ix_["D"]] = Float64(1.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eSQR1", iet_)
        loaset(elftv,it,1,"A")
        loaset(elftv,it,2,"R")
        loaset(elftv,it,3,"D")
        it,iet_,_ = s2mpj_ii( "eSQR2", iet_)
        loaset(elftv,it,1,"D")
        loaset(elftv,it,2,"R")
        it,iet_,_ = s2mpj_ii( "eSC", iet_)
        loaset(elftv,it,1,"AZ")
        loaset(elftv,it,2,"BZ")
        loaset(elftv,it,3,"DZ")
        it,iet_,_ = s2mpj_ii( "eDIST", iet_)
        loaset(elftv,it,1,"X")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"P")
        it,iet_,_ = s2mpj_ii( "eDISTX", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"T")
        loaset(elftv,it,3,"A")
        loaset(elftv,it,4,"D")
        loaset(elftp,it,1,"P")
        it,iet_,_ = s2mpj_ii( "eDISTY", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"T")
        loaset(elftv,it,3,"A")
        loaset(elftv,it,4,"D")
        loaset(elftp,it,1,"P")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "OB"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSC")
        arrset(ielftype,ie,iet_["eSC"])
        vname = "A"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="AZ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "R"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="BZ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "D"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="DZ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "R2SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR2")
        arrset(ielftype,ie,iet_["eSQR2"])
        vname = "D"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="D",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "R"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="R",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "R1SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR1")
        arrset(ielftype,ie,iet_["eSQR1"])
        vname = "D"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="D",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "A"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="A",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "R"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="R",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        for I = Int64(v_["1"]):Int64(v_["NP"])
            ename = "XV1"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eDIST")
            arrset(ielftype,ie,iet_["eDIST"])
            vname = "V1"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="P",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["X"*string(I)]))
            ename = "XV2"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eDISTX")
            arrset(ielftype,ie,iet_["eDISTX"])
            vname = "V1"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "A"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="A",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "D"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="D",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "T"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="T",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="P",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["X"*string(I)]))
            ename = "YW1"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eDIST")
            arrset(ielftype,ie,iet_["eDIST"])
            vname = "W1"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="P",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["Y"*string(I)]))
            ename = "YW2"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eDISTY")
            arrset(ielftype,ie,iet_["eDISTY"])
            vname = "W1"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "A"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="A",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "D"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="D",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "T"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="T",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="P",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["Y"*string(I)]))
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["OBJ"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["OB"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        for I = Int64(v_["1"]):Int64(v_["NP"])
            ig = ig_["IS2"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["XV2"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["YW2"*string(I)])
            loaset(pbm.grelw,ig,posel, 1.)
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["R2SQ"])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(-1.0))
            ig = ig_["OS1"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["XV1"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["YW1"*string(I)])
            loaset(pbm.grelw,ig,posel, 1.)
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["R1SQ"])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(-1.0))
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution             0.87189692
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.cupper[1:pb.nle] = zeros(Float64,pb.nle)
        pb.clower[pb.nle+pb.neq+1:pb.m] = zeros(Float64,pb.nge)
        pb.cupper[1:pb.nge] = fill(Inf,pb.nge)
        pbm.A = spzeros(Float64,0,0)
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-OOR2-MY-6-8"
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

    elseif action == "eSQR1"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        Q = EV_[1]*EV_[3]+EV_[2]
        f_   = Q*Q
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 2.0*Q*EV_[3]
            g_[3] = 2.0*Q*EV_[1]
            g_[2] = 2.0*Q
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,1] = 2.0*EV_[3]*EV_[3]
                H_[1,3] = 2.0*(EV_[1]*EV_[3]+Q)
                H_[3,1] = H_[1,3]
                H_[1,2] = 2.0*EV_[3]
                H_[2,1] = H_[1,2]
                H_[3,3] = 2.0*EV_[1]*EV_[1]
                H_[3,2] = 2.0*EV_[1]
                H_[2,3] = H_[3,2]
                H_[2,2] = 2.0
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eSQR2"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,1,2)
        IV_ =  zeros(Float64,1)
        U_[1,1] = U_[1,1]+1
        U_[1,2] = U_[1,2]+1
        IV_[1] = dot(U_[1,:],EV_)
        f_   = IV_[1]*IV_[1]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = IV_[1]+IV_[1]
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 2.0
                H_ = U_'*H_*U_
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eSC"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        A = EV_[3]+EV_[2]
        AB = 1.0
        AD = 1.0
        B = EV_[1]*EV_[3]+EV_[2]
        BA = EV_[3]
        BB = 1.0
        BD = EV_[1]
        BAD = 1.0
        D = EV_[1]*EV_[3]
        DA = EV_[3]
        DD = EV_[1]
        DAD = 1.0
        E = 2.0*A*D
        EA = 2.0*A*DA
        EB = 2.0*AB*D
        ED = 2.0*(A*DD+AD*D)
        EAB = 2.0*AB*DA
        EAD = 2.0*(AD*DA+A*DAD)
        EBD = 2.0*AB*DD
        EDD = 2.0*(AD*DD+AD*DD)
        P = 2.0*B*D
        PA = 2.0*(B*DA+BA*D)
        loc_PB = 2.0*BB*D
        PD = 2.0*(B*DD+BD*D)
        PAA = 2.0*(BA*DA+BA*DA)
        PAB = 2.0*BB*DA
        PAD = 2.0*(BD*DA+B*DAD+BAD*D+BA*DD)
        PBD = 2.0*BB*DD
        PDD = 2.0*(BD*DD+BD*DD)
        F = D*D+B*B-A*A
        FA = 2.0*(D*DA+B*BA)
        FB = 2.0*(B*BB-A*AB)
        FD = 2.0*(D*DD+B*BD-A*AD)
        FAA = 2.0*(DA*DA+BA*BA)
        FAB = 2.0*BB*BA
        FAD = 2.0*(DD*DA+D*DAD+BD*BA+B*BAD)
        FBB = 2.0*(BB*BB-AB*AB)
        FBD = 2.0*(BD*BB-AD*AB)
        FDD = 2.0*(DD*DD+BD*BD-AD*AD)
        G = D*D-B*B+A*A
        GA = 2.0*(D*DA-B*BA)
        GB = 2.0*(-B*BB+A*AB)
        GD = 2.0*(D*DD-B*BD+A*AD)
        GAA = 2.0*(DA*DA-BA*BA)
        GAB = -2.0*BB*BA
        GAD = 2.0*(DD*DA+D*DAD-BD*BA-B*BAD)
        GBB = 2.0*(-BB*BB+AB*AB)
        GBD = 2.0*(-BD*BB+AD*AB)
        GDD = 2.0*(DD*DD-BD*BD+AD*AD)
        H = F/P
        I = FA-H*PA
        J = FB-H*loc_PB
        K = FD-H*PD
        HA = I/P
        HB = J/P
        HD = K/P
        IA = FAA-HA*PA-H*PAA
        IB = FAB-HB*PA-H*PAB
        ID = FAD-HD*PA-H*PAD
        JB = FBB-HB*loc_PB
        JD = FBD-HD*loc_PB-H*PBD
        KD = FDD-HD*PD-H*PDD
        HAA = (IA-HA*PA)/P
        HAB = (IB-HA*loc_PB)/P
        HAD = (ID-HA*PD)/P
        HBB = (JB-HB*loc_PB)/P
        HBD = (JD-HB*PD)/P
        HDD = (KD-HD*PD)/P
        L = -G/E
        M = -GA-L*EA
        N = -GB-L*EB
        O = -GD-L*ED
        LA = M/E
        LB = N/E
        LD = O/E
        MA = -GAA-LA*EA
        MB = -GAB-LB*EA-L*EAB
        MD = -GAD-LD*EA-L*EAD
        NB = -GBB-LB*EB
        ND = -GBD-LD*EB-L*EBD
        OD = -GDD-LD*ED-L*EDD
        LAA = (MA-LA*EA)/E
        LAB = (MB-LA*EB)/E
        LAD = (MD-LA*ED)/E
        LBB = (NB-LB*EB)/E
        LBD = (ND-LB*ED)/E
        LDD = (OD-LD*ED)/E
        C = acos(H)
        CH = -1.0/sqrt(1.0-H*H)
        CHH = CH*H/(1.0-H*H)
        CA = CH*HA
        CB = CH*HB
        CD = CH*HD
        CAA = CHH*HA*HA+CH*HAA
        CAB = CHH*HA*HB+CH*HAB
        CAD = CHH*HA*HD+CH*HAD
        CBB = CHH*HB*HB+CH*HBB
        CBD = CHH*HB*HD+CH*HBD
        CDD = CHH*HD*HD+CH*HDD
        Q = acos(L)
        QL = -1.0/sqrt(1.0-L*L)
        QLL = QL*L/(1.0-L*L)
        QA = QL*LA
        QB = QL*LB
        QD = QL*LD
        QAA = QLL*LA*LA+QL*LAA
        QAB = QLL*LA*LB+QL*LAB
        QAD = QLL*LA*LD+QL*LAD
        QBB = QLL*LB*LB+QL*LBB
        QBD = QLL*LB*LD+QL*LBD
        QDD = QLL*LD*LD+QL*LDD
        R = B*B*C
        RA = 2.0*B*BA*C+B*B*CA
        RB = 2.0*B*BB*C+B*B*CB
        RD = 2.0*B*BD*C+B*B*CD
        RAA = 2.0*(BA*BA*C+B*BA*CA+B*BA*CA)+B*B*CAA
        RAB = 2.0*(BB*BA*C+B*BA*CB+B*BB*CA)+B*B*CAB
        RAD = 2.0*(BD*BA*C+B*BAD*C+B*BA*CD+B*BD*CA)+B*B*CAD
        RBB = 2.0*(BB*BB*C+B*BB*CB+B*BB*CB)+B*B*CBB
        RBD = 2.0*(BD*BB*C+B*BB*CD+B*BD*CB)+B*B*CBD
        RDD = 2.0*(BD*BD*C+B*BD*CD+B*BD*CD)+B*B*CDD
        S = A*A*Q
        SA = A*A*QA
        SB = 2.0*A*AB*Q+A*A*QB
        SD = 2.0*A*AD*Q+A*A*QD
        SAA = A*A*QAA
        SAB = 2.0*A*AB*QA+A*A*QAB
        SAD = 2.0*A*AD*QA+A*A*QAD
        SBB = 2.0*(AB*AB*Q+A*AB*QB+A*AB*QB)+A*A*QBB
        SBD = 2.0*(AD*AB*Q+A*AB*QD+A*AD*QB)+A*A*QBD
        SDD = 2.0*(AD*AD*Q+A*AD*QD+A*AD*QD)+A*A*QDD
        SQ = sin(Q)
        CQ = L
        W = 0.5*E*SQ
        WA = 0.5*(EA*SQ+E*CQ*QA)
        WB = 0.5*(EB*SQ+E*CQ*QB)
        WD = 0.5*(ED*SQ+E*CQ*QD)
        WAA = 0.5*(EA*CQ*QA+EA*CQ*QA-E*SQ*QA*QA+E*CQ*QAA)
        WAB = 0.5*(EAB*SQ+EA*CQ*QB+EB*CQ*QA-E*SQ*QB*QA+E*CQ*QAB)
        WAD = 0.5*(EAD*SQ+EA*CQ*QD+ED*CQ*QA-E*SQ*QD*QA+E*CQ*QAD)
        WBB = 0.5*(EB*CQ*QB+EB*CQ*QB-E*SQ*QB*QB+E*CQ*QBB)
        WBD = 0.5*(EBD*SQ+EB*CQ*QD+ED*CQ*QB-E*SQ*QD*QB+E*CQ*QBD)
        WDD = 0.5*(EDD*SQ+ED*CQ*QD+ED*CQ*QD-E*SQ*QD*QD+E*CQ*QDD)
        V = S-R+W
        VA = SA-RA+WA
        VB = SB-RB+WB
        VD = SD-RD+WD
        VAA = SAA-RAA+WAA
        VAB = SAB-RAB+WAB
        VAD = SAD-RAD+WAD
        VBB = SBB-RBB+WBB
        VBD = SBD-RBD+WBD
        VDD = SDD-RDD+WDD
        f_   = V
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = VA
            g_[2] = VB
            g_[3] = VD
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,1] = VAA
                H_[1,2] = VAB
                H_[2,1] = H_[1,2]
                H_[1,3] = VAD
                H_[3,1] = H_[1,3]
                H_[2,2] = VBB
                H_[2,3] = VBD
                H_[3,2] = H_[2,3]
                H_[3,3] = VDD
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eDIST"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = (EV_[1]-pbm.elpar[iel_][1])^2
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 2.0*(EV_[1]-pbm.elpar[iel_][1])
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

    elseif action == "eDISTY"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        ST = sin(EV_[2])
        CT = cos(EV_[2])
        B = EV_[1]+EV_[3]*EV_[4]*ST-pbm.elpar[iel_][1]
        BA = EV_[4]*ST
        BD = EV_[3]*ST
        BT = EV_[3]*EV_[4]*CT
        f_   = B*B
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 2.0*B
            g_[3] = 2.0*B*BA
            g_[4] = 2.0*B*BD
            g_[2] = 2.0*B*BT
            if nargout>2
                H_ = zeros(Float64,4,4)
                H_[1,1] = 2.0
                H_[1,3] = 2.0*BA
                H_[3,1] = H_[1,3]
                H_[1,4] = 2.0*BD
                H_[4,1] = H_[1,4]
                H_[1,2] = 2.0*BT
                H_[2,1] = H_[1,2]
                H_[3,3] = 2.0*BA*BA
                H_[3,4] = 2.0*(BD*BA+B*ST)
                H_[4,3] = H_[3,4]
                H_[3,2] = 2.0*(BT*BA+B*EV_[4]*CT)
                H_[2,3] = H_[3,2]
                H_[4,4] = 2.0*BD*BD
                H_[4,2] = 2.0*(BT*BD+B*EV_[3]*CT)
                H_[2,4] = H_[4,2]
                H_[2,2] = 2.0*(BT*BT-B*EV_[3]*EV_[4]*ST)
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eDISTX"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        ST = sin(EV_[2])
        CT = cos(EV_[2])
        B = EV_[1]+EV_[3]*EV_[4]*CT-pbm.elpar[iel_][1]
        BA = EV_[4]*CT
        BD = EV_[3]*CT
        BT = -EV_[3]*EV_[4]*ST
        f_   = B*B
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 2.0*B
            g_[3] = 2.0*B*BA
            g_[4] = 2.0*B*BD
            g_[2] = 2.0*B*BT
            if nargout>2
                H_ = zeros(Float64,4,4)
                H_[1,1] = 2.0
                H_[1,3] = 2.0*BA
                H_[3,1] = H_[1,3]
                H_[1,4] = 2.0*BD
                H_[4,1] = H_[1,4]
                H_[1,2] = 2.0*BT
                H_[2,1] = H_[1,2]
                H_[3,3] = 2.0*BA*BA
                H_[3,4] = 2.0*(BD*BA+B*CT)
                H_[4,3] = H_[3,4]
                H_[3,2] = 2.0*(BT*BA-B*EV_[4]*ST)
                H_[2,3] = H_[3,2]
                H_[4,4] = 2.0*BD*BD
                H_[4,2] = 2.0*(BT*BD-B*EV_[3]*ST)
                H_[2,4] = H_[4,2]
                H_[2,2] = 2.0*(BT*BT-B*EV_[3]*EV_[4]*CT)
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

