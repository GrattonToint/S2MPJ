function HS93(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS93
#    *********
# 
#    A transformer design problem.
# 
#    Source: problem 93 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: Nick Gould, August 1991.
# 
#    classification = "C-OOR2-MY-6-2"
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "HS93"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["N"] = 6
        v_["1"] = 1
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
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        ig,ig_,_ = s2mpj_ii("C1",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"C1")
        ig,ig_,_ = s2mpj_ii("C2",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"C2")
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
        pbm.gconst[ig_["C1"]] = Float64(2.07e+0)
        pbm.gconst[ig_["C2"]] = Float64(1.0e+0)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        if haskey(ix_,"X1")
            pb.x0[ix_["X1"]] = Float64(5.54)
        else
            pb.y0[findfirst(x->x==ig_["X1"],pbm.congrps)] = Float64(5.54)
        end
        if haskey(ix_,"X2")
            pb.x0[ix_["X2"]] = Float64(4.4)
        else
            pb.y0[findfirst(x->x==ig_["X2"],pbm.congrps)] = Float64(4.4)
        end
        if haskey(ix_,"X3")
            pb.x0[ix_["X3"]] = Float64(12.02)
        else
            pb.y0[findfirst(x->x==ig_["X3"],pbm.congrps)] = Float64(12.02)
        end
        if haskey(ix_,"X4")
            pb.x0[ix_["X4"]] = Float64(11.82)
        else
            pb.y0[findfirst(x->x==ig_["X4"],pbm.congrps)] = Float64(11.82)
        end
        if haskey(ix_,"X5")
            pb.x0[ix_["X5"]] = Float64(0.702)
        else
            pb.y0[findfirst(x->x==ig_["X5"],pbm.congrps)] = Float64(0.702)
        end
        if haskey(ix_,"X6")
            pb.x0[ix_["X6"]] = Float64(0.852)
        else
            pb.y0[findfirst(x->x==ig_["X6"],pbm.congrps)] = Float64(0.852)
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eOE1", iet_)
        loaset(elftv,it,1,"X1")
        loaset(elftv,it,2,"X2")
        loaset(elftv,it,3,"X3")
        loaset(elftv,it,4,"X4")
        it,iet_,_ = s2mpj_ii( "eOE2", iet_)
        loaset(elftv,it,1,"X1")
        loaset(elftv,it,2,"X2")
        loaset(elftv,it,3,"X3")
        loaset(elftv,it,4,"X4")
        it,iet_,_ = s2mpj_ii( "eOE3", iet_)
        loaset(elftv,it,1,"X1")
        loaset(elftv,it,2,"X2")
        loaset(elftv,it,3,"X3")
        loaset(elftv,it,4,"X4")
        loaset(elftv,it,5,"X5")
        it,iet_,_ = s2mpj_ii( "eOE4", iet_)
        loaset(elftv,it,1,"X1")
        loaset(elftv,it,2,"X2")
        loaset(elftv,it,3,"X3")
        loaset(elftv,it,4,"X4")
        loaset(elftv,it,5,"X6")
        it,iet_,_ = s2mpj_ii( "eC1E1", iet_)
        loaset(elftv,it,1,"X1")
        loaset(elftv,it,2,"X2")
        loaset(elftv,it,3,"X3")
        loaset(elftv,it,4,"X4")
        loaset(elftv,it,5,"X5")
        loaset(elftv,it,6,"X6")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "OE1"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eOE1")
        arrset(ielftype,ie,iet_["eOE1"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X4",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "OE2"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eOE2")
        arrset(ielftype,ie,iet_["eOE2"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X4",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "OE3"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eOE3")
        arrset(ielftype,ie,iet_["eOE3"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X4",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X5",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "OE4"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eOE4")
        arrset(ielftype,ie,iet_["eOE4"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X4",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X6"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X6",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "C1E1"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eC1E1")
        arrset(ielftype,ie,iet_["eC1E1"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X4",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X5",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X6"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X6",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["OBJ"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["OE1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(2.04e-2))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["OE2"])
        loaset(pbm.grelw,ig,posel,Float64(1.87e-2))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["OE3"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(6.07e-2))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["OE4"])
        loaset(pbm.grelw,ig,posel,Float64(4.37e-2))
        ig = ig_["C1"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["C1E1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(1.0e-3))
        ig = ig_["C2"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["OE3"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(6.2e-4))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["OE4"])
        loaset(pbm.grelw,ig,posel,Float64(5.8e-4))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               135.075961
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
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
        pb.pbclass = "C-OOR2-MY-6-2"
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

    elseif action == "eOE1"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,3,4)
        IV_ =  zeros(Float64,3)
        U_[1,1] = U_[1,1]+1
        U_[2,4] = U_[2,4]+1
        U_[3,1] = U_[3,1]+1
        U_[3,2] = U_[3,2]+1
        U_[3,3] = U_[3,3]+1
        IV_[1] = dot(U_[1,:],EV_)
        IV_[2] = dot(U_[2,:],EV_)
        IV_[3] = dot(U_[3,:],EV_)
        f_   = IV_[1]*IV_[2]*IV_[3]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = IV_[2]*IV_[3]
            g_[2] = IV_[1]*IV_[3]
            g_[3] = IV_[1]*IV_[2]
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,2] = IV_[3]
                H_[2,1] = H_[1,2]
                H_[1,3] = IV_[2]
                H_[3,1] = H_[1,3]
                H_[2,3] = IV_[1]
                H_[3,2] = H_[2,3]
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

    elseif action == "eOE2"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,3,4)
        IV_ =  zeros(Float64,3)
        U_[1,2] = U_[1,2]+1
        U_[2,3] = U_[2,3]+1
        U_[3,1] = U_[3,1]+1
        U_[3,2] = U_[3,2]+1.570000e+00
        U_[3,4] = U_[3,4]+1
        IV_[1] = dot(U_[1,:],EV_)
        IV_[2] = dot(U_[2,:],EV_)
        IV_[3] = dot(U_[3,:],EV_)
        f_   = IV_[1]*IV_[2]*IV_[3]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = IV_[2]*IV_[3]
            g_[2] = IV_[1]*IV_[3]
            g_[3] = IV_[1]*IV_[2]
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,2] = IV_[3]
                H_[2,1] = H_[1,2]
                H_[1,3] = IV_[2]
                H_[3,1] = H_[1,3]
                H_[2,3] = IV_[1]
                H_[3,2] = H_[2,3]
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

    elseif action == "eOE3"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,4,5)
        IV_ =  zeros(Float64,4)
        U_[1,1] = U_[1,1]+1
        U_[2,4] = U_[2,4]+1
        U_[3,5] = U_[3,5]+1
        U_[4,1] = U_[4,1]+1
        U_[4,2] = U_[4,2]+1
        U_[4,3] = U_[4,3]+1
        IV_[1] = dot(U_[1,:],EV_)
        IV_[2] = dot(U_[2,:],EV_)
        IV_[3] = dot(U_[3,:],EV_)
        IV_[4] = dot(U_[4,:],EV_)
        f_   = IV_[1]*IV_[2]*(IV_[3]^2)*IV_[4]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = IV_[2]*(IV_[3]^2)*IV_[4]
            g_[2] = IV_[1]*(IV_[3]^2)*IV_[4]
            g_[3] = IV_[1]*IV_[2]*2.0e+0*IV_[3]*IV_[4]
            g_[4] = IV_[1]*IV_[2]*(IV_[3]^2)
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,4,4)
                H_[1,2] = (IV_[3]^2)*IV_[4]
                H_[2,1] = H_[1,2]
                H_[1,3] = IV_[2]*2.0e+0*IV_[3]*IV_[4]
                H_[3,1] = H_[1,3]
                H_[1,4] = IV_[2]*(IV_[3]^2)
                H_[4,1] = H_[1,4]
                H_[2,3] = IV_[1]*2.0e+0*IV_[3]*IV_[4]
                H_[3,2] = H_[2,3]
                H_[2,4] = IV_[1]*(IV_[3]^2)
                H_[4,2] = H_[2,4]
                H_[3,3] = IV_[1]*IV_[2]*2.0e+0*IV_[4]
                H_[3,4] = IV_[1]*IV_[2]*2.0e+0*IV_[3]
                H_[4,3] = H_[3,4]
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

    elseif action == "eOE4"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,4,5)
        IV_ =  zeros(Float64,4)
        U_[1,2] = U_[1,2]+1
        U_[2,3] = U_[2,3]+1
        U_[3,5] = U_[3,5]+1
        U_[4,1] = U_[4,1]+1
        U_[4,2] = U_[4,2]+1.570000e+00
        U_[4,4] = U_[4,4]+1
        IV_[1] = dot(U_[1,:],EV_)
        IV_[2] = dot(U_[2,:],EV_)
        IV_[3] = dot(U_[3,:],EV_)
        IV_[4] = dot(U_[4,:],EV_)
        f_   = IV_[1]*IV_[2]*(IV_[3]^2)*IV_[4]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = IV_[2]*(IV_[3]^2)*IV_[4]
            g_[2] = IV_[1]*(IV_[3]^2)*IV_[4]
            g_[3] = IV_[1]*IV_[2]*2.0e+0*IV_[3]*IV_[4]
            g_[4] = IV_[1]*IV_[2]*(IV_[3]^2)
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,4,4)
                H_[1,2] = (IV_[3]^2)*IV_[4]
                H_[2,1] = H_[1,2]
                H_[1,3] = IV_[2]*2.0e+0*IV_[3]*IV_[4]
                H_[3,1] = H_[1,3]
                H_[1,4] = IV_[2]*(IV_[3]^2)
                H_[4,1] = H_[1,4]
                H_[2,3] = IV_[1]*2.0e+0*IV_[3]*IV_[4]
                H_[3,2] = H_[2,3]
                H_[2,4] = IV_[1]*(IV_[3]^2)
                H_[4,2] = H_[2,4]
                H_[3,3] = IV_[1]*IV_[2]*2.0e+0*IV_[4]
                H_[3,4] = IV_[1]*IV_[2]*2.0e+0*IV_[3]
                H_[4,3] = H_[3,4]
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

    elseif action == "eC1E1"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]*EV_[2]*EV_[3]*EV_[4]*EV_[5]*EV_[6]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[2]*EV_[3]*EV_[4]*EV_[5]*EV_[6]
            g_[2] = EV_[1]*EV_[3]*EV_[4]*EV_[5]*EV_[6]
            g_[3] = EV_[1]*EV_[2]*EV_[4]*EV_[5]*EV_[6]
            g_[4] = EV_[1]*EV_[2]*EV_[3]*EV_[5]*EV_[6]
            g_[5] = EV_[1]*EV_[2]*EV_[3]*EV_[4]*EV_[6]
            g_[6] = EV_[1]*EV_[2]*EV_[3]*EV_[4]*EV_[5]
            if nargout>2
                H_ = zeros(Float64,6,6)
                H_[1,2] = EV_[3]*EV_[4]*EV_[5]*EV_[6]
                H_[2,1] = H_[1,2]
                H_[1,3] = EV_[2]*EV_[4]*EV_[5]*EV_[6]
                H_[3,1] = H_[1,3]
                H_[1,4] = EV_[2]*EV_[3]*EV_[5]*EV_[6]
                H_[4,1] = H_[1,4]
                H_[1,5] = EV_[2]*EV_[3]*EV_[4]*EV_[6]
                H_[5,1] = H_[1,5]
                H_[1,6] = EV_[2]*EV_[3]*EV_[4]*EV_[5]
                H_[6,1] = H_[1,6]
                H_[2,3] = EV_[1]*EV_[4]*EV_[5]*EV_[6]
                H_[3,2] = H_[2,3]
                H_[2,4] = EV_[1]*EV_[3]*EV_[5]*EV_[6]
                H_[4,2] = H_[2,4]
                H_[2,5] = EV_[1]*EV_[3]*EV_[4]*EV_[6]
                H_[5,2] = H_[2,5]
                H_[2,6] = EV_[1]*EV_[3]*EV_[4]*EV_[5]
                H_[6,2] = H_[2,6]
                H_[3,4] = EV_[1]*EV_[2]*EV_[5]*EV_[6]
                H_[4,3] = H_[3,4]
                H_[3,5] = EV_[1]*EV_[2]*EV_[4]*EV_[6]
                H_[5,3] = H_[3,5]
                H_[3,6] = EV_[1]*EV_[2]*EV_[4]*EV_[5]
                H_[6,3] = H_[3,6]
                H_[4,5] = EV_[1]*EV_[2]*EV_[3]*EV_[6]
                H_[5,4] = H_[4,5]
                H_[4,6] = EV_[1]*EV_[2]*EV_[3]*EV_[5]
                H_[6,4] = H_[4,6]
                H_[5,6] = EV_[1]*EV_[2]*EV_[3]*EV_[4]
                H_[6,5] = H_[5,6]
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

