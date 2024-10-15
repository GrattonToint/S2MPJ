function HS99(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS99
#    *********
# 
#    Source: problem 99 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: Ph. Toint, April 1991.
# 
#    classification = "C-OOR2-AN-7-2"
# 
#    Constants
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "HS99"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["1"] = 1
        v_["2"] = 2
        v_["7"] = 7
        v_["8"] = 8
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["7"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        arrset(pbm.gscale,ig,Float64(-1.0))
        ig,ig_,_ = s2mpj_ii("Q8E",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"Q8E")
        ig,ig_,_ = s2mpj_ii("S8E",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"S8E")
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
        pbm.gconst[ig_["Q8E"]] = Float64(100000.0)
        pbm.gconst[ig_["S8E"]] = Float64(1000.0)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xupper = fill(1.58,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(0.5),pb.n)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eR8T", iet_)
        loaset(elftv,it,1,"X1")
        loaset(elftv,it,2,"X2")
        loaset(elftv,it,3,"X3")
        loaset(elftv,it,4,"X4")
        loaset(elftv,it,5,"X5")
        loaset(elftv,it,6,"X6")
        loaset(elftv,it,7,"X7")
        it,iet_,_ = s2mpj_ii( "eQ8T", iet_)
        loaset(elftv,it,1,"X1")
        loaset(elftv,it,2,"X2")
        loaset(elftv,it,3,"X3")
        loaset(elftv,it,4,"X4")
        loaset(elftv,it,5,"X5")
        loaset(elftv,it,6,"X6")
        loaset(elftv,it,7,"X7")
        it,iet_,_ = s2mpj_ii( "eS8T", iet_)
        loaset(elftv,it,1,"X1")
        loaset(elftv,it,2,"X2")
        loaset(elftv,it,3,"X3")
        loaset(elftv,it,4,"X4")
        loaset(elftv,it,5,"X5")
        loaset(elftv,it,6,"X6")
        loaset(elftv,it,7,"X7")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "R8"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eR8T")
        arrset(ielftype,ie,iet_["eR8T"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.58),Float64(0.5))
        posev = findfirst(x->x=="X1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.58),Float64(0.5))
        posev = findfirst(x->x=="X2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.58),Float64(0.5))
        posev = findfirst(x->x=="X3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.58),Float64(0.5))
        posev = findfirst(x->x=="X4",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.58),Float64(0.5))
        posev = findfirst(x->x=="X5",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X6"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.58),Float64(0.5))
        posev = findfirst(x->x=="X6",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X7"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.58),Float64(0.5))
        posev = findfirst(x->x=="X7",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "Q8"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eQ8T")
        arrset(ielftype,ie,iet_["eQ8T"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.58),Float64(0.5))
        posev = findfirst(x->x=="X1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.58),Float64(0.5))
        posev = findfirst(x->x=="X2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.58),Float64(0.5))
        posev = findfirst(x->x=="X3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.58),Float64(0.5))
        posev = findfirst(x->x=="X4",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.58),Float64(0.5))
        posev = findfirst(x->x=="X5",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X6"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.58),Float64(0.5))
        posev = findfirst(x->x=="X6",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X7"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.58),Float64(0.5))
        posev = findfirst(x->x=="X7",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "S8"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eS8T")
        arrset(ielftype,ie,iet_["eS8T"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.58),Float64(0.5))
        posev = findfirst(x->x=="X1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.58),Float64(0.5))
        posev = findfirst(x->x=="X2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.58),Float64(0.5))
        posev = findfirst(x->x=="X3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.58),Float64(0.5))
        posev = findfirst(x->x=="X4",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.58),Float64(0.5))
        posev = findfirst(x->x=="X5",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X6"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.58),Float64(0.5))
        posev = findfirst(x->x=="X6",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X7"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.58),Float64(0.5))
        posev = findfirst(x->x=="X7",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gL2",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["OBJ"]
        arrset(pbm.grftype,ig,"gL2")
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["R8"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["Q8E"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["Q8"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["S8E"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["S8"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               -831079892.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pbm.A = spzeros(Float64,0,0)
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-OOR2-AN-7-2"
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

    elseif action == "e_globs"

        pbm = args[1]
        arrset(pbm.efpar,1,50.0)
        arrset(pbm.efpar,2,50.0)
        arrset(pbm.efpar,3,75.0)
        arrset(pbm.efpar,4,75.0)
        arrset(pbm.efpar,5,75.0)
        arrset(pbm.efpar,6,100.0)
        arrset(pbm.efpar,7,100.0)
        arrset(pbm.efpar,8,25.0)
        arrset(pbm.efpar,9,25.0)
        arrset(pbm.efpar,10,50.0)
        arrset(pbm.efpar,11,50.0)
        arrset(pbm.efpar,12,50.0)
        arrset(pbm.efpar,13,90.0)
        arrset(pbm.efpar,14,90.0)
        arrset(pbm.efpar,15,32.0)
        return pbm

    elseif action == "eR8T"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        R2 = pbm.efpar[1]*pbm.efpar[8]*cos(EV_[1])
        R3 = pbm.efpar[2]*pbm.efpar[9]*cos(EV_[2])+R2
        R4 = pbm.efpar[3]*pbm.efpar[10]*cos(EV_[3])+R3
        R5 = pbm.efpar[4]*pbm.efpar[11]*cos(EV_[4])+R4
        R6 = pbm.efpar[5]*pbm.efpar[12]*cos(EV_[5])+R5
        R7 = pbm.efpar[6]*pbm.efpar[13]*cos(EV_[6])+R6
        f_   = pbm.efpar[7]*pbm.efpar[14]*cos(EV_[7])+R7
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = -pbm.efpar[1]*pbm.efpar[8]*sin(EV_[1])
            g_[2] = -pbm.efpar[2]*pbm.efpar[9]*sin(EV_[2])
            g_[3] = -pbm.efpar[3]*pbm.efpar[10]*sin(EV_[3])
            g_[4] = -pbm.efpar[4]*pbm.efpar[11]*sin(EV_[4])
            g_[5] = -pbm.efpar[5]*pbm.efpar[12]*sin(EV_[5])
            g_[6] = -pbm.efpar[6]*pbm.efpar[13]*sin(EV_[6])
            g_[7] = -pbm.efpar[7]*pbm.efpar[14]*sin(EV_[7])
            if nargout>2
                H_ = zeros(Float64,7,7)
                H_[1,1] = -pbm.efpar[1]*pbm.efpar[8]*cos(EV_[1])
                H_[2,2] = -pbm.efpar[2]*pbm.efpar[9]*cos(EV_[2])
                H_[3,3] = -pbm.efpar[3]*pbm.efpar[10]*cos(EV_[3])
                H_[4,4] = -pbm.efpar[4]*pbm.efpar[11]*cos(EV_[4])
                H_[5,5] = -pbm.efpar[5]*pbm.efpar[12]*cos(EV_[5])
                H_[6,6] = -pbm.efpar[6]*pbm.efpar[13]*cos(EV_[6])
                H_[7,7] = -pbm.efpar[7]*pbm.efpar[14]*cos(EV_[7])
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eS8T"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        S2 = pbm.efpar[8]*(pbm.efpar[1]*sin(EV_[1])-pbm.efpar[15])
        S3 = pbm.efpar[9]*(pbm.efpar[2]*sin(EV_[2])-pbm.efpar[15])+S2
        S4 = pbm.efpar[10]*(pbm.efpar[3]*sin(EV_[3])-pbm.efpar[15])+S3
        S5 = pbm.efpar[11]*(pbm.efpar[4]*sin(EV_[4])-pbm.efpar[15])+S4
        S6 = pbm.efpar[12]*(pbm.efpar[5]*sin(EV_[5])-pbm.efpar[15])+S5
        S7 = pbm.efpar[13]*(pbm.efpar[6]*sin(EV_[6])-pbm.efpar[15])+S6
        f_   = pbm.efpar[14]*(pbm.efpar[7]*sin(EV_[7])-pbm.efpar[15])+S7
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = pbm.efpar[1]*pbm.efpar[8]*cos(EV_[1])
            g_[2] = pbm.efpar[2]*pbm.efpar[9]*cos(EV_[2])
            g_[3] = pbm.efpar[3]*pbm.efpar[10]*cos(EV_[3])
            g_[4] = pbm.efpar[4]*pbm.efpar[11]*cos(EV_[4])
            g_[5] = pbm.efpar[5]*pbm.efpar[12]*cos(EV_[5])
            g_[6] = pbm.efpar[6]*pbm.efpar[13]*cos(EV_[6])
            g_[7] = pbm.efpar[7]*pbm.efpar[14]*cos(EV_[7])
            if nargout>2
                H_ = zeros(Float64,7,7)
                H_[1,1] = -pbm.efpar[1]*pbm.efpar[8]*sin(EV_[1])
                H_[2,2] = -pbm.efpar[2]*pbm.efpar[9]*sin(EV_[2])
                H_[3,3] = -pbm.efpar[3]*pbm.efpar[10]*sin(EV_[3])
                H_[4,4] = -pbm.efpar[4]*pbm.efpar[11]*sin(EV_[4])
                H_[5,5] = -pbm.efpar[5]*pbm.efpar[12]*sin(EV_[5])
                H_[6,6] = -pbm.efpar[6]*pbm.efpar[13]*sin(EV_[6])
                H_[7,7] = -pbm.efpar[7]*pbm.efpar[14]*sin(EV_[7])
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eQ8T"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        S2 = pbm.efpar[8]*(pbm.efpar[1]*sin(EV_[1])-pbm.efpar[15])
        S3 = pbm.efpar[9]*(pbm.efpar[2]*sin(EV_[2])-pbm.efpar[15])+S2
        S4 = pbm.efpar[10]*(pbm.efpar[3]*sin(EV_[3])-pbm.efpar[15])+S3
        S5 = pbm.efpar[11]*(pbm.efpar[4]*sin(EV_[4])-pbm.efpar[15])+S4
        S6 = pbm.efpar[12]*(pbm.efpar[5]*sin(EV_[5])-pbm.efpar[15])+S5
        S7 = pbm.efpar[13]*(pbm.efpar[6]*sin(EV_[6])-pbm.efpar[15])+S6
        DSD1 = pbm.efpar[1]*pbm.efpar[8]*cos(EV_[1])
        DSD2 = pbm.efpar[2]*pbm.efpar[9]*cos(EV_[2])
        DSD3 = pbm.efpar[3]*pbm.efpar[10]*cos(EV_[3])
        DSD4 = pbm.efpar[4]*pbm.efpar[11]*cos(EV_[4])
        DSD5 = pbm.efpar[5]*pbm.efpar[12]*cos(EV_[5])
        DSD6 = pbm.efpar[6]*pbm.efpar[13]*cos(EV_[6])
        DSD7 = pbm.efpar[7]*pbm.efpar[14]*cos(EV_[7])
        D2SD1 = -pbm.efpar[1]*pbm.efpar[8]*sin(EV_[1])
        D2SD2 = -pbm.efpar[2]*pbm.efpar[9]*sin(EV_[2])
        D2SD3 = -pbm.efpar[3]*pbm.efpar[10]*sin(EV_[3])
        D2SD4 = -pbm.efpar[4]*pbm.efpar[11]*sin(EV_[4])
        D2SD5 = -pbm.efpar[5]*pbm.efpar[12]*sin(EV_[5])
        D2SD6 = -pbm.efpar[6]*pbm.efpar[13]*sin(EV_[6])
        D2SD7 = -pbm.efpar[7]*pbm.efpar[14]*sin(EV_[7])
        Q2 = 0.5*pbm.efpar[8]*pbm.efpar[8]*(pbm.efpar[1]*sin(EV_[1])-pbm.efpar[15])
        Q3 = (0.5*pbm.efpar[9]*pbm.efpar[9]*(pbm.efpar[2]*sin(EV_[2])-pbm.efpar[15])+
             pbm.efpar[9]*S2+Q2)
        Q4 = (0.5*pbm.efpar[10]*pbm.efpar[10]*(pbm.efpar[3]*sin(EV_[3])-pbm.efpar[15])+
             pbm.efpar[10]*S3+Q3)
        Q5 = (0.5*pbm.efpar[11]*pbm.efpar[11]*(pbm.efpar[4]*sin(EV_[4])-pbm.efpar[15])+
             pbm.efpar[11]*S4+Q4)
        Q6 = (0.5*pbm.efpar[12]*pbm.efpar[12]*(pbm.efpar[5]*sin(EV_[5])-pbm.efpar[15])+
             pbm.efpar[12]*S5+Q5)
        Q7 = (0.5*pbm.efpar[13]*pbm.efpar[13]*(pbm.efpar[6]*sin(EV_[6])-pbm.efpar[15])+
             pbm.efpar[13]*S6+Q6)
        f_    = (
              0.5*pbm.efpar[14]*pbm.efpar[14]*(pbm.efpar[7]*sin(EV_[7])-pbm.efpar[15])+pbm.efpar[14]*S7+Q7)
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1]  = (
                  0.5*pbm.efpar[8]*pbm.efpar[8]*pbm.efpar[1]*cos(EV_[1])+(pbm.efpar[14]+pbm.efpar[13]+pbm.efpar[12]+pbm.efpar[11]+pbm.efpar[10]+pbm.efpar[9])*DSD1)
            g_[2]  = (
                  0.5*pbm.efpar[9]*pbm.efpar[9]*pbm.efpar[2]*cos(EV_[2])+(pbm.efpar[14]+pbm.efpar[13]+pbm.efpar[12]+pbm.efpar[11]+pbm.efpar[10])*DSD2)
            g_[3] = (0.5*pbm.efpar[10]*pbm.efpar[10]*pbm.efpar[3]*cos(EV_[3])+
                 (pbm.efpar[14]+pbm.efpar[13]+pbm.efpar[12]+pbm.efpar[11])*DSD3)
            g_[4] = (0.5*pbm.efpar[11]*pbm.efpar[11]*pbm.efpar[4]*cos(EV_[4])+
                 (pbm.efpar[14]+pbm.efpar[13]+pbm.efpar[12])*DSD4)
            g_[5] = (0.5*pbm.efpar[12]*pbm.efpar[12]*pbm.efpar[5]*cos(EV_[5])+
                 (pbm.efpar[14]+pbm.efpar[13])*DSD5)
            g_[6] = (0.5*pbm.efpar[13]*pbm.efpar[13]*pbm.efpar[6]*cos(EV_[6])+
                 pbm.efpar[14]*DSD6)
            g_[7] = 0.5*pbm.efpar[14]*pbm.efpar[14]*pbm.efpar[7]*cos(EV_[7])
            if nargout>2
                H_ = zeros(Float64,7,7)
                H_[1,1] = (-0.5*pbm.efpar[8]*pbm.efpar[8]*pbm.efpar[1]*sin(EV_[1])+
                     (pbm.efpar[14]+pbm.efpar[13]+pbm.efpar[12]+pbm.efpar[11]+pbm.efpar[10]+pbm.efpar[9])*D2SD1)
                H_[2,2] = (-0.5*pbm.efpar[9]*pbm.efpar[9]*pbm.efpar[2]*sin(EV_[2])+
                     (pbm.efpar[14]+pbm.efpar[13]+pbm.efpar[12]+pbm.efpar[11]+pbm.efpar[10])*D2SD2)
                H_[3,3] = (-0.5*pbm.efpar[10]*pbm.efpar[10]*pbm.efpar[3]*sin(EV_[3])+
                     (pbm.efpar[14]+pbm.efpar[13]+pbm.efpar[12]+pbm.efpar[11])*D2SD3)
                H_[4,4] = (-0.5*pbm.efpar[11]*pbm.efpar[11]*pbm.efpar[4]*sin(EV_[4])+
                     (pbm.efpar[14]+pbm.efpar[13]+pbm.efpar[12])*D2SD4)
                H_[5,5] = (-0.5*pbm.efpar[12]*pbm.efpar[12]*pbm.efpar[5]*sin(EV_[5])+
                     (pbm.efpar[14]+pbm.efpar[13])*D2SD5)
                H_[6,6] = (-0.5*pbm.efpar[13]*pbm.efpar[13]*pbm.efpar[6]*sin(EV_[6])+
                     pbm.efpar[14]*D2SD6)
                H_[7,7] = -0.5*pbm.efpar[14]*pbm.efpar[14]*pbm.efpar[7]*sin(EV_[7])
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

    elseif action == "gL2"

        GVAR_   = args[1]
        igr_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_= GVAR_*GVAR_
        if nargout>1
            g_ = GVAR_+GVAR_
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_ = 2.0
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
            pbm.has_globs = [15,0]
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

