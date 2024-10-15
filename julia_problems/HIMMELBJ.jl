function HIMMELBJ(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HIMMELBJ
#    *********
# 
#    An chemical equilibrium problem by A.P. Jones.
#    It has a nonlinear objective and linear constraints
# 
#    Source: problem 6 in
#    D.H. Himmelblau,
#    "Applied nonlinear programming",
#    McGraw-Hill, New-York, 1972.
# 
#    SIF input: Ph. Toint, March 1991.
# 
#    classification = "C-OLR2-MY-45-14"
# 
#    Number of variable sets
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "HIMMELBJ"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["NSETS"] = 7
        v_["NS1"] = 4.0
        v_["NS2"] = 13.0
        v_["NS3"] = 18.0
        v_["NS4"] = 3.0
        v_["NS5"] = 3.0
        v_["NS6"] = 2.0
        v_["NS7"] = 2.0
        v_["NEQ"] = 16
        v_["1"] = 1
        v_["2"] = 2
        v_["3"] = 3
        v_["4"] = 4
        v_["13"] = 13
        v_["18"] = 18
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for K = Int64(v_["1"]):Int64(v_["NSETS"])
            v_["RNSK"] = v_["NS"*string(K)]
            v_["NSK"] = trunc(Int,v_["RNSK"])
            for J = Int64(v_["1"]):Int64(v_["NSK"])
                iv,ix_,_ = s2mpj_ii("X"*string(J)*","*string(K),ix_)
                arrset(pb.xnames,iv,"X"*string(J)*","*string(K))
            end
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X2,1"]
        pbm.A[ig,iv] += Float64(-7.69)
        iv = ix_["X3,1"]
        pbm.A[ig,iv] += Float64(-11.52)
        iv = ix_["X4,1"]
        pbm.A[ig,iv] += Float64(-36.60)
        iv = ix_["X1,2"]
        pbm.A[ig,iv] += Float64(-10.94)
        iv = ix_["X8,2"]
        pbm.A[ig,iv] += Float64(2.5966)
        iv = ix_["X9,2"]
        pbm.A[ig,iv] += Float64(-39.39)
        iv = ix_["X10,2"]
        pbm.A[ig,iv] += Float64(-21.35)
        iv = ix_["X11,2"]
        pbm.A[ig,iv] += Float64(-32.84)
        iv = ix_["X12,2"]
        pbm.A[ig,iv] += Float64(6.26)
        iv = ix_["X1,3"]
        pbm.A[ig,iv] += Float64(10.45)
        iv = ix_["X3,3"]
        pbm.A[ig,iv] += Float64(-0.5)
        iv = ix_["X7,3"]
        pbm.A[ig,iv] += Float64(2.2435)
        iv = ix_["X9,3"]
        pbm.A[ig,iv] += Float64(-39.39)
        iv = ix_["X10,3"]
        pbm.A[ig,iv] += Float64(-21.49)
        iv = ix_["X11,3"]
        pbm.A[ig,iv] += Float64(-32.84)
        iv = ix_["X12,3"]
        pbm.A[ig,iv] += Float64(6.12)
        iv = ix_["X15,3"]
        pbm.A[ig,iv] += Float64(-1.9028)
        iv = ix_["X16,3"]
        pbm.A[ig,iv] += Float64(-2.8889)
        iv = ix_["X17,3"]
        pbm.A[ig,iv] += Float64(-3.3622)
        iv = ix_["X18,3"]
        pbm.A[ig,iv] += Float64(-7.4854)
        iv = ix_["X1,4"]
        pbm.A[ig,iv] += Float64(-15.639)
        iv = ix_["X3,4"]
        pbm.A[ig,iv] += Float64(21.81)
        iv = ix_["X1,5"]
        pbm.A[ig,iv] += Float64(-16.79)
        iv = ix_["X3,5"]
        pbm.A[ig,iv] += Float64(18.9779)
        iv = ix_["X2,6"]
        pbm.A[ig,iv] += Float64(11.959)
        iv = ix_["X2,7"]
        pbm.A[ig,iv] += Float64(12.899)
        ig,ig_,_ = s2mpj_ii("C1",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C1")
        iv = ix_["X1,1"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X1,2"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X1,3"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X15,3"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X16,3"]
        pbm.A[ig,iv] += Float64(2.0)
        iv = ix_["X17,3"]
        pbm.A[ig,iv] += Float64(3.0)
        iv = ix_["X18,3"]
        pbm.A[ig,iv] += Float64(4.0)
        ig,ig_,_ = s2mpj_ii("C2",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C2")
        iv = ix_["X2,1"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X2,2"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X10,2"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X11,2"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X12,2"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X2,3"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X10,3"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X11,3"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X12,3"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X2,6"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X2,7"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("C3",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C3")
        iv = ix_["X3,1"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X3,2"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X3,3"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("C4",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C4")
        iv = ix_["X4,1"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X4,2"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X9,2"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X11,2"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X12,2"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X4,3"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X9,3"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X11,3"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X12,3"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X1,4"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X3,4"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X1,5"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X3,5"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X2,6"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X2,7"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("C5",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C5")
        iv = ix_["X4,1"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X5,2"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X9,2"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X10,2"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X11,2"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X12,2"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X5,3"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X9,3"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X10,3"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X11,3"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X12,3"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("C6",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C6")
        iv = ix_["X6,2"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X6,3"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("C7",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C7")
        iv = ix_["X7,2"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X7,3"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("C8",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C8")
        iv = ix_["X8,2"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X8,3"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("C11",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C11")
        iv = ix_["X14,3"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X15,3"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X16,3"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X17,3"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X18,3"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("C12",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C12")
        iv = ix_["X4,2"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X5,2"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X6,2"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X7,2"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X8,2"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X10,2"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X12,2"]
        pbm.A[ig,iv] += Float64(-2.0)
        iv = ix_["X13,2"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X14,3"]
        pbm.A[ig,iv] += Float64(-4.0)
        iv = ix_["X15,3"]
        pbm.A[ig,iv] += Float64(-3.0)
        iv = ix_["X16,3"]
        pbm.A[ig,iv] += Float64(-2.0)
        iv = ix_["X17,3"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("C13",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C13")
        iv = ix_["X15,3"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X16,3"]
        pbm.A[ig,iv] += Float64(-2.0)
        iv = ix_["X17,3"]
        pbm.A[ig,iv] += Float64(-3.0)
        iv = ix_["X18,3"]
        pbm.A[ig,iv] += Float64(-4.0)
        iv = ix_["X1,4"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X2,4"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X3,4"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("C14",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C14")
        iv = ix_["X1,5"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X2,5"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X3,5"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("C15",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C15")
        iv = ix_["X3,4"]
        pbm.A[ig,iv] += Float64(-4.0)
        iv = ix_["X1,6"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X2,6"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("C16",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C16")
        iv = ix_["X3,5"]
        pbm.A[ig,iv] += Float64(-4.0)
        iv = ix_["X1,7"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X2,7"]
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
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        pbm.gconst[ig_["C1"]] = Float64(0.652981)
        pbm.gconst[ig_["C2"]] = Float64(0.281941)
        pbm.gconst[ig_["C3"]] = Float64(3.705233)
        pbm.gconst[ig_["C4"]] = Float64(47.00022)
        pbm.gconst[ig_["C5"]] = Float64(47.02972)
        pbm.gconst[ig_["C6"]] = Float64(0.08005)
        pbm.gconst[ig_["C7"]] = Float64(0.08813)
        pbm.gconst[ig_["C8"]] = Float64(0.04829)
        pbm.gconst[ig_["C11"]] = Float64(0.0022725)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = fill(1.e-12,pb.n)
        pb.xupper = fill(Inf,pb.n)
        pb.xlower[ix_["X13,2"]] = 0.0155
        pb.xupper[ix_["X13,2"]] = 0.0155
        pb.xlower[ix_["X13,3"]] = 0.0211275
        pb.xupper[ix_["X13,3"]] = 0.0211275
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(0.1),pb.n)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eXLOGX", iet_)
        loaset(elftv,it,1,"X")
        it,iet_,_ = s2mpj_ii( "eXLOGX2", iet_)
        loaset(elftv,it,1,"Y1")
        loaset(elftv,it,2,"Y2")
        loaset(elftv,it,3,"X")
        it,iet_,_ = s2mpj_ii( "eXLOGX3", iet_)
        loaset(elftv,it,1,"Y1")
        loaset(elftv,it,2,"Y2")
        loaset(elftv,it,3,"Y3")
        loaset(elftv,it,4,"X")
        it,iet_,_ = s2mpj_ii( "eXLOGX4", iet_)
        loaset(elftv,it,1,"Y1")
        loaset(elftv,it,2,"Y2")
        loaset(elftv,it,3,"Y3")
        loaset(elftv,it,4,"Y4")
        loaset(elftv,it,5,"X")
        it,iet_,_ = s2mpj_ii( "eXLOGX13", iet_)
        loaset(elftv,it,1,"Y1")
        loaset(elftv,it,2,"Y2")
        loaset(elftv,it,3,"Y3")
        loaset(elftv,it,4,"Y4")
        loaset(elftv,it,5,"Y5")
        loaset(elftv,it,6,"Y6")
        loaset(elftv,it,7,"Y7")
        loaset(elftv,it,8,"Y8")
        loaset(elftv,it,9,"Y9")
        loaset(elftv,it,10,"Y10")
        loaset(elftv,it,11,"Y11")
        loaset(elftv,it,12,"Y12")
        loaset(elftv,it,13,"Y13")
        loaset(elftv,it,14,"X")
        it,iet_,_ = s2mpj_ii( "eXLOGX18", iet_)
        loaset(elftv,it,1,"Y1")
        loaset(elftv,it,2,"Y2")
        loaset(elftv,it,3,"Y3")
        loaset(elftv,it,4,"Y4")
        loaset(elftv,it,5,"Y5")
        loaset(elftv,it,6,"Y6")
        loaset(elftv,it,7,"Y7")
        loaset(elftv,it,8,"Y8")
        loaset(elftv,it,9,"Y9")
        loaset(elftv,it,10,"Y10")
        loaset(elftv,it,11,"Y11")
        loaset(elftv,it,12,"Y12")
        loaset(elftv,it,13,"Y13")
        loaset(elftv,it,14,"Y14")
        loaset(elftv,it,15,"Y15")
        loaset(elftv,it,16,"Y16")
        loaset(elftv,it,17,"Y17")
        loaset(elftv,it,18,"Y18")
        loaset(elftv,it,19,"X")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for K = Int64(v_["1"]):Int64(v_["NSETS"])
            v_["RNSK"] = v_["NS"*string(K)]
            v_["NSK"] = trunc(Int,v_["RNSK"])
            for J = Int64(v_["1"]):Int64(v_["NSK"])
                ename = "A"*string(J)*","*string(K)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eXLOGX")
                arrset(ielftype,ie,iet_["eXLOGX"])
                vname = "X"*string(J)*","*string(K)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
                posev = findfirst(x->x=="X",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
        end
        v_["K"] = 1
        for J = Int64(v_["1"]):Int64(v_["4"])
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eXLOGX4")
            arrset(ielftype,ie,iet_["eXLOGX4"])
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X"*string(J)*","*string(Int64(v_["K"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X1,1"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X2,1"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X3,1"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y3",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X4,1"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y4",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        v_["K"] = 2
        for J = Int64(v_["1"]):Int64(v_["13"])
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eXLOGX13")
            arrset(ielftype,ie,iet_["eXLOGX13"])
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X"*string(J)*","*string(Int64(v_["K"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X1,2"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X2,2"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X3,2"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y3",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X4,2"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y4",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X5,2"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y5",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X6,2"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y6",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X7,2"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y7",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X8,2"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y8",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X9,2"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y9",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X10,2"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y10",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X11,2"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y11",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X12,2"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y12",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X13,2"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y13",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        v_["K"] = 3
        for J = Int64(v_["1"]):Int64(v_["18"])
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eXLOGX18")
            arrset(ielftype,ie,iet_["eXLOGX18"])
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X"*string(J)*","*string(Int64(v_["K"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X1,3"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X2,3"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X3,3"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y3",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X4,3"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y4",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X5,3"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y5",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X6,3"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y6",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X7,3"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y7",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X8,3"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y8",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X9,3"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y9",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X10,3"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y10",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X11,3"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y11",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X12,3"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y12",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X13,3"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y13",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X14,3"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y14",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X15,3"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y15",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X16,3"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y16",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X17,3"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y17",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X18,3"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y18",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        v_["K"] = 4
        for J = Int64(v_["1"]):Int64(v_["3"])
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eXLOGX3")
            arrset(ielftype,ie,iet_["eXLOGX3"])
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X"*string(J)*","*string(Int64(v_["K"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X1,4"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X2,4"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X3,4"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y3",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        v_["K"] = 5
        for J = Int64(v_["1"]):Int64(v_["3"])
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eXLOGX3")
            arrset(ielftype,ie,iet_["eXLOGX3"])
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X"*string(J)*","*string(Int64(v_["K"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X1,5"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X2,5"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X3,5"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y3",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        v_["K"] = 6
        for J = Int64(v_["1"]):Int64(v_["2"])
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eXLOGX2")
            arrset(ielftype,ie,iet_["eXLOGX2"])
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X"*string(J)*","*string(Int64(v_["K"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X1,6"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X2,6"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        v_["K"] = 7
        for J = Int64(v_["1"]):Int64(v_["2"])
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eXLOGX2")
            arrset(ielftype,ie,iet_["eXLOGX2"])
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X"*string(J)*","*string(Int64(v_["K"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X1,7"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(J)*","*string(Int64(v_["K"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X2,7"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,Float64(0.1))
            posev = findfirst(x->x=="Y2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for K = Int64(v_["1"]):Int64(v_["NSETS"])
            v_["RNSK"] = v_["NS"*string(K)]
            v_["NSK"] = trunc(Int,v_["RNSK"])
            for J = Int64(v_["1"]):Int64(v_["NSK"])
                ig = ig_["OBJ"]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A"*string(J)*","*string(K)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,1.)
                posel = posel+1
                loaset(pbm.grelt,ig,posel,ie_["B"*string(J)*","*string(K)])
                loaset(pbm.grelw,ig,posel,Float64(-1.0))
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               -1910.344724
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
        pb.pbclass = "C-OLR2-MY-45-14"
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

    elseif action == "eXLOGX"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        LOGX = log(EV_[1])
        f_   = EV_[1]*LOGX
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = LOGX+1.0
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 1.0/EV_[1]
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eXLOGX2"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,2,3)
        IV_ =  zeros(Float64,2)
        U_[2,1] = U_[2,1]+1
        U_[2,2] = U_[2,2]+1
        U_[1,3] = U_[1,3]+1
        IV_[1] = dot(U_[1,:],EV_)
        IV_[2] = dot(U_[2,:],EV_)
        LOGX = log(IV_[2])
        f_   = IV_[1]*LOGX
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = LOGX
            g_[2] = IV_[1]/IV_[2]
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = 1.0/IV_[2]
                H_[2,1] = H_[1,2]
                H_[2,2] = -IV_[1]/IV_[2]^2
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

    elseif action == "eXLOGX3"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,2,4)
        IV_ =  zeros(Float64,2)
        U_[2,1] = U_[2,1]+1
        U_[2,2] = U_[2,2]+1
        U_[2,3] = U_[2,3]+1
        U_[1,4] = U_[1,4]+1
        IV_[1] = dot(U_[1,:],EV_)
        IV_[2] = dot(U_[2,:],EV_)
        LOGX = log(IV_[2])
        f_   = IV_[1]*LOGX
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = LOGX
            g_[2] = IV_[1]/IV_[2]
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = 1.0/IV_[2]
                H_[2,1] = H_[1,2]
                H_[2,2] = -IV_[1]/IV_[2]^2
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

    elseif action == "eXLOGX4"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,2,5)
        IV_ =  zeros(Float64,2)
        U_[2,1] = U_[2,1]+1
        U_[2,2] = U_[2,2]+1
        U_[2,3] = U_[2,3]+1
        U_[2,4] = U_[2,4]+1
        U_[1,5] = U_[1,5]+1
        IV_[1] = dot(U_[1,:],EV_)
        IV_[2] = dot(U_[2,:],EV_)
        LOGX = log(IV_[2])
        f_   = IV_[1]*LOGX
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = LOGX
            g_[2] = IV_[1]/IV_[2]
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = 1.0/IV_[2]
                H_[2,1] = H_[1,2]
                H_[2,2] = -IV_[1]/IV_[2]^2
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

    elseif action == "eXLOGX13"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,2,14)
        IV_ =  zeros(Float64,2)
        U_[2,1] = U_[2,1]+1
        U_[2,2] = U_[2,2]+1
        U_[2,3] = U_[2,3]+1
        U_[2,4] = U_[2,4]+1
        U_[2,5] = U_[2,5]+1
        U_[2,6] = U_[2,6]+1
        U_[2,7] = U_[2,7]+1
        U_[2,8] = U_[2,8]+1
        U_[2,9] = U_[2,9]+1
        U_[2,10] = U_[2,10]+1
        U_[2,11] = U_[2,11]+1
        U_[2,12] = U_[2,12]+1
        U_[2,13] = U_[2,13]+1
        U_[1,14] = U_[1,14]+1
        IV_[1] = dot(U_[1,:],EV_)
        IV_[2] = dot(U_[2,:],EV_)
        LOGX = log(IV_[2])
        f_   = IV_[1]*LOGX
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = LOGX
            g_[2] = IV_[1]/IV_[2]
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = 1.0/IV_[2]
                H_[2,1] = H_[1,2]
                H_[2,2] = -IV_[1]/IV_[2]^2
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

    elseif action == "eXLOGX18"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,2,19)
        IV_ =  zeros(Float64,2)
        U_[2,1] = U_[2,1]+1
        U_[2,2] = U_[2,2]+1
        U_[2,3] = U_[2,3]+1
        U_[2,4] = U_[2,4]+1
        U_[2,5] = U_[2,5]+1
        U_[2,6] = U_[2,6]+1
        U_[2,7] = U_[2,7]+1
        U_[2,8] = U_[2,8]+1
        U_[2,9] = U_[2,9]+1
        U_[2,10] = U_[2,10]+1
        U_[2,11] = U_[2,11]+1
        U_[2,12] = U_[2,12]+1
        U_[2,13] = U_[2,13]+1
        U_[2,14] = U_[2,14]+1
        U_[2,15] = U_[2,15]+1
        U_[2,16] = U_[2,16]+1
        U_[2,17] = U_[2,17]+1
        U_[2,18] = U_[2,18]+1
        U_[1,19] = U_[1,19]+1
        IV_[1] = dot(U_[1,:],EV_)
        IV_[2] = dot(U_[2,:],EV_)
        LOGX = log(IV_[2])
        f_   = IV_[1]*LOGX
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = LOGX
            g_[2] = IV_[1]/IV_[2]
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = 1.0/IV_[2]
                H_[2,1] = H_[1,2]
                H_[2,2] = -IV_[1]/IV_[2]^2
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

