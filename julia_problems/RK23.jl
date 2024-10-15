function RK23(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : RK23
#    *********
# 
#    Find coefficients for an embedded pair of explicit 2nd
#    and 3rd order Runge Kutta Method such that the leading
#    term in the local truncation error is minimized.
# 
#    Source:
#    Similar ideas for 4th and 5th order pairs are discussed in:
#    Hairer, Norsett and Wanner, Solving Ordinary Differential
#    Equations I, Springer 1980, page 158 ff.
# 
#    SIF input: S. Leyffer, January 1997.
# 
#    classification = "C-LOR2-RN-17-11"
# 
# 
#    ... COMPUTED PARAMETERS
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "RK23"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["ONE"] = 1.0
        v_["THREE"] = 3.0
        v_["FOUR"] = 4.0
        v_["SIX"] = 6.0
        v_["ONETHIRD"] = v_["ONE"]/v_["THREE"]
        v_["ONESIXTH"] = v_["ONE"]/v_["SIX"]
        v_["FOURSIXTH"] = v_["FOUR"]/v_["SIX"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        iv,ix_,_ = s2mpj_ii("C2",ix_)
        arrset(pb.xnames,iv,"C2")
        iv,ix_,_ = s2mpj_ii("A21",ix_)
        arrset(pb.xnames,iv,"A21")
        iv,ix_,_ = s2mpj_ii("C3",ix_)
        arrset(pb.xnames,iv,"C3")
        iv,ix_,_ = s2mpj_ii("A31",ix_)
        arrset(pb.xnames,iv,"A31")
        iv,ix_,_ = s2mpj_ii("A32",ix_)
        arrset(pb.xnames,iv,"A32")
        iv,ix_,_ = s2mpj_ii("B1",ix_)
        arrset(pb.xnames,iv,"B1")
        iv,ix_,_ = s2mpj_ii("B2",ix_)
        arrset(pb.xnames,iv,"B2")
        iv,ix_,_ = s2mpj_ii("B3",ix_)
        arrset(pb.xnames,iv,"B3")
        iv,ix_,_ = s2mpj_ii("BB1",ix_)
        arrset(pb.xnames,iv,"BB1")
        iv,ix_,_ = s2mpj_ii("BB2",ix_)
        arrset(pb.xnames,iv,"BB2")
        iv,ix_,_ = s2mpj_ii("BB3",ix_)
        arrset(pb.xnames,iv,"BB3")
        iv,ix_,_ = s2mpj_ii("TP1",ix_)
        arrset(pb.xnames,iv,"TP1")
        iv,ix_,_ = s2mpj_ii("TM1",ix_)
        arrset(pb.xnames,iv,"TM1")
        iv,ix_,_ = s2mpj_ii("TP2",ix_)
        arrset(pb.xnames,iv,"TP2")
        iv,ix_,_ = s2mpj_ii("TM2",ix_)
        arrset(pb.xnames,iv,"TM2")
        iv,ix_,_ = s2mpj_ii("TP3",ix_)
        arrset(pb.xnames,iv,"TP3")
        iv,ix_,_ = s2mpj_ii("TM3",ix_)
        arrset(pb.xnames,iv,"TM3")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["TP1"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["TM1"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["TP2"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["TM2"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["TP3"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["TM3"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("ROWS1",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ROWS1")
        iv = ix_["A21"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["C2"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("ROWS2",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ROWS2")
        iv = ix_["A31"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["A32"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["C3"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("FIRST2",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"FIRST2")
        iv = ix_["B1"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["B2"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["B3"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("FIRST3",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"FIRST3")
        iv = ix_["BB1"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["BB2"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["BB3"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("SECND2",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"SECND2")
        ig,ig_,_ = s2mpj_ii("SECND3",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"SECND3")
        ig,ig_,_ = s2mpj_ii("THIRD31",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"THIRD31")
        ig,ig_,_ = s2mpj_ii("THIRD32",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"THIRD32")
        ig,ig_,_ = s2mpj_ii("ART1",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ART1")
        iv = ix_["TP1"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["TM2"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("ART2",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ART2")
        iv = ix_["TP2"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["TM2"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("ART3",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ART3")
        iv = ix_["TP3"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["TM3"]
        pbm.A[ig,iv] += Float64(-1.0)
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
        pbm.gconst[ig_["ROWS1"]] = Float64(0.0)
        pbm.gconst[ig_["ROWS2"]] = Float64(0.0)
        pbm.gconst[ig_["FIRST2"]] = Float64(1.0)
        pbm.gconst[ig_["FIRST3"]] = Float64(1.0)
        pbm.gconst[ig_["SECND2"]] = Float64(0.5)
        pbm.gconst[ig_["SECND3"]] = Float64(0.5)
        pbm.gconst[ig_["THIRD31"]] = Float64(v_["ONETHIRD"])
        pbm.gconst[ig_["THIRD32"]] = Float64(v_["ONESIXTH"])
        pbm.gconst[ig_["ART1"]] = Float64(1.0)
        pbm.gconst[ig_["ART2"]] = Float64(1.0)
        pbm.gconst[ig_["ART3"]] = Float64(1.0)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        pb.xlower[ix_["TP1"]] = 0.0
        pb.xlower[ix_["TM1"]] = 0.0
        pb.xlower[ix_["TP2"]] = 0.0
        pb.xlower[ix_["TM2"]] = 0.0
        pb.xlower[ix_["TP3"]] = 0.0
        pb.xlower[ix_["TM3"]] = 0.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        if haskey(ix_,"C2")
            pb.x0[ix_["C2"]] = Float64(1.0)
        else
            pb.y0[findfirst(x->x==ig_["C2"],pbm.congrps)] = Float64(1.0)
        end
        if haskey(ix_,"A21")
            pb.x0[ix_["A21"]] = Float64(1.0)
        else
            pb.y0[findfirst(x->x==ig_["A21"],pbm.congrps)] = Float64(1.0)
        end
        if haskey(ix_,"C3")
            pb.x0[ix_["C3"]] = Float64(0.5)
        else
            pb.y0[findfirst(x->x==ig_["C3"],pbm.congrps)] = Float64(0.5)
        end
        if haskey(ix_,"A31")
            pb.x0[ix_["A31"]] = Float64(0.25)
        else
            pb.y0[findfirst(x->x==ig_["A31"],pbm.congrps)] = Float64(0.25)
        end
        if haskey(ix_,"A32")
            pb.x0[ix_["A32"]] = Float64(0.25)
        else
            pb.y0[findfirst(x->x==ig_["A32"],pbm.congrps)] = Float64(0.25)
        end
        if haskey(ix_,"B1")
            pb.x0[ix_["B1"]] = Float64(0.5)
        else
            pb.y0[findfirst(x->x==ig_["B1"],pbm.congrps)] = Float64(0.5)
        end
        if haskey(ix_,"B2")
            pb.x0[ix_["B2"]] = Float64(0.5)
        else
            pb.y0[findfirst(x->x==ig_["B2"],pbm.congrps)] = Float64(0.5)
        end
        if haskey(ix_,"B3")
            pb.x0[ix_["B3"]] = Float64(0.0)
        else
            pb.y0[findfirst(x->x==ig_["B3"],pbm.congrps)] = Float64(0.0)
        end
        pb.x0[ix_["BB1"]] = Float64(v_["ONESIXTH"])
        pb.x0[ix_["BB2"]] = Float64(v_["ONESIXTH"])
        pb.x0[ix_["BB3"]] = Float64(v_["FOURSIXTH"])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "ePROD", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        it,iet_,_ = s2mpj_ii( "ePRODS", iet_)
        loaset(elftv,it,1,"W1")
        loaset(elftv,it,2,"W2")
        it,iet_,_ = s2mpj_ii( "ePRODQ", iet_)
        loaset(elftv,it,1,"X1")
        loaset(elftv,it,2,"X2")
        it,iet_,_ = s2mpj_ii( "eTPROD", iet_)
        loaset(elftv,it,1,"Y1")
        loaset(elftv,it,2,"Y2")
        loaset(elftv,it,3,"Y3")
        it,iet_,_ = s2mpj_ii( "eQPROD", iet_)
        loaset(elftv,it,1,"Z1")
        loaset(elftv,it,2,"Z2")
        loaset(elftv,it,3,"Z3")
        loaset(elftv,it,4,"Z4")
        it,iet_,_ = s2mpj_ii( "eTPRODS", iet_)
        loaset(elftv,it,1,"U1")
        loaset(elftv,it,2,"U2")
        loaset(elftv,it,3,"U3")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "E1"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD")
        arrset(ielftype,ie,iet_["ePROD"])
        vname = "B2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "C2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E2"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD")
        arrset(ielftype,ie,iet_["ePROD"])
        vname = "B3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "C3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E3"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD")
        arrset(ielftype,ie,iet_["ePROD"])
        vname = "BB2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "C2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E4"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD")
        arrset(ielftype,ie,iet_["ePROD"])
        vname = "BB3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "C3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E5"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePRODS")
        arrset(ielftype,ie,iet_["ePRODS"])
        vname = "BB2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="W1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "C2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="W2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E6"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePRODS")
        arrset(ielftype,ie,iet_["ePRODS"])
        vname = "BB3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="W1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "C3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="W2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E7"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eTPROD")
        arrset(ielftype,ie,iet_["eTPROD"])
        vname = "BB3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "A32"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "C2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E8"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePRODQ")
        arrset(ielftype,ie,iet_["ePRODQ"])
        vname = "BB2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "C2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E9"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePRODQ")
        arrset(ielftype,ie,iet_["ePRODQ"])
        vname = "BB3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "C3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E10"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eQPROD")
        arrset(ielftype,ie,iet_["eQPROD"])
        vname = "BB3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "C3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "A32"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "C2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z4",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E11"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eTPRODS")
        arrset(ielftype,ie,iet_["eTPRODS"])
        vname = "BB3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "A32"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "C2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["SECND2"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E2"])
        loaset(pbm.grelw,ig,posel, 1.)
        ig = ig_["SECND3"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E3"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E4"])
        loaset(pbm.grelw,ig,posel, 1.)
        ig = ig_["THIRD31"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E5"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E6"])
        loaset(pbm.grelw,ig,posel, 1.)
        ig = ig_["THIRD32"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E7"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ART1"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E8"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(4.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E9"])
        loaset(pbm.grelw,ig,posel,Float64(4.0))
        ig = ig_["ART2"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E10"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(8.0))
        ig = ig_["ART3"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E11"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(12.0))
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
        pb.pbclass = "C-LOR2-RN-17-11"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2]
        pb.conderlvl  = pbm.conderlvl;
        return pb, pbm


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

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

    elseif action == "ePRODS"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]*(EV_[2]^2.0)
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[2]^2
            g_[2] = EV_[1]*2.0*EV_[2]
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = 2.0*EV_[2]
                H_[2,1] = H_[1,2]
                H_[2,2] = 2.0*EV_[1]
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "ePRODQ"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]*(EV_[2]^3.0)
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[2]^3.0
            g_[2] = EV_[1]*3.0*(EV_[2]^2.0)
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = 3.0*(EV_[2]^2.0)
                H_[2,1] = H_[1,2]
                H_[2,2] = EV_[1]*6.0*EV_[2]
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eTPROD"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]*EV_[2]*EV_[3]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[2]*EV_[3]
            g_[2] = EV_[1]*EV_[3]
            g_[3] = EV_[1]*EV_[2]
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,2] = EV_[3]
                H_[2,1] = H_[1,2]
                H_[1,3] = EV_[2]
                H_[3,1] = H_[1,3]
                H_[2,3] = EV_[1]
                H_[3,2] = H_[2,3]
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eQPROD"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]*EV_[2]*EV_[3]*EV_[4]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[2]*EV_[3]*EV_[4]
            g_[2] = EV_[1]*EV_[3]*EV_[4]
            g_[3] = EV_[1]*EV_[2]*EV_[4]
            g_[4] = EV_[1]*EV_[2]*EV_[3]
            if nargout>2
                H_ = zeros(Float64,4,4)
                H_[1,2] = EV_[3]*EV_[4]
                H_[2,1] = H_[1,2]
                H_[1,3] = EV_[2]*EV_[4]
                H_[3,1] = H_[1,3]
                H_[1,4] = EV_[2]*EV_[3]
                H_[4,1] = H_[1,4]
                H_[2,3] = EV_[1]*EV_[4]
                H_[3,2] = H_[2,3]
                H_[2,4] = EV_[1]*EV_[3]
                H_[4,2] = H_[2,4]
                H_[3,4] = EV_[1]*EV_[2]
                H_[4,3] = H_[3,4]
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eTPRODS"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]*EV_[2]*(EV_[3]^2.0)
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[2]*(EV_[3]^2.0)
            g_[2] = EV_[1]*(EV_[3]^2.0)
            g_[3] = EV_[1]*EV_[2]*2.0*EV_[3]
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,2] = EV_[3]^2.0
                H_[2,1] = H_[1,2]
                H_[1,3] = EV_[2]*2.0*EV_[3]
                H_[3,1] = H_[1,3]
                H_[2,3] = EV_[1]*2.0*EV_[3]
                H_[3,2] = H_[2,3]
                H_[3,3] = EV_[1]*EV_[2]*2.0
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

