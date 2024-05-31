function HS85(action,args...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS85
#    *********
# 
#    The problem is to optimize the net profit of an hypothetical
#    wood-pulp plant. The constraints include the usual material
#    and energy balances as well as several empirical equations.
# 
#    Source: problem 85 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: Nick Gould, September 1991.
#      SAVEs removed December 3rd 2014
# 
#    classification = "OOI2-MN-5-21"
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "HS85"

    if action == "setup"
        pbm          = PBM(name)
        pb           = PB(name)
        pb.sifpbname = "HS85"
        nargin       = length(args)
        pbm.call     = eval( Meta.parse( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["N"] = 5
        v_["A2"] = 17.505
        v_["A3"] = 11.275
        v_["A4"] = 214.228
        v_["A5"] = 7.458
        v_["A6"] = 0.961
        v_["A7"] = 1.612
        v_["A8"] = 0.146
        v_["A9"] = 107.99
        v_["A10"] = 922.693
        v_["A11"] = 926.832
        v_["A12"] = 18.766
        v_["A13"] = 1072.163
        v_["A14"] = 8961.448
        v_["A15"] = 0.063
        v_["A16"] = 71084.33
        v_["A17"] = 2802713.0
        v_["B2"] = 1053.6667
        v_["B3"] = 35.03
        v_["B4"] = 665.585
        v_["B5"] = 584.463
        v_["B6"] = 265.916
        v_["B7"] = 7.046
        v_["B8"] = 0.222
        v_["B9"] = 273.366
        v_["B10"] = 1286.105
        v_["B11"] = 1444.046
        v_["B12"] = 537.141
        v_["B13"] = 3247.039
        v_["B14"] = 26844.086
        v_["B15"] = 0.386
        v_["B16"] = 140000.0
        v_["B17"] = 12146108.0
        v_["1"] = 1
        v_["2"] = 2
        v_["17"] = 17
        v_["19"] = 19
        v_["20"] = 20
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        xscale  = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2x_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2x_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        ig,ig_,_ = s2x_ii("CON0",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"CON0")
        iv = ix_["X2"]
        pbm.A[ig,iv] += 1.5
        iv = ix_["X3"]
        pbm.A[ig,iv] += -1.0
        for I = Int64(v_["1"]):Int64(v_["20"])
            ig,ig_,_ = s2x_ii("CON"*string(I),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"CON"*string(I))
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
        pbm.congrps = findall(x->x!="<>",gtype)
        pb.nob = ngrp-pb.m
        pbm.objgrps = findall(x->x=="<>",gtype)
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        pbm.gconst[ig_["OBJ"]] = 0.1365
        pbm.gconst[ig_["CON1"]] = 213.1
        for I = Int64(v_["2"]):Int64(v_["17"])
            pbm.gconst[ig_["CON"*string(I)]] = v_["A"*string(I)]
        end
        pbm.gconst[ig_["CON19"]] = -21.0
        pbm.gconst[ig_["CON20"]] = 110.6
        #%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%
        grange = Vector{Float64}(undef,ngrp)
        grange[gegrps,1] = fill(Inf,pb.nge)
        for I = Int64(v_["2"]):Int64(v_["17"])
            v_["DIF"] = v_["B"*string(I)]-v_["A"*string(I)]
            arrset(grange,ig_["CON"*string(I)],v_["DIF"])
        end
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower[ix_["X1"]] = 7.044148e+2
        pb.xupper[ix_["X1"]] = 9.063855e+2
        pb.xlower[ix_["X2"]] = 6.86e+1
        pb.xupper[ix_["X2"]] = 2.8888e+2
        pb.xupper[ix_["X3"]] = 1.3475e+2
        pb.xlower[ix_["X4"]] = 1.930e+2
        pb.xupper[ix_["X4"]] = 2.870966e+2
        pb.xlower[ix_["X5"]] = 2.50e+1
        pb.xupper[ix_["X5"]] = 8.41988e+1
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        if haskey(ix_,"X1")
            pb.x0[ix_["X1"]] = 9.0e+2
        else
            pb.y0[findfirst(x->x==ig_["X1"],pbm.congrps)]pb.y0[findfirst(x->x==ig_[9.0e+2],pbm.congrps)]
        end
        if haskey(ix_,"X2")
            pb.x0[ix_["X2"]] = 8.0e+1
        else
            pb.y0[findfirst(x->x==ig_["X2"],pbm.congrps)] = 8.0e+1
        end
        if haskey(ix_,"X3")
            pb.x0[ix_["X3"]] = 1.15e+2
        else
            pb.y0[findfirst(x->x==ig_["X3"],pbm.congrps)]pb.y0[findfirst(x->x==ig_[1.15e+2],pbm.congrps)]
        end
        if haskey(ix_,"X4")
            pb.x0[ix_["X4"]] = 2.67e+2
        else
            pb.y0[findfirst(x->x==ig_["X4"],pbm.congrps)] = 2.67e+2
        end
        if haskey(ix_,"X5")
            pb.x0[ix_["X5"]] = 2.7e+1
        else
            pb.y0[findfirst(x->x==ig_["X5"],pbm.congrps)]pb.y0[findfirst(x->x==ig_[2.7e+1],pbm.congrps)]
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2x_ii( "Y", iet_)
        loaset(elftv,it,1,"X1")
        loaset(elftv,it,2,"X2")
        loaset(elftv,it,3,"X3")
        loaset(elftv,it,4,"X4")
        loaset(elftv,it,5,"X5")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"PI")
        it,iet_,_ = s2x_ii( "C", iet_)
        loaset(elftv,it,1,"X1")
        loaset(elftv,it,2,"X2")
        loaset(elftv,it,3,"X3")
        loaset(elftv,it,4,"X4")
        loaset(elftv,it,5,"X5")
        loaset(elftp,it,1,"PI")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["19"])
            v_["PI"] = I
            v_["PI"] = 0.01+v_["PI"]
            ename = "C"*string(I)
            ie,ie_,_  = s2x_ii(ename,ie_)
            arrset(pbm.elftype,ie,"C")
            arrset(ielftype, ie, iet_["C"])
            vname = "X1"
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X2"
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X3"
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X3",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X4"
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X4",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X5"
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X5",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="PI",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,v_["PI"])
        end
        for I = Int64(v_["1"]):Int64(v_["20"])
            v_["PI"] = I
            v_["PI"] = 0.01+v_["PI"]
            ename = "Y"*string(I)
            ie,ie_,_  = s2x_ii(ename,ie_)
            arrset(pbm.elftype,ie,"Y")
            arrset(ielftype, ie, iet_["Y"])
            vname = "X1"
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X2"
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X3"
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X3",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X4"
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X4",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X5"
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X5",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="PI",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,v_["PI"])
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["OBJ"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["Y17"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,-5.843e-7)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["Y14"])
        loaset(pbm.grelw,ig,posel,1.17e-4)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["Y13"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,2.358e-5)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["Y16"])
        loaset(pbm.grelw,ig,posel,1.502e-6)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["Y12"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,0.0321)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["Y5"])
        loaset(pbm.grelw,ig,posel,0.00423)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["C18"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.0e-4)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["C19"])
        loaset(pbm.grelw,ig,posel,37.48)
        for I = Int64(v_["1"]):Int64(v_["20"])
            ig = ig_["CON"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["Y"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+pb.neq+1:pb.m] = zeros(Float64,pb.nge)
        pb.cupper[1:pb.nge] = grange[gegrps]
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTIONS %%%%%%%
        lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "OOI2-MN-5-21"
        return pb, pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%
