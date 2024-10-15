function DNIEPER(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DNIEPER
#    *********
# 
#    This problem models the planning of systematic use of water resources
#    in the basin of the river Dnieper.
# 
#    Source: p. 139sq in 
#    B.N. Pshenichnyj
#    "The Linearization Method for Constrained Optimization",
#    Springer Verlag, SCM Series 22, Heidelberg, 1994
# 
#    SIF input: Ph. Toint, December 1994.
# 
#    classification = "C-QOR2-MN-61-24"
# 
#    Problem data
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "DNIEPER"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["C1"] = 5.61
        v_["C2"] = 4.68
        v_["C3"] = 1.62
        v_["C4"] = 1.8
        v_["C5"] = 2.13
        v_["C6"] = 2.1
        v_["C7"] = 1.99
        v_["C8"] = 2.02
        v_["C9"] = 2.14
        v_["C10"] = 2.15
        v_["C11"] = 2.36
        v_["C12"] = 2.63
        v_["C13"] = -0.02
        v_["C14"] = -0.01
        v_["C15"] = -0.16
        v_["C16"] = -0.47
        v_["C17"] = -0.75
        v_["C18"] = -0.94
        v_["C19"] = -0.93
        v_["C20"] = -0.99
        v_["C21"] = -0.42
        v_["C22"] = -0.07
        v_["C23"] = 0.04
        v_["C24"] = -0.06
        v_["1"] = 1
        v_["2"] = 2
        v_["4"] = 4
        v_["5"] = 5
        v_["8"] = 8
        v_["9"] = 9
        v_["12"] = 12
        v_["13"] = 13
        v_["14"] = 14
        v_["16"] = 16
        v_["17"] = 17
        v_["20"] = 20
        v_["21"] = 21
        v_["24"] = 24
        v_["25"] = 25
        v_["36"] = 36
        v_["37"] = 37
        v_["48"] = 48
        v_["49"] = 49
        v_["52"] = 52
        v_["53"] = 53
        v_["56"] = 56
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["56"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        iv,ix_,_ = s2mpj_ii("X0F",ix_)
        arrset(pb.xnames,iv,"X0F")
        iv,ix_,_ = s2mpj_ii("X24F",ix_)
        arrset(pb.xnames,iv,"X24F")
        iv,ix_,_ = s2mpj_ii("X12F",ix_)
        arrset(pb.xnames,iv,"X12F")
        iv,ix_,_ = s2mpj_ii("X36F",ix_)
        arrset(pb.xnames,iv,"X36F")
        iv,ix_,_ = s2mpj_ii("AC",ix_)
        arrset(pb.xnames,iv,"AC")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["1"]):Int64(v_["12"])
            v_["I0"] = 12+I
            v_["I1"] = 24+I
            v_["I2"] = 36+I
            ig,ig_,_ = s2mpj_ii("OBJ",ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X"*string(Int64(v_["I1"]))]
            pbm.A[ig,iv] += Float64(19.95)
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += Float64(0.07656)
            iv = ix_["X"*string(Int64(v_["I2"]))]
            pbm.A[ig,iv] += Float64(-24.89)
            iv = ix_["X"*string(Int64(v_["I0"]))]
            pbm.A[ig,iv] += Float64(-0.7135)
        end
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        arrset(pbm.gscale,ig,Float64(-1.0))
        for I = Int64(v_["1"]):Int64(v_["4"])
            v_["I0"] = 24+I
            ig,ig_,_ = s2mpj_ii("CC"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"CC"*string(I))
            iv = ix_["X"*string(Int64(v_["I0"]))]
            pbm.A[ig,iv] += Float64(-2.68)
        end
        for I = Int64(v_["5"]):Int64(v_["8"])
            v_["I0"] = 24+I
            v_["I1"] = 44+I
            ig,ig_,_ = s2mpj_ii("CC"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"CC"*string(I))
            iv = ix_["X"*string(Int64(v_["I0"]))]
            pbm.A[ig,iv] += Float64(-2.68)
            iv = ix_["X"*string(Int64(v_["I1"]))]
            pbm.A[ig,iv] += Float64(-2.68)
        end
        for I = Int64(v_["9"]):Int64(v_["12"])
            v_["I0"] = 24+I
            ig,ig_,_ = s2mpj_ii("CC"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"CC"*string(I))
            iv = ix_["X"*string(Int64(v_["I0"]))]
            pbm.A[ig,iv] += Float64(-2.68)
        end
        for I = Int64(v_["13"]):Int64(v_["16"])
            v_["I0"] = 12+I
            v_["I1"] = 24+I
            ig,ig_,_ = s2mpj_ii("CC"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"CC"*string(I))
            iv = ix_["X"*string(Int64(v_["I0"]))]
            pbm.A[ig,iv] += Float64(-2.68)
            iv = ix_["X"*string(Int64(v_["I1"]))]
            pbm.A[ig,iv] += Float64(-2.68)
            iv = ix_["AC"]
            pbm.A[ig,iv] += Float64(-1.0)
        end
        for I = Int64(v_["17"]):Int64(v_["20"])
            v_["I0"] = 12+I
            v_["I1"] = 24+I
            v_["I2"] = 36+I
            ig,ig_,_ = s2mpj_ii("CC"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"CC"*string(I))
            iv = ix_["X"*string(Int64(v_["I0"]))]
            pbm.A[ig,iv] += Float64(-2.68)
            iv = ix_["X"*string(Int64(v_["I1"]))]
            pbm.A[ig,iv] += Float64(-2.68)
            iv = ix_["X"*string(Int64(v_["I2"]))]
            pbm.A[ig,iv] += Float64(-2.68)
            iv = ix_["AC"]
            pbm.A[ig,iv] += Float64(-1.0)
        end
        for I = Int64(v_["21"]):Int64(v_["24"])
            v_["I0"] = 12+I
            v_["I1"] = 24+I
            ig,ig_,_ = s2mpj_ii("CC"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"CC"*string(I))
            iv = ix_["X"*string(Int64(v_["I0"]))]
            pbm.A[ig,iv] += Float64(-2.68)
            iv = ix_["X"*string(Int64(v_["I1"]))]
            pbm.A[ig,iv] += Float64(-2.68)
            iv = ix_["AC"]
            pbm.A[ig,iv] += Float64(-1.0)
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
        pbm.gconst[ig_["OBJ"]] = Float64(-112.464)
        for I = Int64(v_["1"]):Int64(v_["24"])
            v_["CST"] = -1.0*v_["C"*string(I)]
            pbm.gconst[ig_["CC"*string(I)]] = Float64(v_["CST"])
        end
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        for I = Int64(v_["1"]):Int64(v_["12"])
            pb.xlower[ix_["X"*string(I)]] = 51.2
            pb.xupper[ix_["X"*string(I)]] = 51.4
        end
        for I = Int64(v_["13"]):Int64(v_["24"])
            pb.xlower[ix_["X"*string(I)]] = 15.0
            pb.xupper[ix_["X"*string(I)]] = 16.1
        end
        for I = Int64(v_["25"]):Int64(v_["36"])
            pb.xlower[ix_["X"*string(I)]] = 0.4
            pb.xupper[ix_["X"*string(I)]] = 4.6
        end
        for I = Int64(v_["37"]):Int64(v_["48"])
            pb.xlower[ix_["X"*string(I)]] = 0.5
            pb.xupper[ix_["X"*string(I)]] = 4.8
        end
        for I = Int64(v_["49"]):Int64(v_["56"])
            pb.xlower[ix_["X"*string(I)]] = 0.0
            pb.xupper[ix_["X"*string(I)]] = 0.7
        end
        pb.xlower[ix_["X0F"]] = 50.82
        pb.xupper[ix_["X0F"]] = 50.82
        pb.xlower[ix_["X24F"]] = 2.0
        pb.xupper[ix_["X24F"]] = 2.0
        pb.xlower[ix_["X12F"]] = 15.5
        pb.xupper[ix_["X12F"]] = 15.5
        pb.xlower[ix_["X36F"]] = 2.3
        pb.xupper[ix_["X36F"]] = 2.3
        pb.xlower[ix_["AC"]] = -Inf
        pb.xupper[ix_["AC"]] = +Inf
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        for I = Int64(v_["1"]):Int64(v_["12"])
            pb.x0[ix_["X"*string(I)]] = Float64(51.35)
        end
        for I = Int64(v_["13"]):Int64(v_["24"])
            pb.x0[ix_["X"*string(I)]] = Float64(15.5)
        end
        for I = Int64(v_["25"]):Int64(v_["36"])
            pb.x0[ix_["X"*string(I)]] = Float64(2.5)
        end
        for I = Int64(v_["37"]):Int64(v_["48"])
            pb.x0[ix_["X"*string(I)]] = Float64(2.6)
        end
        for I = Int64(v_["49"]):Int64(v_["56"])
            pb.x0[ix_["X"*string(I)]] = Float64(0.3)
        end
        pb.x0[ix_["X0F"]] = Float64(50.82)
        pb.x0[ix_["X24F"]] = Float64(2.0)
        pb.x0[ix_["X12F"]] = Float64(15.5)
        pb.x0[ix_["X36F"]] = Float64(2.3)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eWJ", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        it,iet_,_ = s2mpj_ii( "eWK", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        it,iet_,_ = s2mpj_ii( "en2PR", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        it,iet_,_ = s2mpj_ii( "eSQ", iet_)
        loaset(elftv,it,1,"X")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["12"])
            v_["I1"] = 12+I
            v_["I2"] = 36+I
            ename = "E"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
            vname = "X"*string(Int64(v_["I1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["I2"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        ename = "ACSQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        vname = "AC"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        for I = Int64(v_["1"]):Int64(v_["12"])
            v_["I0"] = 24+I
            ename = "W1"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eWJ")
            arrset(ielftype,ie,iet_["eWJ"])
            vname = "X"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["I0"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        ename = "W2"*string(Int64(v_["1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eWJ")
        arrset(ielftype,ie,iet_["eWJ"])
        ename = "W2"*string(Int64(v_["1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X0F"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "W2"*string(Int64(v_["1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X24F"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        for I = Int64(v_["2"]):Int64(v_["12"])
            v_["I0"] = 23+I
            v_["I1"] = -1+I
            ename = "W2"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eWJ")
            arrset(ielftype,ie,iet_["eWJ"])
            vname = "X"*string(Int64(v_["I1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["I0"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        for I = Int64(v_["13"]):Int64(v_["24"])
            v_["I1"] = 24+I
            ename = "W1"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eWK")
            arrset(ielftype,ie,iet_["eWK"])
            vname = "X"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["I1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        ename = "W2"*string(Int64(v_["13"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eWK")
        arrset(ielftype,ie,iet_["eWK"])
        ename = "W2"*string(Int64(v_["13"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X12F"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "W2"*string(Int64(v_["13"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X36F"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        for I = Int64(v_["14"]):Int64(v_["24"])
            v_["I0"] = 23+I
            v_["I1"] = -1+I
            ename = "W2"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eWK")
            arrset(ielftype,ie,iet_["eWK"])
            vname = "X"*string(Int64(v_["I1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["I0"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["12"])
            ig = ig_["OBJ"]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["E"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(2.155))
        end
        ig = ig_["OBJ"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["ACSQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-2000.0))
        for I = Int64(v_["1"]):Int64(v_["24"])
            ig = ig_["CC"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["W1"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["W2"*string(I)])
            loaset(pbm.grelw,ig,posel,Float64(-1.0))
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLUTION             1.87439D+04
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
        pb.pbclass = "C-QOR2-MN-61-24"
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

    elseif action == "eWJ"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        A1 = 34.547
        A2 = -0.55878
        A3 = 8.05339
        A4 = -0.02252
        A5 = -0.29316
        A6 = -0.013521
        A7 = 0.00042
        A8 = 0.00267
        A9 = 0.000281
        A10 = 0.0000032
        f_   = (A1+A2*EV_[1]+A3*EV_[2]+A4*EV_[1]^2+A5*EV_[1]*EV_[2]+A6*EV_[2]^2+
             A7*EV_[1]^3+A8*EV_[1]^2*EV_[2]+A9*EV_[1]*EV_[2]^2+A10*EV_[2]^3)
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = (A2+2.0*A4*EV_[1]+A5*EV_[2]+3.0*A7*EV_[1]^2+2.0*A8*EV_[1]*EV_[2]+
                 A9*EV_[2]^2)
            g_[2] = (A3+A5*EV_[1]+2.0*A6*EV_[2]+A8*EV_[1]^2+2.0*A9*EV_[1]*EV_[2]+
                 3.0*A10*EV_[2]^2)
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = 2.0*A4+6.0*A7*EV_[1]+2.0*A8*EV_[2]
                H_[1,2] = A5+2.0*A8*EV_[1]+2.0*A9*EV_[2]
                H_[2,1] = H_[1,2]
                H_[2,2] = 2.0*A6+2.0*A9*EV_[1]+6.0*A10*EV_[2]
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eWK"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        A1 = 20.923
        A2 = -4.22088
        A3 = 1.42061
        A4 = -0.41040
        A5 = -0.15082
        A7 = -0.00826
        A8 = 0.00404
        A9 = 0.000168
        A10 = -0.000038
        f_   = (A1+A2*EV_[1]+A3*EV_[2]+A4*EV_[1]^2+A5*EV_[1]*EV_[2]+A7*EV_[1]^3+
             A8*EV_[1]^2*EV_[2]+A9*EV_[1]*EV_[2]^2+A10*EV_[2]^3)
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = (A2+2.0*A4*EV_[1]+A5*EV_[2]+3.0*A7*EV_[1]^2+2.0*A8*EV_[1]*EV_[2]+
                 A9*EV_[2]^2)
            g_[2] = A3+A5*EV_[1]+A8*EV_[1]^2+2.0*A9*EV_[1]*EV_[2]+3.0*A10*EV_[2]^2
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = 2.0*A4+6.0*A7*EV_[1]+2.0*A8*EV_[2]
                H_[1,2] = A5+2.0*A8*EV_[1]+2.0*A9*EV_[2]
                H_[2,1] = H_[1,2]
                H_[2,2] = 2.0*A9*EV_[1]+6.0*A10*EV_[2]
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "en2PR"

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

