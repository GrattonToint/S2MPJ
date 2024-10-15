function HS70(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS70
#    *********
# 
#    This problem arises in water flow routing.
# 
#    Source: problem 70 incorrectly stated in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: Nick Gould, August 1991, modified May 2024
# 
#    classification = "C-SQR2-MN-4-1"
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "HS70"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["N"] = 4
        v_["1"] = 1
        v_["2"] = 2
        v_["19"] = 19
        v_["C1"] = 0.1
        v_["C2"] = 1.0
        v_["C3"] = 2.0
        v_["C4"] = 3.0
        v_["C5"] = 4.0
        v_["C6"] = 5.0
        v_["C7"] = 6.0
        v_["C8"] = 7.0
        v_["C9"] = 8.0
        v_["C10"] = 9.0
        v_["C11"] = 10.0
        v_["C12"] = 11.0
        v_["C13"] = 12.0
        v_["C14"] = 13.0
        v_["C15"] = 14.0
        v_["C16"] = 15.0
        v_["C17"] = 16.0
        v_["C18"] = 17.0
        v_["C19"] = 18.0
        v_["Y1"] = 0.00189
        v_["Y2"] = 0.1038
        v_["Y3"] = 0.268
        v_["Y4"] = 0.506
        v_["Y5"] = 0.577
        v_["Y6"] = 0.604
        v_["Y7"] = 0.725
        v_["Y8"] = 0.898
        v_["Y9"] = 0.947
        v_["Y10"] = 0.845
        v_["Y11"] = 0.702
        v_["Y12"] = 0.528
        v_["Y13"] = 0.385
        v_["Y14"] = 0.257
        v_["Y15"] = 0.159
        v_["Y16"] = 0.0869
        v_["Y17"] = 0.0453
        v_["Y18"] = 0.01509
        v_["Y19"] = 0.00189
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
        for I = Int64(v_["1"]):Int64(v_["19"])
            ig,ig_,_ = s2mpj_ii("OBJ"*string(I),ig_)
            arrset(gtype,ig,"<>")
        end
        ig,ig_,_ = s2mpj_ii("C1",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"C1")
        iv = ix_["X3"]
        pbm.A[ig,iv] += Float64(1.0e+0)
        iv = ix_["X4"]
        pbm.A[ig,iv] += Float64(1.0e+0)
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
        for I = Int64(v_["1"]):Int64(v_["19"])
            pbm.gconst[ig_["OBJ"*string(I)]] = Float64(v_["Y"*string(I)])
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = fill(0.00001,pb.n)
        pb.xupper = fill(100.0,pb.n)
        pb.xupper[ix_["X3"]] = 1.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        if haskey(ix_,"X1")
            pb.x0[ix_["X1"]] = Float64(2.0)
        else
            pb.y0[findfirst(x->x==ig_["X1"],pbm.congrps)] = Float64(2.0)
        end
        if haskey(ix_,"X2")
            pb.x0[ix_["X2"]] = Float64(4.0)
        else
            pb.y0[findfirst(x->x==ig_["X2"],pbm.congrps)] = Float64(4.0)
        end
        if haskey(ix_,"X3")
            pb.x0[ix_["X3"]] = Float64(0.04)
        else
            pb.y0[findfirst(x->x==ig_["X3"],pbm.congrps)] = Float64(0.04)
        end
        if haskey(ix_,"X4")
            pb.x0[ix_["X4"]] = Float64(2.0)
        else
            pb.y0[findfirst(x->x==ig_["X4"],pbm.congrps)] = Float64(2.0)
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eY1", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        loaset(elftv,it,3,"V3")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"C")
        it,iet_,_ = s2mpj_ii( "eY2", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        loaset(elftv,it,3,"V3")
        loaset(elftp,it,1,"C")
        it,iet_,_ = s2mpj_ii( "ePROD", iet_)
        loaset(elftv,it,1,"X3")
        loaset(elftv,it,2,"X4")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["19"])
            ename = "Y"*string(I)*","*string(Int64(v_["1"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eY1")
            arrset(ielftype,ie,iet_["eY1"])
            ename = "Y"*string(I)*","*string(Int64(v_["1"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X2"
            iv,ix_,pb  = (
                  s2mpj_nlx(vname,ix_,pb,1,Float64(0.00001),Float64(100.0),nothing))
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "Y"*string(I)*","*string(Int64(v_["1"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X3"
            iv,ix_,pb  = (
                  s2mpj_nlx(vname,ix_,pb,1,Float64(0.00001),Float64(100.0),nothing))
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "Y"*string(I)*","*string(Int64(v_["1"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X4"
            iv,ix_,pb  = (
                  s2mpj_nlx(vname,ix_,pb,1,Float64(0.00001),Float64(100.0),nothing))
            posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "Y"*string(I)*","*string(Int64(v_["1"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            posep = findfirst(x->x=="C",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["C"*string(I)]))
            ename = "Y"*string(I)*","*string(Int64(v_["2"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eY2")
            arrset(ielftype,ie,iet_["eY2"])
            ename = "Y"*string(I)*","*string(Int64(v_["2"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X1"
            iv,ix_,pb  = (
                  s2mpj_nlx(vname,ix_,pb,1,Float64(0.00001),Float64(100.0),nothing))
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "Y"*string(I)*","*string(Int64(v_["2"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X3"
            iv,ix_,pb  = (
                  s2mpj_nlx(vname,ix_,pb,1,Float64(0.00001),Float64(100.0),nothing))
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "Y"*string(I)*","*string(Int64(v_["2"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X4"
            iv,ix_,pb  = (
                  s2mpj_nlx(vname,ix_,pb,1,Float64(0.00001),Float64(100.0),nothing))
            posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "Y"*string(I)*","*string(Int64(v_["2"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            posep = findfirst(x->x=="C",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["C"*string(I)]))
        end
        ename = "C1"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD")
        arrset(ielftype,ie,iet_["ePROD"])
        vname = "X3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.00001),Float64(100.0),nothing))
        posev = findfirst(x->x=="X3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X4"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.00001),Float64(100.0),nothing))
        posev = findfirst(x->x=="X4",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gSQR",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["19"])
            ig = ig_["OBJ"*string(I)]
            arrset(pbm.grftype,ig,"gSQR")
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["Y"*string(I)*","*string(Int64(v_["1"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["Y"*string(I)*","*string(Int64(v_["2"]))])
            loaset(pbm.grelw,ig,posel, 1.)
        end
        ig = ig_["C1"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["C1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               0.007498464
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+pb.neq+1:pb.m] = zeros(Float64,pb.nge)
        pb.cupper[1:pb.nge] = fill(Inf,pb.nge)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-SQR2-MN-4-1"
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
        arrset(pbm.efpar,1,sqrt(1.0e+0/6.2832e+0))
        return pbm

    elseif action == "eY1"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        B = EV_[2]+EV_[3]*(1.0e+0-EV_[2])
        CI = pbm.elpar[iel_][1]/7.658e+0
        P0 = 1.0e+0+1.0e+0/(1.2e+1*EV_[1])
        P0V1 = -1.0e+0/(1.2e+1*EV_[1]^2)
        P0V1V1 = 2.0e+0/(1.2e+1*EV_[1]^3)
        P1 = 1.0e+0/P0
        P1V1 = -P0V1/(P0^2)
        P1V1V1 = (2.0e+0*P0V1^2/P0-P0V1V1)/(P0^2)
        P2 = EV_[2]
        P3 = B^EV_[1]
        LOGB = log(B)
        P3V1 = P3*LOGB
        P3V2 = EV_[1]*(1.0e+0-EV_[3])*B^(EV_[1]-1.0e+0)
        P3V3 = EV_[1]*(1.0e+0-EV_[2])*B^(EV_[1]-1.0e+0)
        P3V1V1 = P3V1*LOGB
        P3V1V2 = P3V2*LOGB+P3*(1.0e+0-EV_[3])/B
        P3V1V3 = P3V3*LOGB+P3*(1.0e+0-EV_[2])/B
        P3V2V2 = EV_[1]*(EV_[1]-1.0e+0)*(1.0e+0-EV_[3])^2*B^(EV_[1]-1.0e+0)
        P3V2V3  = (
              -EV_[1]*B^(EV_[1]-1.0e+0)+EV_[1]*(EV_[1]-1.0e+0)*(1.0e+0-EV_[2])*(1.0e+0-EV_[3])*B^(EV_[1]-2.0e+0))
        P3V3V3 = EV_[1]*(EV_[1]-1.0e+0)*(1.0e+0-EV_[2])^2*B^(EV_[1]-2.0e+0)
        P4 = pbm.efpar[1]*sqrt(EV_[1])
        P4V1 = 5.0e-1*pbm.efpar[1]*sqrt(1.0e+0/EV_[1])
        P4V1V1 = -2.5e-1*pbm.efpar[1]*sqrt(1.0e+0/EV_[1]^3)
        C5 = CI^(-1.0e+0)
        P5 = C5*CI^EV_[1]
        P5V1 = P5*log(CI)
        P5V1V1 = P5V1*log(CI)
        P6 = exp(EV_[1]*(1.0e+0-CI*B))
        P6V1 = P6*(1.0e+0-CI*B)
        P6V2 = -P6*EV_[1]*CI*(1.0e+0-EV_[3])
        P6V3 = -P6*EV_[1]*CI*(1.0e+0-EV_[2])
        P6V1V1 = P6*(1.0e+0-CI*B)^2
        P6V1V2 = P6V2*(1.0e+0-CI*B)-P6*CI*(1.0e+0-EV_[3])
        P6V1V3 = P6V3*(1.0e+0-CI*B)-P6*CI*(1.0e+0-EV_[2])
        P6V2V2 = -P6V2*EV_[1]*CI*(1.0e+0-EV_[3])
        P6V2V3 = -P6V3*EV_[1]*CI*(1.0e+0-EV_[3])+P6*EV_[1]*CI
        P6V3V3 = -P6V3*EV_[1]*CI*(1.0e+0-EV_[2])
        f_   = P1*P2*P3*P4*P5*P6
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = (P1V1*P2*P3*P4*P5*P6+P1*P2*P3V1*P4*P5*P6+P1*P2*P3*P4V1*P5*P6+
                 P1*P2*P3*P4*P5V1*P6+P1*P2*P3*P4*P5*P6V1)
            g_[2] = P1*P3*P4*P5*P6+P1*P2*P3V2*P4*P5*P6+P1*P2*P3*P4*P5*P6V2
            g_[3] = P1*P2*P3V3*P4*P5*P6+P1*P2*P3*P4*P5*P6V3
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,1] = (P1V1V1*P2*P3*P4*P5*P6+P1*P2*P3V1V1*P4*P5*P6+P1*P2*P3*P4V1V1*P5*P6+
                     P1*P2*P3*P4*P5V1V1*P6+P1*P2*P3*P4*P5*P6V1V1+2.0e+0*(P1V1*P2*P3V1*P4*P5*P6+P1V1*P2*P3*P4V1*P5*P6+P1V1*P2*P3*P4*P5V1*P6+P1V1*P2*P3*P4*P5*P6V1+P1*P2*P3V1*P4V1*P5*P6+P1*P2*P3V1*P4*P5V1*P6+P1*P2*P3V1*P4*P5*P6V1+P1*P2*P3*P4V1*P5V1*P6+P1*P2*P3*P4V1*P5*P6V1+P1*P2*P3*P4*P5V1*P6V1))
                H_[1,2]  = (
                      P1V1*(P3*P4*P5*P6+P2*P3V2*P4*P5*P6+P2*P3*P4*P5*P6V2)+P1*(P3V1*P4*P5*P6+(P3*P4V1*P5*P6+P3*P4*P5V1*P6+P3*P4*P5*P6V1)+P2*(P3V1V2*P4*P5*P6+P3V1*P4*P5*P6V2+P3V2*P4V1*P5*P6+P3V2*P4*P5V1*P6+P3V2*P4*P5*P6V1+P3*(P4V1*P5*P6V2+P4*P5V1*P6V2+P4*P5*P6V1V2))))
                H_[2,1] = H_[1,2]
                H_[1,3]  = (
                      P2*(P1V1*P3V3*P4*P5*P6+P1V1*P3*P4*P5*P6V3+P1*(P3V1V3*P4*P5*P6+P3V1*P4*P5*P6V3+P3V3*P4V1*P5*P6+P3V3*P4*P5V1*P6+P3V3*P4*P5*P6V1+P3*(P4*P5V1*P6V3+P4V1*P5*P6V3+P4*P5*P6V1V3))))
                H_[3,1] = H_[1,3]
                H_[2,2]  = (
                      P1*P4*P5*(P2*P3*P6V2V2+P2*P3V2V2*P6+2.0e+0*(P2*P3V2*P6V2+P3*P6V2+P3V2*P6)))
                H_[2,3]  = (
                      P1*P4*P5*(P2*P3V2V3*P6+P2*P3*P6V2V3+P3V3*P6+P3*P6V3+P2*P3V2*P6V3+P2*P3V3*P6V2))
                H_[3,2] = H_[2,3]
                H_[3,3] = P1*P2*P4*P5*(P3V3V3*P6+P3*P6V3V3+2.0e+0*P3V3*P6V3)
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eY2"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        B = EV_[2]+EV_[3]*(1.0e+0-EV_[2])
        CI = pbm.elpar[iel_][1]/7.658e+0
        P0 = 1.0e+0+1.0e+0/(1.2e+1*EV_[1])
        P0V1 = -1.0e+0/(1.2e+1*EV_[1]^2)
        P0V1V1 = 2.0e+0/(1.2e+1*EV_[1]^3)
        P1 = 1.0e+0/P0
        P1V1 = -P0V1/(P0^2)
        P1V1V1 = (2.0e+0*P0V1^2/P0-P0V1V1)/(P0^2)
        P2 = 1.0e+0-EV_[2]
        P3 = (B/EV_[3])^EV_[1]
        LOGB = log(B/EV_[3])
        P3V1 = P3*LOGB
        P3V2 = EV_[1]*(-1.0e+0+1.0e+0/EV_[3])*(B/EV_[3])^(EV_[1]-1.0e+0)
        P3V3 = -EV_[1]*(EV_[2]/EV_[3]^2)*(B/EV_[3])^(EV_[1]-1.0e+0)
        P3V1V1 = P3V1*LOGB
        P3V1V2 = P3V2*LOGB+P3*EV_[3]*(-1.0e+0+1.0e+0/EV_[3])/B
        P3V1V3 = P3V3*LOGB-P3*EV_[2]/(B*EV_[3])
        P3V2V2  = (
              EV_[1]*(EV_[1]-1.0e+0)*(-1.0e+0+1.0e+0/EV_[3])^2*(B/EV_[3])^(EV_[1]-2.0e+0))
        P3V2V3  = (
              EV_[1]*(-1.0e+0/EV_[3]^2)*(B/EV_[3])^(EV_[1]-1.0e+0)+EV_[1]*(EV_[1]-1.0e+0)*(-1.0e+0+1.0e+0/EV_[3])*(-EV_[2]/EV_[3]^2)*(B/EV_[3])^(EV_[1]-2.0e+0))
        P3V3V3 = (2.0e+0*EV_[1]*(EV_[2]/EV_[3]^3)*(B/EV_[3])^(EV_[1]-1.0e+0)+
             EV_[1]*(EV_[1]-1.0e+0)*(EV_[2]/EV_[3]^2)^2*(B/EV_[3])^(EV_[1]-2.0e+0))
        P4 = pbm.efpar[1]*sqrt(EV_[1])
        P4V1 = 5.0e-1*pbm.efpar[1]*sqrt(1.0e+0/EV_[1])
        P4V1V1 = -2.5e-1*pbm.efpar[1]*sqrt(1.0e+0/EV_[1]^3)
        C5 = CI^(-1.0e+0)
        P5 = C5*CI^EV_[1]
        P5V1 = P5*log(CI)
        P5V1V1 = P5V1*log(CI)
        P6 = exp(EV_[1]*(1.0e+0-CI*B/EV_[3]))
        P6V1 = P6*(1.0e+0-CI*B/EV_[3])
        P6V2 = -P6*EV_[1]*CI*(1.0e+0-EV_[3])/EV_[3]
        P6V3 = P6*EV_[1]*CI*EV_[2]/EV_[3]^2
        P6V1V1 = P6*(1.0e+0-CI*B/EV_[3])^2
        P6V1V2 = P6V2*(1.0e+0-CI*B/EV_[3])-P6*CI*(-1.0e+0+1.0e+0/EV_[3])
        P6V1V3 = P6V3*(1.0e+0-CI*B/EV_[3])+P6*CI*EV_[2]/EV_[3]^2
        P6V2V2 = -P6V2*EV_[1]*CI*(1.0e+0-EV_[3])/EV_[3]
        P6V2V3 = -P6V3*EV_[1]*CI*(1.0e+0-EV_[3])/EV_[3]+P6*EV_[1]*CI/EV_[3]^2
        P6V3V3 = P6V3*EV_[1]*CI*EV_[2]/EV_[3]^2-2.0e+0*P6*EV_[1]*CI*EV_[2]/EV_[3]^3
        f_   = P1*P2*P3*P4*P5*P6
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = (P1V1*P2*P3*P4*P5*P6+P1*P2*P3V1*P4*P5*P6+P1*P2*P3*P4V1*P5*P6+
                 P1*P2*P3*P4*P5V1*P6+P1*P2*P3*P4*P5*P6V1)
            g_[2] = -P1*P3*P4*P5*P6+P1*P2*P3V2*P4*P5*P6+P1*P2*P3*P4*P5*P6V2
            g_[3] = P1*P2*P3V3*P4*P5*P6+P1*P2*P3*P4*P5*P6V3
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,1] = (P1V1V1*P2*P3*P4*P5*P6+P1*P2*P3V1V1*P4*P5*P6+P1*P2*P3*P4V1V1*P5*P6+
                     P1*P2*P3*P4*P5V1V1*P6+P1*P2*P3*P4*P5*P6V1V1+2.0e+0*(P1V1*P2*P3V1*P4*P5*P6+P1V1*P2*P3*P4V1*P5*P6+P1V1*P2*P3*P4*P5V1*P6+P1V1*P2*P3*P4*P5*P6V1+P1*P2*P3V1*P4V1*P5*P6+P1*P2*P3V1*P4*P5V1*P6+P1*P2*P3V1*P4*P5*P6V1+P1*P2*P3*P4V1*P5V1*P6+P1*P2*P3*P4V1*P5*P6V1+P1*P2*P3*P4*P5V1*P6V1))
                H_[1,2]  = (
                      P1V1*(-P3*P4*P5*P6+P2*P3V2*P4*P5*P6+P2*P3*P4*P5*P6V2)+P1*(-P3V1*P4*P5*P6-(P3*P4V1*P5*P6+P3*P4*P5V1*P6+P3*P4*P5*P6V1)+P2*(P3V1V2*P4*P5*P6+P3V1*P4*P5*P6V2+P3V2*P4V1*P5*P6+P3V2*P4*P5V1*P6+P3V2*P4*P5*P6V1+P3*(P4V1*P5*P6V2+P4*P5V1*P6V2+P4*P5*P6V1V2))))
                H_[2,1] = H_[1,2]
                H_[1,3]  = (
                      P2*(P1V1*P3V3*P4*P5*P6+P1V1*P3*P4*P5*P6V3+P1*(P3V1V3*P4*P5*P6+P3V1*P4*P5*P6V3+P3V3*P4V1*P5*P6+P3V3*P4*P5V1*P6+P3V3*P4*P5*P6V1+P3*(P4*P5V1*P6V3+P4V1*P5*P6V3+P4*P5*P6V1V3))))
                H_[3,1] = H_[1,3]
                H_[2,2]  = (
                      P1*P4*P5*(P2*P3*P6V2V2+P2*P3V2V2*P6+2.0e+0*(P2*P3V2*P6V2-P3*P6V2-P3V2*P6)))
                H_[2,3]  = (
                      P1*P4*P5*(P2*P3V2V3*P6+P2*P3*P6V2V3-P3V3*P6-P3*P6V3+P2*P3V2*P6V3+P2*P3V3*P6V2))
                H_[3,2] = H_[2,3]
                H_[3,3] = P1*P2*P4*P5*(P3V3V3*P6+P3*P6V3V3+2.0e+0*P3V3*P6V3)
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
                H_[1,2] = 1.0e+0
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

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    elseif action == "gSQR"

        GVAR_   = args[1]
        igr_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_= GVAR_*GVAR_
        if nargout>1
            g_ = GVAR_+GVAR_
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_ = 2.0e+0
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
            pbm.has_globs = [1,0]
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

