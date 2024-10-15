function LIN(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : LIN
#    *********
# 
#    A non-convex global optimization chemical equilibrium problem from the 
#    thesis of W.J. Lin.
#    It has a nonlinear objective and linear constraints.
# 
#    Source: illustrative example (section 4.6) in
#    C.M. McDonald and C.A. Floudas, "Global optimization for the phase 
#    and chemical equilibrium problem: application to the NRTL equation",
#    Computers & Chemical Engineering, (submitted), 1994.
# 
#    SIF input: Marcel Mongeau, 9 February 1994.
# 
#    classification = "C-OLR2-AY-4-2"
# 
#    PARAMETERS likely to be changed for different problems:
# 
#    Number of variable sets (# of phases)
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "LIN"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["P"] = 2
        v_["C"] = 2
        v_["1"] = 1
        v_["2"] = 2
        v_["TAU"*string(Int64(v_["1"]))*","*string(Int64(v_["2"]))] = 3.00498
        v_["TAU"*string(Int64(v_["2"]))*","*string(Int64(v_["1"]))] = 4.69071
        v_["ALF"*string(Int64(v_["1"]))*","*string(Int64(v_["2"]))] = 0.391965
        v_["ALF"*string(Int64(v_["2"]))*","*string(Int64(v_["1"]))] = 0.391965
        v_["INIT"*string(Int64(v_["1"]))] = 0.5
        v_["INIT"*string(Int64(v_["2"]))] = 0.5
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["C"])
            for K = Int64(v_["1"]):Int64(v_["P"])
                iv,ix_,_ = s2mpj_ii("X"*string(I)*","*string(K),ix_)
                arrset(pb.xnames,iv,"X"*string(I)*","*string(K))
            end
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        for I = Int64(v_["1"]):Int64(v_["C"])
            for K = Int64(v_["1"]):Int64(v_["P"])
                ig,ig_,_ = s2mpj_ii("MB"*string(I),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"MB"*string(I))
                iv = ix_["X"*string(I)*","*string(K)]
                pbm.A[ig,iv] += Float64(1.0)
            end
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
        for I = Int64(v_["1"]):Int64(v_["C"])
            pbm.gconst[ig_["MB"*string(I)]] = Float64(v_["INIT"*string(I)])
        end
        v_["ZERO"] = 0.0
        v_["ONE"] = 1.0
        for I = Int64(v_["1"]):Int64(v_["C"])
            v_["ALF"*string(I)*","*string(I)] = v_["ZERO"]
            v_["TAU"*string(I)*","*string(I)] = v_["ZERO"]
        end
        for I = Int64(v_["1"]):Int64(v_["C"])
            for J = Int64(v_["1"]):Int64(v_["C"])
                v_["MALF"] = -1.0*v_["ALF"*string(I)*","*string(J)]
                v_["PROD"] = v_["MALF"]*v_["TAU"*string(I)*","*string(J)]
                v_["G"*string(I)*","*string(J)] = exp(v_["PROD"])
            end
        end
        for I = Int64(v_["1"]):Int64(v_["C"])
            v_["G"*string(I)*","*string(I)] = v_["ONE"]
        end
        for I = Int64(v_["1"]):Int64(v_["C"])
            for J = Int64(v_["1"]):Int64(v_["C"])
                v_["M"*string(I)*","*string(J)]  = (
                      v_["G"*string(I)*","*string(J)]*v_["TAU"*string(I)*","*string(J)])
            end
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = fill(1.e-12,pb.n)
        pb.xupper = fill(Inf,pb.n)
        pb.xupper[ix_["X"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))]]  = (
              v_["INIT1"])
        pb.xupper[ix_["X"*string(Int64(v_["1"]))*","*string(Int64(v_["2"]))]]  = (
              v_["INIT1"])
        pb.xupper[ix_["X"*string(Int64(v_["2"]))*","*string(Int64(v_["1"]))]]  = (
              v_["INIT2"])
        pb.xupper[ix_["X"*string(Int64(v_["2"]))*","*string(Int64(v_["2"]))]]  = (
              v_["INIT2"])
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        pb.x0[ix_["X"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))]]  = (
              Float64(0.5))
        pb.x0[ix_["X"*string(Int64(v_["1"]))*","*string(Int64(v_["2"]))]]  = (
              Float64(0.0))
        pb.x0[ix_["X"*string(Int64(v_["2"]))*","*string(Int64(v_["1"]))]]  = (
              Float64(0.0))
        pb.x0[ix_["X"*string(Int64(v_["2"]))*","*string(Int64(v_["2"]))]]  = (
              Float64(0.5))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eXTAUG1", iet_)
        loaset(elftv,it,1,"Y1")
        loaset(elftv,it,2,"Y2")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"G11")
        loaset(elftp,it,2,"G12")
        loaset(elftp,it,3,"G21")
        loaset(elftp,it,4,"G22")
        loaset(elftp,it,5,"M11")
        loaset(elftp,it,6,"M12")
        loaset(elftp,it,7,"M21")
        loaset(elftp,it,8,"M22")
        it,iet_,_ = s2mpj_ii( "eXTAUG2", iet_)
        loaset(elftv,it,1,"Y1")
        loaset(elftv,it,2,"Y2")
        loaset(elftp,it,1,"G11")
        loaset(elftp,it,2,"G12")
        loaset(elftp,it,3,"G21")
        loaset(elftp,it,4,"G22")
        loaset(elftp,it,5,"M11")
        loaset(elftp,it,6,"M12")
        loaset(elftp,it,7,"M21")
        loaset(elftp,it,8,"M22")
        it,iet_,_ = s2mpj_ii( "eXLOGX", iet_)
        loaset(elftv,it,1,"X")
        it,iet_,_ = s2mpj_ii( "eXLOGXC", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y1")
        loaset(elftv,it,3,"Y2")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for K = Int64(v_["1"]):Int64(v_["P"])
            ename = "A"*string(Int64(v_["1"]))*","*string(K)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eXTAUG1")
            arrset(ielftype,ie,iet_["eXTAUG1"])
            ename = "A"*string(Int64(v_["1"]))*","*string(K)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X"*string(Int64(v_["1"]))*","*string(K)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,nothing)
            posev = findfirst(x->x=="Y1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "A"*string(Int64(v_["1"]))*","*string(K)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X"*string(Int64(v_["2"]))*","*string(K)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,nothing)
            posev = findfirst(x->x=="Y2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "A"*string(Int64(v_["1"]))*","*string(K)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            posep = findfirst(x->x=="G11",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["G"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))]))
            ename = "A"*string(Int64(v_["1"]))*","*string(K)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            posep = findfirst(x->x=="G12",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["G"*string(Int64(v_["1"]))*","*string(Int64(v_["2"]))]))
            ename = "A"*string(Int64(v_["1"]))*","*string(K)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            posep = findfirst(x->x=="G21",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["G"*string(Int64(v_["2"]))*","*string(Int64(v_["1"]))]))
            ename = "A"*string(Int64(v_["1"]))*","*string(K)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            posep = findfirst(x->x=="G22",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["G"*string(Int64(v_["2"]))*","*string(Int64(v_["2"]))]))
            ename = "A"*string(Int64(v_["1"]))*","*string(K)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            posep = findfirst(x->x=="M11",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["M"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))]))
            ename = "A"*string(Int64(v_["1"]))*","*string(K)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            posep = findfirst(x->x=="M12",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["M"*string(Int64(v_["1"]))*","*string(Int64(v_["2"]))]))
            ename = "A"*string(Int64(v_["1"]))*","*string(K)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            posep = findfirst(x->x=="M21",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["M"*string(Int64(v_["2"]))*","*string(Int64(v_["1"]))]))
            ename = "A"*string(Int64(v_["1"]))*","*string(K)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            posep = findfirst(x->x=="M22",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["M"*string(Int64(v_["2"]))*","*string(Int64(v_["2"]))]))
            ename = "A"*string(Int64(v_["2"]))*","*string(K)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eXTAUG2")
            arrset(ielftype,ie,iet_["eXTAUG2"])
            ename = "A"*string(Int64(v_["2"]))*","*string(K)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X"*string(Int64(v_["1"]))*","*string(K)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,nothing)
            posev = findfirst(x->x=="Y1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "A"*string(Int64(v_["2"]))*","*string(K)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X"*string(Int64(v_["2"]))*","*string(K)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,nothing)
            posev = findfirst(x->x=="Y2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "A"*string(Int64(v_["2"]))*","*string(K)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            posep = findfirst(x->x=="G11",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["G"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))]))
            ename = "A"*string(Int64(v_["2"]))*","*string(K)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            posep = findfirst(x->x=="G12",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["G"*string(Int64(v_["1"]))*","*string(Int64(v_["2"]))]))
            ename = "A"*string(Int64(v_["2"]))*","*string(K)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            posep = findfirst(x->x=="G21",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["G"*string(Int64(v_["2"]))*","*string(Int64(v_["1"]))]))
            ename = "A"*string(Int64(v_["2"]))*","*string(K)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            posep = findfirst(x->x=="G22",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["G"*string(Int64(v_["2"]))*","*string(Int64(v_["2"]))]))
            ename = "A"*string(Int64(v_["2"]))*","*string(K)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            posep = findfirst(x->x=="M11",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["M"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))]))
            ename = "A"*string(Int64(v_["2"]))*","*string(K)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            posep = findfirst(x->x=="M12",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["M"*string(Int64(v_["1"]))*","*string(Int64(v_["2"]))]))
            ename = "A"*string(Int64(v_["2"]))*","*string(K)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            posep = findfirst(x->x=="M21",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["M"*string(Int64(v_["2"]))*","*string(Int64(v_["1"]))]))
            ename = "A"*string(Int64(v_["2"]))*","*string(K)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            posep = findfirst(x->x=="M22",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["M"*string(Int64(v_["2"]))*","*string(Int64(v_["2"]))]))
        end
        for K = Int64(v_["1"]):Int64(v_["P"])
            for I = Int64(v_["1"]):Int64(v_["C"])
                ename = "B"*string(I)*","*string(K)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eXLOGX")
                arrset(ielftype,ie,iet_["eXLOGX"])
                vname = "X"*string(I)*","*string(K)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,nothing)
                posev = findfirst(x->x=="X",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
        end
        for K = Int64(v_["1"]):Int64(v_["P"])
            for I = Int64(v_["1"]):Int64(v_["C"])
                ename = "C"*string(I)*","*string(K)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eXLOGXC")
                arrset(ielftype,ie,iet_["eXLOGXC"])
                vname = "X"*string(I)*","*string(K)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,nothing)
                posev = findfirst(x->x=="X",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(Int64(v_["1"]))*","*string(K)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,nothing)
                posev = findfirst(x->x=="Y1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(Int64(v_["2"]))*","*string(K)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.e-12),nothing,nothing)
                posev = findfirst(x->x=="Y2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for K = Int64(v_["1"]):Int64(v_["P"])
            ig = ig_["OBJ"]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["A"*string(Int64(v_["1"]))*","*string(K)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["A"*string(Int64(v_["2"]))*","*string(K)])
            loaset(pbm.grelw,ig,posel, 1.)
            for I = Int64(v_["1"]):Int64(v_["C"])
                ig = ig_["OBJ"]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B"*string(I)*","*string(K)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,1.)
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["C"*string(I)*","*string(K)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(-1.0))
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# Global minimum:        -0.02020 : 
#  XV LIN       X(1,1)    0.00071
#  XV LIN       X(1,2)    0.49929
#  XV LIN       X(2,1)    0.15588
#  XV LIN       X(2,2)    0.34412
# local minimum:         -0.01961 :
#  XV LIN       X(1,1)    0.00213
#  XV LIN       X(1,2)    0.49787
#  XV LIN       X(2,1)    0.46547
#  XV LIN       X(2,2)    0.03453
# local maximum:         -0.01730 :
#  XV LIN       X(1,1)    0.00173
#  XV LIN       X(1,2)    0.49827
#  XV LIN       X(2,1)    0.37544
#  XV LIN       X(2,2)    0.12456
# LO SOLTN               -0.02020
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
        pb.pbclass = "C-OLR2-AY-4-2"
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

    elseif action == "eXTAUG1"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        INSUM1 = EV_[1]*pbm.elpar[iel_][1]+EV_[2]*pbm.elpar[iel_][3]
        INSUM2 = EV_[1]*pbm.elpar[iel_][2]+EV_[2]*pbm.elpar[iel_][4]
        RATIO1 = pbm.elpar[iel_][5]/INSUM1
        RATIO2 = pbm.elpar[iel_][6]/INSUM2
        TERM1 = EV_[1]*RATIO1
        TERM2 = EV_[2]*RATIO2
        SUM = TERM1+TERM2
        SQ1 = TERM1/INSUM1
        SQ2 = TERM2/INSUM2
        SQ11 = SQ1*pbm.elpar[iel_][1]
        SQ12 = SQ2*pbm.elpar[iel_][2]
        SQ21 = SQ1*pbm.elpar[iel_][3]
        SQ22 = SQ2*pbm.elpar[iel_][4]
        TRI1 = RATIO1-SQ11-SQ12
        TRI2 = RATIO2-SQ21-SQ22
        CUB1 = SQ11/INSUM1
        CUB2 = SQ12/INSUM2
        CUB11 = CUB1*pbm.elpar[iel_][1]
        CUB12 = CUB2*pbm.elpar[iel_][2]
        CUBM21 = CUB1*pbm.elpar[iel_][3]
        CUBM22 = CUB2*pbm.elpar[iel_][4]
        CUB21 = SQ21*pbm.elpar[iel_][3]/INSUM1
        CUB22 = SQ22*pbm.elpar[iel_][4]/INSUM2
        H1 = RATIO2-SQ22-2*SQ21
        H2 = pbm.elpar[iel_][6]*pbm.elpar[iel_][2]/INSUM2^2
        H3 = SQ1*pbm.elpar[iel_][3]^2/INSUM1
        H4 = pbm.elpar[iel_][6]*pbm.elpar[iel_][4]/INSUM2^2
        H5 = SQ2*pbm.elpar[iel_][4]^2/INSUM2
        f_   = EV_[1]*SUM
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = SUM+EV_[1]*TRI1
            g_[2] = EV_[1]*TRI2
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = 2*(TRI1+EV_[1]*(-SQ11+CUB11+CUB12))
                H_[1,2] = H1+EV_[1]*(-H2+2*(CUBM21+CUBM22))
                H_[2,1] = H_[1,2]
                H_[2,2] = 2*EV_[1]*(H3-H4+H5)
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eXTAUG2"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        INSUM1 = EV_[1]*pbm.elpar[iel_][1]+EV_[2]*pbm.elpar[iel_][3]
        INSUM2 = EV_[1]*pbm.elpar[iel_][2]+EV_[2]*pbm.elpar[iel_][4]
        RATIO1 = pbm.elpar[iel_][7]/INSUM1
        RATIO2 = pbm.elpar[iel_][8]/INSUM2
        TERM1 = EV_[1]*RATIO1
        TERM2 = EV_[2]*RATIO2
        SUM = TERM1+TERM2
        SQ1 = TERM1/INSUM1
        SQ2 = TERM2/INSUM2
        SQ11 = SQ1*pbm.elpar[iel_][1]
        SQ12 = SQ2*pbm.elpar[iel_][2]
        SQ21 = SQ1*pbm.elpar[iel_][3]
        SQ22 = SQ2*pbm.elpar[iel_][4]
        TRI1 = RATIO1-SQ11-SQ12
        TRI2 = RATIO2-SQ21-SQ22
        CUB1 = SQ11/INSUM1
        CUB2 = SQ12/INSUM2
        CUB11 = CUB1*pbm.elpar[iel_][1]
        CUB12 = CUB2*pbm.elpar[iel_][2]
        CUBM21 = CUB1*pbm.elpar[iel_][3]
        CUBM22 = CUB2*pbm.elpar[iel_][4]
        CUB21 = SQ21*pbm.elpar[iel_][3]/INSUM1
        CUB22 = SQ22*pbm.elpar[iel_][4]/INSUM2
        H1 = RATIO1-SQ11-2*SQ12
        H2 = pbm.elpar[iel_][7]*pbm.elpar[iel_][3]/INSUM1^2
        H3 = SQ1*pbm.elpar[iel_][1]^2/INSUM1
        H4 = pbm.elpar[iel_][7]*pbm.elpar[iel_][1]/INSUM1^2
        H5 = SQ2*pbm.elpar[iel_][2]^2/INSUM2
        f_   = EV_[2]*SUM
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[2]*TRI1
            g_[2] = SUM+EV_[2]*TRI2
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = 2*EV_[2]*(H3-H4+H5)
                H_[1,2] = H1+EV_[2]*(-H2+2*(CUBM21+CUBM22))
                H_[2,1] = H_[1,2]
                H_[2,2] = 2*(TRI2+EV_[2]*(-SQ22+CUB21+CUB22))
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

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

    elseif action == "eXLOGXC"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,2,3)
        IV_ =  zeros(Float64,2)
        U_[2,2] = U_[2,2]+1
        U_[2,3] = U_[2,3]+1
        U_[1,1] = U_[1,1]+1
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

