function BATCH(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : BATCH
#    *********
# 
#    Source: Optimal Design of Multiproduct Batch Plant
#    G.R. Kocis & I.E. Grossmann,
#    "Global OPtimization of Nonconvex Mixed Integer Nonlinear Programmming
#     (MINLP) problems in Process Synthesis", Indust. Engng. Chem. Res.,
#    No. 27, pp 1407--1421, 1988.
# 
#    SIF input: S. Leyffer, October 1997
# 
#    classification = "C-OOR2-AN-46-73"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 6 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "BATCH"

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
        v_["M"] = 6
        v_["N"] = 5
        v_["NU"] = 4
        v_["LOGNU"] = log(4.0)
        v_["VL"] = log(300.0)
        v_["VU"] = log(3000.0)
        v_["H"] = 6000.0
        v_["TLO1"] = 0.729961
        v_["TLO2"] = 0.530628
        v_["TLO3"] = 1.09024
        v_["TLO4"] = -0.133531
        v_["TLO5"] = 0.0487901
        v_["TUP1"] = 2.11626
        v_["TUP2"] = 1.91626
        v_["TUP3"] = 2.47654
        v_["TUP4"] = 1.25276
        v_["TUP5"] = 1.43508
        v_["BLO1"] = 4.45966
        v_["BLO2"] = 3.74950
        v_["BLO3"] = 4.49144
        v_["BLO4"] = 3.14988
        v_["BLO5"] = 3.04452
        v_["BUP1"] = 397.747
        v_["BUP2"] = 882.353
        v_["BUP3"] = 833.333
        v_["BUP4"] = 638.298
        v_["BUP5"] = 666.667
        v_["Q1"] = 250000.0
        v_["Q2"] = 150000.0
        v_["Q3"] = 180000.0
        v_["Q4"] = 160000.0
        v_["Q5"] = 120000.0
        v_["LOGI1"] = log(1.0)
        v_["LOGI2"] = log(2.0)
        v_["LOGI3"] = log(3.0)
        v_["LOGI4"] = log(4.0)
        v_["S1,1"] = log(7.9)
        v_["S2,1"] = log(0.7)
        v_["S3,1"] = log(0.7)
        v_["S4,1"] = log(4.7)
        v_["S5,1"] = log(1.2)
        v_["S1,2"] = log(2.0)
        v_["S2,2"] = log(0.8)
        v_["S3,2"] = log(2.6)
        v_["S4,2"] = log(2.3)
        v_["S5,2"] = log(3.6)
        v_["S1,3"] = log(5.2)
        v_["S2,3"] = log(0.9)
        v_["S3,3"] = log(1.6)
        v_["S4,3"] = log(1.6)
        v_["S5,3"] = log(2.4)
        v_["S1,4"] = log(4.9)
        v_["S2,4"] = log(3.4)
        v_["S3,4"] = log(3.6)
        v_["S4,4"] = log(2.7)
        v_["S5,4"] = log(4.5)
        v_["S1,5"] = log(6.1)
        v_["S2,5"] = log(2.1)
        v_["S3,5"] = log(3.2)
        v_["S4,5"] = log(1.2)
        v_["S5,5"] = log(1.6)
        v_["S1,6"] = log(4.2)
        v_["S2,6"] = log(2.5)
        v_["S3,6"] = log(2.9)
        v_["S4,6"] = log(2.5)
        v_["S5,6"] = log(2.1)
        v_["T1,1"] = log(6.4)
        v_["T2,1"] = log(6.8)
        v_["T3,1"] = log(1.0)
        v_["T4,1"] = log(3.2)
        v_["T5,1"] = log(2.1)
        v_["T1,2"] = log(4.7)
        v_["T2,2"] = log(6.4)
        v_["T3,2"] = log(6.3)
        v_["T4,2"] = log(3.0)
        v_["T5,2"] = log(2.5)
        v_["T1,3"] = log(8.3)
        v_["T2,3"] = log(6.5)
        v_["T3,3"] = log(5.4)
        v_["T4,3"] = log(3.5)
        v_["T5,3"] = log(4.2)
        v_["T1,4"] = log(3.9)
        v_["T2,4"] = log(4.4)
        v_["T3,4"] = log(11.9)
        v_["T4,4"] = log(3.3)
        v_["T5,4"] = log(3.6)
        v_["T1,5"] = log(2.1)
        v_["T2,5"] = log(2.3)
        v_["T3,5"] = log(5.7)
        v_["T4,5"] = log(2.8)
        v_["T5,5"] = log(3.7)
        v_["T1,6"] = log(1.2)
        v_["T2,6"] = log(3.2)
        v_["T3,6"] = log(6.2)
        v_["T4,6"] = log(3.4)
        v_["T5,6"] = log(2.2)
        for J = Int64(v_["1"]):Int64(v_["M"])
            v_["ALPHA"*string(J)] = 250.0
            v_["BETA"*string(J)] = 0.6
        end
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for J = Int64(v_["1"]):Int64(v_["M"])
            iv,ix_,_ = s2mpj_ii("N"*string(J),ix_)
            arrset(pb.xnames,iv,"N"*string(J))
        end
        for J = Int64(v_["1"]):Int64(v_["M"])
            iv,ix_,_ = s2mpj_ii("V"*string(J),ix_)
            arrset(pb.xnames,iv,"V"*string(J))
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("B"*string(I),ix_)
            arrset(pb.xnames,iv,"B"*string(I))
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("TL"*string(I),ix_)
            arrset(pb.xnames,iv,"TL"*string(I))
        end
        for J = Int64(v_["1"]):Int64(v_["M"])
            for K = Int64(v_["1"]):Int64(v_["NU"])
                iv,ix_,_ = s2mpj_ii("Y"*string(K)*","*string(J),ix_)
                arrset(pb.xnames,iv,"Y"*string(K)*","*string(J))
            end
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("COST",ig_)
        arrset(gtype,ig,"<>")
        for I = Int64(v_["1"]):Int64(v_["N"])
            for J = Int64(v_["1"]):Int64(v_["M"])
                ig,ig_,_ = s2mpj_ii("VOL"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,">=")
                arrset(pb.cnames,ig,"VOL"*string(I)*","*string(J))
                iv = ix_["V"*string(J)]
                pbm.A[ig,iv] += Float64(1.0)
                iv = ix_["B"*string(I)]
                pbm.A[ig,iv] += Float64(-1.0)
            end
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            for J = Int64(v_["1"]):Int64(v_["M"])
                ig,ig_,_ = s2mpj_ii("CYCL"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,">=")
                arrset(pb.cnames,ig,"CYCL"*string(I)*","*string(J))
                iv = ix_["N"*string(J)]
                pbm.A[ig,iv] += Float64(1.0)
                iv = ix_["TL"*string(I)]
                pbm.A[ig,iv] += Float64(1.0)
            end
        end
        ig,ig_,_ = s2mpj_ii("HORIZON",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"HORIZON")
        for J = Int64(v_["1"]):Int64(v_["M"])
            for K = Int64(v_["1"]):Int64(v_["NU"])
                ig,ig_,_ = s2mpj_ii("NPAR"*string(J),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"NPAR"*string(J))
                iv = ix_["Y"*string(K)*","*string(J)]
                pbm.A[ig,iv] += Float64(v_["LOGI"*string(K)])
            end
            ig,ig_,_ = s2mpj_ii("NPAR"*string(J),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"NPAR"*string(J))
            iv = ix_["N"*string(J)]
            pbm.A[ig,iv] += Float64(-1.0)
        end
        for J = Int64(v_["1"]):Int64(v_["M"])
            for K = Int64(v_["1"]):Int64(v_["NU"])
                ig,ig_,_ = s2mpj_ii("SOS1"*string(J),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"SOS1"*string(J))
                iv = ix_["Y"*string(K)*","*string(J)]
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
        for I = Int64(v_["1"]):Int64(v_["N"])
            for J = Int64(v_["1"]):Int64(v_["M"])
                pbm.gconst[ig_["VOL"*string(I)*","*string(J)]]  = (
                      Float64(v_["S"*string(I)*","*string(J)]))
            end
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            for J = Int64(v_["1"]):Int64(v_["M"])
                pbm.gconst[ig_["CYCL"*string(I)*","*string(J)]]  = (
                      Float64(v_["T"*string(I)*","*string(J)]))
            end
        end
        pbm.gconst[ig_["HORIZON"]] = Float64(v_["H"])
        for J = Int64(v_["1"]):Int64(v_["M"])
            pbm.gconst[ig_["SOS1"*string(J)]] = Float64(1.0)
        end
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        for J = Int64(v_["1"]):Int64(v_["M"])
            pb.xupper[ix_["N"*string(J)]] = v_["LOGNU"]
            pb.xlower[ix_["V"*string(J)]] = v_["VL"]
            pb.xupper[ix_["V"*string(J)]] = v_["VU"]
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            pb.xlower[ix_["B"*string(I)]] = v_["BLO"*string(I)]
            pb.xupper[ix_["B"*string(I)]] = v_["BUP"*string(I)]
            pb.xlower[ix_["TL"*string(I)]] = v_["TLO"*string(I)]
            pb.xupper[ix_["TL"*string(I)]] = v_["TUP"*string(I)]
        end
        for J = Int64(v_["1"]):Int64(v_["M"])
            for K = Int64(v_["1"]):Int64(v_["NU"])
                pb.xupper[ix_["Y"*string(K)*","*string(J)]] = 1.0
            end
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eEXPXAY", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"A")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for J = Int64(v_["1"]):Int64(v_["M"])
            ename = "EXPO"*string(J)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eEXPXAY")
            arrset(ielftype,ie,iet_["eEXPXAY"])
            vname = "N"*string(J)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "V"*string(J)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="A",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["BETA"*string(J)]))
        end
        for I = Int64(v_["1"]):Int64(v_["M"])
            ename = "EXPC"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eEXPXAY")
            arrset(ielftype,ie,iet_["eEXPXAY"])
            vname = "TL"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "B"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="A",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(-1.0))
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for J = Int64(v_["1"]):Int64(v_["M"])
            ig = ig_["COST"]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["EXPO"*string(J)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["ALPHA"*string(J)]))
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            ig = ig_["HORIZON"]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["EXPC"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["Q"*string(I)]))
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.cupper[1:pb.nle] = zeros(Float64,pb.nle)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.clower[pb.nle+pb.neq+1:pb.m] = zeros(Float64,pb.nge)
        pb.cupper[1:pb.nge] = fill(Inf,pb.nge)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-OOR2-AN-46-73"
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

    elseif action == "eEXPXAY"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        FVALUE = exp(EV_[1]+pbm.elpar[iel_][1]*EV_[2])
        GYVALU = pbm.elpar[iel_][1]*FVALUE
        f_   = FVALUE
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = FVALUE
            g_[2] = GYVALU
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = FVALUE
                H_[1,2] = GYVALU
                H_[2,1] = H_[1,2]
                H_[2,2] = pbm.elpar[iel_][1]*GYVALU
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

