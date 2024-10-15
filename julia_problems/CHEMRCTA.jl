function CHEMRCTA(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : CHEMRCTA
#    *********
# 
#    The tubular chemical reactor model problem by Poore, using a
#    finite difference approximation to the steady state solutions.
# 
#    Source: Problem 8, eqs (8.6)--(8.9) in
#    J.J. More',
#    "A collection of nonlinear model problems"
#    Proceedings of the AMS-SIAM Summer seminar on the Computational
#    Solution of Nonlinear Systems of Equations, Colorado, 1988.
#    Argonne National Laboratory MCS-P60-0289, 1989.
# 
#    SIF input: Ph. Toint, Dec 1989.
#               minor correction by Ph. Shott, Jan 1995 and F Ruediger, Mar 1997.
# 
#    classification = "C-NOR2-MN-V-V"
# 
#    The axial coordinate interval is [0,1]
# 
#    Number of discretized point for the interval [0,1].
#    The number of variables is 2N.
# 
#       Alternative values for the SIF file parameters:
# IE N                   5              $-PARAMETER n = 10
# IE N                   25             $-PARAMETER n = 50
# IE N                   50             $-PARAMETER n = 100
# IE N                   250            $-PARAMETER n = 500    original value
# IE N                   500            $-PARAMETER n = 1000
# IE N                   2500           $-PARAMETER n = 5000
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "CHEMRCTA"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        if nargin<1
            v_["N"] = Int64(5);  #  SIF file default value
        else
            v_["N"] = Int64(args[1]);
        end
        if nargin<2
            v_["PEM"] = Float64(1.0);  #  SIF file default value
        else
            v_["PEM"] = Float64(args[2]);
        end
        if nargin<3
            v_["PEH"] = Float64(5.0);  #  SIF file default value
        else
            v_["PEH"] = Float64(args[3]);
        end
        if nargin<4
            v_["D"] = Float64(0.135);  #  SIF file default value
        else
            v_["D"] = Float64(args[4]);
        end
        if nargin<5
            v_["B"] = Float64(0.5);  #  SIF file default value
        else
            v_["B"] = Float64(args[5]);
        end
        if nargin<6
            v_["BETA"] = Float64(2.0);  #  SIF file default value
        else
            v_["BETA"] = Float64(args[6]);
        end
        if nargin<7
            v_["GAMMA"] = Float64(25.0);  #  SIF file default value
        else
            v_["GAMMA"] = Float64(args[7]);
        end
        v_["1"] = 1
        v_["2"] = 2
        v_["1.0"] = 1.0
        v_["N-1"] = -1+v_["N"]
        v_["1/H"] = Float64(v_["N-1"])
        v_["-1/H"] = -1.0*v_["1/H"]
        v_["H"] = v_["1.0"]/v_["1/H"]
        v_["1/H2"] = v_["1/H"]*v_["1/H"]
        v_["-D"] = -1.0*v_["D"]
        v_["1/PEM"] = v_["1.0"]/v_["PEM"]
        v_["1/H2PEM"] = v_["1/PEM"]*v_["1/H2"]
        v_["-1/H2PM"] = -1.0*v_["1/H2PEM"]
        v_["HPEM"] = v_["PEM"]*v_["H"]
        v_["-HPEM"] = -1.0*v_["HPEM"]
        v_["-2/H2PM"] = v_["-1/H2PM"]+v_["-1/H2PM"]
        v_["CU1"] = 1.0*v_["-HPEM"]
        v_["CUI-1"] = v_["1/H2PEM"]+v_["1/H"]
        v_["CUI"] = v_["-2/H2PM"]+v_["-1/H"]
        v_["BD"] = v_["B"]*v_["D"]
        v_["-BETA"] = -1.0*v_["BETA"]
        v_["1/PEH"] = v_["1.0"]/v_["PEH"]
        v_["1/H2PEH"] = v_["1/PEH"]*v_["1/H2"]
        v_["-1/H2PH"] = -1.0*v_["1/H2PEH"]
        v_["HPEH"] = v_["PEH"]*v_["H"]
        v_["-HPEH"] = -1.0*v_["HPEH"]
        v_["-2/H2PH"] = v_["-1/H2PH"]+v_["-1/H2PH"]
        v_["CT1"] = 1.0*v_["-HPEH"]
        v_["CTI-1"] = v_["1/H2PEH"]+v_["1/H"]
        v_["CTI"] = v_["-2/H2PH"]+v_["-1/H"]
        v_["CTI"] = v_["CTI"]+v_["-BETA"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("T"*string(I),ix_)
            arrset(pb.xnames,iv,"T"*string(I))
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("U"*string(I),ix_)
            arrset(pb.xnames,iv,"U"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("GU"*string(Int64(v_["1"])),ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"GU"*string(Int64(v_["1"])))
        iv = ix_["U"*string(Int64(v_["1"]))]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("GU"*string(Int64(v_["1"])),ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"GU"*string(Int64(v_["1"])))
        iv = ix_["U"*string(Int64(v_["2"]))]
        pbm.A[ig,iv] += Float64(v_["CU1"])
        ig,ig_,_ = s2mpj_ii("GT"*string(Int64(v_["1"])),ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"GT"*string(Int64(v_["1"])))
        iv = ix_["T"*string(Int64(v_["1"]))]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("GT"*string(Int64(v_["1"])),ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"GT"*string(Int64(v_["1"])))
        iv = ix_["T"*string(Int64(v_["2"]))]
        pbm.A[ig,iv] += Float64(v_["CT1"])
        for I = Int64(v_["2"]):Int64(v_["N-1"])
            v_["I-1"] = -1+I
            v_["I+1"] = 1+I
            ig,ig_,_ = s2mpj_ii("GU"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"GU"*string(I))
            iv = ix_["U"*string(Int64(v_["I-1"]))]
            pbm.A[ig,iv] += Float64(v_["CUI-1"])
            iv = ix_["U"*string(I)]
            pbm.A[ig,iv] += Float64(v_["CUI"])
            iv = ix_["U"*string(Int64(v_["I+1"]))]
            pbm.A[ig,iv] += Float64(v_["1/H2PEM"])
            ig,ig_,_ = s2mpj_ii("GT"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"GT"*string(I))
            iv = ix_["T"*string(I)]
            pbm.A[ig,iv] += Float64(v_["BETA"])
            iv = ix_["T"*string(Int64(v_["I-1"]))]
            pbm.A[ig,iv] += Float64(v_["CTI-1"])
            iv = ix_["T"*string(I)]
            pbm.A[ig,iv] += Float64(v_["CTI"])
            iv = ix_["T"*string(Int64(v_["I+1"]))]
            pbm.A[ig,iv] += Float64(v_["1/H2PEH"])
        end
        ig,ig_,_ = s2mpj_ii("GU"*string(Int64(v_["N"])),ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"GU"*string(Int64(v_["N"])))
        iv = ix_["U"*string(Int64(v_["N-1"]))]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("GU"*string(Int64(v_["N"])),ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"GU"*string(Int64(v_["N"])))
        iv = ix_["U"*string(Int64(v_["N"]))]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("GT"*string(Int64(v_["N"])),ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"GT"*string(Int64(v_["N"])))
        iv = ix_["T"*string(Int64(v_["N-1"]))]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("GT"*string(Int64(v_["N"])),ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"GT"*string(Int64(v_["N"])))
        iv = ix_["T"*string(Int64(v_["N"]))]
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
        pbm.gconst[ig_["GU"*string(Int64(v_["1"]))]] = Float64(v_["-HPEM"])
        pbm.gconst[ig_["GT"*string(Int64(v_["1"]))]] = Float64(v_["-HPEH"])
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        for I = Int64(v_["1"]):Int64(v_["N"])
            pb.xlower[ix_["T"*string(I)]] = 0.0000001
        end
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(1.0),pb.n)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eREAC", iet_)
        loaset(elftv,it,1,"U")
        loaset(elftv,it,2,"T")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"G")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["2"]):Int64(v_["N-1"])
            ename = "EU"*string(I)
            ie,ie_,newelt = s2mpj_ii(ename,ie_)
            if newelt > 0
                arrset(pbm.elftype,ie,"eREAC")
                arrset(ielftype,ie,iet_["eREAC"])
            end
            vname = "U"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="U",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "T"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="T",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="G",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["GAMMA"]))
            ename = "ET"*string(I)
            ie,ie_,newelt = s2mpj_ii(ename,ie_)
            if newelt > 0
                arrset(pbm.elftype,ie,"eREAC")
                arrset(ielftype,ie,iet_["eREAC"])
            end
            vname = "U"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="U",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "T"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="T",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="G",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["GAMMA"]))
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["2"]):Int64(v_["N-1"])
            ig = ig_["GU"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["EU"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["-D"]))
            ig = ig_["GT"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["ET"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["BD"]))
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               0.0
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
        pb.pbclass = "C-NOR2-MN-V-V"
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

    elseif action == "eREAC"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        DADT = pbm.elpar[iel_][1]/(EV_[2]*EV_[2])
        D2ADT2 = -2.0*DADT/EV_[2]
        EX = exp(pbm.elpar[iel_][1]-pbm.elpar[iel_][1]/EV_[2])
        UEX = EX*EV_[1]
        f_   = UEX
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EX
            g_[2] = UEX*DADT
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = EX*DADT
                H_[2,1] = H_[1,2]
                H_[2,2] = UEX*(DADT*DADT+D2ADT2)
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

