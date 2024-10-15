function EXPFITA(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : EXPFITA
#    *********
# 
#    One sided rational approximation to the exponential function, as
#    described by Powell.
# 
#    Source:
#    M.J.D. Powell,
#    "A tolerant algorithm for linearly constrained optimization
#    calculations"'
#    Mathematical Programming 45(3), pp.561--562, 1989.
# 
#    SDIF input: Ph. Toint and N. Gould, May 1990.
# 
#    classification = "C-OLR2-AN-5-22"
# 
#    Number of fitting points
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "EXPFITA"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["R"] = 11
        v_["1"] = 1
        v_["5.0"] = 5.0
        v_["R-1"] = -1+v_["R"]
        v_["RR-1"] = Float64(v_["R-1"])
        v_["5/R-1"] = v_["5.0"]/v_["RR-1"]
        for I = Int64(v_["1"]):Int64(v_["R"])
            v_["I-1"] = -1+I
            v_["RI-1"] = Float64(v_["I-1"])
            v_["T"*string(I)] = v_["RI-1"]*v_["5/R-1"]
            v_["ET"*string(I)] = exp(v_["T"*string(I)])
        end
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        iv,ix_,_ = s2mpj_ii("P0",ix_)
        arrset(pb.xnames,iv,"P0")
        iv,ix_,_ = s2mpj_ii("P1",ix_)
        arrset(pb.xnames,iv,"P1")
        iv,ix_,_ = s2mpj_ii("P2",ix_)
        arrset(pb.xnames,iv,"P2")
        iv,ix_,_ = s2mpj_ii("Q1",ix_)
        arrset(pb.xnames,iv,"Q1")
        iv,ix_,_ = s2mpj_ii("Q2",ix_)
        arrset(pb.xnames,iv,"Q2")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        for I = Int64(v_["1"]):Int64(v_["R"])
            v_["TM5"] = -5.0+v_["T"*string(I)]
            v_["TM5SQ"] = v_["TM5"]*v_["TM5"]
            v_["QC1"] = v_["TM5"]*v_["ET"*string(I)]
            v_["QC2"] = v_["TM5SQ"]*v_["ET"*string(I)]
            v_["-QC1"] = -1.0*v_["QC1"]
            v_["-QC2"] = -1.0*v_["QC2"]
            v_["2T"] = v_["T"*string(I)]*v_["T"*string(I)]
            ig,ig_,_ = s2mpj_ii("C"*string(I),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"C"*string(I))
            iv = ix_["P0"]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["P1"]
            pbm.A[ig,iv] += Float64(v_["T"*string(I)])
            iv = ix_["P2"]
            pbm.A[ig,iv] += Float64(v_["2T"])
            iv = ix_["Q1"]
            pbm.A[ig,iv] += Float64(v_["-QC1"])
            iv = ix_["Q2"]
            pbm.A[ig,iv] += Float64(v_["-QC2"])
            ig,ig_,_ = s2mpj_ii("B"*string(I),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"B"*string(I))
            iv = ix_["Q1"]
            pbm.A[ig,iv] += Float64(v_["TM5"])
            iv = ix_["Q2"]
            pbm.A[ig,iv] += Float64(v_["TM5SQ"])
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
        for I = Int64(v_["1"]):Int64(v_["R"])
            pbm.gconst[ig_["C"*string(I)]] = Float64(v_["ET"*string(I)])
            pbm.gconst[ig_["B"*string(I)]] = Float64(-0.99999)
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        if haskey(ix_,"P0")
            pb.x0[ix_["P0"]] = Float64(1.0)
        else
            pb.y0[findfirst(x->x==ig_["P0"],pbm.congrps)] = Float64(1.0)
        end
        if haskey(ix_,"P1")
            pb.x0[ix_["P1"]] = Float64(1.0)
        else
            pb.y0[findfirst(x->x==ig_["P1"],pbm.congrps)] = Float64(1.0)
        end
        if haskey(ix_,"P2")
            pb.x0[ix_["P2"]] = Float64(6.0)
        else
            pb.y0[findfirst(x->x==ig_["P2"],pbm.congrps)] = Float64(6.0)
        end
        if haskey(ix_,"Q1")
            pb.x0[ix_["Q1"]] = Float64(0.0)
        else
            pb.y0[findfirst(x->x==ig_["Q1"],pbm.congrps)] = Float64(0.0)
        end
        if haskey(ix_,"Q2")
            pb.x0[ix_["Q2"]] = Float64(0.0)
        else
            pb.y0[findfirst(x->x==ig_["Q2"],pbm.congrps)] = Float64(0.0)
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eFIT", iet_)
        loaset(elftv,it,1,"P0")
        loaset(elftv,it,2,"P1")
        loaset(elftv,it,3,"P2")
        loaset(elftv,it,4,"Q1")
        loaset(elftv,it,5,"Q2")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"T")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["R"])
            ename = "F"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eFIT")
            arrset(ielftype,ie,iet_["eFIT"])
            vname = "P0"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="P0",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "P1"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="P1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "P2"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="P2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "Q1"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="Q1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "Q2"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="Q2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="T",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["T"*string(I)]))
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["R"])
            ig = ig_["OBJ"]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["F"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
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
        pb.pbclass = "C-OLR2-AN-5-22"
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

    elseif action == "eFIT"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        TM5 = pbm.elpar[iel_][1]-5.0
        TM5SQ = TM5*TM5
        T2 = pbm.elpar[iel_][1]*pbm.elpar[iel_][1]
        ET = exp(pbm.elpar[iel_][1])
        QT = 1.0+EV_[4]*TM5+EV_[5]*TM5SQ
        ETQT = ET*QT
        ETQT2 = ETQT*QT
        ETQT3 = ETQT2*QT
        PT = EV_[1]+EV_[2]*pbm.elpar[iel_][1]+EV_[3]*T2
        F = PT/ETQT-1.0
        TWOF = F+F
        DFDP0 = 1.0/ETQT
        DFDP1 = pbm.elpar[iel_][1]/ETQT
        DFDP2 = T2/ETQT
        DFDQ1 = -PT*TM5/ETQT2
        DFDQ2 = -PT*TM5SQ/ETQT2
        D2P0Q1 = -TM5/ETQT2
        D2P0Q2 = -TM5SQ/ETQT2
        D2P1Q1 = -pbm.elpar[iel_][1]*TM5/ETQT2
        D2P1Q2 = -pbm.elpar[iel_][1]*TM5SQ/ETQT2
        D2P2Q1 = -T2*TM5/ETQT2
        D2P2Q2 = -T2*TM5SQ/ETQT2
        D2Q1Q1 = 2.0*PT*TM5SQ/ETQT3
        D2Q1Q2 = 2.0*PT*TM5SQ*TM5/ETQT3
        D2Q2Q2 = 2.0*PT*TM5SQ*TM5SQ/ETQT3
        f_   = F*F
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = TWOF*DFDP0
            g_[2] = TWOF*DFDP1
            g_[3] = TWOF*DFDP2
            g_[4] = TWOF*DFDQ1
            g_[5] = TWOF*DFDQ2
            if nargout>2
                H_ = zeros(Float64,5,5)
                H_[1,1] = 2.0*DFDP0*DFDP0
                H_[1,2] = 2.0*DFDP0*DFDP1
                H_[2,1] = H_[1,2]
                H_[1,3] = 2.0*DFDP0*DFDP2
                H_[3,1] = H_[1,3]
                H_[2,2] = 2.0*DFDP1*DFDP1
                H_[2,3] = 2.0*DFDP1*DFDP2
                H_[3,2] = H_[2,3]
                H_[3,3] = 2.0*DFDP2*DFDP2
                H_[1,4] = TWOF*D2P0Q1+2.0*DFDP0*DFDQ1
                H_[4,1] = H_[1,4]
                H_[1,5] = TWOF*D2P0Q2+2.0*DFDP0*DFDQ2
                H_[5,1] = H_[1,5]
                H_[2,4] = TWOF*D2P1Q1+2.0*DFDP1*DFDQ1
                H_[4,2] = H_[2,4]
                H_[2,5] = TWOF*D2P1Q2+2.0*DFDP1*DFDQ2
                H_[5,2] = H_[2,5]
                H_[3,4] = TWOF*D2P2Q1+2.0*DFDP2*DFDQ1
                H_[4,3] = H_[3,4]
                H_[3,5] = TWOF*D2P2Q2+2.0*DFDP2*DFDQ2
                H_[5,3] = H_[3,5]
                H_[4,4] = TWOF*D2Q1Q1+2.0*DFDQ1*DFDQ1
                H_[4,5] = TWOF*D2Q1Q2+2.0*DFDQ1*DFDQ2
                H_[5,4] = H_[4,5]
                H_[5,5] = TWOF*D2Q2Q2+2.0*DFDQ2*DFDQ2
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

