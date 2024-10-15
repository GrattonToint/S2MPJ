function HIMMELBI(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HIMMELBI
#    *********
# 
#    An unpleasant weapon assignment problem by Bracken and McCormick.
# 
#    The real problem has integer variables.
#    Also, the sign of ci have been reversed in order to have a
#    meaningful constraints on the total number of weapons (a fully
#    desirable situation).
# 
#    Source: problem 23 in
#    D.H. Himmelblau,
#    "Applied nonlinear programming",
#    McGraw-Hill, New-York, 1972.
# 
#    SIF input: Ph. Toint, March 1991.
#               minor correction by Ph. Shott, Jan 1995.
# 
#    classification = "C-OLR2-MN-100-12"
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "HIMMELBI"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["N"] = 100
        v_["NT"] = 20
        v_["A1,1"] = 1.0
        v_["A2,1"] = 0.84
        v_["A3,1"] = 0.96
        v_["A4,1"] = 1.0
        v_["A5,1"] = 0.92
        v_["A1,2"] = 0.95
        v_["A2,2"] = 0.83
        v_["A3,2"] = 0.95
        v_["A4,2"] = 1.0
        v_["A5,2"] = 0.94
        v_["A1,3"] = 1.0
        v_["A2,3"] = 0.85
        v_["A3,3"] = 0.95
        v_["A4,3"] = 1.0
        v_["A5,3"] = 0.92
        v_["A1,4"] = 1.0
        v_["A2,4"] = 0.84
        v_["A3,4"] = 0.96
        v_["A4,4"] = 1.0
        v_["A5,4"] = 0.95
        v_["A1,5"] = 1.0
        v_["A2,5"] = 0.85
        v_["A3,5"] = 0.96
        v_["A4,5"] = 1.0
        v_["A5,5"] = 0.95
        v_["A1,6"] = 0.85
        v_["A2,6"] = 0.81
        v_["A3,6"] = 0.90
        v_["A4,6"] = 1.0
        v_["A5,6"] = 0.98
        v_["A1,7"] = 0.90
        v_["A2,7"] = 0.81
        v_["A3,7"] = 0.92
        v_["A4,7"] = 1.0
        v_["A5,7"] = 0.98
        v_["A1,8"] = 0.85
        v_["A2,8"] = 0.82
        v_["A3,8"] = 0.91
        v_["A4,8"] = 1.0
        v_["A5,8"] = 1.0
        v_["A1,9"] = 0.80
        v_["A2,9"] = 0.80
        v_["A3,9"] = 0.92
        v_["A4,9"] = 1.0
        v_["A5,9"] = 1.0
        v_["A1,10"] = 1.0
        v_["A2,10"] = 0.86
        v_["A3,10"] = 0.95
        v_["A4,10"] = 0.96
        v_["A5,10"] = 0.90
        v_["A1,11"] = 1.0
        v_["A2,11"] = 1.0
        v_["A3,11"] = 0.99
        v_["A4,11"] = 0.91
        v_["A5,11"] = 0.95
        v_["A1,12"] = 1.0
        v_["A2,12"] = 0.98
        v_["A3,12"] = 0.98
        v_["A4,12"] = 0.92
        v_["A5,12"] = 0.96
        v_["A1,13"] = 1.0
        v_["A2,13"] = 1.0
        v_["A3,13"] = 0.99
        v_["A4,13"] = 0.91
        v_["A5,13"] = 0.91
        v_["A1,14"] = 1.0
        v_["A2,14"] = 0.88
        v_["A3,14"] = 0.98
        v_["A4,14"] = 0.92
        v_["A5,14"] = 0.98
        v_["A1,15"] = 1.0
        v_["A2,15"] = 0.87
        v_["A3,15"] = 0.97
        v_["A4,15"] = 0.98
        v_["A5,15"] = 0.99
        v_["A1,16"] = 1.0
        v_["A2,16"] = 0.88
        v_["A3,16"] = 0.98
        v_["A4,16"] = 0.93
        v_["A5,16"] = 0.99
        v_["A1,17"] = 1.0
        v_["A2,17"] = 0.85
        v_["A3,17"] = 0.95
        v_["A4,17"] = 1.0
        v_["A5,17"] = 1.0
        v_["A1,18"] = 0.95
        v_["A2,18"] = 0.84
        v_["A3,18"] = 0.92
        v_["A4,18"] = 1.0
        v_["A5,18"] = 1.0
        v_["A1,19"] = 1.0
        v_["A2,19"] = 0.85
        v_["A3,19"] = 0.93
        v_["A4,19"] = 1.0
        v_["A5,19"] = 1.0
        v_["A1,20"] = 1.0
        v_["A2,20"] = 0.85
        v_["A3,20"] = 0.92
        v_["A4,20"] = 1.0
        v_["A5,20"] = 1.0
        v_["B1"] = 30.0
        v_["B6"] = 100.0
        v_["B10"] = 40.0
        v_["B14"] = 50.0
        v_["B15"] = 70.0
        v_["B16"] = 35.0
        v_["B20"] = 10.0
        v_["U1"] = 60.0
        v_["U2"] = 50.0
        v_["U3"] = 50.0
        v_["U4"] = 75.0
        v_["U5"] = 40.0
        v_["U6"] = 60.0
        v_["U7"] = 35.0
        v_["U8"] = 30.0
        v_["U9"] = 25.0
        v_["U10"] = 150.0
        v_["U11"] = 30.0
        v_["U12"] = 45.0
        v_["U13"] = 125.0
        v_["U14"] = 200.0
        v_["U15"] = 200.0
        v_["U16"] = 130.0
        v_["U17"] = 100.0
        v_["U18"] = 100.0
        v_["U19"] = 100.0
        v_["U20"] = 150.0
        v_["C1"] = 200.0
        v_["C2"] = 100.0
        v_["C3"] = 300.0
        v_["C4"] = 150.0
        v_["C5"] = 250.0
        v_["1"] = 1
        v_["2"] = 2
        v_["3"] = 3
        v_["4"] = 4
        v_["5"] = 5
        v_["NW"] = 0.0
        for I = Int64(v_["1"]):Int64(v_["5"])
            v_["NW"] = v_["NW"]+v_["C"*string(I)]
        end
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for J = Int64(v_["1"]):Int64(v_["NT"])
            for I = Int64(v_["1"]):Int64(v_["5"])
                iv,ix_,_ = s2mpj_ii("X"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"X"*string(I)*","*string(J))
            end
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for J = Int64(v_["1"]):Int64(v_["NT"])
            ig,ig_,_ = s2mpj_ii("P"*string(J),ig_)
            arrset(gtype,ig,"<>")
            v_["1/UJ"] = 1.0/v_["U"*string(J)]
            arrset(pbm.gscale,ig,Float64(v_["1/UJ"]))
        end
        v_["J"] = 1
        for I = Int64(v_["1"]):Int64(v_["5"])
            ig,ig_,_ = s2mpj_ii("CB"*string(Int64(v_["J"])),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"CB"*string(Int64(v_["J"])))
            iv = ix_["X"*string(I)*","*string(Int64(v_["J"]))]
            pbm.A[ig,iv] += Float64(1.0)
        end
        v_["J"] = 6
        for I = Int64(v_["1"]):Int64(v_["5"])
            ig,ig_,_ = s2mpj_ii("CB"*string(Int64(v_["J"])),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"CB"*string(Int64(v_["J"])))
            iv = ix_["X"*string(I)*","*string(Int64(v_["J"]))]
            pbm.A[ig,iv] += Float64(1.0)
        end
        v_["J"] = 10
        for I = Int64(v_["1"]):Int64(v_["5"])
            ig,ig_,_ = s2mpj_ii("CB"*string(Int64(v_["J"])),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"CB"*string(Int64(v_["J"])))
            iv = ix_["X"*string(I)*","*string(Int64(v_["J"]))]
            pbm.A[ig,iv] += Float64(1.0)
        end
        v_["J"] = 14
        for I = Int64(v_["1"]):Int64(v_["5"])
            ig,ig_,_ = s2mpj_ii("CB"*string(Int64(v_["J"])),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"CB"*string(Int64(v_["J"])))
            iv = ix_["X"*string(I)*","*string(Int64(v_["J"]))]
            pbm.A[ig,iv] += Float64(1.0)
        end
        v_["J"] = 15
        for I = Int64(v_["1"]):Int64(v_["5"])
            ig,ig_,_ = s2mpj_ii("CB"*string(Int64(v_["J"])),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"CB"*string(Int64(v_["J"])))
            iv = ix_["X"*string(I)*","*string(Int64(v_["J"]))]
            pbm.A[ig,iv] += Float64(1.0)
        end
        v_["J"] = 16
        for I = Int64(v_["1"]):Int64(v_["5"])
            ig,ig_,_ = s2mpj_ii("CB"*string(Int64(v_["J"])),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"CB"*string(Int64(v_["J"])))
            iv = ix_["X"*string(I)*","*string(Int64(v_["J"]))]
            pbm.A[ig,iv] += Float64(1.0)
        end
        v_["J"] = 20
        for I = Int64(v_["1"]):Int64(v_["5"])
            ig,ig_,_ = s2mpj_ii("CB"*string(Int64(v_["J"])),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"CB"*string(Int64(v_["J"])))
            iv = ix_["X"*string(I)*","*string(Int64(v_["J"]))]
            pbm.A[ig,iv] += Float64(1.0)
        end
        for I = Int64(v_["1"]):Int64(v_["5"])
            for J = Int64(v_["1"]):Int64(v_["NT"])
                ig,ig_,_ = s2mpj_ii("CC"*string(I),ig_)
                arrset(gtype,ig,"<=")
                arrset(pb.cnames,ig,"CC"*string(I))
                iv = ix_["X"*string(I)*","*string(J)]
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
        for J = Int64(v_["1"]):Int64(v_["NT"])
            pbm.gconst[ig_["P"*string(J)]] = Float64(1.0)
        end
        v_["J"] = 1
        pbm.gconst[ig_["CB"*string(Int64(v_["J"]))]]  = (
              Float64(v_["B"*string(Int64(v_["J"]))]))
        v_["J"] = 6
        pbm.gconst[ig_["CB"*string(Int64(v_["J"]))]]  = (
              Float64(v_["B"*string(Int64(v_["J"]))]))
        v_["J"] = 10
        pbm.gconst[ig_["CB"*string(Int64(v_["J"]))]]  = (
              Float64(v_["B"*string(Int64(v_["J"]))]))
        v_["J"] = 14
        pbm.gconst[ig_["CB"*string(Int64(v_["J"]))]]  = (
              Float64(v_["B"*string(Int64(v_["J"]))]))
        v_["J"] = 15
        pbm.gconst[ig_["CB"*string(Int64(v_["J"]))]]  = (
              Float64(v_["B"*string(Int64(v_["J"]))]))
        v_["J"] = 16
        pbm.gconst[ig_["CB"*string(Int64(v_["J"]))]]  = (
              Float64(v_["B"*string(Int64(v_["J"]))]))
        v_["J"] = 20
        pbm.gconst[ig_["CB"*string(Int64(v_["J"]))]]  = (
              Float64(v_["B"*string(Int64(v_["J"]))]))
        for I = Int64(v_["1"]):Int64(v_["5"])
            pbm.gconst[ig_["CC"*string(I)]] = Float64(v_["C"*string(I)])
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xupper = fill(v_["NW"],pb.n)
        pb.xlower = zeros(Float64,pb.n)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "en5PEXP", iet_)
        loaset(elftv,it,1,"Y1")
        loaset(elftv,it,2,"Y2")
        loaset(elftv,it,3,"Y3")
        loaset(elftv,it,4,"Y4")
        loaset(elftv,it,5,"Y5")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"A1")
        loaset(elftp,it,2,"A2")
        loaset(elftp,it,3,"A3")
        loaset(elftp,it,4,"A4")
        loaset(elftp,it,5,"A5")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for J = Int64(v_["1"]):Int64(v_["NT"])
            ename = "PP"*string(J)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"en5PEXP")
            arrset(ielftype,ie,iet_["en5PEXP"])
            vname = "X"*string(Int64(v_["1"]))*","*string(J)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(v_["NW"]),nothing)
            posev = findfirst(x->x=="Y1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["2"]))*","*string(J)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(v_["NW"]),nothing)
            posev = findfirst(x->x=="Y2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["3"]))*","*string(J)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(v_["NW"]),nothing)
            posev = findfirst(x->x=="Y3",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["4"]))*","*string(J)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(v_["NW"]),nothing)
            posev = findfirst(x->x=="Y4",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["5"]))*","*string(J)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(v_["NW"]),nothing)
            posev = findfirst(x->x=="Y5",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="A1",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["A"*string(Int64(v_["1"]))*","*string(J)]))
            posep = findfirst(x->x=="A2",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["A"*string(Int64(v_["2"]))*","*string(J)]))
            posep = findfirst(x->x=="A3",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["A"*string(Int64(v_["3"]))*","*string(J)]))
            posep = findfirst(x->x=="A4",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["A"*string(Int64(v_["4"]))*","*string(J)]))
            posep = findfirst(x->x=="A5",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["A"*string(Int64(v_["5"]))*","*string(J)]))
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for J = Int64(v_["1"]):Int64(v_["NT"])
            ig = ig_["P"*string(J)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["PP"*string(J)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               -1735.56958
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.cupper[1:pb.nle] = zeros(Float64,pb.nle)
        pb.clower[pb.nle+pb.neq+1:pb.m] = zeros(Float64,pb.nge)
        pb.cupper[1:pb.nge] = fill(Inf,pb.nge)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-OLR2-MN-100-12"
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

    elseif action == "en5PEXP"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        LA1 = log(pbm.elpar[iel_][1])
        LA2 = log(pbm.elpar[iel_][2])
        LA3 = log(pbm.elpar[iel_][3])
        LA4 = log(pbm.elpar[iel_][4])
        LA5 = log(pbm.elpar[iel_][5])
        FF  = (
              pbm.elpar[iel_][1]^EV_[1]*pbm.elpar[iel_][2]^EV_[2]*pbm.elpar[iel_][3]^EV_[3]*pbm.elpar[iel_][4]^EV_[4]*pbm.elpar[iel_][5]^EV_[5])
        f_   = FF
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = LA1*FF
            g_[2] = LA2*FF
            g_[3] = LA3*FF
            g_[4] = LA4*FF
            g_[5] = LA5*FF
            if nargout>2
                H_ = zeros(Float64,5,5)
                H_[1,1] = LA1*LA1*FF
                H_[1,2] = LA1*LA2*FF
                H_[2,1] = H_[1,2]
                H_[1,3] = LA1*LA3*FF
                H_[3,1] = H_[1,3]
                H_[1,4] = LA1*LA4*FF
                H_[4,1] = H_[1,4]
                H_[1,5] = LA1*LA5*FF
                H_[5,1] = H_[1,5]
                H_[2,2] = LA2*LA2*FF
                H_[2,3] = LA2*LA3*FF
                H_[3,2] = H_[2,3]
                H_[2,4] = LA2*LA4*FF
                H_[4,2] = H_[2,4]
                H_[2,5] = LA2*LA5*FF
                H_[5,2] = H_[2,5]
                H_[3,3] = LA3*LA3*FF
                H_[3,4] = LA3*LA4*FF
                H_[4,3] = H_[3,4]
                H_[3,5] = LA3*LA5*FF
                H_[5,3] = H_[3,5]
                H_[4,4] = LA4*LA4*FF
                H_[4,5] = LA4*LA5*FF
                H_[5,4] = H_[4,5]
                H_[5,5] = LA5*LA5*FF
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

