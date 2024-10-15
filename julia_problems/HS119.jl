function HS119(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS119
#    *********
# 
#    Source: problem 119 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    Original Source: problem 7 in
#    A.R. Colville
#    "A comparative study on nonlinear programming"
#    IBM Scientific Center Report 320-2949, New York, 1968.
# 
#    SIF input: A.R. Conn, March 1991.
# 
#    classification = "C-OLR2-AN-16-8"
# 
#    Set useful parameters
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "HS119"

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
        v_["2"] = 2
        v_["3"] = 3
        v_["4"] = 4
        v_["5"] = 5
        v_["6"] = 6
        v_["7"] = 7
        v_["8"] = 8
        v_["9"] = 9
        v_["10"] = 10
        v_["11"] = 11
        v_["12"] = 12
        v_["13"] = 13
        v_["14"] = 14
        v_["15"] = 15
        v_["16"] = 16
        for I = Int64(v_["1"]):Int64(v_["16"])
            for J = Int64(v_["1"]):Int64(v_["16"])
                v_["A"*string(I)*","*string(J)] = 0.0
            end
        end
        for I = Int64(v_["1"]):Int64(v_["8"])
            for J = Int64(v_["1"]):Int64(v_["16"])
                v_["B"*string(I)*","*string(J)] = 0.0
            end
        end
        for I = Int64(v_["1"]):Int64(v_["16"])
            v_["A"*string(I)*","*string(I)] = 1.0
        end
        v_["A"*string(Int64(v_["1"]))*","*string(Int64(v_["4"]))] = 1.0
        v_["A"*string(Int64(v_["1"]))*","*string(Int64(v_["7"]))] = 1.0
        v_["A"*string(Int64(v_["1"]))*","*string(Int64(v_["8"]))] = 1.0
        v_["A"*string(Int64(v_["1"]))*","*string(Int64(v_["16"]))] = 1.0
        v_["A"*string(Int64(v_["2"]))*","*string(Int64(v_["3"]))] = 1.0
        v_["A"*string(Int64(v_["2"]))*","*string(Int64(v_["7"]))] = 1.0
        v_["A"*string(Int64(v_["2"]))*","*string(Int64(v_["10"]))] = 1.0
        v_["A"*string(Int64(v_["3"]))*","*string(Int64(v_["7"]))] = 1.0
        v_["A"*string(Int64(v_["3"]))*","*string(Int64(v_["9"]))] = 1.0
        v_["A"*string(Int64(v_["3"]))*","*string(Int64(v_["10"]))] = 1.0
        v_["A"*string(Int64(v_["3"]))*","*string(Int64(v_["14"]))] = 1.0
        v_["A"*string(Int64(v_["4"]))*","*string(Int64(v_["7"]))] = 1.0
        v_["A"*string(Int64(v_["4"]))*","*string(Int64(v_["11"]))] = 1.0
        v_["A"*string(Int64(v_["4"]))*","*string(Int64(v_["15"]))] = 1.0
        v_["A"*string(Int64(v_["5"]))*","*string(Int64(v_["6"]))] = 1.0
        v_["A"*string(Int64(v_["5"]))*","*string(Int64(v_["10"]))] = 1.0
        v_["A"*string(Int64(v_["5"]))*","*string(Int64(v_["12"]))] = 1.0
        v_["A"*string(Int64(v_["5"]))*","*string(Int64(v_["16"]))] = 1.0
        v_["A"*string(Int64(v_["6"]))*","*string(Int64(v_["8"]))] = 1.0
        v_["A"*string(Int64(v_["6"]))*","*string(Int64(v_["15"]))] = 1.0
        v_["A"*string(Int64(v_["7"]))*","*string(Int64(v_["11"]))] = 1.0
        v_["A"*string(Int64(v_["7"]))*","*string(Int64(v_["13"]))] = 1.0
        v_["A"*string(Int64(v_["8"]))*","*string(Int64(v_["10"]))] = 1.0
        v_["A"*string(Int64(v_["8"]))*","*string(Int64(v_["15"]))] = 1.0
        v_["A"*string(Int64(v_["9"]))*","*string(Int64(v_["12"]))] = 1.0
        v_["A"*string(Int64(v_["9"]))*","*string(Int64(v_["16"]))] = 1.0
        v_["A"*string(Int64(v_["10"]))*","*string(Int64(v_["14"]))] = 1.0
        v_["A"*string(Int64(v_["11"]))*","*string(Int64(v_["13"]))] = 1.0
        v_["A"*string(Int64(v_["12"]))*","*string(Int64(v_["14"]))] = 1.0
        v_["A"*string(Int64(v_["13"]))*","*string(Int64(v_["14"]))] = 1.0
        v_["B"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))] = 0.22
        v_["B"*string(Int64(v_["2"]))*","*string(Int64(v_["1"]))] = -1.46
        v_["B"*string(Int64(v_["3"]))*","*string(Int64(v_["1"]))] = 1.29
        v_["B"*string(Int64(v_["4"]))*","*string(Int64(v_["1"]))] = -1.10
        v_["B"*string(Int64(v_["7"]))*","*string(Int64(v_["1"]))] = 1.12
        v_["B"*string(Int64(v_["1"]))*","*string(Int64(v_["2"]))] = 0.20
        v_["B"*string(Int64(v_["3"]))*","*string(Int64(v_["2"]))] = -0.89
        v_["B"*string(Int64(v_["4"]))*","*string(Int64(v_["2"]))] = -1.06
        v_["B"*string(Int64(v_["6"]))*","*string(Int64(v_["2"]))] = -1.72
        v_["B"*string(Int64(v_["8"]))*","*string(Int64(v_["2"]))] = 0.45
        v_["B"*string(Int64(v_["1"]))*","*string(Int64(v_["3"]))] = 0.19
        v_["B"*string(Int64(v_["2"]))*","*string(Int64(v_["3"]))] = -1.30
        v_["B"*string(Int64(v_["4"]))*","*string(Int64(v_["3"]))] = 0.95
        v_["B"*string(Int64(v_["6"]))*","*string(Int64(v_["3"]))] = -0.33
        v_["B"*string(Int64(v_["8"]))*","*string(Int64(v_["3"]))] = 0.26
        v_["B"*string(Int64(v_["1"]))*","*string(Int64(v_["4"]))] = 0.25
        v_["B"*string(Int64(v_["2"]))*","*string(Int64(v_["4"]))] = 1.82
        v_["B"*string(Int64(v_["4"]))*","*string(Int64(v_["4"]))] = -0.54
        v_["B"*string(Int64(v_["5"]))*","*string(Int64(v_["4"]))] = -1.43
        v_["B"*string(Int64(v_["7"]))*","*string(Int64(v_["4"]))] = 0.31
        v_["B"*string(Int64(v_["8"]))*","*string(Int64(v_["4"]))] = -1.10
        v_["B"*string(Int64(v_["1"]))*","*string(Int64(v_["5"]))] = 0.15
        v_["B"*string(Int64(v_["2"]))*","*string(Int64(v_["5"]))] = -1.15
        v_["B"*string(Int64(v_["3"]))*","*string(Int64(v_["5"]))] = -1.16
        v_["B"*string(Int64(v_["5"]))*","*string(Int64(v_["5"]))] = 1.51
        v_["B"*string(Int64(v_["6"]))*","*string(Int64(v_["5"]))] = 1.62
        v_["B"*string(Int64(v_["8"]))*","*string(Int64(v_["5"]))] = 0.58
        v_["B"*string(Int64(v_["1"]))*","*string(Int64(v_["6"]))] = 0.11
        v_["B"*string(Int64(v_["3"]))*","*string(Int64(v_["6"]))] = -0.96
        v_["B"*string(Int64(v_["4"]))*","*string(Int64(v_["6"]))] = -1.78
        v_["B"*string(Int64(v_["5"]))*","*string(Int64(v_["6"]))] = 0.59
        v_["B"*string(Int64(v_["6"]))*","*string(Int64(v_["6"]))] = 1.24
        v_["B"*string(Int64(v_["1"]))*","*string(Int64(v_["7"]))] = 0.12
        v_["B"*string(Int64(v_["2"]))*","*string(Int64(v_["7"]))] = 0.80
        v_["B"*string(Int64(v_["4"]))*","*string(Int64(v_["7"]))] = -0.41
        v_["B"*string(Int64(v_["5"]))*","*string(Int64(v_["7"]))] = -0.33
        v_["B"*string(Int64(v_["6"]))*","*string(Int64(v_["7"]))] = 0.21
        v_["B"*string(Int64(v_["7"]))*","*string(Int64(v_["7"]))] = 1.12
        v_["B"*string(Int64(v_["8"]))*","*string(Int64(v_["7"]))] = -1.03
        v_["B"*string(Int64(v_["1"]))*","*string(Int64(v_["8"]))] = 0.13
        v_["B"*string(Int64(v_["3"]))*","*string(Int64(v_["8"]))] = -0.49
        v_["B"*string(Int64(v_["5"]))*","*string(Int64(v_["8"]))] = -0.43
        v_["B"*string(Int64(v_["6"]))*","*string(Int64(v_["8"]))] = -0.26
        v_["B"*string(Int64(v_["8"]))*","*string(Int64(v_["8"]))] = 0.10
        v_["B"*string(Int64(v_["1"]))*","*string(Int64(v_["9"]))] = 1.00
        v_["B"*string(Int64(v_["7"]))*","*string(Int64(v_["9"]))] = -0.36
        v_["B"*string(Int64(v_["2"]))*","*string(Int64(v_["10"]))] = 1.00
        v_["B"*string(Int64(v_["3"]))*","*string(Int64(v_["11"]))] = 1.00
        v_["B"*string(Int64(v_["4"]))*","*string(Int64(v_["12"]))] = 1.00
        v_["B"*string(Int64(v_["5"]))*","*string(Int64(v_["13"]))] = 1.00
        v_["B"*string(Int64(v_["6"]))*","*string(Int64(v_["14"]))] = 1.00
        v_["B"*string(Int64(v_["7"]))*","*string(Int64(v_["15"]))] = 1.00
        v_["B"*string(Int64(v_["8"]))*","*string(Int64(v_["16"]))] = 1.00
        v_["C"*string(Int64(v_["1"]))] = 2.5
        v_["C"*string(Int64(v_["2"]))] = 1.1
        v_["C"*string(Int64(v_["3"]))] = -3.1
        v_["C"*string(Int64(v_["4"]))] = -3.5
        v_["C"*string(Int64(v_["5"]))] = 1.3
        v_["C"*string(Int64(v_["6"]))] = 2.1
        v_["C"*string(Int64(v_["7"]))] = 2.3
        v_["C"*string(Int64(v_["8"]))] = -1.5
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["16"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["1"]):Int64(v_["16"])
            ig,ig_,_ = s2mpj_ii("OG"*string(I),ig_)
            arrset(gtype,ig,"<>")
        end
        for I = Int64(v_["1"]):Int64(v_["8"])
            for J = Int64(v_["1"]):Int64(v_["16"])
                ig,ig_,_ = s2mpj_ii("G"*string(I),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"G"*string(I))
                iv = ix_["X"*string(J)]
                pbm.A[ig,iv] += Float64(v_["B"*string(I)*","*string(J)])
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
        for I = Int64(v_["1"]):Int64(v_["8"])
            pbm.gconst[ig_["G"*string(I)]] = Float64(v_["C"*string(I)])
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xupper = fill(5.0,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(10.0),pb.n)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "ePROD", iet_)
        loaset(elftv,it,1,"U1")
        loaset(elftv,it,2,"U2")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"AIJ")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["16"])
            for J = Int64(v_["1"]):Int64(v_["16"])
                ename = "S"*string(I)*","*string(J)
                ie,ie_,newelt = s2mpj_ii(ename,ie_)
                if newelt > 0
                    arrset(pbm.elftype,ie,"ePROD")
                    arrset(ielftype,ie,iet_["ePROD"])
                end
                vname = "X"*string(I)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(5.0),Float64(10.0))
                posev = findfirst(x->x=="U1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(5.0),Float64(10.0))
                posev = findfirst(x->x=="U2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                posep = findfirst(x->x=="AIJ",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["A"*string(I)*","*string(J)]))
            end
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["16"])
            for J = Int64(v_["1"]):Int64(v_["16"])
                ig = ig_["OG"*string(I)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["S"*string(I)*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,1.)
            end
        end
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
        pb.pbclass = "C-OLR2-AN-16-8"
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
        TU1P1 = 2.0*EV_[1]+1
        TU2P1 = 2.0*EV_[2]+1
        FIRST = EV_[1]^2+EV_[1]+1.0
        SECOND = EV_[2]^2+EV_[2]+1.0
        f_   = pbm.elpar[iel_][1]*FIRST*SECOND
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = pbm.elpar[iel_][1]*TU1P1*SECOND
            g_[2] = pbm.elpar[iel_][1]*TU2P1*FIRST
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = pbm.elpar[iel_][1]*2.0*SECOND
                H_[1,2] = pbm.elpar[iel_][1]*TU1P1*TU2P1
                H_[2,1] = H_[1,2]
                H_[2,2] = pbm.elpar[iel_][1]*2.0*FIRST
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

