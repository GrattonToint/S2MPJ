function HS117(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS117
#    *********
# 
#    Source: problem 117 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: Nick Gould, August 1991.
# 
#    classification = "C-OQR2-AN-15-5"
# 
#    Number of constraints
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "HS117"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["N"] = 5
        v_["M"] = 10
        v_["M+N"] = v_["M"]+v_["N"]
        v_["1"] = 1
        v_["E1"] = -15.0
        v_["E2"] = -27.0
        v_["E3"] = -36.0
        v_["E4"] = -18.0
        v_["E5"] = -12.0
        v_["C1,1"] = 30.0
        v_["C2,1"] = -20.0
        v_["C3,1"] = -10.0
        v_["C4,1"] = 32.0
        v_["C5,1"] = -10.0
        v_["C1,2"] = -20.0
        v_["C2,2"] = 39.0
        v_["C3,2"] = -6.0
        v_["C4,2"] = -31.0
        v_["C5,2"] = 32.0
        v_["C1,3"] = -10.0
        v_["C2,3"] = -6.0
        v_["C3,3"] = 10.0
        v_["C4,3"] = -6.0
        v_["C5,3"] = -10.0
        v_["C1,4"] = 32.0
        v_["C2,4"] = -31.0
        v_["C3,4"] = -6.0
        v_["C4,4"] = 39.0
        v_["C5,4"] = -20.0
        v_["C1,5"] = -10.0
        v_["C2,5"] = 32.0
        v_["C3,5"] = -10.0
        v_["C4,5"] = -20.0
        v_["C5,5"] = 30.0
        v_["D1"] = 4.0
        v_["D2"] = 8.0
        v_["D3"] = 10.0
        v_["D4"] = 6.0
        v_["D5"] = 2.0
        v_["A1,1"] = -16.0
        v_["A2,1"] = 0.0
        v_["A3,1"] = -3.5
        v_["A4,1"] = 0.0
        v_["A5,1"] = 0.0
        v_["A6,1"] = 2.0
        v_["A7,1"] = -1.0
        v_["A8,1"] = -1.0
        v_["A9,1"] = 1.0
        v_["A10,1"] = 1.0
        v_["A1,2"] = 2.0
        v_["A2,2"] = -2.0
        v_["A3,2"] = 0.0
        v_["A4,2"] = -2.0
        v_["A5,2"] = -9.0
        v_["A6,2"] = 0.0
        v_["A7,2"] = -1.0
        v_["A8,2"] = -2.0
        v_["A9,2"] = 2.0
        v_["A10,2"] = 1.0
        v_["A1,3"] = 0.0
        v_["A2,3"] = 0.0
        v_["A3,3"] = 2.0
        v_["A4,3"] = 0.0
        v_["A5,3"] = -2.0
        v_["A6,3"] = -4.0
        v_["A7,3"] = -1.0
        v_["A8,3"] = -3.0
        v_["A9,3"] = 3.0
        v_["A10,3"] = 1.0
        v_["A1,4"] = 1.0
        v_["A2,4"] = 4.0
        v_["A3,4"] = 0.0
        v_["A4,4"] = -4.0
        v_["A5,4"] = 1.0
        v_["A6,4"] = 0.0
        v_["A7,4"] = -1.0
        v_["A8,4"] = -2.0
        v_["A9,4"] = 4.0
        v_["A10,4"] = 1.0
        v_["A1,5"] = 0.0
        v_["A2,5"] = 2.0
        v_["A3,5"] = 0.0
        v_["A4,5"] = -1.0
        v_["A5,5"] = -2.8
        v_["A6,5"] = 0.0
        v_["A7,5"] = -1.0
        v_["A8,5"] = -1.0
        v_["A9,5"] = 5.0
        v_["A10,5"] = 1.0
        v_["B1"] = -40.0
        v_["B2"] = -2.0
        v_["B3"] = -0.25
        v_["B4"] = -4.0
        v_["B5"] = -4.0
        v_["B6"] = -1.0
        v_["B7"] = -40.0
        v_["B8"] = -60.0
        v_["B9"] = 5.0
        v_["B10"] = 1.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["M+N"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for J = Int64(v_["1"]):Int64(v_["M"])
            v_["-BJ"] = -1.0*v_["B"*string(J)]
            ig,ig_,_ = s2mpj_ii("OBJ",ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X"*string(J)]
            pbm.A[ig,iv] += Float64(v_["-BJ"])
        end
        for J = Int64(v_["1"]):Int64(v_["N"])
            for K = Int64(v_["1"]):Int64(v_["N"])
                v_["M+K"] = v_["M"]+K
                v_["2CKJ"] = 2.0*v_["C"*string(K)*","*string(J)]
                ig,ig_,_ = s2mpj_ii("C"*string(J),ig_)
                arrset(gtype,ig,">=")
                arrset(pb.cnames,ig,"C"*string(J))
                iv = ix_["X"*string(Int64(v_["M+K"]))]
                pbm.A[ig,iv] += Float64(v_["2CKJ"])
            end
            for K = Int64(v_["1"]):Int64(v_["M"])
                v_["-AKJ"] = -1.0*v_["A"*string(K)*","*string(J)]
                ig,ig_,_ = s2mpj_ii("C"*string(J),ig_)
                arrset(gtype,ig,">=")
                arrset(pb.cnames,ig,"C"*string(J))
                iv = ix_["X"*string(K)]
                pbm.A[ig,iv] += Float64(v_["-AKJ"])
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
        for J = Int64(v_["1"]):Int64(v_["N"])
            v_["-EJ"] = -1.0*v_["E"*string(J)]
            pbm.gconst[ig_["C"*string(J)]] = Float64(v_["-EJ"])
        end
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(0.001),pb.n)
        pb.y0 = fill(Float64(0.001),pb.m)
        if haskey(ix_,"X7")
            pb.x0[ix_["X7"]] = Float64(60.0)
        else
            pb.y0[findfirst(x->x==ig_["X7"],pbm.congrps)] = Float64(60.0)
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eSQUARE", iet_)
        loaset(elftv,it,1,"XJ")
        it,iet_,_ = s2mpj_ii( "eCUBE", iet_)
        loaset(elftv,it,1,"XJ")
        it,iet_,_ = s2mpj_ii( "ePROD", iet_)
        loaset(elftv,it,1,"XI")
        loaset(elftv,it,2,"XJ")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for J = Int64(v_["1"]):Int64(v_["N"])
            v_["M+J"] = v_["M"]+J
            ename = "2D"*string(J)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eCUBE")
            arrset(ielftype,ie,iet_["eCUBE"])
            vname = "X"*string(Int64(v_["M+J"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.001))
            posev = findfirst(x->x=="XJ",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "3D"*string(J)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSQUARE")
            arrset(ielftype,ie,iet_["eSQUARE"])
            vname = "X"*string(Int64(v_["M+J"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.001))
            posev = findfirst(x->x=="XJ",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            for K = Int64(v_["1"]):Int64(v_["N"])
                v_["M+K"] = v_["M"]+K
                ename = "C"*string(K)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"ePROD")
                arrset(ielftype,ie,iet_["ePROD"])
                vname = "X"*string(Int64(v_["M+K"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.001))
                posev = findfirst(x->x=="XI",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(Int64(v_["M+J"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.001))
                posev = findfirst(x->x=="XJ",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for J = Int64(v_["1"]):Int64(v_["N"])
            v_["2DJ"] = 2.0*v_["D"*string(J)]
            ig = ig_["OBJ"]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["2D"*string(J)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["2DJ"]))
            v_["3DJ"] = 3.0*v_["D"*string(J)]
            ig = ig_["C"*string(J)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["3D"*string(J)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["3DJ"]))
            for K = Int64(v_["1"]):Int64(v_["N"])
                ig = ig_["OBJ"]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["C"*string(K)*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["C"*string(K)*","*string(J)]))
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               32.34867897
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
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
        pb.pbclass = "C-OQR2-AN-15-5"
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

    elseif action == "eSQUARE"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]*EV_[1]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 2.0e+0*EV_[1]
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 2.0e+0
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eCUBE"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]*EV_[1]*EV_[1]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 3.0e+0*EV_[1]*EV_[1]
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 6.0e+0*EV_[1]
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

