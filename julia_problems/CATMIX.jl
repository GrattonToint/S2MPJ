function CATMIX(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : CATMIX
#    *********
# 
#    Determine the optimal mixing policy of two catalysts along the
#    length of a tubular plug flow reactor involving several reactions
# 
#    This is problem 14 in the COPS (Version 2) collection of 
#    E. Dolan and J. More'
#    see "Benchmarking Optimization Software with COPS"
#    Argonne National Labs Technical Report ANL/MCS-246 (2000)
# 
#    SIF input: Nick Gould, November 2000
# 
#    classification = "C-OOR2-AN-V-V"
# 
#    The number of subintervals
# 
# 
#       Alternative values for the SIF file parameters:
# IE NH                  100            $-PARAMETER
# IE NH                  200            $-PARAMETER
# IE NH                  400            $-PARAMETER
# IE NH                  800            $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "CATMIX"

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
            v_["NH"] = Int64(10);  #  SIF file default value
        else
            v_["NH"] = Int64(args[1]);
        end
# IE NH                  3000           $-PARAMETER
# IE NH                  30000          $-PARAMETER
        v_["TF"] = 1.0
        v_["X1u0"] = 1.0
        v_["X2u0"] = 0.0
        v_["ALPHA"] = 0.0
        v_["RNH"] = Float64(v_["NH"])
        v_["H"] = v_["TF"]/v_["RNH"]
        v_["0"] = 0
        v_["1"] = 1
        v_["NH-1"] = -1+v_["NH"]
        v_["ALPHAH"] = v_["ALPHA"]*v_["H"]
        v_["H/2"] = 0.5*v_["H"]
        v_["-H/2"] = -1.0*v_["H/2"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["0"]):Int64(v_["NH"])
            iv,ix_,_ = s2mpj_ii("U"*string(I),ix_)
            arrset(pb.xnames,iv,"U"*string(I))
            iv,ix_,_ = s2mpj_ii("X1"*string(I),ix_)
            arrset(pb.xnames,iv,"X1"*string(I))
            iv,ix_,_ = s2mpj_ii("X2"*string(I),ix_)
            arrset(pb.xnames,iv,"X2"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X1"*string(Int64(v_["NH"]))]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X2"*string(Int64(v_["NH"]))]
        pbm.A[ig,iv] += Float64(1.0)
        for I = Int64(v_["0"]):Int64(v_["NH-1"])
            v_["I+1"] = 1+I
            ig,ig_,_ = s2mpj_ii("ODE1"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"ODE1"*string(I))
            iv = ix_["X1"*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["X1"*string(Int64(v_["I+1"]))]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("ODE2"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"ODE2"*string(I))
            iv = ix_["X2"*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["X2"*string(Int64(v_["I+1"]))]
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
        pbm.gconst[ig_["OBJ"]] = Float64(1.0)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        for I = Int64(v_["0"]):Int64(v_["NH"])
            pb.xlower[ix_["U"*string(I)]] = 0.0
            pb.xupper[ix_["U"*string(I)]] = 1.0
        end
        pb.xlower[ix_["X1"*string(Int64(v_["0"]))]] = v_["X1u0"]
        pb.xupper[ix_["X1"*string(Int64(v_["0"]))]] = v_["X1u0"]
        pb.xlower[ix_["X2"*string(Int64(v_["0"]))]] = v_["X2u0"]
        pb.xupper[ix_["X2"*string(Int64(v_["0"]))]] = v_["X2u0"]
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        for I = Int64(v_["0"]):Int64(v_["NH"])
            if haskey(ix_,"U"*string(I))
                pb.x0[ix_["U"*string(I)]] = Float64(0.0)
            else
                pb.y0[findfirst(x->x==ig_["U"*string(I)],pbm.congrps)] = Float64(0.0)
            end
            if haskey(ix_,"X1"*string(I))
                pb.x0[ix_["X1"*string(I)]] = Float64(1.0)
            else
                pb.y0[findfirst(x->x==ig_["X1"*string(I)],pbm.congrps)] = Float64(1.0)
            end
            if haskey(ix_,"X2"*string(I))
                pb.x0[ix_["X2"*string(I)]] = Float64(0.0)
            else
                pb.y0[findfirst(x->x==ig_["X2"*string(I)],pbm.congrps)] = Float64(0.0)
            end
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eDIFSQ", iet_)
        loaset(elftv,it,1,"U1")
        loaset(elftv,it,2,"U2")
        it,iet_,_ = s2mpj_ii( "eP1", iet_)
        loaset(elftv,it,1,"U")
        loaset(elftv,it,2,"X1")
        loaset(elftv,it,3,"X2")
        it,iet_,_ = s2mpj_ii( "eP2", iet_)
        loaset(elftv,it,1,"U")
        loaset(elftv,it,2,"X")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["0"]):Int64(v_["NH-1"])
            v_["I+1"] = 1+I
            ename = "O"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eDIFSQ")
            arrset(ielftype,ie,iet_["eDIFSQ"])
            vname = "U"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="U1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "U"*string(Int64(v_["I+1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="U2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        for I = Int64(v_["0"]):Int64(v_["NH"])
            ename = "P1"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eP1")
            arrset(ielftype,ie,iet_["eP1"])
            vname = "U"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="U",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X1"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X2"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "P2"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eP2")
            arrset(ielftype,ie,iet_["eP2"])
            vname = "U"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="U",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X2"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["0"]):Int64(v_["NH-1"])
            ig = ig_["OBJ"]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["O"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["ALPHAH"]))
            ig = ig_["ODE1"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["P1"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["H/2"]))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["P1"*string(Int64(v_["I+1"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["H/2"]))
            ig = ig_["ODE2"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["P1"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["-H/2"]))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["P1"*string(Int64(v_["I+1"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["-H/2"]))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["P2"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["H/2"]))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["P2"*string(Int64(v_["I+1"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["H/2"]))
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLUTION             -4.7748D-02   $ (NH=100)
# LO SOLUTION             -4.8016D-02   $ (NH=200)
# LO SOLUTION             -4.7862D-02   $ (NH=400)
# LO SOLUTION             -4.7185D-02   $ (NH=800)
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
        pb.pbclass = "C-OOR2-AN-V-V"
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

    elseif action == "eDIFSQ"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,1,2)
        IV_ =  zeros(Float64,1)
        U_[1,1] = U_[1,1]+1
        U_[1,2] = U_[1,2]-1
        IV_[1] = dot(U_[1,:],EV_)
        f_   = IV_[1]*IV_[1]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = IV_[1]+IV_[1]
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 2.0
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

    elseif action == "eP1"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,2,3)
        IV_ =  zeros(Float64,2)
        U_[1,1] = U_[1,1]+1
        U_[2,2] = U_[2,2]-1
        U_[2,3] = U_[2,3]+10
        IV_[1] = dot(U_[1,:],EV_)
        IV_[2] = dot(U_[2,:],EV_)
        f_   = IV_[1]*IV_[2]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = IV_[2]
            g_[2] = IV_[1]
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = 1.0
                H_[2,1] = H_[1,2]
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

    elseif action == "eP2"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = (EV_[1]-1.0)*EV_[2]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[2]
            g_[2] = EV_[1]-1.0
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

