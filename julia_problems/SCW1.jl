function SCW1(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : SCW1
#    *********
# 
#    Source: a discretization of an infinite-demsional problem proposed 
#    by Simon Chandler-Wilde (U. Reading):
# 
#    Given a function u in C[0,2 pi] with ||u||_infty <= 1, find the 
#    supremum of c^2(u) + s^2(u), where
#      c(u) = int_0^2 pi cos(t)u(t) dt and
#      s(u) = int_0^2 pi sin(t)u(t) dt      
# 
#    The discretized version ignores the required continuity, and 
#    posits a piecewise constant solution that oscilates between
#    plus and minus one. The anticipated solution is -16.
# 
#    SIF input: Nick Gould, July 2020
# 
#    classification = "C-SLR2-MN-V-V"
# 
#    Number of internal knots
# 
#       Alternative values for the SIF file parameters:
# IE K                   1              $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "SCW1"

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
            v_["K"] = Int64(7);  #  SIF file default value
        else
            v_["K"] = Int64(args[1]);
        end
# IE K                   3              $-PARAMETER
# IE K                   7              $-PARAMETER     original value
# IE K                   10             $-PARAMETER
# IE K                   100            $-PARAMETER
        v_["0"] = 0
        v_["1"] = 1
        v_["2"] = 2
        v_["K+1"] = 1+v_["K"]
        v_["RK+1"] = Float64(v_["K+1"])
        v_["PI/4"] = atan(1.0)
        v_["PI"] = 4.0*v_["PI/4"]
        v_["2PI"] = 2.0*v_["PI"]
        v_["2PI/K+1"] = v_["2PI"]/v_["RK+1"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["0"]):Int64(v_["K+1"])
            iv,ix_,_ = s2mpj_ii("T"*string(I),ix_)
            arrset(pb.xnames,iv,"T"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("S",ig_)
        arrset(gtype,ig,"<>")
        ig,ig_,_ = s2mpj_ii("C",ig_)
        arrset(gtype,ig,"<>")
        for I = Int64(v_["0"]):Int64(v_["K"])
            v_["I+1"] = 1+I
            ig,ig_,_ = s2mpj_ii("CON"*string(I),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"CON"*string(I))
            iv = ix_["T"*string(Int64(v_["I+1"]))]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["T"*string(I)]
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
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower[ix_["T"*string(Int64(v_["0"]))]] = 0.0
        pb.xupper[ix_["T"*string(Int64(v_["0"]))]] = 0.0
        for I = Int64(v_["1"]):Int64(v_["K"])
            pb.xlower[ix_["T"*string(I)]] = 0.0
            pb.xupper[ix_["T"*string(I)]] = v_["2PI"]
        end
        pb.xlower[ix_["T"*string(Int64(v_["K+1"]))]] = v_["2PI"]
        pb.xupper[ix_["T"*string(Int64(v_["K+1"]))]] = v_["2PI"]
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        if haskey(ix_,"T"*string(Int64(v_["1"])))
            pb.x0[ix_["T"*string(Int64(v_["1"]))]] = Float64(1.0)
        else
            pb.y0[findfirst(x->x==ig_["T"*string(Int64(v_["1"]))],pbm.congrps)]  = (
                  Float64(1.0))
        end
        for I = Int64(v_["2"]):Int64(v_["K"])
            v_["RI"] = Float64(I)
            v_["START"] = v_["RI"]*v_["2PI/K+1"]
            if haskey(ix_,"T"*string(I))
                pb.x0[ix_["T"*string(I)]] = Float64(0.0)
            else
                pb.y0[findfirst(x->x==ig_["T"*string(I)],pbm.congrps)] = Float64(0.0)
            end
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eSINT", iet_)
        loaset(elftv,it,1,"T")
        it,iet_,_ = s2mpj_ii( "eCOST", iet_)
        loaset(elftv,it,1,"T")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["0"]):Int64(v_["K+1"])
            ename = "S"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSINT")
            arrset(ielftype,ie,iet_["eSINT"])
            vname = "T"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="T",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "C"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eCOST")
            arrset(ielftype,ie,iet_["eCOST"])
            vname = "T"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="T",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gMAXSQ",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["C"]
        arrset(pbm.grftype,ig,"gMAXSQ")
        for I = Int64(v_["0"]):Int64(v_["2"]):Int64(v_["K"])
            v_["I+1"] = 1+I
            ig = ig_["C"]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["S"*string(Int64(v_["I+1"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(1.0))
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["S"*string(I)])
            loaset(pbm.grelw,ig,posel,Float64(-1.0))
        end
        for I = Int64(v_["1"]):Int64(v_["2"]):Int64(v_["K"])
            v_["I+1"] = 1+I
            ig = ig_["C"]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["S"*string(Int64(v_["I+1"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(-1.0))
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["S"*string(I)])
            loaset(pbm.grelw,ig,posel,Float64(1.0))
        end
        ig = ig_["S"]
        arrset(pbm.grftype,ig,"gMAXSQ")
        for I = Int64(v_["0"]):Int64(v_["2"]):Int64(v_["K"])
            v_["I+1"] = 1+I
            ig = ig_["S"]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["C"*string(Int64(v_["I+1"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(-1.0))
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["C"*string(I)])
            loaset(pbm.grelw,ig,posel,Float64(1.0))
        end
        for I = Int64(v_["1"]):Int64(v_["2"]):Int64(v_["K"])
            v_["I+1"] = 1+I
            ig = ig_["S"]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["C"*string(Int64(v_["I+1"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(1.0))
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["C"*string(I)])
            loaset(pbm.grelw,ig,posel,Float64(-1.0))
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
# LO SCW                 0.0
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
        pb.pbclass = "C-SLR2-MN-V-V"
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

    elseif action == "eSINT"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        S = sin(EV_[1])
        f_   = S
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = cos(EV_[1])
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = -S
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eCOST"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        C = cos(EV_[1])
        f_   = C
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = -sin(EV_[1])
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = -C
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

    elseif action == "gMAXSQ"

        GVAR_   = args[1]
        igr_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_= -GVAR_*GVAR_
        if nargout>1
            g_ = -GVAR_-GVAR_
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_ = -2.0e+0
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

