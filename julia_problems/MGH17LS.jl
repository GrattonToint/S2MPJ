function MGH17LS(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : MGH17LS
#    *********
# 
#    NIST Data fitting problem MGH17.
# 
#    Fit: y = b1 + b2*exp[-x*b4] + b3*exp[-x*b5] + e
# 
#    Source:  Problem from the NIST nonlinear regression test set
#      http://www.itl.nist.gov/div898/strd/nls/nls_main.shtml
# 
#    Reference: Osborne, M. R. (1972).  
#     Some aspects of nonlinear least squares calculations.  
#     In Numerical Methods for Nonlinear Optimization, Lootsma (Ed).  
#     New York, NY:  Academic Press, pp. 171-189.
# 
#    SIF input: Nick Gould and Tyrone Rees, Oct 2015
# 
#    classification = "C-SUR2-MN-5-0"
# 
#    Number of data values
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "MGH17LS"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["M"] = 33
        v_["N"] = 5
        v_["1"] = 1
        v_["X1"] = 0.0E+0
        v_["X2"] = 1.0E+1
        v_["X3"] = 2.0E+1
        v_["X4"] = 3.0E+1
        v_["X5"] = 4.0E+1
        v_["X6"] = 5.0E+1
        v_["X7"] = 6.0E+1
        v_["X8"] = 7.0E+1
        v_["X9"] = 8.0E+1
        v_["X10"] = 9.0E+1
        v_["X11"] = 1.0E+2
        v_["X12"] = 1.1E+2
        v_["X13"] = 1.2E+2
        v_["X14"] = 1.3E+2
        v_["X15"] = 1.4E+2
        v_["X16"] = 1.5E+2
        v_["X17"] = 1.6E+2
        v_["X18"] = 1.7E+2
        v_["X19"] = 1.8E+2
        v_["X20"] = 1.9E+2
        v_["X21"] = 2.0E+2
        v_["X22"] = 2.1E+2
        v_["X23"] = 2.2E+2
        v_["X24"] = 2.3E+2
        v_["X25"] = 2.4E+2
        v_["X26"] = 2.5E+2
        v_["X27"] = 2.6E+2
        v_["X28"] = 2.7E+2
        v_["X29"] = 2.8E+2
        v_["X30"] = 2.9E+2
        v_["X31"] = 3.0E+2
        v_["X32"] = 3.1E+2
        v_["X33"] = 3.2E+2
        v_["Y1"] = 8.44E-1
        v_["Y2"] = 9.08E-1
        v_["Y3"] = 9.32E-1
        v_["Y4"] = 9.36E-1
        v_["Y5"] = 9.25E-1
        v_["Y6"] = 9.08E-1
        v_["Y7"] = 8.81E-1
        v_["Y8"] = 8.50E-1
        v_["Y9"] = 8.18E-1
        v_["Y10"] = 7.84E-1
        v_["Y11"] = 7.51E-1
        v_["Y12"] = 7.18E-1
        v_["Y13"] = 6.85E-1
        v_["Y14"] = 6.58E-1
        v_["Y15"] = 6.28E-1
        v_["Y16"] = 6.03E-1
        v_["Y17"] = 5.80E-1
        v_["Y18"] = 5.58E-1
        v_["Y19"] = 5.38E-1
        v_["Y20"] = 5.22E-1
        v_["Y21"] = 5.06E-1
        v_["Y22"] = 4.90E-1
        v_["Y23"] = 4.78E-1
        v_["Y24"] = 4.67E-1
        v_["Y25"] = 4.57E-1
        v_["Y26"] = 4.48E-1
        v_["Y27"] = 4.38E-1
        v_["Y28"] = 4.31E-1
        v_["Y29"] = 4.24E-1
        v_["Y30"] = 4.20E-1
        v_["Y31"] = 4.14E-1
        v_["Y32"] = 4.11E-1
        v_["Y33"] = 4.06E-1
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("B"*string(I),ix_)
            arrset(pb.xnames,iv,"B"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["1"]):Int64(v_["M"])
            ig,ig_,_ = s2mpj_ii("F"*string(I),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["B1"]
            pbm.A[ig,iv] += Float64(1.0)
        end
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        for I = Int64(v_["1"]):Int64(v_["M"])
            pbm.gconst[ig_["F"*string(I)]] = Float64(v_["Y"*string(I)])
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.x0[ix_["B1"]] = Float64(50.0)
        pb.x0[ix_["B2"]] = Float64(150.0)
        pb.x0[ix_["B3"]] = Float64(-100.0)
        pb.x0[ix_["B4"]] = Float64(1.0)
        pb.x0[ix_["B5"]] = Float64(2.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eE2", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"X")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["M"])
            ename = "EA"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eE2")
            arrset(ielftype,ie,iet_["eE2"])
            vname = "B2"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "B4"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="X",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["X"*string(I)]))
            ename = "EB"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eE2")
            arrset(ielftype,ie,iet_["eE2"])
            vname = "B3"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "B5"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="X",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["X"*string(I)]))
        end
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gL2",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for ig in 1:ngrp
            arrset(pbm.grftype,ig,"gL2")
        end
        for I = Int64(v_["1"]):Int64(v_["M"])
            ig = ig_["F"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["EA"*string(I)])
            loaset(pbm.grelw,ig,posel,1.)
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["EB"*string(I)])
            loaset(pbm.grelw,ig,posel,1.)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
        pb.objlower = 0.0
#    Solution
# LO SOLTN               
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-SUR2-MN-5-0"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        return pb, pbm

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "eE2"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        E = exp(-EV_[2]*pbm.elpar[iel_][1])
        f_   = EV_[1]*E
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = E
            g_[2] = -EV_[1]*pbm.elpar[iel_][1]*E
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = -pbm.elpar[iel_][1]*E
                H_[2,1] = H_[1,2]
                H_[2,2] = EV_[1]*pbm.elpar[iel_][1]*pbm.elpar[iel_][1]*E
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

    elseif action == "gL2"

        GVAR_   = args[1]
        igr_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_= GVAR_*GVAR_
        if nargout>1
            g_ = GVAR_+GVAR_
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_ = 2.0
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

