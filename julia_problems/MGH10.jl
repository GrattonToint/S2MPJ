function MGH10(action::String,args::Union{Any}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : MGH10
#    *********
# 
#    NIST Data fitting problem MGH10 given as an inconsistent set of
#    nonlinear equations.
# 
#    Fit: y = b1 * exp[b2/(x+b3)] + e
# 
#    Source:  Problem from the NIST nonlinear regression test set
#      http://www.itl.nist.gov/div898/strd/nls/nls_main.shtml
# 
#    Reference: Meyer, R. R. (1970).  
#      Theoretical and computational aspects of nonlinear 
#      regression.  In Nonlinear Programming, Rosen, 
#      Mangasarian and Ritter (Eds).  
#      New York, NY: Academic Press, pp. 465-486.
# 
#    SIF input: Nick Gould and Tyrone Rees, Oct 2015
#               correction by S. Gratton & Ph. Toint, May 2024
# 
#    classification = "C-CNOR2-MN-3-16"
# 
#    Number of data values
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 21 VI 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "MGH10"
    if ( !isdefined(@__MODULE__, :s2mpj_ii) )
        error( "Please include(\"s2mpjlib.jl\") using \"s2mpjlib.jl\" from the S2MPJ distribution before calling MGH10.")
    end

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["M"] = 16
        v_["N"] = 3
        v_["1"] = 1
        v_["X1"] = 5.00E+01
        v_["X2"] = 5.50E+01
        v_["X3"] = 6.00E+01
        v_["X4"] = 6.50E+01
        v_["X5"] = 7.00E+01
        v_["X6"] = 7.50E+01
        v_["X7"] = 8.00E+01
        v_["X8"] = 8.50E+01
        v_["X9"] = 9.00E+01
        v_["X10"] = 9.50E+01
        v_["X11"] = 1.00E+02
        v_["X12"] = 1.05E+02
        v_["X13"] = 1.10E+02
        v_["X14"] = 1.15E+02
        v_["X15"] = 1.20E+02
        v_["X16"] = 1.25E+02
        v_["Y1"] = 3.478E+04
        v_["Y2"] = 2.861E+04
        v_["Y3"] = 2.365E+04
        v_["Y4"] = 1.963E+04
        v_["Y5"] = 1.637E+04
        v_["Y6"] = 1.372E+04
        v_["Y7"] = 1.154E+04
        v_["Y8"] = 9.744E+03
        v_["Y9"] = 8.261E+03
        v_["Y10"] = 7.030E+03
        v_["Y11"] = 6.005E+03
        v_["Y12"] = 5.147E+03
        v_["Y13"] = 4.427E+03
        v_["Y14"] = 3.820E+03
        v_["Y15"] = 3.307E+03
        v_["Y16"] = 2.872E+03
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        irA   = Int64[]
        icA   = Int64[]
        valA  = Float64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("B"*string(I),ix_)
            arrset(pb.xnames,iv,"B"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype = String[]
        for I = Int64(v_["1"]):Int64(v_["M"])
            ig,ig_,_ = s2mpj_ii("F"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"F"*string(I))
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
        for I = Int64(v_["1"]):Int64(v_["M"])
            pbm.gconst[ig_["F"*string(I)]] = Float64(v_["Y"*string(I)])
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        if haskey(ix_,"B1")
            pb.x0[ix_["B1"]] = Float64(2.0)
        else
            pb.y0[findfirst(x->x==ig_["B1"],pbm.congrps)] = Float64(2.0)
        end
        if haskey(ix_,"B2")
            pb.x0[ix_["B2"]] = Float64(400000.0)
        else
            pb.y0[findfirst(x->x==ig_["B2"],pbm.congrps)] = Float64(400000.0)
        end
        if haskey(ix_,"B3")
            pb.x0[ix_["B3"]] = Float64(25000.0)
        else
            pb.y0[findfirst(x->x==ig_["B3"],pbm.congrps)] = Float64(25000.0)
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eE12", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        loaset(elftv,it,3,"V3")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"X")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["M"])
            ename = "E"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eE12")
            arrset(ielftype,ie,iet_["eE12"])
            vname = "B1"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "B2"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "B3"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="X",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["X"*string(I)]))
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["M"])
            ig = ig_["F"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["E"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
        pb.objlower = 0.0
#    Solution
# LO SOLTN               
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-CNOR2-MN-3-16"
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

    elseif action == "eE12"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        V3PX = EV_[3]+pbm.elpar[iel_][1]
        V3PX2 = V3PX*V3PX
        V3PX3 = V3PX*V3PX2
        E = exp(EV_[2]/V3PX)
        F = EV_[1]*E
        f_   = F
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = E
            g_[2] = F/V3PX
            g_[3] = -F*EV_[2]/V3PX2
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,2] = E/V3PX
                H_[2,1] = H_[1,2]
                H_[1,3] = -EV_[2]*E/V3PX2
                H_[3,1] = H_[1,3]
                H_[2,2] = F/V3PX2
                H_[2,3] = -F*(1.0/V3PX2+EV_[2]/V3PX3)
                H_[3,2] = H_[2,3]
                H_[3,3] = F*EV_[2]*(EV_[2]/V3PX^4+2.0/V3PX3)
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

