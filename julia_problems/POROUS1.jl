function POROUS1(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : POROUS1
#    *********
# 
#    The problem is to solve the porous medium equation on the unit square.
#    The equation is
# 
#        \Delta ( u^2 ) + d \frac{\partial}{\partial x_1}( u^3 ) + f = 0
# 
#    within the domain.  The boundary condition are that u = 1 on the bottom
#    and left sides and u = 0 on the top and right sides.  Discretization is
#    using the usual central differences. The function f is a point source of
#    maginitude 50 at the lower left grid point.  The initial approximation
#    is a discretization of 1 - x_1 x_2.
# 
#    Source: example 3.2.4 in
#    S. Eisenstat and H. Walker,
#    "Choosing the forcing terms in an inexact Newton method"
#    Report 6/94/75, Dept of Maths, Utah State University, 1994.
# 
#    SIF input: Ph. Toint, July 1994. Corrected November 2002.
# 
#    classification = "C-NOR2-MN-V-V"
# 
#    P is the number of points in one side of the unit square.
#    There are P*P variables.
# 
#       Alternative values for the SIF file parameters:
# IE P                   32             $-PARAMETER      original value
# IE P                   64             $-PARAMETER 
# IE P                   72             $-PARAMETER 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "POROUS1"

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
            v_["P"] = Int64(5);  #  SIF file default value
        else
            v_["P"] = Int64(args[1]);
        end
        if nargin<2
            v_["D"] = Float64(50.0);  #  SIF file default value
        else
            v_["D"] = Float64(args[2]);
        end
        v_["1"] = 1
        v_["2"] = 2
        v_["P-1"] = -1+v_["P"]
        v_["RP-1"] = Float64(v_["P-1"])
        v_["H"] = 1.0/v_["RP-1"]
        v_["H2"] = v_["H"]*v_["H"]
        v_["1/H2"] = 1.0/v_["H2"]
        v_["2H"] = 2.0*v_["H"]
        v_["D/2H"] = v_["D"]/v_["2H"]
        v_["-D/2H"] = -1.0*v_["D/2H"]
        v_["-4/H2"] = -4.0*v_["1/H2"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for J = Int64(v_["1"]):Int64(v_["P"])
            for I = Int64(v_["1"]):Int64(v_["P"])
                iv,ix_,_ = s2mpj_ii("U"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"U"*string(I)*","*string(J))
            end
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["2"]):Int64(v_["P-1"])
            for J = Int64(v_["2"]):Int64(v_["P-1"])
                ig,ig_,_ = s2mpj_ii("G"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"G"*string(I)*","*string(J))
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
        pbm.gconst[ig_["G"*string(Int64(v_["P-1"]))*","*string(Int64(v_["P-1"]))]]  = (
              Float64(-50.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        for J = Int64(v_["1"]):Int64(v_["P"])
            pb.xlower[ix_["U"*string(Int64(v_["1"]))*","*string(J)]] = 1.0
            pb.xupper[ix_["U"*string(Int64(v_["1"]))*","*string(J)]] = 1.0
            pb.xlower[ix_["U"*string(Int64(v_["P"]))*","*string(J)]] = 0.0
            pb.xupper[ix_["U"*string(Int64(v_["P"]))*","*string(J)]] = 0.0
        end
        for I = Int64(v_["2"]):Int64(v_["P-1"])
            pb.xlower[ix_["U"*string(I)*","*string(Int64(v_["P"]))]] = 1.0
            pb.xupper[ix_["U"*string(I)*","*string(Int64(v_["P"]))]] = 1.0
            pb.xlower[ix_["U"*string(I)*","*string(Int64(v_["1"]))]] = 0.0
            pb.xupper[ix_["U"*string(I)*","*string(Int64(v_["1"]))]] = 0.0
        end
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        for I = Int64(v_["1"]):Int64(v_["P"])
            for J = Int64(v_["1"]):Int64(v_["P"])
                v_["RI"] = Float64(I)
                v_["RJ"] = Float64(J)
                v_["I-1"] = -1.0+v_["RI"]
                v_["J-1"] = -1.0+v_["RJ"]
                v_["X1"] = v_["I-1"]*v_["H"]
                v_["X2"] = v_["J-1"]*v_["H"]
                v_["X1X2"] = v_["X1"]*v_["X2"]
                v_["MX1X2"] = -1.0*v_["X1X2"]
                v_["UIJ"] = 1.0+v_["MX1X2"]
                pb.x0[ix_["U"*string(I)*","*string(J)]] = Float64(v_["UIJ"])
            end
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eSQ", iet_)
        loaset(elftv,it,1,"U")
        it,iet_,_ = s2mpj_ii( "eCB", iet_)
        loaset(elftv,it,1,"U")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["P"])
            for J = Int64(v_["1"]):Int64(v_["P"])
                ename = "US"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eSQ")
                arrset(ielftype,ie,iet_["eSQ"])
                vname = "U"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="U",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                ename = "UC"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eCB")
                arrset(ielftype,ie,iet_["eCB"])
                vname = "U"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="U",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["2"]):Int64(v_["P-1"])
            v_["I-1"] = -1+I
            v_["I+1"] = 1+I
            for J = Int64(v_["2"]):Int64(v_["P-1"])
                v_["J-1"] = -1+J
                v_["J+1"] = 1+J
                ig = ig_["G"*string(I)*","*string(J)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["US"*string(Int64(v_["I+1"]))*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["1/H2"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["US"*string(Int64(v_["I-1"]))*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["1/H2"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["US"*string(I)*","*string(Int64(v_["J-1"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["1/H2"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["US"*string(I)*","*string(Int64(v_["J+1"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["1/H2"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["US"*string(I)*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-4/H2"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["UC"*string(Int64(v_["I+1"]))*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["D/2H"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["UC"*string(Int64(v_["I-1"]))*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-D/2H"]))
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN               0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pbm.A = spzeros(Float64,0,0)
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

    elseif action == "eSQ"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]*EV_[1]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[1]+EV_[1]
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 2.0
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eCB"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]*EV_[1]*EV_[1]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 3.0*EV_[1]*EV_[1]
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 6.0*EV_[1]
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

