function HANGING(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HANGING
#    *********
# 
#    A catenary problem in 3 dimensions.  A rectangular grid is hung from its
#    4 corners under gravity.  The problem is to determine the resulting shape.
# 
#    Source:  
#    an example in a talk by Nesterova and Vial, LLN, 1994.
# 
#    SIF input: Ph. Toint, November 1994.
# 
#    classification = "C-LQR2-AY-V-V"
# 
#    dimension of the grid
# 
#       Alternative values for the SIF file parameters:
# IE NX                  3              $-PARAMETER n = 27
# IE NY                  3              $-PARAMETER
# 
# IE NX                  5              $-PARAMETER n = 90
# IE NY                  6              $-PARAMETER
# 
# IE NX                  10             $-PARAMETER n = 300  original value
# IE NY                  10             $-PARAMETER
# 
# IE NX                  20             $-PARAMETER n = 1800
# IE NY                  30             $-PARAMETER
# 
# IE NX                  40             $-PARAMETER n = 3600
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "HANGING"

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
            v_["NX"] = Int64(3);  #  SIF file default value
        else
            v_["NX"] = Int64(args[1]);
        end
# IE NY                  30             $-PARAMETER
        if nargin<2
            v_["NY"] = Int64(3);  #  SIF file default value
        else
            v_["NY"] = Int64(args[2]);
        end
        v_["LX"] = 1.8
        v_["LY"] = 1.8
        v_["1"] = 1
        v_["NX-1"] = -1+v_["NX"]
        v_["NY-1"] = -1+v_["NY"]
        v_["LX2"] = v_["LX"]*v_["LX"]
        v_["LY2"] = v_["LY"]*v_["LY"]
        v_["RNX"] = Float64(v_["NX"])
        v_["RNY"] = Float64(v_["NY"])
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["NX"])
            for J = Int64(v_["1"]):Int64(v_["NY"])
                iv,ix_,_ = s2mpj_ii("X"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"X"*string(I)*","*string(J))
                iv,ix_,_ = s2mpj_ii("Y"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"Y"*string(I)*","*string(J))
                iv,ix_,_ = s2mpj_ii("Z"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"Z"*string(I)*","*string(J))
            end
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["1"]):Int64(v_["NX"])
            for J = Int64(v_["1"]):Int64(v_["NY"])
                ig,ig_,_ = s2mpj_ii("OBJ",ig_)
                arrset(gtype,ig,"<>")
                iv = ix_["Z"*string(I)*","*string(J)]
                pbm.A[ig,iv] += Float64(1.0)
            end
        end
        for I = Int64(v_["1"]):Int64(v_["NX"])
            for J = Int64(v_["1"]):Int64(v_["NY-1"])
                ig,ig_,_ = s2mpj_ii("RC"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"<=")
                arrset(pb.cnames,ig,"RC"*string(I)*","*string(J))
            end
        end
        for I = Int64(v_["1"]):Int64(v_["NX-1"])
            for J = Int64(v_["1"]):Int64(v_["NY"])
                ig,ig_,_ = s2mpj_ii("DC"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"<=")
                arrset(pb.cnames,ig,"DC"*string(I)*","*string(J))
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
        for I = Int64(v_["1"]):Int64(v_["NX"])
            for J = Int64(v_["1"]):Int64(v_["NY-1"])
                pbm.gconst[ig_["RC"*string(I)*","*string(J)]] = Float64(v_["LX2"])
            end
        end
        for I = Int64(v_["1"]):Int64(v_["NX-1"])
            for J = Int64(v_["1"]):Int64(v_["NY"])
                pbm.gconst[ig_["DC"*string(I)*","*string(J)]] = Float64(v_["LY2"])
            end
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        pb.xlower[ix_["X"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))]] = 0.0
        pb.xupper[ix_["X"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))]] = 0.0
        pb.xlower[ix_["Y"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))]] = 0.0
        pb.xupper[ix_["Y"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))]] = 0.0
        pb.xlower[ix_["Z"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))]] = 0.0
        pb.xupper[ix_["Z"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))]] = 0.0
        pb.xlower[ix_["X"*string(Int64(v_["NX"]))*","*string(Int64(v_["1"]))]]  = (
              v_["RNX"])
        pb.xupper[ix_["X"*string(Int64(v_["NX"]))*","*string(Int64(v_["1"]))]]  = (
              v_["RNX"])
        pb.xlower[ix_["Y"*string(Int64(v_["NX"]))*","*string(Int64(v_["1"]))]] = 0.0
        pb.xupper[ix_["Y"*string(Int64(v_["NX"]))*","*string(Int64(v_["1"]))]] = 0.0
        pb.xlower[ix_["Z"*string(Int64(v_["NX"]))*","*string(Int64(v_["1"]))]] = 0.0
        pb.xupper[ix_["Z"*string(Int64(v_["NX"]))*","*string(Int64(v_["1"]))]] = 0.0
        pb.xlower[ix_["X"*string(Int64(v_["1"]))*","*string(Int64(v_["NY"]))]] = 0.0
        pb.xupper[ix_["X"*string(Int64(v_["1"]))*","*string(Int64(v_["NY"]))]] = 0.0
        pb.xlower[ix_["Y"*string(Int64(v_["1"]))*","*string(Int64(v_["NY"]))]]  = (
              v_["RNY"])
        pb.xupper[ix_["Y"*string(Int64(v_["1"]))*","*string(Int64(v_["NY"]))]]  = (
              v_["RNY"])
        pb.xlower[ix_["Z"*string(Int64(v_["1"]))*","*string(Int64(v_["NY"]))]] = 0.0
        pb.xupper[ix_["Z"*string(Int64(v_["1"]))*","*string(Int64(v_["NY"]))]] = 0.0
        pb.xlower[ix_["X"*string(Int64(v_["NX"]))*","*string(Int64(v_["NY"]))]]  = (
              v_["RNX"])
        pb.xupper[ix_["X"*string(Int64(v_["NX"]))*","*string(Int64(v_["NY"]))]]  = (
              v_["RNX"])
        pb.xlower[ix_["Y"*string(Int64(v_["NX"]))*","*string(Int64(v_["NY"]))]]  = (
              v_["RNY"])
        pb.xupper[ix_["Y"*string(Int64(v_["NX"]))*","*string(Int64(v_["NY"]))]]  = (
              v_["RNY"])
        pb.xlower[ix_["Z"*string(Int64(v_["NX"]))*","*string(Int64(v_["NY"]))]]  = (
              0.0)
        pb.xupper[ix_["Z"*string(Int64(v_["NX"]))*","*string(Int64(v_["NY"]))]]  = (
              0.0)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        for I = Int64(v_["1"]):Int64(v_["NX"])
            v_["I-1"] = -1+I
            v_["RI-1"] = Float64(v_["I-1"])
            for J = Int64(v_["1"]):Int64(v_["NY"])
                v_["J-1"] = -1+J
                v_["RJ-1"] = Float64(v_["J-1"])
                pb.x0[ix_["X"*string(I)*","*string(J)]] = Float64(v_["RI-1"])
                pb.x0[ix_["Y"*string(I)*","*string(J)]] = Float64(v_["RJ-1"])
            end
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eISQ", iet_)
        loaset(elftv,it,1,"XX")
        loaset(elftv,it,2,"YY")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for J = Int64(v_["1"]):Int64(v_["NY-1"])
            v_["J+1"] = 1+J
            for I = Int64(v_["1"]):Int64(v_["NX"])
                ename = "RX"*string(I)*","*string(J)
                ie,ie_,newelt = s2mpj_ii(ename,ie_)
                if newelt > 0
                    arrset(pbm.elftype,ie,"eISQ")
                    arrset(ielftype,ie,iet_["eISQ"])
                end
                vname = "X"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="XX",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(I)*","*string(Int64(v_["J+1"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="YY",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                ename = "RY"*string(I)*","*string(J)
                ie,ie_,newelt = s2mpj_ii(ename,ie_)
                if newelt > 0
                    arrset(pbm.elftype,ie,"eISQ")
                    arrset(ielftype,ie,iet_["eISQ"])
                end
                vname = "Y"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="XX",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "Y"*string(I)*","*string(Int64(v_["J+1"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="YY",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                ename = "RZ"*string(I)*","*string(J)
                ie,ie_,newelt = s2mpj_ii(ename,ie_)
                if newelt > 0
                    arrset(pbm.elftype,ie,"eISQ")
                    arrset(ielftype,ie,iet_["eISQ"])
                end
                vname = "Z"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="XX",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "Z"*string(I)*","*string(Int64(v_["J+1"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="YY",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
        end
        for I = Int64(v_["1"]):Int64(v_["NX-1"])
            v_["I+1"] = 1+I
            for J = Int64(v_["1"]):Int64(v_["NY"])
                ename = "DX"*string(I)*","*string(J)
                ie,ie_,newelt = s2mpj_ii(ename,ie_)
                if newelt > 0
                    arrset(pbm.elftype,ie,"eISQ")
                    arrset(ielftype,ie,iet_["eISQ"])
                end
                vname = "X"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="XX",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(Int64(v_["I+1"]))*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="YY",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                ename = "DY"*string(I)*","*string(J)
                ie,ie_,newelt = s2mpj_ii(ename,ie_)
                if newelt > 0
                    arrset(pbm.elftype,ie,"eISQ")
                    arrset(ielftype,ie,iet_["eISQ"])
                end
                vname = "Y"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="XX",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "Y"*string(Int64(v_["I+1"]))*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="YY",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                ename = "DZ"*string(I)*","*string(J)
                ie,ie_,newelt = s2mpj_ii(ename,ie_)
                if newelt > 0
                    arrset(pbm.elftype,ie,"eISQ")
                    arrset(ielftype,ie,iet_["eISQ"])
                end
                vname = "Z"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="XX",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "Z"*string(Int64(v_["I+1"]))*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="YY",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["NX"])
            for J = Int64(v_["1"]):Int64(v_["NY-1"])
                ig = ig_["RC"*string(I)*","*string(J)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["RX"*string(I)*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,1.)
                posel = posel+1
                loaset(pbm.grelt,ig,posel,ie_["RY"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel, 1.)
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["RZ"*string(I)*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,1.)
            end
        end
        for I = Int64(v_["1"]):Int64(v_["NX-1"])
            for J = Int64(v_["1"]):Int64(v_["NY"])
                ig = ig_["DC"*string(I)*","*string(J)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["DX"*string(I)*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,1.)
                posel = posel+1
                loaset(pbm.grelt,ig,posel,ie_["DY"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel, 1.)
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["DZ"*string(I)*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,1.)
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN(3,3)          -6.1184107487
# LO SOLTN(5,6)          -77.260229515
# LO SOLTN(10,10)        -620.17603242
# LO SOLTN(20,30)        -1025.4292887
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.cupper[1:pb.nle] = zeros(Float64,pb.nle)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-LQR2-AY-V-V"
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

    elseif action == "eISQ"

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

