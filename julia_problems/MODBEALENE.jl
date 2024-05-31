function MODBEALENE(action,args...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : MODBEALENE
#    *********
#    A variation on Beale's problem in 2 variables
#    This is a nonlinear equation variant of MODBEALE
# 
#    Source: An adaptation by Ph. Toint of Problem 5 in
#    J.J. More', B.S. Garbow and K.E. Hillstrom,
#    "Testing Unconstrained Optimization Software",
#    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
# 
#    See also Buckley#89.
#    SIF input: Ph. Toint, Mar 2003.
#               Nick Gould (nonlinear equation version), Jan 2019
# 
#    classification = "NOR2-AN-V-V"
# 
#    The number of variables is  2 * N/2
# 
#       Alternative values for the SIF file parameters:
# IE N/2                 1              $-PARAMETER     original value
# IE N/2                 2              $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "MODBEALENE"

    if action == "setup"
        pbm          = PBM(name)
        pb           = PB(name)
        pb.sifpbname = "MODBEALENE"
        nargin       = length(args)
        pbm.call     = eval( Meta.parse( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        if nargin<1
            v_["N/2"] = Int64(5);  #  SIF file default value
        else
            v_["N/2"] = Int64(args[1]);
        end
# IE N/2                 100            $-PARAMETER
# IE N/2                 1000           $-PARAMETER
# IE N/2                 10000          $-PARAMETER
        if nargin<2
            v_["ALPHA"] = Float64(50.0);  #  SIF file default value
        else
            v_["ALPHA"] = Float64(args[2]);
        end
        v_["1"] = 1
        v_["N"] = v_["N/2"]+v_["N/2"]
        v_["N/2-1"] = -1+v_["N/2"]
        v_["ALPHINV"] = 1.0/v_["ALPHA"]
        v_["RALPHINV"] = sqrt(v_["ALPHINV"])
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        xscale  = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for J = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2x_ii("X"*string(J),ix_)
            arrset(pb.xnames,iv,"X"*string(J))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["1"]):Int64(v_["N/2-1"])
            v_["I-1"] = -1+I
            v_["2I-1"] = v_["I-1"]+v_["I-1"]
            v_["J"] = 1+v_["2I-1"]
            v_["J+1"] = 1+v_["J"]
            v_["J+2"] = 2+v_["J"]
            ig,ig_,_ = s2x_ii("BA"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"BA"*string(I))
            ig,ig_,_ = s2x_ii("BB"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"BB"*string(I))
            ig,ig_,_ = s2x_ii("BC"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"BC"*string(I))
            ig,ig_,_ = s2x_ii("L"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"L"*string(I))
            iv = ix_["X"*string(Int64(v_["J+1"]))]
            pbm.A[ig,iv] += Float64(6.0)
            iv = ix_["X"*string(Int64(v_["J+2"]))]
            pbm.A[ig,iv] += Float64(-1.0)
            arrset(pbm.gscale,ig,Float64(v_["RALPHINV"]))
        end
        ig,ig_,_ = s2x_ii("BA"*string(Int64(v_["N/2"])),ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"BA"*string(Int64(v_["N/2"])))
        ig,ig_,_ = s2x_ii("BB"*string(Int64(v_["N/2"])),ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"BB"*string(Int64(v_["N/2"])))
        ig,ig_,_ = s2x_ii("BC"*string(Int64(v_["N/2"])),ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"BC"*string(Int64(v_["N/2"])))
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
        pbm.congrps = findall(x->x!="<>",gtype)
        pb.nob = ngrp-pb.m
        pbm.objgrps = findall(x->x=="<>",gtype)
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        for I = Int64(v_["1"]):Int64(v_["N/2"])
            pbm.gconst[ig_["BA"*string(I)]] = Float64(1.5)
            pbm.gconst[ig_["BB"*string(I)]] = Float64(2.25)
            pbm.gconst[ig_["BC"*string(I)]] = Float64(2.625)
        end
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(1.0),pb.n)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2x_ii( "ePRODB", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"POW")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["N/2"])
            v_["I-1"] = -1+I
            v_["2I-1"] = v_["I-1"]+v_["I-1"]
            v_["J"] = 1+v_["2I-1"]
            v_["J+1"] = 1+v_["J"]
            ename = "AE"*string(I)
            ie,ie_,newelt = s2x_ii(ename,ie_)
            if newelt > 0
                arrset(pbm.elftype,ie,"ePRODB")
                arrset(ielftype,ie,iet_["ePRODB"])
            end
            vname = "X"*string(Int64(v_["J"]))
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,1.0)
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["J+1"]))
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,1.0)
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="POW",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(1.0))
            ename = "BE"*string(I)
            ie,ie_,newelt = s2x_ii(ename,ie_)
            if newelt > 0
                arrset(pbm.elftype,ie,"ePRODB")
                arrset(ielftype,ie,iet_["ePRODB"])
            end
            vname = "X"*string(Int64(v_["J"]))
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,1.0)
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["J+1"]))
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,1.0)
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="POW",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(2.0))
            ename = "CE"*string(I)
            ie,ie_,newelt = s2x_ii(ename,ie_)
            if newelt > 0
                arrset(pbm.elftype,ie,"ePRODB")
                arrset(ielftype,ie,iet_["ePRODB"])
            end
            vname = "X"*string(Int64(v_["J"]))
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,1.0)
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["J+1"]))
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,1.0)
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="POW",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(3.0))
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N/2"])
            ig = ig_["BA"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["AE"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
            ig = ig_["BB"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["BE"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
            ig = ig_["BC"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["CE"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
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
        lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "NOR2-AN-V-V"
        return pb, pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "ePRODB"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        T = 1.0-EV_[2]^pbm.elpar[iel_][1]
        POWM1 = pbm.elpar[iel_][1]-1.0
        W = -pbm.elpar[iel_][1]*EV_[2]^POWM1
        f_   = EV_[1]*T
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = T
            g_[2] = EV_[1]*W
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = 0.0
                H_[1,2] = W
                H_[2,1] = H_[1,2]
                H_[2,2] = -EV_[1]*pbm.elpar[iel_][1]*POWM1*EV_[2]^(pbm.elpar[iel_][1]-2.0)
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

    elseif action in  ["fx","fgx","fgHx","cx","cJx","cJHx","cIx","cIJx","cIJHx","cIJxv","fHxv","cJxv","Lxy","Lgxy","LgHxy","LIxy","LIgxy","LIgHxy","LHxyv","LIHxyv"]

        pbm = args[1]
        if pbm.name == name
            pbm.has_globs = [0,0]
            return s2x_eval(action,args...)
        else
            println("ERROR: please run "*name*" with action = setup")
            return ntuple(i->undef,args[end])
        end

    else
        println("ERROR: unknown action "*action*" requested from "*name*"%s.jl")
        return ntuple(i->undef,args[end])
    end

end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
