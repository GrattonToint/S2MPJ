function BLOCKQP3(action,args...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : BLOCKQP3
#    *********
# 
#    A non-convex quadratic program with some structure.
# 
#    The objective function is of the form
#       sum (i=1,n) (i/n) ( x_i y_i ) + sum(j=1,b) z_j^2
#    (cf BLOCKQP1.SIF where the factor (i/n) is missing).
#    There are n equality constraints of the form
#         x_i + y_i + sum (j=1,b) z_j = b   
#    There is an inequality constraint of the form
#       sum(i=1,n) x_i + y_i + sum(j=1,b) z_j >= b + 1
#    Finally, there are simple bounds
#          1 <= x_i, y_i <= 1    (i=1,n)
#          0 <= z_j <= 2         (j=1,b)
# 
#    SIF input: Nick Gould, June 1994
# 
#    classification = "QLR2-AN-V-V"
# 
#    The number of equality constraints
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "BLOCKQP3"

    if action == "setup"
        pbm          = PBM(name)
        pb           = PB(name)
        pb.sifpbname = "BLOCKQP3"
        nargin       = length(args)
        pbm.call     = eval( Meta.parse( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        if nargin<1
            v_["N"] = 10;  #  SIF file default value
        else
            v_["N"] = args[1];
        end
        if nargin<2
            v_["B"] = 5;  #  SIF file default value
        else
            v_["B"] = args[2];
        end
        v_["1"] = 1
        v_["RB"] = v_["B"]
        v_["RB+1"] = 1+v_["RB"]
        v_["RN"] = v_["N"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        xscale  = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2x_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
            iv,ix_,_ = s2x_ii("Y"*string(I),ix_)
            arrset(pb.xnames,iv,"Y"*string(I))
        end
        for J = Int64(v_["1"]):Int64(v_["B"])
            iv,ix_,_ = s2x_ii("Z"*string(J),ix_)
            arrset(pb.xnames,iv,"Z"*string(J))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2x_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        for I = Int64(v_["1"]):Int64(v_["N"])
            ig,ig_,_ = s2x_ii(I,ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,I)
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += 1.0
            iv = ix_["Y"*string(I)]
            pbm.A[ig,iv] += 1.0
            ig,ig_,_ = s2x_ii("E"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"E"*string(I))
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += 1.0
            iv = ix_["Y"*string(I)]
            pbm.A[ig,iv] += 1.0
        end
        for J = Int64(v_["1"]):Int64(v_["B"])
            ig,ig_,_ = s2x_ii("I",ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"I")
            iv = ix_["Z"*string(J)]
            pbm.A[ig,iv] += 1.0
            for I = Int64(v_["1"]):Int64(v_["N"])
                ig,ig_,_ = s2x_ii("E"*string(I),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"E"*string(I))
                iv = ix_["Z"*string(J)]
                pbm.A[ig,iv] += 1.0
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
        pbm.congrps = findall(x->x!="<>",gtype)
        pb.nob = ngrp-pb.m
        pbm.objgrps = findall(x->x=="<>",gtype)
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        for I = Int64(v_["1"]):Int64(v_["N"])
            pbm.gconst[ig_["E"*string(I)]] = v_["RB"]
        end
        pbm.gconst[ig_["I"]] = v_["RB+1"]
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        for I = Int64(v_["1"]):Int64(v_["N"])
            pb.xlower[ix_["X"*string(I)]] = -1.0
            pb.xupper[ix_["X"*string(I)]] = 1.0
            pb.xlower[ix_["Y"*string(I)]] = -1.0
            pb.xupper[ix_["Y"*string(I)]] = 1.0
        end
        for J = Int64(v_["1"]):Int64(v_["B"])
            pb.xlower[ix_["Z"*string(J)]] = 0.0
            pb.xupper[ix_["Z"*string(J)]] = 2.0
        end
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(0.5,pb.n,)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2x_ii( "SQ", iet_)
        loaset(elftv,it,1,"Z")
        it,iet_,_ = s2x_ii( "PROD", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["N"])
            ename = "P"*string(I)
            ie,ie_,_  = s2x_ii(ename,ie_)
            arrset(pbm.elftype,ie,"PROD")
            arrset(ielftype, ie, iet_["PROD"])
            vname = "X"*string(I)
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,0.5)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "Y"*string(I)
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,0.5)
            posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        for I = Int64(v_["1"]):Int64(v_["B"])
            ename = "S"*string(I)
            ie,ie_,_  = s2x_ii(ename,ie_)
            arrset(pbm.elftype,ie,"SQ")
            arrset(ielftype, ie, iet_["SQ"])
            vname = "Z"*string(I)
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,0.5)
            posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            v_["RI"] = I
            v_["I/N"] = v_["RI"]/v_["RN"]
            ig = ig_["OBJ"]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["P"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,v_["I/N"])
        end
        for J = Int64(v_["1"]):Int64(v_["B"])
            ig = ig_["OBJ"]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["S"*string(J)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.clower[pb.nle+pb.neq+1:pb.m] = zeros(Float64,pb.nge)
        pb.cupper[1:pb.nge] = fill(Inf,pb.nge)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTIONS %%%%%%%
        lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "QLR2-AN-V-V"
        return pb, pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "SQ"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = 0.5*EV_[1]*EV_[1]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[1]
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 1.0
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "PROD"

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

