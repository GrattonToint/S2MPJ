function RAYBENDL(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    A ray bending problem.  A ray across a inhomogeneous 2D medium is
#    represented by a piecewise linear curve whose knots can be chosen.  
#    The problem is then to optimize the positions of these knots in order 
#    to obtain a ray path corresponding to the minimum travel time from 
#    source to receiver,  according to Fermat principle.
# 
#    The problem becomes harder and harder when the dimesnion increases
#    because the knots are getting closer and closer and the objective
#    has a nondifferentiable kink when two knots coincide.  The difficulty
#    is less apparent when exact second derivatives are not used.
# 
#    Source: a test example in
#    T.J. Moser, G. Nolet and R. Snieder,
#    "Ray bending revisited",
#    Bulletin of the Seism. Society of America 21(1).
# 
#    SIF input: Ph Toint, Dec 1991.
# 
#    classification = "C-OXR2-MY-V-0"
# 
#    number of  knots  ( >= 4 )
#    ( n = 2( NKNOTS - 1 ) ) 
# 
#       Alternative values for the SIF file parameters:
# IE NKNOTS              4              $-PARAMETER n = 6
# IE NKNOTS              11             $-PARAMETER n = 20
# IE NKNOTS              21             $-PARAMETER n = 40     original value
# IE NKNOTS              32             $-PARAMETER n = 62
# IE NKNOTS              64             $-PARAMETER n = 126
# IE NKNOTS              512            $-PARAMETER n = 1022
# IE NKNOTS              1024           $-PARAMETER n = 2046
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "RAYBENDL"

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
            v_["NKNOTS"] = Int64(4);  #  SIF file default value
        else
            v_["NKNOTS"] = Int64(args[1]);
        end
        v_["XSRC"] = 0.0
        v_["ZSRC"] = 0.0
        v_["XRCV"] = 100.0
        v_["ZRCV"] = 100.0
        v_["NK-1"] = -1+v_["NKNOTS"]
        v_["NK-2"] = -2+v_["NKNOTS"]
        v_["0"] = 0
        v_["1"] = 1
        v_["2"] = 2
        v_["3"] = 3
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["0"]):Int64(v_["NKNOTS"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
            iv,ix_,_ = s2mpj_ii("Z"*string(I),ix_)
            arrset(pb.xnames,iv,"Z"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["1"]):Int64(v_["NKNOTS"])
            ig,ig_,_ = s2mpj_ii("TIME"*string(I),ig_)
            arrset(gtype,ig,"<>")
            arrset(pbm.gscale,ig,Float64(2.0))
        end
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        pb.xlower[ix_["X"*string(Int64(v_["0"]))]] = v_["XSRC"]
        pb.xupper[ix_["X"*string(Int64(v_["0"]))]] = v_["XSRC"]
        pb.xlower[ix_["Z"*string(Int64(v_["0"]))]] = v_["ZSRC"]
        pb.xupper[ix_["Z"*string(Int64(v_["0"]))]] = v_["ZSRC"]
        pb.xlower[ix_["X"*string(Int64(v_["NKNOTS"]))]] = v_["XRCV"]
        pb.xupper[ix_["X"*string(Int64(v_["NKNOTS"]))]] = v_["XRCV"]
        pb.xlower[ix_["Z"*string(Int64(v_["NKNOTS"]))]] = v_["ZRCV"]
        pb.xupper[ix_["Z"*string(Int64(v_["NKNOTS"]))]] = v_["ZRCV"]
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        v_["XRANGE"] = v_["XRCV"]-v_["XSRC"]
        v_["ZRANGE"] = v_["ZRCV"]-v_["ZSRC"]
        v_["RKNOTS"] = Float64(v_["NKNOTS"])
        for I = Int64(v_["0"]):Int64(v_["NKNOTS"])
            v_["REALI"] = Float64(I)
            v_["FRAC"] = v_["REALI"]/v_["RKNOTS"]
            v_["XINCR"] = v_["FRAC"]*v_["XRANGE"]
            v_["ZINCR"] = v_["FRAC"]*v_["ZRANGE"]
            v_["XC"] = v_["XSRC"]+v_["XINCR"]
            v_["ZC"] = v_["ZSRC"]+v_["ZINCR"]
            pb.x0[ix_["X"*string(I)]] = Float64(v_["XC"])
            pb.x0[ix_["Z"*string(I)]] = Float64(v_["ZC"])
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eTT", iet_)
        loaset(elftv,it,1,"X1")
        loaset(elftv,it,2,"X2")
        loaset(elftv,it,3,"Z1")
        loaset(elftv,it,4,"Z2")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["NKNOTS"])
            v_["I-1"] = -1+I
            ename = "T"*string(I)
            ie,ie_,newelt = s2mpj_ii(ename,ie_)
            if newelt > 0
                arrset(pbm.elftype,ie,"eTT")
                arrset(ielftype,ie,iet_["eTT"])
            end
            vname = "X"*string(Int64(v_["I-1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "Z"*string(Int64(v_["I-1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="Z1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "Z"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="Z2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["NKNOTS"])
            ig = ig_["TIME"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["T"*string(I)])
            loaset(pbm.grelw,ig,posel,1.)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#   Solution of the continuous problem
# LO RAYBENDL            96.2424
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.A = spzeros(Float64,0,0)
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-OXR2-MY-V-0"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        return pb, pbm

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "e_globs"

        pbm = args[1]
        arrset(pbm.efpar,1,0.01)
        return pbm

    elseif action == "eTT"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,3,4)
        IV_ =  zeros(Float64,3)
        U_[3,1] = U_[3,1]-1
        U_[3,2] = U_[3,2]+1
        U_[1,3] = U_[1,3]+1
        U_[2,4] = U_[2,4]+1
        IV_[1] = dot(U_[1,:],EV_)
        IV_[2] = dot(U_[2,:],EV_)
        IV_[3] = dot(U_[3,:],EV_)
        C0 = 1.0+pbm.efpar[1]*IV_[1]
        C1 = 1.0+pbm.efpar[1]*IV_[2]
        DCDZ = pbm.efpar[1]
        V = 1.0/C1+1.0/C0
        VDZ0 = -DCDZ/(C0*C0)
        VDZ1 = -DCDZ/(C1*C1)
        VDZ0Z0 = 2.0*DCDZ*DCDZ/C0^3
        VDZ1Z1 = 2.0*DCDZ*DCDZ/C1^3
        DZ1 = IV_[2]-IV_[1]
        R = sqrt(IV_[3]*IV_[3]+DZ1*DZ1)
        RDX = IV_[3]/R
        RDZ1 = DZ1/R
        RDZ0 = -RDZ1
        RDXDX = (1.0-IV_[3]*IV_[3]/(R*R))/R
        RDXZ1 = -IV_[3]*DZ1/R^3
        RDXZ0 = -RDXZ1
        RDZ1Z1 = (1.0-DZ1*DZ1/(R*R))/R
        RDZ0Z0 = RDZ1Z1
        RDZ0Z1 = -RDZ1Z1
        f_   = V*R
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[3] = V*RDX
            g_[1] = V*RDZ0+VDZ0*R
            g_[2] = V*RDZ1+VDZ1*R
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[3,3] = V*RDXDX
                H_[3,1] = VDZ0*RDX+V*RDXZ0
                H_[1,3] = H_[3,1]
                H_[3,2] = VDZ1*RDX+V*RDXZ1
                H_[2,3] = H_[3,2]
                H_[1,1] = V*RDZ0Z0+VDZ0Z0*R+2.0*VDZ0*RDZ0
                H_[1,2] = V*RDZ0Z1+VDZ1*RDZ0+VDZ0*RDZ1
                H_[2,1] = H_[1,2]
                H_[2,2] = V*RDZ1Z1+VDZ1Z1*R+2.0*VDZ1*RDZ1
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
            pbm.has_globs = [1,0]
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

