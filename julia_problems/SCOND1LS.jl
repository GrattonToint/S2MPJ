function SCOND1LS(action,args...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : SCOND1LS
#    *********
# 
#    The semiconductor problem by Rheinboldt, using a finite difference
#    approximation.
#    This is the least-squares version of problem SEMICON1.
# 
#    Source: problem 10 in
#    J.J. More',
#    "A collection of nonlinear model problems"
#    Proceedings of the AMS-SIAM Summer seminar on the Computational
#    Solution of Nonlinear Systems of Equations, Colorado, 1988.
#    Argonne National Laboratory MCS-P60-0289, 1989.
# 
#    SIF input: Ph. Toint, Dec 1989.
# 
#    classification = "SBR2-AN-V-V"
# 
#    N  = Number of discretized point inside the interval [a, b]
#    LN = Index of the last negative discretization point
#         (the interest is in the negative part)
# 
#       Alternative values for the SIF file parameters:
# IE N                   10             $-PARAMETER     original value
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "SCOND1LS"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        self.call    = eval( Meta.parse( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        if nargin<1
            v_["N"] = Int64(10);  #  SIF file default value
        else
            v_["N"] = Int64(args[1]);
        end
        if nargin<2
            v_["LN"] = Int64(9);  #  SIF file default value
        else
            v_["LN"] = Int64(args[2]);
        end
# IE N                   50             $-PARAMETER
# IE LN                  45             $-PARAMETER
# IE N                   100            $-PARAMETER
# IE LN                  90             $-PARAMETER
# IE N                   500            $-PARAMETER
# IE LN                  450            $-PARAMETER
# IE N                   1000           $-PARAMETER
# IE LN                  900            $-PARAMETER
# IE N                   5000           $-PARAMETER
# IE LN                  4500           $-PARAMETER
        if nargin<3
            v_["LAMBDA"] = Float64(1.0);  #  SIF file default value
        else
            v_["LAMBDA"] = Float64(args[3]);
        end
        v_["A"] = -0.00009
        v_["B"] = 0.00001
        v_["UA"] = 0.0
        v_["UB"] = 700.0
        v_["CA"] = 1.0e12
        v_["CB"] = 1.0e13
        v_["BETA"] = 40.0
        v_["LN+1"] = 1+v_["LN"]
        v_["N+1"] = 1+v_["N"]
        v_["-A"] = -1.0*v_["A"]
        v_["B-A"] = v_["B"]+v_["-A"]
        v_["RN+1"] = Float64(v_["N+1"])
        v_["TMP"] = 1.0/v_["RN+1"]
        v_["H"] = v_["B-A"]*v_["TMP"]
        v_["H2"] = v_["H"]*v_["H"]
        v_["LB"] = v_["LAMBDA"]*v_["BETA"]
        v_["H2CA"] = v_["H2"]*v_["CA"]
        v_["H2CB"] = v_["H2"]*v_["CB"]
        v_["LH2CA"] = v_["LAMBDA"]*v_["H2CA"]
        v_["LH2CB"] = v_["LAMBDA"]*v_["H2CB"]
        v_["LUA"] = v_["LAMBDA"]*v_["UA"]
        v_["LUB"] = v_["LAMBDA"]*v_["UB"]
        v_["ULW"] = -5.0+v_["LUA"]
        v_["UUP"] = 5.0+v_["LUB"]
        v_["-LB"] = -1.0*v_["LB"]
        v_["-LUB"] = -1.0*v_["LUB"]
        v_["-LH2CB"] = -1.0*v_["LH2CB"]
        v_["0"] = 0
        v_["1"] = 1
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["0"]):Int64(v_["N+1"])
            iv,ix_,_ = s2mpj_ii("U"*string(I),ix_)
            arrset(pb.xnames,iv,"U"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            v_["I+1"] = 1+I
            v_["I-1"] = -1+I
            ig,ig_,_ = s2mpj_ii("G"*string(I),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["U"*string(Int64(v_["I-1"]))]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["U"*string(I)]
            pbm.A[ig,iv] += Float64(-2.0)
            iv = ix_["U"*string(Int64(v_["I+1"]))]
            pbm.A[ig,iv] += Float64(1.0)
        end
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        for I = Int64(v_["1"]):Int64(v_["LN"])
            pbm.gconst[ig_["G"*string(I)]] = Float64(v_["LH2CA"])
        end
        for I = Int64(v_["LN+1"]):Int64(v_["N"])
            pbm.gconst[ig_["G"*string(I)]] = Float64(v_["-LH2CB"])
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xupper = fill(v_["UUP"],pb.n)
        pb.xlower = fill(v_["ULW"],pb.n)
        pb.xlower[ix_["U"*string(Int64(v_["0"]))]] = v_["LUA"]
        pb.xupper[ix_["U"*string(Int64(v_["0"]))]] = v_["LUA"]
        pb.xlower[ix_["U"*string(Int64(v_["N+1"]))]] = v_["LUB"]
        pb.xupper[ix_["U"*string(Int64(v_["N+1"]))]] = v_["LUB"]
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(0.0),pb.n)
        pb.x0[ix_["U"*string(Int64(v_["0"]))]] = Float64(v_["LUA"])
        pb.x0[ix_["U"*string(Int64(v_["N+1"]))]] = Float64(v_["LUB"])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eWE1", iet_)
        loaset(elftv,it,1,"X")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"LAC")
        loaset(elftp,it,2,"LAB")
        loaset(elftp,it,3,"LU")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["N"])
            ename = "EA"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eWE1")
            arrset(ielftype, ie, iet_["eWE1"])
            vname = "U"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,v_["ULW"],v_["UUP"],0.0)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="LAC",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["LH2CA"]))
            posep = findfirst(x->x=="LAB",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["-LB"]))
            posep = findfirst(x->x=="LU",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["LUA"]))
            ename = "EB"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eWE1")
            arrset(ielftype, ie, iet_["eWE1"])
            vname = "U"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,v_["ULW"],v_["UUP"],0.0)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="LAC",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["-LH2CB"]))
            posep = findfirst(x->x=="LAB",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["LB"]))
            posep = findfirst(x->x=="LU",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["LUB"]))
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
        for I = Int64(v_["1"]):Int64(v_["N"])
            ig = ig_["G"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["EA"*string(I)])
            loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["EB"*string(I)])
            loaset(pbm.grelw,ig,posel, 1.)
        end
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "SBR2-AN-V-V"
        return pb, pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "eWE1"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        FVAL  = (
              pbm.elpar[iel_][1]*exp(pbm.elpar[iel_][2]*(EV_[1]-pbm.elpar[iel_][3])))
        f_   = FVAL
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = pbm.elpar[iel_][2]*FVAL
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = pbm.elpar[iel_][2]*pbm.elpar[iel_][2]*FVAL
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

    elseif action in  ["fx","fgx","fgHx","cx","cJx","cJHx","cIx","cIJx","cIJHx","cIJxv","fHxv","cJxv","Lxy","Lgxy","LgHxy","LIxy","LIgxy","LIgHxy","LHxyv","LIHxyv"]

        pbm = args[1]
        if pbm.name == name
            pbm.has_globs = [0,0]
            return s2mpj_eval(action,args...)
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

