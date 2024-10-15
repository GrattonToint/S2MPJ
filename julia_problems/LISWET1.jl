function LISWET1(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : LISWET1
#    *********
# 
#    A k-convex approximation problem posed as a 
#    convex quadratic problem, with variable dimensions.
# 
#    Formulation:
#    -----------
# 
#                 n+k             2
#    minimize 1/2 sum ( x  - c  )
#                 i=1    i    i
# 
#    subject to
# 
#                  k              k-i
#                 sum ( k ) ( -1 )    x     > 0
#                 i=0 ( i )            j+i  = 
# 
#    where c  = g( t ) + small perturbation, t  = (i-1)/(n+k-1)
#           i       i                         i 
# 
#    Case 1: g(t) = sqrt(t)
# 
#    NB. Perturbations are not random as Li and Swetits's 
#        random number generator is undefined.
# 
#    Source:
#    W. Li and J. Swetits,
#    "A Newton method for convex regression, data smoothing and
#    quadratic programming with bounded constraints",
#    SIAM J. Optimization 3 (3) pp 466-488, 1993.
# 
#    SIF input: Nick Gould, August 1994.
# 
#    classification = "C-QLR2-AN-V-V"
# 
#       Alternative values for the SIF file parameters:
# IE N                   100            $-PARAMETER 103 variables original value 
# IE K                   3              $-PARAMETER original value
# 
# IE N                   100            $-PARAMETER 104 variables    
# IE K                   4              $-PARAMETER
# 
# IE N                   100            $-PARAMETER 105 variables    
# IE K                   5              $-PARAMETER
# 
# IE N                   100            $-PARAMETER 106 variables    
# IE K                   6              $-PARAMETER
# 
# IE N                   400            $-PARAMETER 402 variables    
# IE K                   2              $-PARAMETER
# 
# IE N                   400            $-PARAMETER 403 variables    
# IE K                   3              $-PARAMETER
# 
# IE N                   2000           $-PARAMETER 2001 variables    
# IE K                   1              $-PARAMETER
# 
# IE N                   2000           $-PARAMETER 2002 variables    
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "LISWET1"

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
            v_["N"] = Int64(50);  #  SIF file default value
        else
            v_["N"] = Int64(args[1]);
        end
# IE K                   2              $-PARAMETER
        if nargin<2
            v_["K"] = Int64(3);  #  SIF file default value
        else
            v_["K"] = Int64(args[2]);
        end
# IE N                   10000          $-PARAMETER 10001 variables    
# IE K                   1              $-PARAMETER
# IE N                   10000          $-PARAMETER 10002 variables    
# IE K                   2              $-PARAMETER
        v_["0"] = 0
        v_["1"] = 1
        v_["ONE"] = 1.0
        v_["HALF"] = 0.5
        v_["N+K"] = v_["N"]+v_["K"]
        v_["N+K-1"] = -1+v_["N+K"]
        v_["RN+K-1"] = Float64(v_["N+K-1"])
        v_["CONST"] = 0.0
        v_["B"*string(Int64(v_["0"]))] = v_["ONE"]
        for I = Int64(v_["1"]):Int64(v_["K"])
            v_["I-1"] = -1+I
            v_["RI"] = Float64(I)
            v_["B"*string(I)] = v_["B"*string(Int64(v_["I-1"]))]*v_["RI"]
        end
        v_["C"*string(Int64(v_["0"]))] = v_["ONE"]
        v_["PLUSMINUS"] = v_["ONE"]
        for I = Int64(v_["1"]):Int64(v_["K"])
            v_["K-I"] = v_["K"]-I
            v_["PLUSMINUS"] = -1.0*v_["PLUSMINUS"]
            v_["C"*string(I)] = v_["B"*string(Int64(v_["K"]))]/v_["B"*string(I)]
            v_["C"*string(I)] = v_["C"*string(I)]/v_["B"*string(Int64(v_["K-I"]))]
            v_["C"*string(I)] = v_["C"*string(I)]*v_["PLUSMINUS"]
        end
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N+K"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["1"]):Int64(v_["N+K"])
            v_["I-1"] = -1+I
            v_["RI"] = Float64(I)
            v_["RI-1"] = Float64(v_["I-1"])
            v_["TI"] = v_["RI-1"]/v_["RN+K-1"]
            v_["GT"] = sqrt(v_["TI"])
            v_["RANDOM"] = sin(v_["RI"])
            v_["RANDOM"] = 0.1*v_["RANDOM"]
            v_["CI"] = v_["GT"]+v_["RANDOM"]
            v_["-CI"] = -1.0*v_["CI"]
            v_["-CI*CI"] = v_["-CI"]*v_["CI"]
            v_["CONST"] = v_["CONST"]+v_["-CI*CI"]
            ig,ig_,_ = s2mpj_ii("OBJ",ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += Float64(v_["-CI"])
        end
        for J = Int64(v_["1"]):Int64(v_["N"])
            v_["J+K"] = J+v_["K"]
            for I = Int64(v_["0"]):Int64(v_["K"])
                v_["J+K-I"] = v_["J+K"]-I
                ig,ig_,_ = s2mpj_ii("CON"*string(J),ig_)
                arrset(gtype,ig,">=")
                arrset(pb.cnames,ig,"CON"*string(J))
                iv = ix_["X"*string(Int64(v_["J+K-I"]))]
                pbm.A[ig,iv] += Float64(v_["C"*string(I)])
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
        v_["CONST"] = v_["HALF"]*v_["CONST"]
        pbm.gconst[ig_["OBJ"]] = Float64(v_["CONST"])
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eSQ", iet_)
        loaset(elftv,it,1,"X")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["N+K"])
            ename = "XSQ"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSQ")
            arrset(ielftype,ie,iet_["eSQ"])
            vname = "X"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N+K"])
            ig = ig_["OBJ"]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["XSQ"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               
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
        pb.pbclass = "C-QLR2-AN-V-V"
        pb.x0          = zeros(Float64,pb.n)
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
        f_   = 5.0e-1*EV_[1]*EV_[1]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[1]
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 1.0e+0
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

