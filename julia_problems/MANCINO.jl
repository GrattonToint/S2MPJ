function MANCINO(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : MANCINO
#    *********
# 
#    Mancino's function with variable dimension.
# 
#    Source:
#    E. Spedicato,
#    "Computational experience with quasi-Newton algorithms for
#    minimization problems of moderate size",
#    Report N-175, CISE, Milano, 1975.
# 
#    See also Buckley #51 (p. 72), Schittkowski #391 (for N = 30)
# 
#    SIF input: Ph. Toint, Dec 1989.
#               correction by Ph. Shott, January, 1995.
#               correction by S. Gratton & Ph. Toint, May 2024
# 
#    classification = "C-SUR2-AN-V-0"
# 
#    The definitions
#      s_{i,j} = \sin \log v_{i,j}   and s_{i,j} = \cos \log v_{i,j}
#    have been used.  It seems that the additional exponent ALPHA
#    in Buckley is a typo.
# 
#    Number of variables
# 
#       Alternative values for the SIF file parameters:
# IE N                   10             $-PARAMETER
# IE N                   20             $-PARAMETER
# IE N                   30             $-PARAMETER Schittkowski #391
# IE N                   50             $-PARAMETER
# IE N                   100            $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "MANCINO"

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
            v_["N"] = Int64(10);  #  SIF file default value
        else
            v_["N"] = Int64(args[1]);
        end
        if nargin<2
            v_["ALPHA"] = Int64(5);  #  SIF file default value
        else
            v_["ALPHA"] = Int64(args[2]);
        end
        if nargin<3
            v_["BETA"] = Float64(14.0);  #  SIF file default value
        else
            v_["BETA"] = Float64(args[3]);
        end
        if nargin<4
            v_["GAMMA"] = Int64(3);  #  SIF file default value
        else
            v_["GAMMA"] = Int64(args[4]);
        end
        v_["RALPHA"] = Float64(v_["ALPHA"])
        v_["RN"] = Float64(v_["N"])
        v_["N-1"] = -1+v_["N"]
        v_["RN-1"] = Float64(v_["N-1"])
        v_["N-1SQ"] = v_["RN-1"]*v_["RN-1"]
        v_["BETAN"] = v_["BETA"]*v_["RN"]
        v_["BETAN2"] = v_["BETAN"]*v_["BETAN"]
        v_["AL+1"] = 1.0+v_["RALPHA"]
        v_["A1SQ"] = v_["AL+1"]*v_["AL+1"]
        v_["F0"] = v_["A1SQ"]*v_["N-1SQ"]
        v_["F1"] = -1.0*v_["F0"]
        v_["F2"] = v_["BETAN2"]+v_["F1"]
        v_["F3"] = 1.0/v_["F2"]
        v_["F4"] = v_["BETAN"]*v_["F3"]
        v_["A"] = -1.0*v_["F4"]
        v_["-N/2"] = -0.5*v_["RN"]
        v_["1"] = 1
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            ig,ig_,_ = s2mpj_ii("G"*string(I),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += Float64(v_["BETAN"])
        end
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        for I = Int64(v_["1"]):Int64(v_["N"])
            v_["RI"] = Float64(I)
            v_["I-N/2"] = v_["RI"]+v_["-N/2"]
            v_["CI"] = 1.0
            for J = Int64(v_["1"]):Int64(v_["GAMMA"])
                v_["CI"] = v_["CI"]*v_["I-N/2"]
            end
            pbm.gconst[ig_["G"*string(I)]] = Float64(v_["CI"])
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        for I = Int64(v_["1"]):Int64(v_["N"])
            v_["I-1"] = -1+I
            v_["RI"] = Float64(I)
            v_["H"] = 0.0
            for J = Int64(v_["1"]):Int64(v_["I-1"])
                v_["RJ"] = Float64(J)
                v_["1/J"] = 1.0/v_["RJ"]
                v_["I/J"] = v_["RI"]*v_["1/J"]
                v_["SQI/J"] = sqrt(v_["I/J"])
                v_["LIJ"] = log(v_["SQI/J"])
                v_["SIJ"] = sin(v_["LIJ"])
                v_["CIJ"] = cos(v_["LIJ"])
                v_["SA"] = 1.0
                v_["CA"] = 1.0
                for K = Int64(v_["1"]):Int64(v_["ALPHA"])
                    v_["SA"] = v_["SA"]*v_["SIJ"]
                    v_["CA"] = v_["CA"]*v_["CIJ"]
                end
                v_["SCA"] = v_["SA"]+v_["CA"]
                v_["HIJ"] = v_["SQI/J"]*v_["SCA"]
                v_["H"] = v_["H"]+v_["HIJ"]
            end
            v_["I+1"] = 1+I
            for J = Int64(v_["I+1"]):Int64(v_["N"])
                v_["RJ"] = Float64(J)
                v_["1/J"] = 1.0/v_["RJ"]
                v_["I/J"] = v_["RI"]*v_["1/J"]
                v_["SQI/J"] = sqrt(v_["I/J"])
                v_["LIJ"] = log(v_["SQI/J"])
                v_["SIJ"] = sin(v_["LIJ"])
                v_["CIJ"] = cos(v_["LIJ"])
                v_["SA"] = 1.0
                v_["CA"] = 1.0
                for K = Int64(v_["1"]):Int64(v_["ALPHA"])
                    v_["SA"] = v_["SA"]*v_["SIJ"]
                    v_["CA"] = v_["CA"]*v_["CIJ"]
                end
                v_["SCA"] = v_["SA"]+v_["CA"]
                v_["HIJ"] = v_["SQI/J"]*v_["SCA"]
                v_["H"] = v_["H"]+v_["HIJ"]
            end
            v_["I-N/2"] = v_["RI"]+v_["-N/2"]
            v_["CI"] = 1.0
            for J = Int64(v_["1"]):Int64(v_["GAMMA"])
                v_["CI"] = v_["CI"]*v_["I-N/2"]
            end
            v_["TMP"] = v_["H"]+v_["CI"]
            v_["XI0"] = v_["TMP"]*v_["A"]
            pb.x0[ix_["X"*string(I)]] = Float64(v_["XI0"])
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eMANC", iet_)
        loaset(elftv,it,1,"X")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"II")
        loaset(elftp,it,2,"JJ")
        loaset(elftp,it,3,"AL")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["N"])
            v_["RI"] = Float64(I)
            v_["I-1"] = -1+I
            for J = Int64(v_["1"]):Int64(v_["I-1"])
                v_["RJ"] = Float64(J)
                ename = "E"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eMANC")
                arrset(ielftype,ie,iet_["eMANC"])
                vname = "X"*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="X",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                posep = findfirst(x->x=="II",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["RI"]))
                posep = findfirst(x->x=="JJ",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["RJ"]))
                posep = findfirst(x->x=="AL",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["RALPHA"]))
            end
            v_["I+1"] = 1+I
            for J = Int64(v_["I+1"]):Int64(v_["N"])
                v_["RJ"] = Float64(J)
                ename = "E"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eMANC")
                arrset(ielftype,ie,iet_["eMANC"])
                vname = "X"*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="X",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                posep = findfirst(x->x=="II",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["RI"]))
                posep = findfirst(x->x=="JJ",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["RJ"]))
                posep = findfirst(x->x=="AL",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["RALPHA"]))
            end
        end
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gL2",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            v_["I-1"] = -1+I
            ig = ig_["G"*string(I)]
            arrset(pbm.grftype,ig,"gL2")
            for J = Int64(v_["1"]):Int64(v_["I-1"])
                ig = ig_["G"*string(I)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["E"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel,1.)
            end
            v_["I+1"] = 1+I
            for J = Int64(v_["I+1"]):Int64(v_["N"])
                ig = ig_["G"*string(I)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["E"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel,1.)
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
        pb.objlower = 0.0
#    Solution
# LO SOLTN               0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-SUR2-AN-V-0"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        return pb, pbm

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "eMANC"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        IAL = pbm.elpar[iel_][3]
        IA1 = IAL-1
        A2 = pbm.elpar[iel_][3]-2.0
        IA2 = IAL-2
        IA3 = IAL-3
        INVIJ = EV_[1]*EV_[1]+pbm.elpar[iel_][1]/pbm.elpar[iel_][2]
        VIJ = sqrt(INVIJ)
        V2 = VIJ*VIJ
        DVIJ = EV_[1]/VIJ
        LIJ = log(VIJ)
        SIJ = sin(LIJ)
        CIJ = cos(LIJ)
        DSDX = CIJ*DVIJ/VIJ
        DCDX = -SIJ*DVIJ/VIJ
        SUMAL = SIJ^IAL+CIJ^IAL
        DSUMAL = pbm.elpar[iel_][3]*(DSDX*SIJ^IA1+DCDX*CIJ^IA1)
        SCIJ = SIJ*CIJ
        DSCIJ = SIJ*DCDX+DSDX*CIJ
        SAL = SIJ^IA2-CIJ^IA2
        DSAL = A2*(DSDX*SIJ^IA3-DCDX*CIJ^IA3)
        B = SUMAL+pbm.elpar[iel_][3]*SCIJ*SAL
        DBDX = DSUMAL+pbm.elpar[iel_][3]*(DSCIJ*SAL+SCIJ*DSAL)
        f_   = VIJ*SUMAL
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[1]*B/VIJ
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = (B+EV_[1]*DBDX)/VIJ-B*EV_[1]*DVIJ/V2
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

