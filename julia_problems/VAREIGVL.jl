function VAREIGVL(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : VAREIGVL
#    *********
# 
#    The variational eigenvalue by Auchmuty.
#    This problems features a banded matrix of bandwidth 2M+1 = 9.
# 
#    This problem has N least-squares groups, each having a linear part
#    only and N nonlinear elements,
#    plus a least q-th power group having N nonlinear elements.
# 
#    Source: problem 1 in
#    J.J. More',
#    "A collection of nonlinear model problems"
#    Proceedings of the AMS-SIAM Summer seminar on the Computational
#    Solution of Nonlinear Systems of Equations, Colorado, 1988.
#    Argonne National Laboratory MCS-P60-0289, 1989.
# 
#    SIF input: Ph. Toint, Dec 1989.
#               correction by Ph. Shott, January, 1995
#               and Nick Gould, December, 2019, May 2024
# 
#    classification = "C-OUR2-AN-V-0"
# 
#    Number of variables -1 (variable)
# 
#       Alternative values for the SIF file parameters:
# IE N                   19             $-PARAMETER
# IE N                   49             $-PARAMETER     original value
# IE N                   99             $-PARAMETER
# IE N                   499            $-PARAMETER
# IE N                   999            $-PARAMETER
# IE N                   4999           $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "VAREIGVL"

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
            v_["N"] = Int64(19);  #  SIF file default value
        else
            v_["N"] = Int64(args[1]);
        end
# IE M                   4              $-PARAMETER  .le. N
# IE M                   5              $-PARAMETER  .le. N
        if nargin<2
            v_["M"] = Int64(6);  #  SIF file default value
        else
            v_["M"] = Int64(args[2]);
        end
        if nargin<3
            v_["Q"] = Float64(1.5);  #  SIF file default value
        else
            v_["Q"] = Float64(args[3]);
        end
        v_["1"] = 1
        v_["-1.0"] = -1.0
        v_["N+1"] = 1+v_["N"]
        v_["-M"] = -1*v_["M"]
        v_["M+1"] = 1+v_["M"]
        v_["N-M"] = v_["N"]+v_["-M"]
        v_["N-M+1"] = 1+v_["N-M"]
        v_["N2"] = v_["N"]*v_["N"]
        v_["RN2"] = Float64(v_["N2"])
        v_["-1/N2"] = v_["-1.0"]/v_["RN2"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        iv,ix_,_ = s2mpj_ii("MU",ix_)
        arrset(pb.xnames,iv,"MU")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["1"]):Int64(v_["M"])
            v_["RI"] = Float64(I)
            v_["-I"] = -1.0*v_["RI"]
            v_["I+M"] = I+v_["M"]
            for J = Int64(v_["1"]):Int64(v_["I+M"])
                v_["RJ"] = Float64(J)
                v_["IJ"] = v_["RI"]*v_["RJ"]
                v_["SIJ"] = sin(v_["IJ"])
                v_["J-I"] = v_["RJ"]+v_["-I"]
                v_["J-ISQ"] = v_["J-I"]*v_["J-I"]
                v_["ARG"] = v_["J-ISQ"]*v_["-1/N2"]
                v_["EX"] = exp(v_["ARG"])
                v_["AIJ"] = v_["SIJ"]*v_["EX"]
                ig,ig_,_ = s2mpj_ii("G"*string(I),ig_)
                arrset(gtype,ig,"<>")
                iv = ix_["X"*string(J)]
                pbm.A[ig,iv] += Float64(v_["AIJ"])
            end
        end
        for I = Int64(v_["M+1"]):Int64(v_["N-M"])
            v_["RI"] = Float64(I)
            v_["-I"] = -1.0*v_["RI"]
            v_["I-M"] = I+v_["-M"]
            v_["I+M"] = I+v_["M"]
            for J = Int64(v_["I-M"]):Int64(v_["I+M"])
                v_["RJ"] = Float64(J)
                v_["IJ"] = v_["RI"]*v_["RJ"]
                v_["SIJ"] = sin(v_["IJ"])
                v_["J-I"] = v_["RJ"]+v_["-I"]
                v_["J-ISQ"] = v_["J-I"]*v_["J-I"]
                v_["ARG"] = v_["J-ISQ"]*v_["-1/N2"]
                v_["EX"] = exp(v_["ARG"])
                v_["AIJ"] = v_["SIJ"]*v_["EX"]
                ig,ig_,_ = s2mpj_ii("G"*string(I),ig_)
                arrset(gtype,ig,"<>")
                iv = ix_["X"*string(J)]
                pbm.A[ig,iv] += Float64(v_["AIJ"])
            end
        end
        for I = Int64(v_["N-M+1"]):Int64(v_["N"])
            v_["RI"] = Float64(I)
            v_["-I"] = -1.0*v_["RI"]
            v_["I-M"] = I+v_["-M"]
            for J = Int64(v_["I-M"]):Int64(v_["N"])
                v_["RJ"] = Float64(J)
                v_["IJ"] = v_["RI"]*v_["RJ"]
                v_["SIJ"] = sin(v_["IJ"])
                v_["J-I"] = v_["RJ"]+v_["-I"]
                v_["J-ISQ"] = v_["J-I"]*v_["J-I"]
                v_["ARG"] = v_["J-ISQ"]*v_["-1/N2"]
                v_["EX"] = exp(v_["ARG"])
                v_["AIJ"] = v_["SIJ"]*v_["EX"]
                ig,ig_,_ = s2mpj_ii("G"*string(I),ig_)
                arrset(gtype,ig,"<>")
                iv = ix_["X"*string(J)]
                pbm.A[ig,iv] += Float64(v_["AIJ"])
            end
        end
        ig,ig_,_ = s2mpj_ii("G"*string(Int64(v_["N+1"])),ig_)
        arrset(gtype,ig,"<>")
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(1.0),pb.n)
        pb.x0[ix_["MU"]] = Float64(0.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "en2PR", iet_)
        loaset(elftv,it,1,"M")
        loaset(elftv,it,2,"X")
        it,iet_,_ = s2mpj_ii( "eSQ", iet_)
        loaset(elftv,it,1,"X")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["N"])
            ename = "P"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
            vname = "MU"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="M",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "S"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSQ")
            arrset(ielftype,ie,iet_["eSQ"])
            vname = "X"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gLQ",igt_)
        it,igt_,_ = s2mpj_ii("gLQ",igt_)
        grftp = Vector{Vector{String}}()
        loaset(grftp,it,1,"POWER")
        it,igt_,_ = s2mpj_ii("gLQ2",igt_)
        it,igt_,_ = s2mpj_ii("gLQ2",igt_)
        loaset(grftp,it,1,"POWER")
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            ig = ig_["G"*string(I)]
            arrset(pbm.grftype,ig,"gLQ")
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["P"*string(I)])
            loaset(pbm.grelw,ig,posel,Float64(-1.0))
            posgp = findfirst(x->x=="POWER",grftp[igt_[pbm.grftype[ig]]])
            loaset(pbm.grpar,ig,posgp,Float64(2.0))
        end
        ig = ig_["G"*string(Int64(v_["N+1"]))]
        arrset(pbm.grftype,ig,"gLQ2")
        for I = Int64(v_["1"]):Int64(v_["N"])
            ig = ig_["G"*string(Int64(v_["N+1"]))]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["S"*string(I)])
            loaset(pbm.grelw,ig,posel,1.)
        end
        ig = ig_["G"*string(Int64(v_["N+1"]))]
        posgp = findfirst(x->x=="POWER",grftp[igt_[pbm.grftype[ig]]])
        loaset(pbm.grpar,ig,posgp,Float64(v_["Q"]))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN               0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-OUR2-AN-V-0"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        return pb, pbm

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "en2PR"

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

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    elseif action == "gLQ"

        GVAR_   = args[1]
        igr_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        IPOWER = pbm.grpar[igr_][1]
        PM1 = IPOWER-1
        f_= GVAR_^IPOWER/pbm.grpar[igr_][1]
        if nargout>1
            g_ = GVAR_^PM1
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_ = PM1*GVAR_^(IPOWER-2)
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "gLQ2"

        GVAR_   = args[1]
        igr_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_= GVAR_^pbm.grpar[igr_][1]/pbm.grpar[igr_][1]
        if nargout>1
            g_ = GVAR_^(pbm.grpar[igr_][1]-1.0e0)
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_ = (pbm.grpar[igr_][1]-1.0e0)*GVAR_^(pbm.grpar[igr_][1]-2.0e0)
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

