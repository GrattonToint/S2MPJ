function DRCAV1LQ(action,args...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DRCAV1LQ
#    *********
# 
#    This system of nonlinear equations models the stream function
#    corresponding to an incompressible fluid flow in a driven cavity 
#    (after elimination of the vorticity). The system is solved in the
#    least-squares sense.  The nonlinear system formulation is problem DRCAVTY1.
# 
#    The problem is nonconvex.
#    It differs from the problems DRCAV2LQ and DRCAV3LQ by the value 
#    chosen for the Reynolds number.
# 
#    Source:  
#    P.N. Brown and Y. Saad, 
#    "Hybrid Krylov Methods for Nonlinear Systems of Equations",
#    SIAM J. Sci. Stat. Comput. 11, pp. 450-481, 1990.
#    The boundary conditions have been set according to
#    I.E. Kaporin and O. Axelsson,
#    "On a class of nonlinear equation solvers based on the residual norm
#    reduction over a sequence of affine subspaces",
#    SIAM J, Sci. Comput. 16(1), 1995.
# 
#    SIF input: Ph. Toint, Jan 1995.
# 
#    classification = "OXR2-MY-V-V"
# 
#    Discretization mesh: n = (M+3)**2 - fixed variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "DRCAV1LQ"

    if action == "setup"
        pbm          = PBM(name)
        pb           = PB(name)
        pb.sifpbname = "DRCAV1LQ"
        nargin       = length(args)
        pbm.call     = eval( Meta.parse( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        if nargin<1
            v_["M"] = Int64(10);  #  SIF file default value
        else
            v_["M"] = Int64(args[1]);
        end
#       Alternative values for the SIF file parameters:
# IE M                   31             $-PARAMETER  n =   961
# IE M                   63             $-PARAMETER  n =  3969
# IE M                   100            $-PARAMETER  n = 10000
        v_["M+2"] = 2+v_["M"]
        v_["RM+2"] = Float64(v_["M+2"])
        v_["H"] = 1.0/v_["RM+2"]
        if nargin<2
            v_["RE"] = Float64(500.0);  #  SIF file default value
        else
            v_["RE"] = Float64(args[2]);
        end
        v_["-1"] = -1
        v_["0"] = 0
        v_["1"] = 1
        v_["M+1"] = 1+v_["M"]
        v_["H/2"] = 0.5*v_["H"]
        v_["-H/2"] = -1.0*v_["H/2"]
        v_["RE/4"] = 0.25*v_["RE"]
        v_["-RE/4"] = -1.0*v_["RE/4"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        xscale  = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["-1"]):Int64(v_["M+2"])
            for J = Int64(v_["-1"]):Int64(v_["M+2"])
                iv,ix_,_ = s2x_ii("Y"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"Y"*string(I)*","*string(J))
            end
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["1"]):Int64(v_["M"])
            v_["I-2"] = -2+I
            v_["I-1"] = -1+I
            v_["I+1"] = 1+I
            v_["I+2"] = 2+I
            for J = Int64(v_["1"]):Int64(v_["M"])
                v_["J-2"] = -2+J
                v_["J-1"] = -1+J
                v_["J+1"] = 1+J
                v_["J+2"] = 2+J
                ig,ig_,_ = s2x_ii("E"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"<>")
                iv = ix_["Y"*string(I)*","*string(J)]
                pbm.A[ig,iv] += Float64(20.0)
                iv = ix_["Y"*string(Int64(v_["I-1"]))*","*string(J)]
                pbm.A[ig,iv] += Float64(-8.0)
                iv = ix_["Y"*string(Int64(v_["I+1"]))*","*string(J)]
                pbm.A[ig,iv] += Float64(-8.0)
                iv = ix_["Y"*string(I)*","*string(Int64(v_["J-1"]))]
                pbm.A[ig,iv] += Float64(-8.0)
                iv = ix_["Y"*string(I)*","*string(Int64(v_["J+1"]))]
                pbm.A[ig,iv] += Float64(-8.0)
                iv = ix_["Y"*string(Int64(v_["I-1"]))*","*string(Int64(v_["J+1"]))]
                pbm.A[ig,iv] += Float64(2.0)
                iv = ix_["Y"*string(Int64(v_["I+1"]))*","*string(Int64(v_["J-1"]))]
                pbm.A[ig,iv] += Float64(2.0)
                iv = ix_["Y"*string(Int64(v_["I-1"]))*","*string(Int64(v_["J-1"]))]
                pbm.A[ig,iv] += Float64(2.0)
                iv = ix_["Y"*string(Int64(v_["I+1"]))*","*string(Int64(v_["J+1"]))]
                pbm.A[ig,iv] += Float64(2.0)
                iv = ix_["Y"*string(Int64(v_["I-2"]))*","*string(J)]
                pbm.A[ig,iv] += Float64(1.0)
                iv = ix_["Y"*string(Int64(v_["I+2"]))*","*string(J)]
                pbm.A[ig,iv] += Float64(1.0)
                iv = ix_["Y"*string(I)*","*string(Int64(v_["J-2"]))]
                pbm.A[ig,iv] += Float64(1.0)
                iv = ix_["Y"*string(I)*","*string(Int64(v_["J+2"]))]
                pbm.A[ig,iv] += Float64(1.0)
            end
        end
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        for J = Int64(v_["-1"]):Int64(v_["M+2"])
            pb.xlower[ix_["Y"*string(Int64(v_["-1"]))*","*string(J)]] = 0.0
            pb.xupper[ix_["Y"*string(Int64(v_["-1"]))*","*string(J)]] = 0.0
            pb.xlower[ix_["Y"*string(Int64(v_["0"]))*","*string(J)]] = 0.0
            pb.xupper[ix_["Y"*string(Int64(v_["0"]))*","*string(J)]] = 0.0
        end
        for I = Int64(v_["1"]):Int64(v_["M"])
            pb.xlower[ix_["Y"*string(I)*","*string(Int64(v_["-1"]))]] = 0.0
            pb.xupper[ix_["Y"*string(I)*","*string(Int64(v_["-1"]))]] = 0.0
            pb.xlower[ix_["Y"*string(I)*","*string(Int64(v_["0"]))]] = 0.0
            pb.xupper[ix_["Y"*string(I)*","*string(Int64(v_["0"]))]] = 0.0
        end
        for I = Int64(v_["1"]):Int64(v_["M"])
            pb.xlower[ix_["Y"*string(I)*","*string(Int64(v_["M+1"]))]] = 0.0
            pb.xupper[ix_["Y"*string(I)*","*string(Int64(v_["M+1"]))]] = 0.0
            pb.xlower[ix_["Y"*string(I)*","*string(Int64(v_["M+2"]))]] = 0.0
            pb.xupper[ix_["Y"*string(I)*","*string(Int64(v_["M+2"]))]] = 0.0
        end
        for J = Int64(v_["-1"]):Int64(v_["M+2"])
            pb.xlower[ix_["Y"*string(Int64(v_["M+1"]))*","*string(J)]] = v_["-H/2"]
            pb.xupper[ix_["Y"*string(Int64(v_["M+1"]))*","*string(J)]] = v_["-H/2"]
            pb.xlower[ix_["Y"*string(Int64(v_["M+2"]))*","*string(J)]] = v_["H/2"]
            pb.xupper[ix_["Y"*string(Int64(v_["M+2"]))*","*string(J)]] = v_["H/2"]
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2x_ii( "eIPR", iet_)
        loaset(elftv,it,1,"A1")
        loaset(elftv,it,2,"A2")
        loaset(elftv,it,3,"B1")
        loaset(elftv,it,4,"B2")
        loaset(elftv,it,5,"B3")
        loaset(elftv,it,6,"B4")
        loaset(elftv,it,7,"B5")
        loaset(elftv,it,8,"B6")
        loaset(elftv,it,9,"B7")
        loaset(elftv,it,10,"B8")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["M"])
            v_["I-2"] = -2+I
            v_["I-1"] = -1+I
            v_["I+1"] = 1+I
            v_["I+2"] = 2+I
            for J = Int64(v_["1"]):Int64(v_["M"])
                v_["J-2"] = -2+J
                v_["J-1"] = -1+J
                v_["J+1"] = 1+J
                v_["J+2"] = 2+J
                ename = "X"*string(I)*","*string(J)
                ie,ie_,newelt = s2x_ii(ename,ie_)
                if newelt > 0
                    arrset(pbm.elftype,ie,"eIPR")
                    arrset(ielftype,ie,iet_["eIPR"])
                end
                vname = "Y"*string(I)*","*string(Int64(v_["J+1"]))
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="A1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "Y"*string(I)*","*string(Int64(v_["J-1"]))
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="A2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "Y"*string(Int64(v_["I-2"]))*","*string(J)
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="B1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "Y"*string(Int64(v_["I-1"]))*","*string(Int64(v_["J-1"]))
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="B2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "Y"*string(Int64(v_["I-1"]))*","*string(Int64(v_["J+1"]))
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="B3",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "Y"*string(Int64(v_["I-1"]))*","*string(J)
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="B4",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "Y"*string(Int64(v_["I+1"]))*","*string(J)
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="B5",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "Y"*string(Int64(v_["I+1"]))*","*string(Int64(v_["J-1"]))
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="B6",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "Y"*string(Int64(v_["I+1"]))*","*string(Int64(v_["J+1"]))
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="B7",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "Y"*string(Int64(v_["I+2"]))*","*string(J)
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="B8",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                ename = "Z"*string(I)*","*string(J)
                ie,ie_,newelt = s2x_ii(ename,ie_)
                if newelt > 0
                    arrset(pbm.elftype,ie,"eIPR")
                    arrset(ielftype,ie,iet_["eIPR"])
                end
                vname = "Y"*string(Int64(v_["I+1"]))*","*string(J)
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="A1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "Y"*string(Int64(v_["I-1"]))*","*string(J)
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="A2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "Y"*string(I)*","*string(Int64(v_["J-2"]))
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="B1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "Y"*string(Int64(v_["I-1"]))*","*string(Int64(v_["J-1"]))
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="B2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "Y"*string(Int64(v_["I+1"]))*","*string(Int64(v_["J-1"]))
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="B3",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "Y"*string(I)*","*string(Int64(v_["J-1"]))
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="B4",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "Y"*string(I)*","*string(Int64(v_["J+1"]))
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="B5",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "Y"*string(Int64(v_["I-1"]))*","*string(Int64(v_["J+1"]))
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="B6",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "Y"*string(Int64(v_["I+1"]))*","*string(Int64(v_["J+1"]))
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="B7",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "Y"*string(I)*","*string(Int64(v_["J+2"]))
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="B8",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
        end
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2x_ii("gL2",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for ig in 1:ngrp
            arrset(pbm.grftype,ig,"gL2")
        end
        for I = Int64(v_["1"]):Int64(v_["M"])
            for J = Int64(v_["1"]):Int64(v_["M"])
                ig = ig_["E"*string(I)*","*string(J)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["X"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel,Float64(v_["RE/4"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["Z"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel,Float64(v_["-RE/4"]))
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "OXR2-MY-V-V"
        pb.x0          = zeros(Float64,pb.n)
        return pb, pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "eIPR"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,2,10)
        IV_ =  zeros(Float64,2)
        U_[1,1] = U_[1,1]+1
        U_[1,2] = U_[1,2]-1
        U_[2,3] = U_[2,3]+1
        U_[2,4] = U_[2,4]+1
        U_[2,5] = U_[2,5]+1
        U_[2,6] = U_[2,6]-4
        U_[2,7] = U_[2,7]+4
        U_[2,8] = U_[2,8]-1
        U_[2,9] = U_[2,9]-1
        U_[2,10] = U_[2,10]-1
        IV_[1] = dot(U_[1,:],EV_)
        IV_[2] = dot(U_[2,:],EV_)
        f_   = IV_[1]*IV_[2]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = IV_[2]
            g_[2] = IV_[1]
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = 1.0
                H_[2,1] = H_[1,2]
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

