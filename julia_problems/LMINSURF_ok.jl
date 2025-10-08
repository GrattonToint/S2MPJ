function LMINSURF(action,args...)#SG
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem: LMINSURF
#    *********
# 
#    The linear minimum surface problem.
# 
#    The problem comes from the discretization of the minimum surface
#    problem on the unit square: given a set of boundary conditions on
#    the four sides of the square, one must find the surface which
#    meets these boundary conditions and is of minimum area.
# 
#    The unit square is discretized into (p-1)**2 little squares. The
#    heights of the considered surface above the corners of these little
#    squares are the problem variables,  There are p**2 of them.
#    Given these heights, the area above a little square is
#    approximated by the
#      S(i,j) = sqrt( 1 + 0.5(p-1)**2 ( a(i,j)**2 + b(i,j)**2 ) ) / (p-1)**2
#    where
#      a(i,j) = x(i,j) - x(i+1,j+1)
#    and
#      b(i,j) = x(i+1,j) - x(i,j+1)
# 
#    In the Linear Mininum Surface, the boundary conditions are given
#    as the heights of a given plane above the square boundaries.  This
#    plane is specified by its height above the (0,0) point (H00 below),
#    and its slopes along the first and second coordinate
#    directions in the plane (these slopes are denoted SLOPEJ and SLOPEI below).
# 
#    Source:
#    A Griewank and Ph. Toint,
#    "Partitioned variable metric updates for large structured
#    optimization problems",
#    Numerische Mathematik 39:429-448, 1982.
# 
#    SIF input: Ph. Toint, Dec 1989.
# 
#    classification = "OXR2-MY-V-0"
# 
#    P is the number of points in one side of the unit square
# 
# IE P                   4              $-PARAMETER n = 16     original value
# IE P                   7              $-PARAMETER n = 49
# IE P                   8              $-PARAMETER n = 64
# IE P                   11             $-PARAMETER n = 121
# IE P                   31             $-PARAMETER n = 961
# IE P                   32             $-PARAMETER n = 1024
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "LMINSURF"

    if action == "setup"
        pbm          = PBM(name)
        pb.sifpbname = "LMINSURF"
        nargin       = length(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        if nargin<=1
            v_["P"] = 75;  #%  SIF file default value  SG comment avec # et pas %
        else
            v_["P"] = args[1]; # acces avec [ et pas (
        end
        v_["H00"] = 1.0
        v_["SLOPEJ"] = 4.0
        v_["SLOPEI"] = 8.0
        v_["TWOP"] = v_["P"]+v_["P"]
        v_["P-1"] = -1+v_["P"]
        v_["PP-1"] = v_["P"]*v_["P-1"]
        v_["RP-1"] = v_["P-1"]
        v_["INVP-1"] = 1.0/v_["RP-1"]
        v_["RP-1SQ"] = v_["INVP-1"]*v_["INVP-1"]
        v_["SCALE"] = 1.0/v_["RP-1SQ"]
        v_["SQP-1"] = v_["RP-1"]*v_["RP-1"]
        v_["PARAM"] = 0.5*v_["SQP-1"]
        v_["1"] = 1
        v_["2"] = 2
        v_["STON"] = v_["INVP-1"]*v_["SLOPEI"]
        v_["WTOE"] = v_["INVP-1"]*v_["SLOPEJ"]
        v_["H01"] = v_["H00"]+v_["SLOPEJ"]
        v_["H10"] = v_["H00"]+v_["SLOPEI"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        for J = int(v_["1"]):int(v_["P"]) # pas de :
            for I = int(v_["1"]):int(v_["P"])
                iv,ix_,_ = s2x_ii("X"*string(I)*","*string(J),ix_)
                update_list(pb.xnames,iv,"X"*string(I)*","*string(J))
            end
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        for I = int(v_["1"]):int(v_["P-1"])
            for J = int(v_["1"]):int(v_["P-1"])
                ig,ig_,_ = s2x_ii("S"*string(I)*","*string(J),ig_)
                update_list(pbm.grnames,ig,"S"*string(I)*","*string(J))
                update_list(gtype,ig,"<>")
                update_list(pbm.gscale,ig,v_["SCALE"])
            end
        end
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%%  CONSTANTS %%%%%%%%%%%%%%%%%%%
        pbm.gconst = fill(-1.0,ngrp)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        for J = int(v_["1"]):int(v_["P"])
            v_["J-1"] = -1+J
            v_["RJ-1"] = v_["J-1"]
            v_["TH"] = v_["RJ-1"]*v_["WTOE"]
            v_["TL"] = v_["TH"]+v_["H00"]
            v_["TU"] = v_["TH"]+v_["H10"]
            pb.xlower[ix_["X"*string(v_["1"])*","*string(J)]] = v_["TL"]
            pb.xupper[ix_["X"*string(v_["1"])*","*string(J)]] = v_["TL"]
            pb.xlower[ix_["X"*string(v_["P"])*","*string(J)]] = v_["TU"]
            pb.xupper[ix_["X"*string(v_["P"])*","*string(J)]] = v_["TU"]
        end
        for I = int(v_["2"]):int(v_["P-1"])
            v_["I-1"] = -1+I
            v_["RI-1"] = v_["I-1"]
            v_["TV"] = v_["RI-1"]*v_["STON"]
            v_["TR"] = v_["TV"]+v_["H00"]
            v_["TL"] = v_["TV"]+v_["H01"]
            pb.xlower[ix_["X"*string(I)*","*string(v_["P"])]] = v_["TL"]
            pb.xupper[ix_["X"*string(I)*","*string(v_["P"])]] = v_["TL"]
            pb.xlower[ix_["X"*string(I)*","*string(v_["1"])]] = v_["TR"]
            pb.xupper[ix_["X"*string(I)*","*string(v_["1"])]] = v_["TR"]
        end
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(0.0,pb.n,)
        for J = int(v_["1"]):int(v_["P"])
            v_["J-1"] = -1+J
            v_["RJ-1"] = v_["J-1"]
            v_["TH"] = v_["RJ-1"]*v_["WTOE"]
            v_["TL"] = v_["TH"]+v_["H00"]
            v_["TU"] = v_["TH"]+v_["H10"]
            pb.x0[ix_["X"*string(v_["1"])*","*string(J)]] = v_["TL"]
            pb.x0[ix_["X"*string(v_["P"])*","*string(J)]] = v_["TU"]
        end
        for I = int(v_["2"]):int(v_["P-1"])
            v_["I-1"] = -1+I
            v_["RI-1"] = v_["I-1"]
            v_["TV"] = v_["RI-1"]*v_["STON"]
            v_["TR"] = v_["TV"]+v_["H00"]
            v_["TL"] = v_["TV"]+v_["H01"]
            pb.x0[ix_["X"*string(I)*","*string(v_["P"])]] = v_["TL"]
            pb.x0[ix_["X"*string(I)*","*string(v_["1"])]] = v_["TR"]
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = Dict{String,Int}()
        it,iet_,_ = s2x_ii( "ISQ", iet_)
        update_list(elftv[it],0,"V1")  # SG "V1" et pas 'V1'
        update_list(elftv[it],2,"V2")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = Dict{String,Int}()
        for I = int(v_["1"]):int(v_["P-1"])
            v_["I+1"] = 1+I
            for J = int(v_["1"]):int(v_["P-1"])
                v_["J+1"] = 1+J
                ename = "A"*string(I)*","*string(J)
                ie,ie_,_  = s2x_ii(ename,ie_)
                update_list(pbm.enames,ie,ename)
                update_list(pbm.elftype,ie,"ISQ")
                ielftype = iet_["ISQ"]; # SG et pas 'ISQ'
                elv = elftv[ielftype];
                vname = "X"*string(I)*","*string(J)
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,None,None,0.0)
                posev = find(x->x=="V1",elv)
                update_list(pbm.elvar[ie],posev,iv)
                vname = "X"*string(v_["I+1"])*","*string(v_["J+1"])
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,None,None,0.0)
                posev = find(x->x=="V2",elv)
                update_list(pbm.elvar[ie],posev,iv)
                ename = "B"*string(I)*","*string(J)
                ie,ie_,_  = s2x_ii(ename,ie_)
                update_list(pbm.enames,ie,ename)
                update_list(pbm.elftype,ie,"ISQ")
                ielftype = iet_["ISQ"]; # SG et pas 'ISQ'
                elv = elftv[ielftype];
                vname = "X"*string(v_["I+1"])*","*string(J)
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,None,None,0.0)
                posev = find(x->x=="V1",elv)
                update_list(pbm.elvar[ie],posev,iv)
                vname = "X"*string(I)*","*string(v_["J+1"])
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,None,None,0.0)
                posev = find(x->x=="V2",elv)
                update_list(pbm.elvar[ie],posev,iv)
            end
        end
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2x_ii("SQROOT",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = [Int[] for _ in 1:ngrp ]
        for I = int(v_["1"]):int(v_["P-1"])
            for J = int(v_["1"]):int(v_["P-1"])
                ig = ig_["S"*string(I)*","*string(J)]
                update_list(pbm.grftype,ig,"SQROOT")
                posel = length(pbm.grelt[ig])+1
                pbm.grelt  = (
                      update_list(pbm.grelt[ig],posel,ie_["A"*string(I)*","*string(J)]))
                update_list(pbm.grelw[ig],posel,v_["PARAM"])
                posel = length(pbm.grelt[ig])+1
                pbm.grelt  = (
                      update_list(pbm.grelt[ig],posel,ie_["B"*string(I)*","*string(J)]))
                update_list(pbm.grelw[ig],posel,v_["PARAM"])
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%% RETURN VALUES FROM THE SETUP ACTIONS %%%%%%%
        pb.pbclass = "OXR2-MY-V-0"
        return pb, pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "ISQ"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        U_ = spzeros(Float64,1,2)
        IV_ =  zeros(Float64,1)
        U_[1,1] = U_[1,1]+1
        U_[1,2] = U_[1,2]-1
        IV_[1] = U_[1,:]*EV_
        f_   = IV_[1] * IV_[1]
        if nargout>1
            try
                dim = length(IV_)
            catch
                dim = length(EV_)
            end
            g_ = zeros(Float64,dim)
            g_[0] = IV_[1] + IV_[1]
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[0,0] = 2.0
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

    elseif action == "SQROOT"

        GVAR_   = args[1]
        igr_    = args[2]
        nargout = args[3]
        pbm.g_SQRAL = SQRT(GVAR_)
        f_= pbm.g_SQRAL
        if nargout>1
            g_ = 0.5D0 / pbm.g_SQRAL
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_ = -0.25D0 / ( pbm.g_SQRAL * GVAR_ )
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

        if pbm.name == name
            pbm.has_globs = [0,0]
            return s2x_eval(action,args...)
        else
            println("ERROR: please run "*pb.name*" with action = setup")
#            return ntuple(1->undef,args[end]) # chez moi cela ne passe pas.
            return
        end

    else
        println("ERROR: unknown action " * action * " requested from " * pb.name * ".jl")
#       return ntuple(1->undef,args[end]) # SG ci-dessus
        return 
    end

end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

