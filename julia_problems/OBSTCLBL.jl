function OBSTCLBL(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : OBSTCLBL
#    *********
# 
#    A quadratic obstacle problem by Dembo and Tulowitzki
# 
#    The problem comes from the obstacle problem on a rectangle.
#    The rectangle is discretized into (px-1)(py-1) little rectangles. The
#    heights of the considered surface above the corners of these little
#    rectangles are the problem variables,  There are px*py of them.
# 
#    Source:
#    R. Dembo and U. Tulowitzki,
#    "On the minimization of quadratic functions subject to box
#    constraints",
#    WP 71, Yale University (new Haven, USA), 1983.
# 
#    See also More 1989 (Problem B, Starting point L)
# 
#    SIF input: Ph. Toint, Dec 1989.
# 
#    classification = "C-QBR2-AY-V-0"
# 
#    PX is the number of points along the X side of the rectangle
#    PY is the number of points along the Y side of the rectangle
# 
#       Alternative values for the SIF file parameters:
# IE PX                  4              $-PARAMETER n = 16
# IE PY                  4              $-PARAMETER
# 
# IE PX                  10             $-PARAMETER n = 100     original value
# IE PY                  10             $-PARAMETER             original value
# 
# IE PX                  23             $-PARAMETER n = 529
# IE PY                  23             $-PARAMETER
# 
# IE PX                  32             $-PARAMETER n = 1024
# IE PY                  32             $-PARAMETER
# 
# IE PX                  75             $-PARAMETER n = 5625
# IE PY                  75             $-PARAMETER
# 
# IE PX                  100            $-PARAMETER n = 10000
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "OBSTCLBL"

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
            v_["PX"] = Int64(5);  #  SIF file default value
        else
            v_["PX"] = Int64(args[1]);
        end
# IE PY                  100            $-PARAMETER
        if nargin<2
            v_["PY"] = Int64(20);  #  SIF file default value
        else
            v_["PY"] = Int64(args[2]);
        end
# IE PX                  125            $-PARAMETER n = 15625
# IE PY                  125            $-PARAMETER
        if nargin<3
            v_["C"] = Float64(1.0);  #  SIF file default value
        else
            v_["C"] = Float64(args[3]);
        end
        v_["PX-1"] = -1+v_["PX"]
        v_["RPX-1"] = Float64(v_["PX-1"])
        v_["HX"] = 1.0/v_["RPX-1"]
        v_["1/HX"] = 1.0/v_["HX"]
        v_["PY-1"] = -1+v_["PY"]
        v_["RPY-1"] = Float64(v_["PY-1"])
        v_["HY"] = 1.0/v_["RPY-1"]
        v_["1/HY"] = 1.0/v_["HY"]
        v_["HXHY"] = v_["HX"]*v_["HY"]
        v_["HX/HY"] = v_["HX"]*v_["1/HY"]
        v_["HY/HX"] = v_["HY"]*v_["1/HX"]
        v_["HY/4HX"] = 0.25*v_["HY/HX"]
        v_["HX/4HY"] = 0.25*v_["HX/HY"]
        v_["C0"] = v_["HXHY"]*v_["C"]
        v_["LC"] = -1.0*v_["C0"]
        v_["1"] = 1
        v_["2"] = 2
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for J = Int64(v_["1"]):Int64(v_["PX"])
            for I = Int64(v_["1"]):Int64(v_["PY"])
                iv,ix_,_ = s2mpj_ii("X"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"X"*string(I)*","*string(J))
            end
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["2"]):Int64(v_["PY-1"])
            for J = Int64(v_["2"]):Int64(v_["PX-1"])
                ig,ig_,_ = s2mpj_ii("G"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"<>")
                iv = ix_["X"*string(I)*","*string(J)]
                pbm.A[ig,iv] += Float64(v_["LC"])
            end
        end
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        for J = Int64(v_["1"]):Int64(v_["PX"])
            pb.xlower[ix_["X"*string(Int64(v_["1"]))*","*string(J)]] = 0.0
            pb.xupper[ix_["X"*string(Int64(v_["1"]))*","*string(J)]] = 0.0
            pb.xlower[ix_["X"*string(Int64(v_["PY"]))*","*string(J)]] = 0.0
            pb.xupper[ix_["X"*string(Int64(v_["PY"]))*","*string(J)]] = 0.0
        end
        for I = Int64(v_["2"]):Int64(v_["PY-1"])
            pb.xlower[ix_["X"*string(I)*","*string(Int64(v_["PX"]))]] = 0.0
            pb.xupper[ix_["X"*string(I)*","*string(Int64(v_["PX"]))]] = 0.0
            pb.xlower[ix_["X"*string(I)*","*string(Int64(v_["1"]))]] = 0.0
            pb.xupper[ix_["X"*string(I)*","*string(Int64(v_["1"]))]] = 0.0
        end
        for I = Int64(v_["2"]):Int64(v_["PY-1"])
            v_["I-1"] = -1+I
            v_["RI-1"] = Float64(v_["I-1"])
            v_["XI1"] = v_["RI-1"]*v_["HY"]
            v_["3XI1"] = 9.2*v_["XI1"]
            v_["SXI1"] = sin(v_["3XI1"])
            for J = Int64(v_["2"]):Int64(v_["PX-1"])
                v_["J-1"] = -1+J
                v_["RJ-1"] = Float64(v_["J-1"])
                v_["XI2"] = v_["RJ-1"]*v_["HX"]
                v_["3XI2"] = 9.3*v_["XI2"]
                v_["SXI2"] = sin(v_["3XI2"])
                v_["L1"] = v_["SXI1"]*v_["SXI2"]
                v_["L2"] = v_["L1"]*v_["L1"]
                v_["LOW"] = v_["L2"]*v_["L1"]
                pb.xlower[ix_["X"*string(I)*","*string(J)]] = v_["LOW"]
            end
        end
        for I = Int64(v_["2"]):Int64(v_["PY-1"])
            v_["I-1"] = -1+I
            v_["RI-1"] = Float64(v_["I-1"])
            v_["XI1"] = v_["RI-1"]*v_["HY"]
            v_["3XI1"] = 9.2*v_["XI1"]
            v_["SXI1"] = sin(v_["3XI1"])
            for J = Int64(v_["2"]):Int64(v_["PX-1"])
                v_["J-1"] = -1+J
                v_["RJ-1"] = Float64(v_["J-1"])
                v_["XI2"] = v_["RJ-1"]*v_["HX"]
                v_["3XI2"] = 9.3*v_["XI2"]
                v_["SXI2"] = sin(v_["3XI2"])
                v_["L1"] = v_["SXI1"]*v_["SXI2"]
                v_["L2"] = v_["L1"]*v_["L1"]
                v_["UPP"] = 0.02+v_["L2"]
                pb.xupper[ix_["X"*string(I)*","*string(J)]] = v_["UPP"]
            end
        end
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        for J = Int64(v_["1"]):Int64(v_["PX"])
            pb.x0[ix_["X"*string(Int64(v_["1"]))*","*string(J)]] = Float64(0.0)
            pb.x0[ix_["X"*string(Int64(v_["PY"]))*","*string(J)]] = Float64(0.0)
        end
        for I = Int64(v_["2"]):Int64(v_["PY-1"])
            pb.x0[ix_["X"*string(I)*","*string(Int64(v_["PX"]))]] = Float64(0.0)
            pb.x0[ix_["X"*string(I)*","*string(Int64(v_["1"]))]] = Float64(0.0)
        end
        for I = Int64(v_["2"]):Int64(v_["PY-1"])
            v_["I-1"] = -1+I
            v_["RI-1"] = Float64(v_["I-1"])
            v_["XI1"] = v_["RI-1"]*v_["HY"]
            v_["3XI1"] = 9.2*v_["XI1"]
            v_["SXI1"] = sin(v_["3XI1"])
            for J = Int64(v_["2"]):Int64(v_["PX-1"])
                v_["J-1"] = -1+J
                v_["RJ-1"] = Float64(v_["J-1"])
                v_["XI2"] = v_["RJ-1"]*v_["HX"]
                v_["3XI2"] = 9.3*v_["XI2"]
                v_["SXI2"] = sin(v_["3XI2"])
                v_["L1"] = v_["SXI1"]*v_["SXI2"]
                v_["L2"] = v_["L1"]*v_["L1"]
                v_["LOW"] = v_["L2"]*v_["L1"]
                pb.x0[ix_["X"*string(I)*","*string(J)]] = Float64(v_["LOW"])
            end
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eISQ", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["2"]):Int64(v_["PY-1"])
            v_["I-1"] = -1+I
            v_["I+1"] = 1+I
            for J = Int64(v_["2"]):Int64(v_["PX-1"])
                v_["J-1"] = -1+J
                v_["J+1"] = 1+J
                ename = "A"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eISQ")
                arrset(ielftype,ie,iet_["eISQ"])
                vname = "X"*string(Int64(v_["I+1"]))*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                ename = "B"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eISQ")
                arrset(ielftype,ie,iet_["eISQ"])
                vname = "X"*string(I)*","*string(Int64(v_["J+1"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                ename = "C"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eISQ")
                arrset(ielftype,ie,iet_["eISQ"])
                vname = "X"*string(Int64(v_["I-1"]))*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                ename = "D"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eISQ")
                arrset(ielftype,ie,iet_["eISQ"])
                vname = "X"*string(I)*","*string(Int64(v_["J-1"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["2"]):Int64(v_["PY-1"])
            for J = Int64(v_["2"]):Int64(v_["PX-1"])
                ig = ig_["G"*string(I)*","*string(J)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel,Float64(v_["HY/4HX"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel,Float64(v_["HX/4HY"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["C"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel,Float64(v_["HY/4HX"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["D"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel,Float64(v_["HX/4HY"]))
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN(4)            -0.0081108
# LO SOLTN(10)           2.87503823
# LO SOLTN(23)           6.51932527
# LO SOLTN(32)           6.88708670
# LO SOLTN(75)           ???
# LO SOLTN(100)          ???
# LO SOLTN(125)          ???
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-QBR2-AY-V-0"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
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

