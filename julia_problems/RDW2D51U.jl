function RDW2D51U(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : RDW2D51U
#    *********
# 
#    A finite-element approximation to the distributed optimal control problem
# 
#       min 1/2||u-v||_L2^2 + beta ||f||_L2^2
# 
#    subject to - nabla^2 u = f
# 
#    where v is given on and within the boundary of a unit [0,1] box in 
#    2 dimensions, and u = v on its boundary. The discretization uses 
#    quadrilateral elememts.
# 
#    The problem is stated as a quadratic program
# 
#    Source:  example 5.1 in 
#     T. Rees, H. S. Dollar and A. J. Wathen
#     "Optimal solvers for PDE-constrained optimization"
#     SIAM J. Sci. Comp. (to appear) 2009
# 
#    SIF input: Nick Gould, May 2009
#               correction by S. Gratton & Ph. Toint, May 2024
# 
#    classification = "C-QLR2-AN-V-V"
# 
#    Number of nodes in each direction (a power of 2)
# 
#       Alternative values for the SIF file parameters:
# IE N                   2             $-PARAMETER
# IE N                   4             $-PARAMETER
# IE N                   8             $-PARAMETER
# IE N                   16            $-PARAMETER
# IE N                   32            $-PARAMETER
# IE N                   64            $-PARAMETER
# IE N                   128           $-PARAMETER
# IE N                   256           $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "RDW2D51U"

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
            v_["N"] = Int64(4);  #  SIF file default value
        else
            v_["N"] = Int64(args[1]);
        end
# IE N                   512           $-PARAMETER
# IE N                   1024          $-PARAMETER
# IE N                   2048          $-PARAMETER
# IE N                   4096          $-PARAMETER
# IE N                   8192          $-PARAMETER
# IE N                   16384         $-PARAMETER
        if nargin<2
            v_["BETA"] = Float64(0.01);  #  SIF file default value
        else
            v_["BETA"] = Float64(args[2]);
        end
        v_["ZERO"] = 0.0
        v_["ONE"] = 1.0
        v_["SIX"] = 6.0
        v_["THIRTYSIX"] = 36.0
        v_["0"] = 0
        v_["1"] = 1
        v_["2"] = 2
        v_["N-"] = -1+v_["N"]
        v_["N-1"] = -1+v_["N"]
        v_["N-2"] = -2+v_["N"]
        v_["N/2"] = trunc(Int,(v_["N"]/v_["2"]))
        v_["N/2+1"] = v_["N/2"]+v_["1"]
        v_["RN"] = Float64(v_["N"])
        v_["H"] = v_["ONE"]/v_["RN"]
        v_["H**2"] = v_["H"]*v_["H"]
        v_["H**2/36"] = v_["H**2"]/v_["THIRTYSIX"]
        v_["-H**2/36"] = -1.0*v_["H**2/36"]
        v_["2BETA"] = 2.0*v_["BETA"]
        v_["2BH**2/36"] = v_["2BETA"]*v_["H**2/36"]
        v_["1/6"] = v_["ONE"]/v_["SIX"]
        for I = Int64(v_["0"]):Int64(v_["N/2"])
            v_["RI"] = Float64(I)
            v_["2RI"] = 2.0*v_["RI"]
            v_["2RIH"] = v_["2RI"]*v_["H"]
            v_["2RIH-1"] = v_["2RIH"]-v_["ONE"]
            v_["2RIH-1S"] = v_["2RIH-1"]*v_["2RIH-1"]
            for J = Int64(v_["0"]):Int64(v_["N/2"])
                v_["RJ"] = Float64(J)
                v_["2RJ"] = 2.0*v_["RJ"]
                v_["2RJH"] = v_["2RJ"]*v_["H"]
                v_["2RJH-1"] = v_["2RJH"]-v_["ONE"]
                v_["2RJH-1S"] = v_["2RJH-1"]*v_["2RJH-1"]
                v_["V"] = v_["2RIH-1S"]*v_["2RJH-1S"]
                v_["V"*string(I)*","*string(J)] = v_["V"]
            end
        end
        for I = Int64(v_["N/2+1"]):Int64(v_["N"])
            for J = Int64(v_["N/2+1"]):Int64(v_["N"])
                v_["V"*string(I)*","*string(J)] = v_["ZERO"]
            end
        end
        for I = Int64(v_["0"]):Int64(v_["N/2"])
            for J = Int64(v_["N/2+1"]):Int64(v_["N"])
                v_["V"*string(I)*","*string(J)] = v_["ZERO"]
            end
        end
        for I = Int64(v_["N/2+1"]):Int64(v_["N"])
            for J = Int64(v_["0"]):Int64(v_["N/2"])
                v_["V"*string(I)*","*string(J)] = v_["ZERO"]
            end
        end
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["0"]):Int64(v_["N"])
            for J = Int64(v_["0"]):Int64(v_["N"])
                iv,ix_,_ = s2mpj_ii("F"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"F"*string(I)*","*string(J))
            end
        end
        for I = Int64(v_["0"]):Int64(v_["N"])
            for J = Int64(v_["0"]):Int64(v_["N"])
                iv,ix_,_ = s2mpj_ii("U"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"U"*string(I)*","*string(J))
            end
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        for I = Int64(v_["1"]):Int64(v_["N-1"])
            for J = Int64(v_["1"]):Int64(v_["N-1"])
                ig,ig_,_ = s2mpj_ii("L"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"L"*string(I)*","*string(J))
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
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        pb.xlower[ix_["U"*string(Int64(v_["0"]))*","*string(Int64(v_["0"]))]]  = (
              v_["V"*string(Int64(v_["0"]))*","*string(Int64(v_["0"]))])
        pb.xupper[ix_["U"*string(Int64(v_["0"]))*","*string(Int64(v_["0"]))]]  = (
              v_["V"*string(Int64(v_["0"]))*","*string(Int64(v_["0"]))])
        pb.xlower[ix_["U"*string(Int64(v_["N"]))*","*string(Int64(v_["0"]))]]  = (
              v_["V"*string(Int64(v_["N"]))*","*string(Int64(v_["0"]))])
        pb.xupper[ix_["U"*string(Int64(v_["N"]))*","*string(Int64(v_["0"]))]]  = (
              v_["V"*string(Int64(v_["N"]))*","*string(Int64(v_["0"]))])
        pb.xlower[ix_["U"*string(Int64(v_["0"]))*","*string(Int64(v_["N"]))]]  = (
              v_["V"*string(Int64(v_["0"]))*","*string(Int64(v_["N"]))])
        pb.xupper[ix_["U"*string(Int64(v_["0"]))*","*string(Int64(v_["N"]))]]  = (
              v_["V"*string(Int64(v_["0"]))*","*string(Int64(v_["N"]))])
        pb.xlower[ix_["U"*string(Int64(v_["N"]))*","*string(Int64(v_["N"]))]]  = (
              v_["V"*string(Int64(v_["N"]))*","*string(Int64(v_["N"]))])
        pb.xupper[ix_["U"*string(Int64(v_["N"]))*","*string(Int64(v_["N"]))]]  = (
              v_["V"*string(Int64(v_["N"]))*","*string(Int64(v_["N"]))])
        for I = Int64(v_["1"]):Int64(v_["N-1"])
            pb.xlower[ix_["U"*string(Int64(v_["0"]))*","*string(I)]]  = (
                  v_["V"*string(Int64(v_["0"]))*","*string(I)])
            pb.xupper[ix_["U"*string(Int64(v_["0"]))*","*string(I)]]  = (
                  v_["V"*string(Int64(v_["0"]))*","*string(I)])
            pb.xlower[ix_["U"*string(Int64(v_["N"]))*","*string(I)]]  = (
                  v_["V"*string(Int64(v_["N"]))*","*string(I)])
            pb.xupper[ix_["U"*string(Int64(v_["N"]))*","*string(I)]]  = (
                  v_["V"*string(Int64(v_["N"]))*","*string(I)])
            pb.xlower[ix_["U"*string(I)*","*string(Int64(v_["0"]))]]  = (
                  v_["V"*string(I)*","*string(Int64(v_["0"]))])
            pb.xupper[ix_["U"*string(I)*","*string(Int64(v_["0"]))]]  = (
                  v_["V"*string(I)*","*string(Int64(v_["0"]))])
            pb.xlower[ix_["U"*string(I)*","*string(Int64(v_["N"]))]]  = (
                  v_["V"*string(I)*","*string(Int64(v_["N"]))])
            pb.xupper[ix_["U"*string(I)*","*string(Int64(v_["N"]))]]  = (
                  v_["V"*string(I)*","*string(Int64(v_["N"]))])
        end
        pb.xlower[ix_["F"*string(Int64(v_["0"]))*","*string(Int64(v_["0"]))]] = 0.0
        pb.xupper[ix_["F"*string(Int64(v_["0"]))*","*string(Int64(v_["0"]))]] = 0.0
        pb.xlower[ix_["F"*string(Int64(v_["N"]))*","*string(Int64(v_["0"]))]] = 0.0
        pb.xupper[ix_["F"*string(Int64(v_["N"]))*","*string(Int64(v_["0"]))]] = 0.0
        pb.xlower[ix_["F"*string(Int64(v_["0"]))*","*string(Int64(v_["N"]))]] = 0.0
        pb.xupper[ix_["F"*string(Int64(v_["0"]))*","*string(Int64(v_["N"]))]] = 0.0
        pb.xlower[ix_["F"*string(Int64(v_["N"]))*","*string(Int64(v_["N"]))]] = 0.0
        pb.xupper[ix_["F"*string(Int64(v_["N"]))*","*string(Int64(v_["N"]))]] = 0.0
        for I = Int64(v_["1"]):Int64(v_["N-1"])
            pb.xlower[ix_["F"*string(Int64(v_["0"]))*","*string(I)]] = 0.0
            pb.xupper[ix_["F"*string(Int64(v_["0"]))*","*string(I)]] = 0.0
            pb.xlower[ix_["F"*string(Int64(v_["N"]))*","*string(I)]] = 0.0
            pb.xupper[ix_["F"*string(Int64(v_["N"]))*","*string(I)]] = 0.0
            pb.xlower[ix_["F"*string(I)*","*string(Int64(v_["0"]))]] = 0.0
            pb.xupper[ix_["F"*string(I)*","*string(Int64(v_["0"]))]] = 0.0
            pb.xlower[ix_["F"*string(I)*","*string(Int64(v_["N"]))]] = 0.0
            pb.xupper[ix_["F"*string(I)*","*string(Int64(v_["N"]))]] = 0.0
        end
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        for I = Int64(v_["0"]):Int64(v_["N"])
            for J = Int64(v_["0"]):Int64(v_["N"])
                v_["I+J"] = I+J
                v_["RI+J"] = Float64(v_["I+J"])
                v_["RI+J/N"] = v_["RI+J"]/v_["RN"]
            end
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eM", iet_)
        loaset(elftv,it,1,"U1")
        loaset(elftv,it,2,"U2")
        loaset(elftv,it,3,"U3")
        loaset(elftv,it,4,"U4")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"V1")
        loaset(elftp,it,2,"V2")
        loaset(elftp,it,3,"V3")
        loaset(elftp,it,4,"V4")
        it,iet_,_ = s2mpj_ii( "eM0", iet_)
        loaset(elftv,it,1,"F1")
        loaset(elftv,it,2,"F2")
        loaset(elftv,it,3,"F3")
        loaset(elftv,it,4,"F4")
        it,iet_,_ = s2mpj_ii( "eA", iet_)
        loaset(elftv,it,1,"U1")
        loaset(elftv,it,2,"U2")
        loaset(elftv,it,3,"U3")
        loaset(elftv,it,4,"U4")
        it,iet_,_ = s2mpj_ii( "eB", iet_)
        loaset(elftv,it,1,"U1")
        loaset(elftv,it,2,"U2")
        loaset(elftv,it,3,"U3")
        loaset(elftv,it,4,"U4")
        it,iet_,_ = s2mpj_ii( "eC", iet_)
        loaset(elftv,it,1,"U1")
        loaset(elftv,it,2,"U2")
        loaset(elftv,it,3,"U3")
        loaset(elftv,it,4,"U4")
        it,iet_,_ = s2mpj_ii( "eD", iet_)
        loaset(elftv,it,1,"U1")
        loaset(elftv,it,2,"U2")
        loaset(elftv,it,3,"U3")
        loaset(elftv,it,4,"U4")
        it,iet_,_ = s2mpj_ii( "eP", iet_)
        loaset(elftv,it,1,"F1")
        loaset(elftv,it,2,"F2")
        loaset(elftv,it,3,"F3")
        loaset(elftv,it,4,"F4")
        it,iet_,_ = s2mpj_ii( "eQ", iet_)
        loaset(elftv,it,1,"F1")
        loaset(elftv,it,2,"F2")
        loaset(elftv,it,3,"F3")
        loaset(elftv,it,4,"F4")
        it,iet_,_ = s2mpj_ii( "eR", iet_)
        loaset(elftv,it,1,"F1")
        loaset(elftv,it,2,"F2")
        loaset(elftv,it,3,"F3")
        loaset(elftv,it,4,"F4")
        it,iet_,_ = s2mpj_ii( "eS", iet_)
        loaset(elftv,it,1,"F1")
        loaset(elftv,it,2,"F2")
        loaset(elftv,it,3,"F3")
        loaset(elftv,it,4,"F4")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["0"]):Int64(v_["N-1"])
            v_["I+"] = I+v_["1"]
            for J = Int64(v_["0"]):Int64(v_["N-1"])
                v_["J+"] = J+v_["1"]
                ename = "E"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eM")
                arrset(ielftype,ie,iet_["eM"])
                vname = "U"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="U1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "U"*string(I)*","*string(Int64(v_["J+"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="U2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "U"*string(Int64(v_["I+"]))*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="U3",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "U"*string(Int64(v_["I+"]))*","*string(Int64(v_["J+"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="U4",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                posep = findfirst(x->x=="V1",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["V"*string(I)*","*string(J)]))
                posep = findfirst(x->x=="V2",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["V"*string(I)*","*string(Int64(v_["J+"]))]))
                posep = findfirst(x->x=="V3",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["V"*string(Int64(v_["I+"]))*","*string(J)]))
                posep = findfirst(x->x=="V4",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["V"*string(Int64(v_["I+"]))*","*string(Int64(v_["J+"]))]))
            end
        end
        for I = Int64(v_["0"]):Int64(v_["N-1"])
            v_["I+"] = I+v_["1"]
            for J = Int64(v_["0"]):Int64(v_["N-1"])
                v_["J+"] = J+v_["1"]
                ename = "F"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eM0")
                arrset(ielftype,ie,iet_["eM0"])
                vname = "F"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="F1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "F"*string(I)*","*string(Int64(v_["J+"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="F2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "F"*string(Int64(v_["I+"]))*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="F3",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "F"*string(Int64(v_["I+"]))*","*string(Int64(v_["J+"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="F4",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
        end
        for I = Int64(v_["0"]):Int64(v_["N-1"])
            v_["I+"] = I+v_["1"]
            for J = Int64(v_["0"]):Int64(v_["N-1"])
                v_["J+"] = J+v_["1"]
                ename = "A"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eA")
                arrset(ielftype,ie,iet_["eA"])
                vname = "U"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="U1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "U"*string(I)*","*string(Int64(v_["J+"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="U2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "U"*string(Int64(v_["I+"]))*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="U3",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "U"*string(Int64(v_["I+"]))*","*string(Int64(v_["J+"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="U4",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                ename = "B"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eB")
                arrset(ielftype,ie,iet_["eB"])
                vname = "U"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="U1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "U"*string(I)*","*string(Int64(v_["J+"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="U2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "U"*string(Int64(v_["I+"]))*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="U3",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "U"*string(Int64(v_["I+"]))*","*string(Int64(v_["J+"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="U4",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                ename = "C"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eC")
                arrset(ielftype,ie,iet_["eC"])
                vname = "U"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="U1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "U"*string(I)*","*string(Int64(v_["J+"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="U2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "U"*string(Int64(v_["I+"]))*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="U3",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "U"*string(Int64(v_["I+"]))*","*string(Int64(v_["J+"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="U4",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                ename = "D"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eD")
                arrset(ielftype,ie,iet_["eD"])
                vname = "U"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="U1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "U"*string(I)*","*string(Int64(v_["J+"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="U2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "U"*string(Int64(v_["I+"]))*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="U3",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "U"*string(Int64(v_["I+"]))*","*string(Int64(v_["J+"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="U4",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
        end
        for I = Int64(v_["0"]):Int64(v_["N-1"])
            v_["I+"] = I+v_["1"]
            for J = Int64(v_["0"]):Int64(v_["N-1"])
                v_["J+"] = J+v_["1"]
                ename = "P"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eP")
                arrset(ielftype,ie,iet_["eP"])
                vname = "F"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="F1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "F"*string(I)*","*string(Int64(v_["J+"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="F2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "F"*string(Int64(v_["I+"]))*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="F3",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "F"*string(Int64(v_["I+"]))*","*string(Int64(v_["J+"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="F4",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                ename = "Q"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eQ")
                arrset(ielftype,ie,iet_["eQ"])
                vname = "F"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="F1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "F"*string(I)*","*string(Int64(v_["J+"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="F2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "F"*string(Int64(v_["I+"]))*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="F3",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "F"*string(Int64(v_["I+"]))*","*string(Int64(v_["J+"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="F4",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                ename = "R"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eR")
                arrset(ielftype,ie,iet_["eR"])
                vname = "F"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="F1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "F"*string(I)*","*string(Int64(v_["J+"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="F2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "F"*string(Int64(v_["I+"]))*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="F3",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "F"*string(Int64(v_["I+"]))*","*string(Int64(v_["J+"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="F4",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                ename = "S"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eS")
                arrset(ielftype,ie,iet_["eS"])
                vname = "F"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="F1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "F"*string(I)*","*string(Int64(v_["J+"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="F2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "F"*string(Int64(v_["I+"]))*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="F3",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "F"*string(Int64(v_["I+"]))*","*string(Int64(v_["J+"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="F4",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["0"]):Int64(v_["N-1"])
            for J = Int64(v_["0"]):Int64(v_["N-1"])
                ig = ig_["OBJ"]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["E"*string(I)*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["H**2/36"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["F"*string(I)*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["2BH**2/36"]))
            end
        end
        for I = Int64(v_["1"]):Int64(v_["N-2"])
            v_["I+"] = I+v_["1"]
            for J = Int64(v_["1"]):Int64(v_["N-2"])
                v_["J+"] = J+v_["1"]
                ig = ig_["L"*string(I)*","*string(J)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A"*string(I)*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["1/6"]))
                ig = ig_["L"*string(I)*","*string(Int64(v_["J+"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B"*string(I)*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["1/6"]))
                ig = ig_["L"*string(Int64(v_["I+"]))*","*string(J)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["C"*string(I)*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["1/6"]))
                ig = ig_["L"*string(Int64(v_["I+"]))*","*string(Int64(v_["J+"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["D"*string(I)*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["1/6"]))
                ig = ig_["L"*string(I)*","*string(J)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["P"*string(I)*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-H**2/36"]))
                ig = ig_["L"*string(I)*","*string(Int64(v_["J+"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["Q"*string(I)*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-H**2/36"]))
                ig = ig_["L"*string(Int64(v_["I+"]))*","*string(J)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["R"*string(I)*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-H**2/36"]))
                ig = ig_["L"*string(Int64(v_["I+"]))*","*string(Int64(v_["J+"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["S"*string(I)*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-H**2/36"]))
            end
        end
        for I = Int64(v_["1"]):Int64(v_["N-2"])
            v_["I+"] = I+v_["1"]
            ig = ig_["L"*string(I)*","*string(Int64(v_["N-"]))]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["A"*string(I)*","*string(Int64(v_["N-"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["1/6"]))
            ig = ig_["L"*string(I)*","*string(Int64(v_["1"]))]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["B"*string(I)*","*string(Int64(v_["0"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["1/6"]))
            ig = ig_["L"*string(Int64(v_["I+"]))*","*string(Int64(v_["N-"]))]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["C"*string(I)*","*string(Int64(v_["N-"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["1/6"]))
            ig = ig_["L"*string(Int64(v_["I+"]))*","*string(Int64(v_["1"]))]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["D"*string(I)*","*string(Int64(v_["0"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["1/6"]))
            ig = ig_["L"*string(I)*","*string(Int64(v_["N-"]))]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["P"*string(I)*","*string(Int64(v_["N-"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["-H**2/36"]))
            ig = ig_["L"*string(I)*","*string(Int64(v_["1"]))]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["Q"*string(I)*","*string(Int64(v_["0"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["-H**2/36"]))
            ig = ig_["L"*string(Int64(v_["I+"]))*","*string(Int64(v_["N-"]))]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["R"*string(I)*","*string(Int64(v_["N-"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["-H**2/36"]))
            ig = ig_["L"*string(Int64(v_["I+"]))*","*string(Int64(v_["1"]))]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["S"*string(I)*","*string(Int64(v_["0"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["-H**2/36"]))
        end
        for J = Int64(v_["1"]):Int64(v_["N-2"])
            v_["J+"] = J+v_["1"]
            ig = ig_["L"*string(Int64(v_["N-"]))*","*string(J)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["A"*string(Int64(v_["N-"]))*","*string(J)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["1/6"]))
            ig = ig_["L"*string(Int64(v_["N-"]))*","*string(Int64(v_["J+"]))]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["B"*string(Int64(v_["N-"]))*","*string(J)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["1/6"]))
            ig = ig_["L"*string(Int64(v_["1"]))*","*string(J)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["C"*string(Int64(v_["0"]))*","*string(J)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["1/6"]))
            ig = ig_["L"*string(Int64(v_["1"]))*","*string(Int64(v_["J+"]))]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["D"*string(Int64(v_["0"]))*","*string(J)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["1/6"]))
            ig = ig_["L"*string(Int64(v_["N-"]))*","*string(J)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["P"*string(Int64(v_["N-"]))*","*string(J)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["-H**2/36"]))
            ig = ig_["L"*string(Int64(v_["N-"]))*","*string(Int64(v_["J+"]))]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["Q"*string(Int64(v_["N-"]))*","*string(J)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["-H**2/36"]))
            ig = ig_["L"*string(Int64(v_["1"]))*","*string(J)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["R"*string(Int64(v_["0"]))*","*string(J)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["-H**2/36"]))
            ig = ig_["L"*string(Int64(v_["1"]))*","*string(Int64(v_["J+"]))]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["S"*string(Int64(v_["0"]))*","*string(J)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["-H**2/36"]))
        end
        ig = ig_["L"*string(Int64(v_["N-"]))*","*string(Int64(v_["N-"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["A"*string(Int64(v_["N-"]))*","*string(Int64(v_["N-"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["1/6"]))
        ig = ig_["L"*string(Int64(v_["N-"]))*","*string(Int64(v_["1"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["B"*string(Int64(v_["N-"]))*","*string(Int64(v_["0"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["1/6"]))
        ig = ig_["L"*string(Int64(v_["1"]))*","*string(Int64(v_["N-"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["C"*string(Int64(v_["0"]))*","*string(Int64(v_["N-"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["1/6"]))
        ig = ig_["L"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["D"*string(Int64(v_["0"]))*","*string(Int64(v_["0"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["1/6"]))
        ig = ig_["L"*string(Int64(v_["N-"]))*","*string(Int64(v_["N-"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["P"*string(Int64(v_["N-"]))*","*string(Int64(v_["N-"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["-H**2/36"]))
        ig = ig_["L"*string(Int64(v_["N-"]))*","*string(Int64(v_["1"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["Q"*string(Int64(v_["N-"]))*","*string(Int64(v_["0"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["-H**2/36"]))
        ig = ig_["L"*string(Int64(v_["1"]))*","*string(Int64(v_["N-"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["R"*string(Int64(v_["0"]))*","*string(Int64(v_["N-"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["-H**2/36"]))
        ig = ig_["L"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["S"*string(Int64(v_["0"]))*","*string(Int64(v_["0"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["-H**2/36"]))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pbm.A = spzeros(Float64,0,0)
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


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "eM"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        UV1 = EV_[1]-pbm.elpar[iel_][1]
        UV2 = EV_[2]-pbm.elpar[iel_][2]
        UV3 = EV_[3]-pbm.elpar[iel_][3]
        UV4 = EV_[4]-pbm.elpar[iel_][4]
        f_   = (2.0*UV1^2+2.0*UV2^2+2.0*UV3^2+2.0*UV4^2+2.0*UV1*UV2+2.0*UV1*UV3+
             UV1*UV4+UV2*UV3+2.0*UV2*UV4+2.0*UV3*UV4)
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 4.0*UV1+2.0*UV2+2.0*UV3+UV4
            g_[2] = 2.0*UV1+4.0*UV2+UV3+2.0*UV4
            g_[3] = 2.0*UV1+UV2+4.0*UV3+2.0*UV4
            g_[4] = UV1+2.0*UV2+2.0*UV3+4.0*UV4
            if nargout>2
                H_ = zeros(Float64,4,4)
                H_[1,1] = 4.0
                H_[1,2] = 2.0
                H_[2,1] = H_[1,2]
                H_[1,3] = 2.0
                H_[3,1] = H_[1,3]
                H_[1,4] = 1.0
                H_[4,1] = H_[1,4]
                H_[2,2] = 4.0
                H_[2,3] = 1.0
                H_[3,2] = H_[2,3]
                H_[2,4] = 2.0
                H_[4,2] = H_[2,4]
                H_[3,3] = 4.0
                H_[3,4] = 2.0
                H_[4,3] = H_[3,4]
                H_[4,4] = 4.0
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eM0"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = (2.0*EV_[1]^2+2.0*EV_[2]^2+2.0*EV_[3]^2+2.0*EV_[4]^2+2.0*EV_[1]*EV_[2]+
             2.0*EV_[1]*EV_[3]+EV_[1]*EV_[4]+EV_[2]*EV_[3]+2.0*EV_[2]*EV_[4]+2.0*EV_[3]*EV_[4])
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 4.0*EV_[1]+2.0*EV_[2]+2.0*EV_[3]+EV_[4]
            g_[2] = 2.0*EV_[1]+4.0*EV_[2]+EV_[3]+2.0*EV_[4]
            g_[3] = 2.0*EV_[1]+EV_[2]+4.0*EV_[3]+2.0*EV_[4]
            g_[4] = EV_[1]+2.0*EV_[2]+2.0*EV_[3]+4.0*EV_[4]
            if nargout>2
                H_ = zeros(Float64,4,4)
                H_[1,1] = 4.0
                H_[1,2] = 2.0
                H_[2,1] = H_[1,2]
                H_[1,3] = 2.0
                H_[3,1] = H_[1,3]
                H_[1,4] = 1.0
                H_[4,1] = H_[1,4]
                H_[2,2] = 4.0
                H_[2,3] = 1.0
                H_[3,2] = H_[2,3]
                H_[2,4] = 2.0
                H_[4,2] = H_[2,4]
                H_[3,3] = 4.0
                H_[3,4] = 2.0
                H_[4,3] = H_[3,4]
                H_[4,4] = 4.0
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eA"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        C1 = 4.0
        C2 = -1.0
        C3 = -1.0
        C4 = -2.0
        f_   = C1*EV_[1]+C2*EV_[2]+C3*EV_[3]+C4*EV_[4]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = C1
            g_[2] = C2
            g_[3] = C3
            g_[4] = C4
            if nargout>2
                H_ = zeros(Float64,4,4)
                H_[1,1] = 0.0
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eB"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        C1 = -1.0
        C2 = 4.0
        C3 = -2.0
        C4 = -1.0
        f_   = C1*EV_[1]+C2*EV_[2]+C3*EV_[3]+C4*EV_[4]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = C1
            g_[2] = C2
            g_[3] = C3
            g_[4] = C4
            if nargout>2
                H_ = zeros(Float64,4,4)
                H_[1,1] = 0.0
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eC"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        C1 = -1.0
        C2 = -2.0
        C3 = 4.0
        C4 = -1.0
        f_   = C1*EV_[1]+C2*EV_[2]+C3*EV_[3]+C4*EV_[4]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = C1
            g_[2] = C2
            g_[3] = C3
            g_[4] = C4
            if nargout>2
                H_ = zeros(Float64,4,4)
                H_[1,1] = 0.0
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eD"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        C1 = -2.0
        C2 = -1.0
        C3 = -1.0
        C4 = 4.0
        f_   = C1*EV_[1]+C2*EV_[2]+C3*EV_[3]+C4*EV_[4]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = C1
            g_[2] = C2
            g_[3] = C3
            g_[4] = C4
            if nargout>2
                H_ = zeros(Float64,4,4)
                H_[1,1] = 0.0
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eP"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        C1 = 4.0
        C2 = 2.0
        C3 = 2.0
        C4 = 1.0
        f_   = C1*EV_[1]+C2*EV_[2]+C3*EV_[3]+C4*EV_[4]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = C1
            g_[2] = C2
            g_[3] = C3
            g_[4] = C4
            if nargout>2
                H_ = zeros(Float64,4,4)
                H_[1,1] = 0.0
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eQ"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        C1 = 2.0
        C2 = 4.0
        C3 = 1.0
        C4 = 2.0
        f_   = C1*EV_[1]+C2*EV_[2]+C3*EV_[3]+C4*EV_[4]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = C1
            g_[2] = C2
            g_[3] = C3
            g_[4] = C4
            if nargout>2
                H_ = zeros(Float64,4,4)
                H_[1,1] = 0.0
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eR"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        C1 = 2.0
        C2 = 1.0
        C3 = 4.0
        C4 = 2.0
        f_   = C1*EV_[1]+C2*EV_[2]+C3*EV_[3]+C4*EV_[4]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = C1
            g_[2] = C2
            g_[3] = C3
            g_[4] = C4
            if nargout>2
                H_ = zeros(Float64,4,4)
                H_[1,1] = 0.0
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eS"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        C1 = 1.0
        C2 = 2.0
        C3 = 2.0
        C4 = 4.0
        f_   = C1*EV_[1]+C2*EV_[2]+C3*EV_[3]+C4*EV_[4]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = C1
            g_[2] = C2
            g_[3] = C3
            g_[4] = C4
            if nargout>2
                H_ = zeros(Float64,4,4)
                H_[1,1] = 0.0
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

