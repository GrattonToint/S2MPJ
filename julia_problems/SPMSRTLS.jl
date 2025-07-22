function SPMSRTLS(action::String,args::Union{Any}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : SPMSRTLS
#    *********
# 
#    Liu and Nocedal tridiagonal matrix square root problem.
# 
#    Source:  problem 151 (p. 93) in
#    A.R. Buckley,
#    "Test functions for unconstrained minimization",
#    TR 1989CS-3, Mathematics, statistics and computing centre,
#    Dalhousie University, Halifax (CDN), 1989.
# 
#    This is a least-squares variant of problem SPMSQRT.
# 
#    SIF input: Ph. Toint, Dec 1989.
# 
#    classification = "C-CSUR2-AN-V-V"
# 
#    M is the dimension of the matrix
#    The number of variables is 3*M-2
# 
#       Alternative values for the SIF file parameters:
# IE M                   10             $-PARAMETER n = 28     original value
# IE M                   34             $-PARAMETER n = 100
# IE M                   167            $-PARAMETER n = 499
# IE M                   334            $-PARAMETER n = 1000
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 22 VII 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "SPMSRTLS"
    if ( !isdefined(@__MODULE__, :s2mpj_ii) )
        error( "Please include(\"s2mpjlib.jl\") using \"s2mpjlib.jl\" from the S2MPJ distribution before calling SPMSRTLS.")
    end

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
            v_["M"] = Int64(1667);  #  SIF file default value
        else
            v_["M"] = Int64(args[1]);
        end
# IE M                   3334           $-PARAMETER n = 10000
        v_["M-1"] = -1+v_["M"]
        v_["M-2"] = -2+v_["M"]
        v_["M-3"] = -3+v_["M"]
        v_["1"] = 1
        v_["2"] = 2
        v_["3"] = 3
        v_["4"] = 4
        v_["B1,1"] = sin(1.0)
        v_["B1,2"] = sin(4.0)
        v_["K"] = 2
        for I = Int64(v_["2"]):Int64(v_["M-1"])
            v_["I-1"] = -1+I
            v_["I+1"] = 1+I
            v_["K"] = 1+v_["K"]
            v_["RK"] = Float64(v_["K"])
            v_["RKSQR"] = v_["RK"]*v_["RK"]
            v_["B"*string(I)*","*string(Int64(v_["I-1"]))] = sin(v_["RKSQR"])
            v_["K"] = 1+v_["K"]
            v_["RK"] = Float64(v_["K"])
            v_["RKSQR"] = v_["RK"]*v_["RK"]
            v_["B"*string(I)*","*string(I)] = sin(v_["RKSQR"])
            v_["K"] = 1+v_["K"]
            v_["RK"] = Float64(v_["K"])
            v_["RKSQR"] = v_["RK"]*v_["RK"]
            v_["B"*string(I)*","*string(Int64(v_["I+1"]))] = sin(v_["RKSQR"])
        end
        v_["K"] = 1+v_["K"]
        v_["RK"] = Float64(v_["K"])
        v_["RKSQR"] = v_["RK"]*v_["RK"]
        v_["B"*string(Int64(v_["M"]))*","*string(Int64(v_["M-1"]))]  = (
              sin(v_["RKSQR"]))
        v_["K"] = 1+v_["K"]
        v_["RK"] = Float64(v_["K"])
        v_["RKSQR"] = v_["RK"]*v_["RK"]
        v_["B"*string(Int64(v_["M"]))*","*string(Int64(v_["M"]))] = sin(v_["RKSQR"])
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        irA   = Int64[]
        icA   = Int64[]
        valA  = Float64[]
        iv,ix_,_  = (
              s2mpj_ii("X"*string(Int64(v_["1"]))*","*string(Int64(v_["1"])),ix_))
        arrset(pb.xnames,iv,"X"*string(Int64(v_["1"]))*","*string(Int64(v_["1"])))
        iv,ix_,_  = (
              s2mpj_ii("X"*string(Int64(v_["1"]))*","*string(Int64(v_["2"])),ix_))
        arrset(pb.xnames,iv,"X"*string(Int64(v_["1"]))*","*string(Int64(v_["2"])))
        for I = Int64(v_["2"]):Int64(v_["M-1"])
            v_["I-1"] = -1+I
            v_["I+1"] = 1+I
            iv,ix_,_ = s2mpj_ii("X"*string(I)*","*string(Int64(v_["I-1"])),ix_)
            arrset(pb.xnames,iv,"X"*string(I)*","*string(Int64(v_["I-1"])))
            iv,ix_,_ = s2mpj_ii("X"*string(I)*","*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I)*","*string(I))
            iv,ix_,_ = s2mpj_ii("X"*string(I)*","*string(Int64(v_["I+1"])),ix_)
            arrset(pb.xnames,iv,"X"*string(I)*","*string(Int64(v_["I+1"])))
        end
        iv,ix_,_  = (
              s2mpj_ii("X"*string(Int64(v_["M"]))*","*string(Int64(v_["M-1"])),ix_))
        arrset(pb.xnames,iv,"X"*string(Int64(v_["M"]))*","*string(Int64(v_["M-1"])))
        iv,ix_,_  = (
              s2mpj_ii("X"*string(Int64(v_["M"]))*","*string(Int64(v_["M"])),ix_))
        arrset(pb.xnames,iv,"X"*string(Int64(v_["M"]))*","*string(Int64(v_["M"])))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype = String[]
        for J = Int64(v_["1"]):Int64(v_["3"])
            ig,ig_,_ = s2mpj_ii("E"*string(Int64(v_["1"]))*","*string(J),ig_)
            arrset(gtype,ig,"<>")
        end
        for J = Int64(v_["1"]):Int64(v_["4"])
            ig,ig_,_ = s2mpj_ii("E"*string(Int64(v_["2"]))*","*string(J),ig_)
            arrset(gtype,ig,"<>")
        end
        for I = Int64(v_["3"]):Int64(v_["M-2"])
            v_["I-2"] = -2+I
            v_["I+2"] = 2+I
            for J = Int64(v_["I-2"]):Int64(v_["I+2"])
                ig,ig_,_ = s2mpj_ii("E"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"<>")
            end
        end
        for J = Int64(v_["M-3"]):Int64(v_["M"])
            ig,ig_,_ = s2mpj_ii("E"*string(Int64(v_["M-1"]))*","*string(J),ig_)
            arrset(gtype,ig,"<>")
        end
        for J = Int64(v_["M-2"]):Int64(v_["M"])
            ig,ig_,_ = s2mpj_ii("E"*string(Int64(v_["M"]))*","*string(J),ig_)
            arrset(gtype,ig,"<>")
        end
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        v_["ENTRY"] = v_["B1,1"]*v_["B1,1"]
        v_["PROD"] = v_["B1,2"]*v_["B2,1"]
        v_["ENTRY"] = v_["ENTRY"]+v_["PROD"]
        pbm.gconst[ig_["E1,1"]] = Float64(v_["ENTRY"])
        for I = Int64(v_["2"]):Int64(v_["M-1"])
            v_["I-1"] = -1+I
            v_["I+1"] = 1+I
            v_["ENTRY"]  = (
                  v_["B"*string(I)*","*string(I)]*v_["B"*string(I)*","*string(I)])
            v_["PROD"]  = (
                  v_["B"*string(Int64(v_["I-1"]))*","*string(I)]*v_["B"*string(I)*","*string(Int64(v_["I-1"]))])
            v_["ENTRY"] = v_["ENTRY"]+v_["PROD"]
            v_["PROD"]  = (
                  v_["B"*string(Int64(v_["I+1"]))*","*string(I)]*v_["B"*string(I)*","*string(Int64(v_["I+1"]))])
            v_["ENTRY"] = v_["ENTRY"]+v_["PROD"]
            pbm.gconst[ig_["E"*string(I)*","*string(I)]] = Float64(v_["ENTRY"])
        end
        v_["ENTRY"]  = (
              v_["B"*string(Int64(v_["M"]))*","*string(Int64(v_["M"]))]*v_["B"*string(Int64(v_["M"]))*","*string(Int64(v_["M"]))])
        v_["PROD"]  = (
              v_["B"*string(Int64(v_["M-1"]))*","*string(Int64(v_["M"]))]*v_["B"*string(Int64(v_["M"]))*","*string(Int64(v_["M-1"]))])
        v_["ENTRY"] = v_["ENTRY"]+v_["PROD"]
        pbm.gconst[ig_["E"*string(Int64(v_["M"]))*","*string(Int64(v_["M"]))]]  = (
              Float64(v_["ENTRY"]))
        for I = Int64(v_["1"]):Int64(v_["M-1"])
            v_["I+1"] = 1+I
            v_["ENTRY"]  = (
                  v_["B"*string(Int64(v_["I+1"]))*","*string(I)]*v_["B"*string(I)*","*string(I)])
            v_["PROD"]  = (
                  v_["B"*string(Int64(v_["I+1"]))*","*string(Int64(v_["I+1"]))]*v_["B"*string(Int64(v_["I+1"]))*","*string(I)])
            v_["ENTRY"] = v_["ENTRY"]+v_["PROD"]
            pbm.gconst[ig_["E"*string(Int64(v_["I+1"]))*","*string(I)]]  = (
                  Float64(v_["ENTRY"]))
        end
        for I = Int64(v_["2"]):Int64(v_["M"])
            v_["I-1"] = -1+I
            v_["ENTRY"]  = (
                  v_["B"*string(Int64(v_["I-1"]))*","*string(I)]*v_["B"*string(I)*","*string(I)])
            v_["PROD"]  = (
                  v_["B"*string(Int64(v_["I-1"]))*","*string(Int64(v_["I-1"]))]*v_["B"*string(Int64(v_["I-1"]))*","*string(I)])
            v_["ENTRY"] = v_["ENTRY"]+v_["PROD"]
            pbm.gconst[ig_["E"*string(Int64(v_["I-1"]))*","*string(I)]]  = (
                  Float64(v_["ENTRY"]))
        end
        for I = Int64(v_["2"]):Int64(v_["M-1"])
            v_["I-1"] = -1+I
            v_["I+1"] = 1+I
            v_["ENTRY"]  = (
                  v_["B"*string(Int64(v_["I+1"]))*","*string(I)]*v_["B"*string(I)*","*string(Int64(v_["I-1"]))])
            pbm.gconst[ig_["E"*string(Int64(v_["I+1"]))*","*string(Int64(v_["I-1"]))]]  = (
                  Float64(v_["ENTRY"]))
        end
        for I = Int64(v_["2"]):Int64(v_["M-1"])
            v_["I-1"] = -1+I
            v_["I+1"] = 1+I
            v_["ENTRY"]  = (
                  v_["B"*string(Int64(v_["I-1"]))*","*string(I)]*v_["B"*string(I)*","*string(Int64(v_["I+1"]))])
            pbm.gconst[ig_["E"*string(Int64(v_["I-1"]))*","*string(Int64(v_["I+1"]))]]  = (
                  Float64(v_["ENTRY"]))
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        v_["PROD"] = 0.2*v_["B1,1"]
        pb.x0[ix_["X1,1"]] = Float64(v_["PROD"])
        v_["PROD"] = 0.2*v_["B1,2"]
        pb.x0[ix_["X1,2"]] = Float64(v_["PROD"])
        for I = Int64(v_["2"]):Int64(v_["M-1"])
            v_["I-1"] = -1+I
            v_["I+1"] = 1+I
            v_["PROD"] = 0.2*v_["B"*string(I)*","*string(Int64(v_["I-1"]))]
            pb.x0[ix_["X"*string(I)*","*string(Int64(v_["I-1"]))]] = Float64(v_["PROD"])
            v_["PROD"] = 0.2*v_["B"*string(I)*","*string(I)]
            pb.x0[ix_["X"*string(I)*","*string(I)]] = Float64(v_["PROD"])
            v_["PROD"] = 0.2*v_["B"*string(I)*","*string(Int64(v_["I+1"]))]
            pb.x0[ix_["X"*string(I)*","*string(Int64(v_["I+1"]))]] = Float64(v_["PROD"])
        end
        v_["PROD"] = 0.2*v_["B"*string(Int64(v_["M"]))*","*string(Int64(v_["M-1"]))]
        pb.x0[ix_["X"*string(Int64(v_["M"]))*","*string(Int64(v_["M-1"]))]]  = (
              Float64(v_["PROD"]))
        v_["PROD"] = 0.2*v_["B"*string(Int64(v_["M"]))*","*string(Int64(v_["M"]))]
        pb.x0[ix_["X"*string(Int64(v_["M"]))*","*string(Int64(v_["M"]))]]  = (
              Float64(v_["PROD"]))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "ePROD2", iet_)
        loaset(elftv,it,1,"VI")
        loaset(elftv,it,2,"VJ")
        it,iet_,_ = s2mpj_ii( "eSQ", iet_)
        loaset(elftv,it,1,"V")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "D"*string(Int64(v_["1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        ename = "D"*string(Int64(v_["1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "G"*string(Int64(v_["1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD2")
        arrset(ielftype,ie,iet_["ePROD2"])
        ename = "G"*string(Int64(v_["1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["1"]))*","*string(Int64(v_["2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VI",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "G"*string(Int64(v_["1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["2"]))*","*string(Int64(v_["1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VJ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "H"*string(Int64(v_["1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD2")
        arrset(ielftype,ie,iet_["ePROD2"])
        ename = "H"*string(Int64(v_["1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VI",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "H"*string(Int64(v_["1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["1"]))*","*string(Int64(v_["2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VJ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "R"*string(Int64(v_["1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD2")
        arrset(ielftype,ie,iet_["ePROD2"])
        ename = "R"*string(Int64(v_["1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["1"]))*","*string(Int64(v_["2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VI",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "R"*string(Int64(v_["1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["2"]))*","*string(Int64(v_["2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VJ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "S"*string(Int64(v_["1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD2")
        arrset(ielftype,ie,iet_["ePROD2"])
        ename = "S"*string(Int64(v_["1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["1"]))*","*string(Int64(v_["2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VI",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "S"*string(Int64(v_["1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["2"]))*","*string(Int64(v_["3"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VJ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "B"*string(Int64(v_["2"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD2")
        arrset(ielftype,ie,iet_["ePROD2"])
        ename = "B"*string(Int64(v_["2"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["2"]))*","*string(Int64(v_["1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VI",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "B"*string(Int64(v_["2"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VJ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "C"*string(Int64(v_["2"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD2")
        arrset(ielftype,ie,iet_["ePROD2"])
        ename = "C"*string(Int64(v_["2"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["2"]))*","*string(Int64(v_["2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VI",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "C"*string(Int64(v_["2"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["2"]))*","*string(Int64(v_["1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VJ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F"*string(Int64(v_["2"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD2")
        arrset(ielftype,ie,iet_["ePROD2"])
        ename = "F"*string(Int64(v_["2"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["2"]))*","*string(Int64(v_["1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VI",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F"*string(Int64(v_["2"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["1"]))*","*string(Int64(v_["2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VJ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "D"*string(Int64(v_["2"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        ename = "D"*string(Int64(v_["2"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["2"]))*","*string(Int64(v_["2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "G"*string(Int64(v_["2"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD2")
        arrset(ielftype,ie,iet_["ePROD2"])
        ename = "G"*string(Int64(v_["2"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["2"]))*","*string(Int64(v_["3"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VI",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "G"*string(Int64(v_["2"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["3"]))*","*string(Int64(v_["2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VJ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "H"*string(Int64(v_["2"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD2")
        arrset(ielftype,ie,iet_["ePROD2"])
        ename = "H"*string(Int64(v_["2"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["2"]))*","*string(Int64(v_["2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VI",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "H"*string(Int64(v_["2"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["2"]))*","*string(Int64(v_["3"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VJ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "R"*string(Int64(v_["2"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD2")
        arrset(ielftype,ie,iet_["ePROD2"])
        ename = "R"*string(Int64(v_["2"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["2"]))*","*string(Int64(v_["3"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VI",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "R"*string(Int64(v_["2"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["3"]))*","*string(Int64(v_["3"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VJ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "S"*string(Int64(v_["2"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD2")
        arrset(ielftype,ie,iet_["ePROD2"])
        ename = "S"*string(Int64(v_["2"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["2"]))*","*string(Int64(v_["3"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VI",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "S"*string(Int64(v_["2"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["3"]))*","*string(Int64(v_["4"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VJ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        for I = Int64(v_["3"]):Int64(v_["M-2"])
            v_["I-2"] = -2+I
            v_["I-1"] = -1+I
            v_["I+1"] = 1+I
            v_["I+2"] = 2+I
            ename = "A"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"ePROD2")
            arrset(ielftype,ie,iet_["ePROD2"])
            vname = "X"*string(I)*","*string(Int64(v_["I-1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
            posev = findfirst(x->x=="VI",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["I-1"]))*","*string(Int64(v_["I-2"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
            posev = findfirst(x->x=="VJ",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"ePROD2")
            arrset(ielftype,ie,iet_["ePROD2"])
            vname = "X"*string(I)*","*string(Int64(v_["I-1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
            posev = findfirst(x->x=="VI",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["I-1"]))*","*string(Int64(v_["I-1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
            posev = findfirst(x->x=="VJ",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "C"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"ePROD2")
            arrset(ielftype,ie,iet_["ePROD2"])
            vname = "X"*string(I)*","*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
            posev = findfirst(x->x=="VI",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(I)*","*string(Int64(v_["I-1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
            posev = findfirst(x->x=="VJ",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "F"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"ePROD2")
            arrset(ielftype,ie,iet_["ePROD2"])
            vname = "X"*string(I)*","*string(Int64(v_["I-1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
            posev = findfirst(x->x=="VI",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["I-1"]))*","*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
            posev = findfirst(x->x=="VJ",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "D"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSQ")
            arrset(ielftype,ie,iet_["eSQ"])
            vname = "X"*string(I)*","*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "G"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"ePROD2")
            arrset(ielftype,ie,iet_["ePROD2"])
            vname = "X"*string(I)*","*string(Int64(v_["I+1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
            posev = findfirst(x->x=="VI",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["I+1"]))*","*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
            posev = findfirst(x->x=="VJ",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "H"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"ePROD2")
            arrset(ielftype,ie,iet_["ePROD2"])
            vname = "X"*string(I)*","*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
            posev = findfirst(x->x=="VI",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(I)*","*string(Int64(v_["I+1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
            posev = findfirst(x->x=="VJ",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "R"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"ePROD2")
            arrset(ielftype,ie,iet_["ePROD2"])
            vname = "X"*string(I)*","*string(Int64(v_["I+1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
            posev = findfirst(x->x=="VI",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["I+1"]))*","*string(Int64(v_["I+1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
            posev = findfirst(x->x=="VJ",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "S"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"ePROD2")
            arrset(ielftype,ie,iet_["ePROD2"])
            vname = "X"*string(I)*","*string(Int64(v_["I+1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
            posev = findfirst(x->x=="VI",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["I+1"]))*","*string(Int64(v_["I+2"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
            posev = findfirst(x->x=="VJ",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        ename = "A"*string(Int64(v_["M-1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD2")
        arrset(ielftype,ie,iet_["ePROD2"])
        ename = "A"*string(Int64(v_["M-1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["M-1"]))*","*string(Int64(v_["M-2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VI",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "A"*string(Int64(v_["M-1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["M-2"]))*","*string(Int64(v_["M-3"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VJ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "B"*string(Int64(v_["M-1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD2")
        arrset(ielftype,ie,iet_["ePROD2"])
        ename = "B"*string(Int64(v_["M-1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["M-1"]))*","*string(Int64(v_["M-2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VI",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "B"*string(Int64(v_["M-1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["M-2"]))*","*string(Int64(v_["M-2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VJ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "C"*string(Int64(v_["M-1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD2")
        arrset(ielftype,ie,iet_["ePROD2"])
        ename = "C"*string(Int64(v_["M-1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["M-1"]))*","*string(Int64(v_["M-1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VI",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "C"*string(Int64(v_["M-1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["M-1"]))*","*string(Int64(v_["M-2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VJ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F"*string(Int64(v_["M-1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD2")
        arrset(ielftype,ie,iet_["ePROD2"])
        ename = "F"*string(Int64(v_["M-1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["M-1"]))*","*string(Int64(v_["M-2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VI",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F"*string(Int64(v_["M-1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["M-2"]))*","*string(Int64(v_["M-1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VJ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "D"*string(Int64(v_["M-1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        ename = "D"*string(Int64(v_["M-1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["M-1"]))*","*string(Int64(v_["M-1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "G"*string(Int64(v_["M-1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD2")
        arrset(ielftype,ie,iet_["ePROD2"])
        ename = "G"*string(Int64(v_["M-1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["M-1"]))*","*string(Int64(v_["M"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VI",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "G"*string(Int64(v_["M-1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["M"]))*","*string(Int64(v_["M-1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VJ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "H"*string(Int64(v_["M-1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD2")
        arrset(ielftype,ie,iet_["ePROD2"])
        ename = "H"*string(Int64(v_["M-1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["M-1"]))*","*string(Int64(v_["M-1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VI",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "H"*string(Int64(v_["M-1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["M-1"]))*","*string(Int64(v_["M"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VJ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "R"*string(Int64(v_["M-1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD2")
        arrset(ielftype,ie,iet_["ePROD2"])
        ename = "R"*string(Int64(v_["M-1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["M-1"]))*","*string(Int64(v_["M"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VI",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "R"*string(Int64(v_["M-1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["M"]))*","*string(Int64(v_["M"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VJ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "A"*string(Int64(v_["M"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD2")
        arrset(ielftype,ie,iet_["ePROD2"])
        ename = "A"*string(Int64(v_["M"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["M"]))*","*string(Int64(v_["M-1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VI",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "A"*string(Int64(v_["M"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["M-1"]))*","*string(Int64(v_["M-2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VJ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "B"*string(Int64(v_["M"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD2")
        arrset(ielftype,ie,iet_["ePROD2"])
        ename = "B"*string(Int64(v_["M"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["M"]))*","*string(Int64(v_["M-1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VI",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "B"*string(Int64(v_["M"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["M-1"]))*","*string(Int64(v_["M-1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VJ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "C"*string(Int64(v_["M"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD2")
        arrset(ielftype,ie,iet_["ePROD2"])
        ename = "C"*string(Int64(v_["M"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["M"]))*","*string(Int64(v_["M"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VI",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "C"*string(Int64(v_["M"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["M"]))*","*string(Int64(v_["M-1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VJ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F"*string(Int64(v_["M"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD2")
        arrset(ielftype,ie,iet_["ePROD2"])
        ename = "F"*string(Int64(v_["M"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["M"]))*","*string(Int64(v_["M-1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VI",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F"*string(Int64(v_["M"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["M-1"]))*","*string(Int64(v_["M"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="VJ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "D"*string(Int64(v_["M"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        ename = "D"*string(Int64(v_["M"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["M"]))*","*string(Int64(v_["M"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
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
        ig = ig_["E"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["D"*string(Int64(v_["1"]))])
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["G"*string(Int64(v_["1"]))])
        loaset(pbm.grelw,ig,posel, 1.)
        ig = ig_["E"*string(Int64(v_["1"]))*","*string(Int64(v_["2"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["H"*string(Int64(v_["1"]))])
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["R"*string(Int64(v_["1"]))])
        loaset(pbm.grelw,ig,posel, 1.)
        ig = ig_["E"*string(Int64(v_["1"]))*","*string(Int64(v_["3"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["S"*string(Int64(v_["1"]))])
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["E"*string(Int64(v_["2"]))*","*string(Int64(v_["1"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["B"*string(Int64(v_["2"]))])
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["C"*string(Int64(v_["2"]))])
        loaset(pbm.grelw,ig,posel, 1.)
        ig = ig_["E"*string(Int64(v_["2"]))*","*string(Int64(v_["2"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F"*string(Int64(v_["2"]))])
        loaset(pbm.grelw,ig,posel,1.)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["D"*string(Int64(v_["2"]))])
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["G"*string(Int64(v_["2"]))])
        loaset(pbm.grelw,ig,posel, 1.)
        ig = ig_["E"*string(Int64(v_["2"]))*","*string(Int64(v_["3"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["H"*string(Int64(v_["2"]))])
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["R"*string(Int64(v_["2"]))])
        loaset(pbm.grelw,ig,posel, 1.)
        ig = ig_["E"*string(Int64(v_["2"]))*","*string(Int64(v_["4"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["S"*string(Int64(v_["2"]))])
        loaset(pbm.grelw,ig,posel,1.)
        for I = Int64(v_["3"]):Int64(v_["M-2"])
            v_["I-2"] = -2+I
            v_["I-1"] = -1+I
            v_["I+1"] = 1+I
            v_["I+2"] = 2+I
            ig = ig_["E"*string(I)*","*string(Int64(v_["I-2"]))]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["A"*string(I)])
            loaset(pbm.grelw,ig,posel,1.)
            ig = ig_["E"*string(I)*","*string(Int64(v_["I-1"]))]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["B"*string(I)])
            loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["C"*string(I)])
            loaset(pbm.grelw,ig,posel, 1.)
            ig = ig_["E"*string(I)*","*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["F"*string(I)])
            loaset(pbm.grelw,ig,posel,1.)
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["D"*string(I)])
            loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["G"*string(I)])
            loaset(pbm.grelw,ig,posel, 1.)
            ig = ig_["E"*string(I)*","*string(Int64(v_["I+1"]))]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["H"*string(I)])
            loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["R"*string(I)])
            loaset(pbm.grelw,ig,posel, 1.)
            ig = ig_["E"*string(I)*","*string(Int64(v_["I+2"]))]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["S"*string(I)])
            loaset(pbm.grelw,ig,posel,1.)
        end
        ig = ig_["E"*string(Int64(v_["M-1"]))*","*string(Int64(v_["M-3"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["A"*string(Int64(v_["M-1"]))])
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["E"*string(Int64(v_["M-1"]))*","*string(Int64(v_["M-2"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["B"*string(Int64(v_["M-1"]))])
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["C"*string(Int64(v_["M-1"]))])
        loaset(pbm.grelw,ig,posel, 1.)
        ig = ig_["E"*string(Int64(v_["M-1"]))*","*string(Int64(v_["M-1"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F"*string(Int64(v_["M-1"]))])
        loaset(pbm.grelw,ig,posel,1.)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["D"*string(Int64(v_["M-1"]))])
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["G"*string(Int64(v_["M-1"]))])
        loaset(pbm.grelw,ig,posel, 1.)
        ig = ig_["E"*string(Int64(v_["M-1"]))*","*string(Int64(v_["M"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["H"*string(Int64(v_["M-1"]))])
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["R"*string(Int64(v_["M-1"]))])
        loaset(pbm.grelw,ig,posel, 1.)
        ig = ig_["E"*string(Int64(v_["M"]))*","*string(Int64(v_["M-2"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["A"*string(Int64(v_["M"]))])
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["E"*string(Int64(v_["M"]))*","*string(Int64(v_["M-1"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["B"*string(Int64(v_["M"]))])
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["C"*string(Int64(v_["M"]))])
        loaset(pbm.grelw,ig,posel, 1.)
        ig = ig_["E"*string(Int64(v_["M"]))*","*string(Int64(v_["M"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F"*string(Int64(v_["M"]))])
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["D"*string(Int64(v_["M"]))])
        loaset(pbm.grelw,ig,posel, 1.)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-CSUR2-AN-V-V"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        return pb, pbm


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "ePROD2"

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

