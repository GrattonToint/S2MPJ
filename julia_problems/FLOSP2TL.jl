function FLOSP2TL(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : FLOSP2TL
#    *********
# 
#    A  two-dimensional base  flow  problem in an inclined enclosure.
# 
#    Temperature constant at y = +/- 1 boundary conditions
#    Low Reynold's number
# 
#    The flow is considered in a square of length 2,  centered on the
#    origin and aligned with the x-y axes. The square is divided into
#    4 n ** 2  sub-squares,  each of  length 1 / n.  The differential
#    equation is replaced by  discrete nonlinear equations at each of 
#    the grid points. 
# 
#    The differential equation relates the vorticity, temperature and
#    a stream function.
#    
#    Source: 
#    J. N. Shadid
#    "Experimental and computational study of the stability
#    of Natural convection flow in an inclined enclosure",
#    Ph. D. Thesis, University of Minnesota, 1989,
#    problem SP2 (pp.128-130), 
# 
#    SIF input: Nick Gould, August 1993.
# 
#    classification = "C-NQR2-MY-V-V"
# 
#    Half the number of discretization intervals
#    Number of variables = 3(2M+1)**2 
# 
#       Alternative values for the SIF file parameters:
# IE M                   1              $-PARAMETER n=27
# IE M                   2              $-PARAMETER n=75
# IE M                   5              $-PARAMETER n=363     original value
# IE M                   8              $-PARAMETER n=867
# IE M                   10             $-PARAMETER n=1323
# IE M                   15             $-PARAMETER n=2883
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "FLOSP2TL"

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
            v_["M"] = Int64(1);  #  SIF file default value
        else
            v_["M"] = Int64(args[1]);
        end
        if nargin<2
            v_["RA"] = Float64(1.0e+3);  #  SIF file default value
        else
            v_["RA"] = Float64(args[2]);
        end
        v_["PI/4"] = atan(1.0)
        v_["PI"] = 4.0*v_["PI/4"]
        v_["AX"] = 1.0
        v_["THETA"] = 0.5*v_["PI"]
        v_["A1"] = 0.0
        v_["A2"] = 1.0
        v_["A3"] = 0.0
        v_["B1"] = 0.0
        v_["B2"] = 1.0
        v_["B3"] = 1.0
        v_["F1"] = 1.0
        v_["F2"] = 0.0
        v_["F3"] = 0.0
        v_["G1"] = 1.0
        v_["G2"] = 0.0
        v_["G3"] = 0.0
        v_["M-1"] = -1+v_["M"]
        v_["-M"] = -1*v_["M"]
        v_["-M+1"] = -1*v_["M-1"]
        v_["1/H"] = Float64(v_["M"])
        v_["-1/H"] = -1.0*v_["1/H"]
        v_["2/H"] = 2.0*v_["1/H"]
        v_["-2/H"] = -2.0*v_["1/H"]
        v_["H"] = 1.0/v_["1/H"]
        v_["H2"] = v_["H"]*v_["H"]
        v_["1/H2"] = v_["1/H"]*v_["1/H"]
        v_["-2/H2"] = -2.0*v_["1/H2"]
        v_["1/2H"] = 0.5*v_["1/H"]
        v_["-1/2H"] = -0.5*v_["1/H"]
        v_["AXX"] = v_["AX"]*v_["AX"]
        v_["SINTHETA"] = sin(v_["THETA"])
        v_["COSTHETA"] = cos(v_["THETA"])
        v_["PI1"] = v_["AX"]*v_["RA"]
        v_["PI1"] = v_["PI1"]*v_["COSTHETA"]
        v_["PI1"] = -0.5*v_["PI1"]
        v_["-PI1"] = -1.0*v_["PI1"]
        v_["PI2"] = v_["AXX"]*v_["RA"]
        v_["PI2"] = v_["PI2"]*v_["SINTHETA"]
        v_["PI2"] = 0.5*v_["PI2"]
        v_["-PI2"] = -1.0*v_["PI2"]
        v_["2A1"] = 2.0*v_["A1"]
        v_["2B1"] = 2.0*v_["B1"]
        v_["2F1"] = 2.0*v_["F1"]
        v_["2G1"] = 2.0*v_["G1"]
        v_["2F1/AX"] = v_["2F1"]/v_["AX"]
        v_["2G1/AX"] = v_["2G1"]/v_["AX"]
        v_["AX/2"] = 0.5*v_["AX"]
        v_["AXX/2"] = 0.5*v_["AXX"]
        v_["AXX/4"] = 0.25*v_["AXX"]
        v_["2AX"] = 2.0*v_["AX"]
        v_["2AXX"] = 2.0*v_["AXX"]
        v_["2/AX"] = 2.0/v_["AX"]
        v_["2/AXH"] = v_["2/H"]/v_["AX"]
        v_["-2/AXH"] = -1.0*v_["2/AXH"]
        v_["PI1/2H"] = v_["PI1"]*v_["1/2H"]
        v_["-PI1/2H"] = v_["PI1"]*v_["-1/2H"]
        v_["PI2/2H"] = v_["PI2"]*v_["1/2H"]
        v_["-PI2/2H"] = v_["PI2"]*v_["-1/2H"]
        v_["2A1/H"] = v_["2A1"]*v_["1/H"]
        v_["-2A1/H"] = v_["2A1"]*v_["-1/H"]
        v_["2B1/H"] = v_["2B1"]*v_["1/H"]
        v_["-2B1/H"] = v_["2B1"]*v_["-1/H"]
        v_["2F1/AXH"] = v_["2F1/AX"]*v_["1/H"]
        v_["-2F1/AXH"] = v_["2F1/AX"]*v_["-1/H"]
        v_["2G1/AXH"] = v_["2G1/AX"]*v_["1/H"]
        v_["-2G1/AXH"] = v_["2G1/AX"]*v_["-1/H"]
        v_["AX/H2"] = v_["AX"]*v_["1/H2"]
        v_["-AX/H2"] = -1.0*v_["AX/H2"]
        v_["AX/4H2"] = 0.25*v_["AX/H2"]
        v_["-AX/4H2"] = -0.25*v_["AX/H2"]
        v_["AXX/H2"] = v_["AXX"]*v_["1/H2"]
        v_["-2AXX/H2"] = -2.0*v_["AXX/H2"]
        v_["1"] = 1
        v_["2"] = 2
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for J = Int64(v_["-M"]):Int64(v_["M"])
            for I = Int64(v_["-M"]):Int64(v_["M"])
                iv,ix_,_ = s2mpj_ii("OM"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"OM"*string(I)*","*string(J))
                iv,ix_,_ = s2mpj_ii("PH"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"PH"*string(I)*","*string(J))
                iv,ix_,_ = s2mpj_ii("PS"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"PS"*string(I)*","*string(J))
            end
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for J = Int64(v_["-M+1"]):Int64(v_["M-1"])
            v_["J+"] = 1+J
            v_["J-"] = -1+J
            for I = Int64(v_["-M+1"]):Int64(v_["M-1"])
                v_["I+"] = 1+I
                v_["I-"] = -1+I
                ig,ig_,_ = s2mpj_ii("S"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"S"*string(I)*","*string(J))
                iv = ix_["OM"*string(I)*","*string(J)]
                pbm.A[ig,iv] += Float64(v_["-2/H2"])
                iv = ix_["OM"*string(Int64(v_["I+"]))*","*string(J)]
                pbm.A[ig,iv] += Float64(v_["1/H2"])
                iv = ix_["OM"*string(Int64(v_["I-"]))*","*string(J)]
                pbm.A[ig,iv] += Float64(v_["1/H2"])
                iv = ix_["OM"*string(I)*","*string(J)]
                pbm.A[ig,iv] += Float64(v_["-2AXX/H2"])
                iv = ix_["OM"*string(I)*","*string(Int64(v_["J+"]))]
                pbm.A[ig,iv] += Float64(v_["AXX/H2"])
                iv = ix_["OM"*string(I)*","*string(Int64(v_["J-"]))]
                pbm.A[ig,iv] += Float64(v_["AXX/H2"])
                iv = ix_["PH"*string(Int64(v_["I+"]))*","*string(J)]
                pbm.A[ig,iv] += Float64(v_["-PI1/2H"])
                iv = ix_["PH"*string(Int64(v_["I-"]))*","*string(J)]
                pbm.A[ig,iv] += Float64(v_["PI1/2H"])
                iv = ix_["PH"*string(I)*","*string(Int64(v_["J+"]))]
                pbm.A[ig,iv] += Float64(v_["-PI2/2H"])
                iv = ix_["PH"*string(I)*","*string(Int64(v_["J-"]))]
                pbm.A[ig,iv] += Float64(v_["PI2/2H"])
                ig,ig_,_ = s2mpj_ii("V"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"V"*string(I)*","*string(J))
                iv = ix_["PS"*string(I)*","*string(J)]
                pbm.A[ig,iv] += Float64(v_["-2/H2"])
                iv = ix_["PS"*string(Int64(v_["I+"]))*","*string(J)]
                pbm.A[ig,iv] += Float64(v_["1/H2"])
                iv = ix_["PS"*string(Int64(v_["I-"]))*","*string(J)]
                pbm.A[ig,iv] += Float64(v_["1/H2"])
                iv = ix_["PS"*string(I)*","*string(J)]
                pbm.A[ig,iv] += Float64(v_["-2AXX/H2"])
                iv = ix_["PS"*string(I)*","*string(Int64(v_["J+"]))]
                pbm.A[ig,iv] += Float64(v_["AXX/H2"])
                iv = ix_["PS"*string(I)*","*string(Int64(v_["J-"]))]
                pbm.A[ig,iv] += Float64(v_["AXX/H2"])
                iv = ix_["OM"*string(I)*","*string(J)]
                pbm.A[ig,iv] += Float64(v_["AXX/4"])
                ig,ig_,_ = s2mpj_ii("E"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"E"*string(I)*","*string(J))
                iv = ix_["PH"*string(I)*","*string(J)]
                pbm.A[ig,iv] += Float64(v_["-2/H2"])
                iv = ix_["PH"*string(Int64(v_["I+"]))*","*string(J)]
                pbm.A[ig,iv] += Float64(v_["1/H2"])
                iv = ix_["PH"*string(Int64(v_["I-"]))*","*string(J)]
                pbm.A[ig,iv] += Float64(v_["1/H2"])
                iv = ix_["PH"*string(I)*","*string(J)]
                pbm.A[ig,iv] += Float64(v_["-2AXX/H2"])
                iv = ix_["PH"*string(I)*","*string(Int64(v_["J+"]))]
                pbm.A[ig,iv] += Float64(v_["AXX/H2"])
                iv = ix_["PH"*string(I)*","*string(Int64(v_["J-"]))]
                pbm.A[ig,iv] += Float64(v_["AXX/H2"])
            end
        end
        for K = Int64(v_["-M"]):Int64(v_["M"])
            ig,ig_,_ = s2mpj_ii("T"*string(K)*","*string(Int64(v_["M"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"T"*string(K)*","*string(Int64(v_["M"])))
            iv = ix_["PH"*string(K)*","*string(Int64(v_["M"]))]
            pbm.A[ig,iv] += Float64(v_["2A1/H"])
            ig,ig_,_ = s2mpj_ii("T"*string(K)*","*string(Int64(v_["M"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"T"*string(K)*","*string(Int64(v_["M"])))
            iv = ix_["PH"*string(K)*","*string(Int64(v_["M-1"]))]
            pbm.A[ig,iv] += Float64(v_["-2A1/H"])
            ig,ig_,_ = s2mpj_ii("T"*string(K)*","*string(Int64(v_["M"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"T"*string(K)*","*string(Int64(v_["M"])))
            iv = ix_["PH"*string(K)*","*string(Int64(v_["M"]))]
            pbm.A[ig,iv] += Float64(v_["A2"])
            ig,ig_,_ = s2mpj_ii("T"*string(K)*","*string(Int64(v_["-M"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"T"*string(K)*","*string(Int64(v_["-M"])))
            iv = ix_["PH"*string(K)*","*string(Int64(v_["-M+1"]))]
            pbm.A[ig,iv] += Float64(v_["2B1/H"])
            ig,ig_,_ = s2mpj_ii("T"*string(K)*","*string(Int64(v_["-M"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"T"*string(K)*","*string(Int64(v_["-M"])))
            iv = ix_["PH"*string(K)*","*string(Int64(v_["-M"]))]
            pbm.A[ig,iv] += Float64(v_["-2B1/H"])
            ig,ig_,_ = s2mpj_ii("T"*string(K)*","*string(Int64(v_["-M"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"T"*string(K)*","*string(Int64(v_["-M"])))
            iv = ix_["PH"*string(K)*","*string(Int64(v_["-M"]))]
            pbm.A[ig,iv] += Float64(v_["B2"])
            ig,ig_,_ = s2mpj_ii("T"*string(Int64(v_["M"]))*","*string(K),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"T"*string(Int64(v_["M"]))*","*string(K))
            iv = ix_["PH"*string(Int64(v_["M"]))*","*string(K)]
            pbm.A[ig,iv] += Float64(v_["2F1/AXH"])
            ig,ig_,_ = s2mpj_ii("T"*string(Int64(v_["M"]))*","*string(K),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"T"*string(Int64(v_["M"]))*","*string(K))
            iv = ix_["PH"*string(Int64(v_["M-1"]))*","*string(K)]
            pbm.A[ig,iv] += Float64(v_["-2F1/AXH"])
            ig,ig_,_ = s2mpj_ii("T"*string(Int64(v_["M"]))*","*string(K),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"T"*string(Int64(v_["M"]))*","*string(K))
            iv = ix_["PH"*string(Int64(v_["M"]))*","*string(K)]
            pbm.A[ig,iv] += Float64(v_["F2"])
            ig,ig_,_ = s2mpj_ii("T"*string(Int64(v_["-M"]))*","*string(K),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"T"*string(Int64(v_["-M"]))*","*string(K))
            iv = ix_["PH"*string(Int64(v_["-M+1"]))*","*string(K)]
            pbm.A[ig,iv] += Float64(v_["2G1/AXH"])
            ig,ig_,_ = s2mpj_ii("T"*string(Int64(v_["-M"]))*","*string(K),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"T"*string(Int64(v_["-M"]))*","*string(K))
            iv = ix_["PH"*string(Int64(v_["-M"]))*","*string(K)]
            pbm.A[ig,iv] += Float64(v_["-2G1/AXH"])
            ig,ig_,_ = s2mpj_ii("T"*string(Int64(v_["-M"]))*","*string(K),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"T"*string(Int64(v_["-M"]))*","*string(K))
            iv = ix_["PH"*string(Int64(v_["-M"]))*","*string(K)]
            pbm.A[ig,iv] += Float64(v_["G2"])
            ig,ig_,_ = s2mpj_ii("V"*string(K)*","*string(Int64(v_["M"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"V"*string(K)*","*string(Int64(v_["M"])))
            iv = ix_["PS"*string(K)*","*string(Int64(v_["M"]))]
            pbm.A[ig,iv] += Float64(v_["-2/H"])
            ig,ig_,_ = s2mpj_ii("V"*string(K)*","*string(Int64(v_["M"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"V"*string(K)*","*string(Int64(v_["M"])))
            iv = ix_["PS"*string(K)*","*string(Int64(v_["M-1"]))]
            pbm.A[ig,iv] += Float64(v_["2/H"])
            ig,ig_,_ = s2mpj_ii("V"*string(K)*","*string(Int64(v_["-M"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"V"*string(K)*","*string(Int64(v_["-M"])))
            iv = ix_["PS"*string(K)*","*string(Int64(v_["-M+1"]))]
            pbm.A[ig,iv] += Float64(v_["2/H"])
            ig,ig_,_ = s2mpj_ii("V"*string(K)*","*string(Int64(v_["-M"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"V"*string(K)*","*string(Int64(v_["-M"])))
            iv = ix_["PS"*string(K)*","*string(Int64(v_["-M"]))]
            pbm.A[ig,iv] += Float64(v_["-2/H"])
            ig,ig_,_ = s2mpj_ii("V"*string(Int64(v_["M"]))*","*string(K),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"V"*string(Int64(v_["M"]))*","*string(K))
            iv = ix_["PS"*string(Int64(v_["M"]))*","*string(K)]
            pbm.A[ig,iv] += Float64(v_["-2/AXH"])
            ig,ig_,_ = s2mpj_ii("V"*string(Int64(v_["M"]))*","*string(K),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"V"*string(Int64(v_["M"]))*","*string(K))
            iv = ix_["PS"*string(Int64(v_["M-1"]))*","*string(K)]
            pbm.A[ig,iv] += Float64(v_["2/AXH"])
            ig,ig_,_ = s2mpj_ii("V"*string(Int64(v_["-M"]))*","*string(K),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"V"*string(Int64(v_["-M"]))*","*string(K))
            iv = ix_["PS"*string(Int64(v_["-M+1"]))*","*string(K)]
            pbm.A[ig,iv] += Float64(v_["2/AXH"])
            ig,ig_,_ = s2mpj_ii("V"*string(Int64(v_["-M"]))*","*string(K),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"V"*string(Int64(v_["-M"]))*","*string(K))
            iv = ix_["PS"*string(Int64(v_["-M"]))*","*string(K)]
            pbm.A[ig,iv] += Float64(v_["-2/AXH"])
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
        for K = Int64(v_["-M"]):Int64(v_["M"])
            pbm.gconst[ig_["T"*string(K)*","*string(Int64(v_["M"]))]]  = (
                  Float64(v_["A3"]))
            pbm.gconst[ig_["T"*string(K)*","*string(Int64(v_["-M"]))]]  = (
                  Float64(v_["B3"]))
            pbm.gconst[ig_["T"*string(Int64(v_["M"]))*","*string(K)]]  = (
                  Float64(v_["F3"]))
            pbm.gconst[ig_["T"*string(Int64(v_["-M"]))*","*string(K)]]  = (
                  Float64(v_["G3"]))
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        for K = Int64(v_["-M"]):Int64(v_["M"])
            pb.xlower[ix_["PS"*string(K)*","*string(Int64(v_["-M"]))]] = 1.0
            pb.xupper[ix_["PS"*string(K)*","*string(Int64(v_["-M"]))]] = 1.0
            pb.xlower[ix_["PS"*string(Int64(v_["-M"]))*","*string(K)]] = 1.0
            pb.xupper[ix_["PS"*string(Int64(v_["-M"]))*","*string(K)]] = 1.0
            pb.xlower[ix_["PS"*string(K)*","*string(Int64(v_["M"]))]] = 1.0
            pb.xupper[ix_["PS"*string(K)*","*string(Int64(v_["M"]))]] = 1.0
            pb.xlower[ix_["PS"*string(Int64(v_["M"]))*","*string(K)]] = 1.0
            pb.xupper[ix_["PS"*string(Int64(v_["M"]))*","*string(K)]] = 1.0
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "ePROD", iet_)
        loaset(elftv,it,1,"PSIM")
        loaset(elftv,it,2,"PSIP")
        loaset(elftv,it,3,"PHIM")
        loaset(elftv,it,4,"PHIP")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for J = Int64(v_["-M+1"]):Int64(v_["M-1"])
            v_["J+"] = 1+J
            v_["J-"] = -1+J
            for I = Int64(v_["-M+1"]):Int64(v_["M-1"])
                v_["I+"] = 1+I
                v_["I-"] = -1+I
                ename = "E"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"ePROD")
                arrset(ielftype,ie,iet_["ePROD"])
                vname = "PS"*string(I)*","*string(Int64(v_["J+"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="PSIP",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "PS"*string(I)*","*string(Int64(v_["J-"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="PSIM",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "PH"*string(Int64(v_["I+"]))*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="PHIP",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "PH"*string(Int64(v_["I-"]))*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="PHIM",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                ename = "F"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"ePROD")
                arrset(ielftype,ie,iet_["ePROD"])
                vname = "PS"*string(Int64(v_["I+"]))*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="PSIP",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "PS"*string(Int64(v_["I-"]))*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="PSIM",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "PH"*string(I)*","*string(Int64(v_["J+"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="PHIP",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "PH"*string(I)*","*string(Int64(v_["J-"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="PHIM",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for J = Int64(v_["-M+1"]):Int64(v_["M-1"])
            for I = Int64(v_["-M+1"]):Int64(v_["M-1"])
                ig = ig_["E"*string(I)*","*string(J)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["E"*string(I)*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-AX/4H2"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["F"*string(I)*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["AX/4H2"]))
            end
        end
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-NQR2-MY-V-V"
        pb.x0          = zeros(Float64,pb.n)
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2]
        pb.conderlvl  = pbm.conderlvl;
        return pb, pbm


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "ePROD"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,2,4)
        IV_ =  zeros(Float64,2)
        U_[1,2] = U_[1,2]+1
        U_[1,1] = U_[1,1]-1
        U_[2,4] = U_[2,4]+1
        U_[2,3] = U_[2,3]-1
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

