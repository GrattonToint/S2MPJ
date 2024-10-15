function DRUGDISE(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DRUGDISE
#    *********
# 
#    This is a variant of the drug displacement problem DRUGDIS where the
#    state equations have been Expanded in term of more intermediate
#    functions, each one of them being less nonlinear.
# 
#    The problem is based on the kinetic model of Aarons and Rowland which
#    simulates the interaction of the two drugs (warfarin and phenylnutazone)
#    in a patient bloodstream.  The state variable are the concentrations of
#    unbound warfarin (w) and phenylbutazone (p).  The problem is to control
#    the rate of injection (u) of the pain-killing phenylbutazone so that both
#    drugs reach a specified steady-state in minimum time and the concentration
#    of warfarin does not rise above a toxicity level.
# 
#    The problem is discretized using the trapeziodal rule.  It is non-convex.
# 
#    Source:
#    H. Maurer and M. Wiegand,
#    "Numerical solution of a drug displacement problem with bounded state
#    variables",
#    Optimal Control Applications and Methods 13, pp. 43-55, 1992.
# 
#    SIF input: Ph. Toint, Nov 1993.
# 
#    classification = "C-LOR2-MY-V-V"
# 
#    Discretization: specify the number of interior points + 1
# 
#       Alternative values for the SIF file parameters:
# IE NI                  10             $-PARAMETER n=63, m=50 
# IE NI                  100            $-PARAMETER n=603, m=500   original value
# IE NI                  100            $-PARAMETER n=6003, m=5000 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "DRUGDISE"

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
            v_["NI"] = Int64(10);  #  SIF file default value
        else
            v_["NI"] = Int64(args[1]);
        end
        if nargin<2
            v_["TOXIC"] = Float64(0.026);  #  SIF file default value
        else
            v_["TOXIC"] = Float64(args[2]);
        end
        if nargin<3
            v_["WSS"] = Float64(0.02);  #  SIF file default value
        else
            v_["WSS"] = Float64(args[3]);
        end
        if nargin<4
            v_["UMAX"] = Float64(8.0);  #  SIF file default value
        else
            v_["UMAX"] = Float64(args[4]);
        end
        if nargin<5
            v_["PSTART"] = Float64(0.0);  #  SIF file default value
        else
            v_["PSTART"] = Float64(args[5]);
        end
        if nargin<6
            v_["PFINAL"] = Float64(2.0);  #  SIF file default value
        else
            v_["PFINAL"] = Float64(args[6]);
        end
        if nargin<7
            v_["Z"] = Float64(46.4);  #  SIF file default value
        else
            v_["Z"] = Float64(args[7]);
        end
        v_["AVP"] = v_["PSTART"]+v_["PFINAL"]
        v_["AVP"] = 0.5*v_["AVP"]
        v_["-Z"] = -1.0*v_["Z"]
        v_["-ZZ"] = v_["Z"]*v_["-Z"]
        v_["NI-1"] = -1+v_["NI"]
        v_["RNI"] = Float64(v_["NI"])
        v_["-1/NI"] = -1.0/v_["RNI"]
        v_["-Z/NI"] = v_["Z"]*v_["-1/NI"]
        v_["0"] = 0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        iv,ix_,_ = s2mpj_ii("TF",ix_)
        arrset(pb.xnames,iv,"TF")
        arrset(pb.xscale,iv,200.0)
        for I = Int64(v_["0"]):Int64(v_["NI"])
            iv,ix_,_ = s2mpj_ii("W"*string(I),ix_)
            arrset(pb.xnames,iv,"W"*string(I))
            arrset(pb.xscale,iv,0.02)
        end
        for I = Int64(v_["0"]):Int64(v_["NI"])
            iv,ix_,_ = s2mpj_ii("P"*string(I),ix_)
            arrset(pb.xnames,iv,"P"*string(I))
        end
        for I = Int64(v_["0"]):Int64(v_["NI-1"])
            iv,ix_,_ = s2mpj_ii("U"*string(I),ix_)
            arrset(pb.xnames,iv,"U"*string(I))
        end
        for I = Int64(v_["0"]):Int64(v_["NI-1"])
            iv,ix_,_ = s2mpj_ii("A"*string(I),ix_)
            arrset(pb.xnames,iv,"A"*string(I))
            arrset(pb.xscale,iv,200.0)
        end
        for I = Int64(v_["0"]):Int64(v_["NI-1"])
            iv,ix_,_ = s2mpj_ii("B"*string(I),ix_)
            arrset(pb.xnames,iv,"B"*string(I))
            arrset(pb.xscale,iv,200.0)
        end
        for I = Int64(v_["0"]):Int64(v_["NI-1"])
            iv,ix_,_ = s2mpj_ii("C"*string(I),ix_)
            arrset(pb.xnames,iv,"C"*string(I))
            arrset(pb.xscale,iv,0.0000001)
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("TFINAL",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["TF"]
        pbm.A[ig,iv] += Float64(1.0)
        arrset(pbm.gscale,ig,Float64(100.0))
        for I = Int64(v_["0"]):Int64(v_["NI-1"])
            v_["I+1"] = 1+I
            ig,ig_,_ = s2mpj_ii("EW"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"EW"*string(I))
            iv = ix_["W"*string(Int64(v_["I+1"]))]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["W"*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            arrset(pbm.gscale,ig,Float64(0.02))
            ig,ig_,_ = s2mpj_ii("EP"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"EP"*string(I))
            iv = ix_["P"*string(Int64(v_["I+1"]))]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["P"*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("EA"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"EA"*string(I))
            iv = ix_["A"*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["P"*string(I)]
            pbm.A[ig,iv] += Float64(v_["-Z"])
            arrset(pbm.gscale,ig,Float64(200.0))
            ig,ig_,_ = s2mpj_ii("EB"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"EB"*string(I))
            iv = ix_["B"*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["W"*string(I)]
            pbm.A[ig,iv] += Float64(v_["-Z"])
            arrset(pbm.gscale,ig,Float64(200.0))
            ig,ig_,_ = s2mpj_ii("EC"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"EC"*string(I))
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
        for I = Int64(v_["0"]):Int64(v_["NI-1"])
            pbm.gconst[ig_["EA"*string(I)]] = Float64(232.0)
            pbm.gconst[ig_["EB"*string(I)]] = Float64(232.0)
        end
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        for I = Int64(v_["0"]):Int64(v_["NI-1"])
            pb.xlower[ix_["C"*string(I)]] = -Inf
            pb.xupper[ix_["C"*string(I)]] = +Inf
        end
        pb.xlower[ix_["TF"]] = 200.0
        for I = Int64(v_["0"]):Int64(v_["NI"])
            pb.xupper[ix_["W"*string(I)]] = v_["TOXIC"]
        end
        for I = Int64(v_["0"]):Int64(v_["NI-1"])
            pb.xupper[ix_["U"*string(I)]] = v_["UMAX"]
        end
        pb.xlower[ix_["W"*string(Int64(v_["0"]))]] = v_["WSS"]
        pb.xupper[ix_["W"*string(Int64(v_["0"]))]] = v_["WSS"]
        pb.xlower[ix_["W"*string(Int64(v_["NI"]))]] = v_["WSS"]
        pb.xupper[ix_["W"*string(Int64(v_["NI"]))]] = v_["WSS"]
        pb.xlower[ix_["P"*string(Int64(v_["0"]))]] = v_["PSTART"]
        pb.xupper[ix_["P"*string(Int64(v_["0"]))]] = v_["PSTART"]
        pb.xlower[ix_["P"*string(Int64(v_["NI"]))]] = v_["PFINAL"]
        pb.xupper[ix_["P"*string(Int64(v_["NI"]))]] = v_["PFINAL"]
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        v_["2W/10"] = 0.2*v_["WSS"]
        v_["2P/10"] = 0.2*v_["AVP"]
        v_["2(W+P)/10"] = v_["2W/10"]+v_["2P/10"]
        v_["D"] = 1.0+v_["2(W+P)/10"]
        v_["DD"] = v_["D"]*v_["D"]
        v_["ZP"] = v_["AVP"]*v_["Z"]
        v_["ZW"] = v_["WSS"]*v_["Z"]
        v_["AA"] = v_["DD"]+v_["ZP"]
        v_["AA"] = 232.0+v_["AA"]
        v_["BB"] = v_["DD"]+v_["ZW"]
        v_["BB"] = 232.0+v_["BB"]
        v_["AB"] = v_["AA"]*v_["BB"]
        v_["WP"] = v_["WSS"]*v_["AVP"]
        v_["-ZZWP"] = v_["WP"]*v_["-ZZ"]
        v_["CD"] = v_["AB"]+v_["-ZZWP"]
        v_["CC"] = v_["DD"]/v_["CD"]
        for I = Int64(v_["0"]):Int64(v_["NI-1"])
            pb.x0[ix_["W"*string(I)]] = Float64(v_["WSS"])
            pb.x0[ix_["P"*string(I)]] = Float64(v_["AVP"])
            pb.x0[ix_["U"*string(I)]] = Float64(v_["UMAX"])
            pb.x0[ix_["A"*string(I)]] = Float64(v_["AA"])
            pb.x0[ix_["B"*string(I)]] = Float64(v_["BB"])
            pb.x0[ix_["C"*string(I)]] = Float64(v_["CC"])
        end
        pb.x0[ix_["TF"]] = Float64(240.0)
        pb.x0[ix_["W"*string(Int64(v_["NI"]))]] = Float64(v_["WSS"])
        pb.x0[ix_["P"*string(Int64(v_["NI"]))]] = Float64(v_["PFINAL"])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "en3S", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        loaset(elftv,it,3,"V3")
        loaset(elftv,it,4,"V4")
        it,iet_,_ = s2mpj_ii( "en3D2", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        loaset(elftv,it,3,"V3")
        loaset(elftv,it,4,"V4")
        loaset(elftv,it,5,"V5")
        it,iet_,_ = s2mpj_ii( "eDSQ", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        it,iet_,_ = s2mpj_ii( "en3PR", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        loaset(elftv,it,3,"V3")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["0"]):Int64(v_["NI-1"])
            ename = "WA"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"en3S")
            arrset(ielftype,ie,iet_["en3S"])
            vname = "TF"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "C"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "A"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "W"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V4",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "WB"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"en3D2")
            arrset(ielftype,ie,iet_["en3D2"])
            vname = "TF"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "C"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "W"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "U"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V4",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "P"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V5",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "PA"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"en3D2")
            arrset(ielftype,ie,iet_["en3D2"])
            vname = "TF"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "C"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "B"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "U"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V4",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "P"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V5",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "PB"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"en3S")
            arrset(ielftype,ie,iet_["en3S"])
            vname = "TF"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "C"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "P"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "W"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V4",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "DD"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eDSQ")
            arrset(ielftype,ie,iet_["eDSQ"])
            vname = "W"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "P"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "CA"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"en3PR")
            arrset(ielftype,ie,iet_["en3PR"])
            vname = "C"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "A"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "B"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "CB"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"en3PR")
            arrset(ielftype,ie,iet_["en3PR"])
            vname = "C"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "P"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "W"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["0"]):Int64(v_["NI-1"])
            ig = ig_["EW"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["WA"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["-1/NI"]))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["WB"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["-Z/NI"]))
            ig = ig_["EP"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["PA"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["-1/NI"]))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["PB"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["-Z/NI"]))
            ig = ig_["EA"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["DD"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(-1.0))
            ig = ig_["EB"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["DD"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(-1.0))
            ig = ig_["EC"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["CA"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(1.0))
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["DD"*string(I)])
            loaset(pbm.grelw,ig,posel,Float64(-1.0))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["CB"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["-ZZ"]))
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 200.0
#    Solution
# LO SOLTN               ????
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
        pb.pbclass = "C-LOR2-MY-V-V"
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

    elseif action == "e_globs"

        pbm = args[1]
        arrset(pbm.efpar,1,0.02)
        return pbm

    elseif action == "en3S"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        WSSMV4 = pbm.efpar[1]-EV_[4]
        f_   = EV_[1]*EV_[2]*EV_[3]*WSSMV4
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[2]*EV_[3]*WSSMV4
            g_[2] = EV_[1]*EV_[3]*WSSMV4
            g_[3] = EV_[1]*EV_[2]*WSSMV4
            g_[4] = -EV_[1]*EV_[2]*EV_[3]
            if nargout>2
                H_ = zeros(Float64,4,4)
                H_[1,2] = EV_[3]*WSSMV4
                H_[2,1] = H_[1,2]
                H_[1,3] = EV_[2]*WSSMV4
                H_[3,1] = H_[1,3]
                H_[1,4] = -EV_[2]*EV_[3]
                H_[4,1] = H_[1,4]
                H_[2,3] = EV_[1]*WSSMV4
                H_[3,2] = H_[2,3]
                H_[2,4] = -EV_[1]*EV_[3]
                H_[4,2] = H_[2,4]
                H_[3,4] = -EV_[1]*EV_[2]
                H_[4,3] = H_[3,4]
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "en3D2"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,4,5)
        IV_ =  zeros(Float64,4)
        U_[1,1] = U_[1,1]+1
        U_[2,2] = U_[2,2]+1
        U_[3,3] = U_[3,3]+1
        U_[4,4] = U_[4,4]+1
        U_[4,5] = U_[4,5]-2
        IV_[1] = dot(U_[1,:],EV_)
        IV_[2] = dot(U_[2,:],EV_)
        IV_[3] = dot(U_[3,:],EV_)
        IV_[4] = dot(U_[4,:],EV_)
        f_   = IV_[1]*IV_[2]*IV_[3]*IV_[4]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = IV_[2]*IV_[3]*IV_[4]
            g_[2] = IV_[1]*IV_[3]*IV_[4]
            g_[3] = IV_[1]*IV_[2]*IV_[4]
            g_[4] = IV_[1]*IV_[2]*IV_[3]
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,4,4)
                H_[1,2] = IV_[3]*IV_[4]
                H_[2,1] = H_[1,2]
                H_[1,3] = IV_[2]*IV_[4]
                H_[3,1] = H_[1,3]
                H_[1,4] = IV_[2]*IV_[3]
                H_[4,1] = H_[1,4]
                H_[2,3] = IV_[1]*IV_[4]
                H_[3,2] = H_[2,3]
                H_[2,4] = IV_[1]*IV_[3]
                H_[4,2] = H_[2,4]
                H_[3,4] = IV_[1]*IV_[2]
                H_[4,3] = H_[3,4]
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

    elseif action == "eDSQ"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,1,2)
        IV_ =  zeros(Float64,1)
        U_[1,1] = U_[1,1]+2.000000e-01
        U_[1,2] = U_[1,2]+2.000000e-01
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

    elseif action == "en3PR"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]*EV_[2]*EV_[3]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[2]*EV_[3]
            g_[2] = EV_[1]*EV_[3]
            g_[3] = EV_[1]*EV_[2]
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,2] = EV_[3]
                H_[2,1] = H_[1,2]
                H_[1,3] = EV_[2]
                H_[3,1] = H_[1,3]
                H_[2,3] = EV_[1]
                H_[3,2] = H_[2,3]
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
            pbm.has_globs = [1,0]
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

