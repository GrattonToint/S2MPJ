function FEEDLOC(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : FEEDLOC
#    *********
# 
#    Feed tray location & determination of optimum number of trays 
#    in a distillation column
# 
#    SIF input: S. Leyffer, October 1997
# 
#    classification = "C-LOR2-AN-90-259"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "FEEDLOC"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["1"] = 1
        v_["2"] = 2
        v_["3"] = 3
        v_["M"] = 2
        v_["NMAX"] = 12
        v_["NMAX-1"] = -1+v_["NMAX"]
        v_["F"] = 100.0
        v_["AL1"] = 1.0
        v_["AL2"] = 5.13435
        v_["XF1"] = 0.80
        v_["XF2"] = 0.20
        v_["SPEC"] = 0.001
        v_["BIGM"] = 1000.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["NMAX"])
            iv,ix_,_ = s2mpj_ii("S"*string(I),ix_)
            arrset(pb.xnames,iv,"S"*string(I))
            iv,ix_,_ = s2mpj_ii("W"*string(I),ix_)
            arrset(pb.xnames,iv,"W"*string(I))
            iv,ix_,_ = s2mpj_ii("Z"*string(I),ix_)
            arrset(pb.xnames,iv,"Z"*string(I))
        end
        iv,ix_,_ = s2mpj_ii("N",ix_)
        arrset(pb.xnames,iv,"N")
        for I = Int64(v_["1"]):Int64(v_["NMAX"])
            for J = Int64(v_["1"]):Int64(v_["M"])
                iv,ix_,_ = s2mpj_ii("X"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"X"*string(I)*","*string(J))
                iv,ix_,_ = s2mpj_ii("Y"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"Y"*string(I)*","*string(J))
            end
        end
        iv,ix_,_ = s2mpj_ii("L",ix_)
        arrset(pb.xnames,iv,"L")
        iv,ix_,_ = s2mpj_ii("V",ix_)
        arrset(pb.xnames,iv,"V")
        iv,ix_,_ = s2mpj_ii("R",ix_)
        arrset(pb.xnames,iv,"R")
        iv,ix_,_ = s2mpj_ii("P1",ix_)
        arrset(pb.xnames,iv,"P1")
        iv,ix_,_ = s2mpj_ii("P2",ix_)
        arrset(pb.xnames,iv,"P2")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["R"]
        pbm.A[ig,iv] += Float64(1.0)
        for I = Int64(v_["1"]):Int64(v_["NMAX"])
            ig,ig_,_ = s2mpj_ii("FENTR",ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"FENTR")
            iv = ix_["W"*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            ig,ig_,_ = s2mpj_ii("NTRAY",ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"NTRAY")
            iv = ix_["S"*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            ig,ig_,_ = s2mpj_ii("NDEF1",ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"NDEF1")
            iv = ix_["Z"*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            v_["RI"] = Float64(I)
            ig,ig_,_ = s2mpj_ii("NDEF2",ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"NDEF2")
            iv = ix_["S"*string(I)]
            pbm.A[ig,iv] += Float64(v_["RI"])
        end
        ig,ig_,_ = s2mpj_ii("NDEF1",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NDEF1")
        iv = ix_["N"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("NDEF2",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NDEF2")
        iv = ix_["N"]
        pbm.A[ig,iv] += Float64(-1.0)
        for I = Int64(v_["1"]):Int64(v_["NMAX-1"])
            v_["I+1"] = 1+I
            ig,ig_,_ = s2mpj_ii("NIL"*string(I),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"NIL"*string(I))
            iv = ix_["Z"*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["Z"*string(Int64(v_["I+1"]))]
            pbm.A[ig,iv] += Float64(-1.0)
        end
        for I = Int64(v_["1"]):Int64(v_["NMAX"])
            v_["RI"] = Float64(I)
            ig,ig_,_ = s2mpj_ii("ENTX",ig_)
            arrset(gtype,ig,"<=")
            arrset(pb.cnames,ig,"ENTX")
            iv = ix_["W"*string(I)]
            pbm.A[ig,iv] += Float64(v_["RI"])
            v_["RI"] = -1.0*v_["RI"]
            iv = ix_["S"*string(I)]
            pbm.A[ig,iv] += Float64(v_["RI"])
        end
        for I = Int64(v_["1"]):Int64(v_["NMAX"])
            ig,ig_,_ = s2mpj_ii("LASTX"*string(I),ig_)
            arrset(gtype,ig,"<=")
            arrset(pb.cnames,ig,"LASTX"*string(I))
            iv = ix_["S"*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["Z"*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
        end
        for I = Int64(v_["1"]):Int64(v_["NMAX"])
            ig,ig_,_ = s2mpj_ii("ZNOT"*string(I),ig_)
            arrset(gtype,ig,"<=")
            arrset(pb.cnames,ig,"ZNOT"*string(I))
            iv = ix_["Z"*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            for K = Int64(I):Int64(v_["NMAX"])
                ig,ig_,_ = s2mpj_ii("ZNOT"*string(I),ig_)
                arrset(gtype,ig,"<=")
                arrset(pb.cnames,ig,"ZNOT"*string(I))
                iv = ix_["S"*string(K)]
                pbm.A[ig,iv] += Float64(-1.0)
            end
        end
        for I = Int64(v_["1"]):Int64(v_["NMAX"])
            ig,ig_,_ = s2mpj_ii("FEEDX"*string(I),ig_)
            arrset(gtype,ig,"<=")
            arrset(pb.cnames,ig,"FEEDX"*string(I))
            iv = ix_["W"*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["Z"*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
        end
        for I = Int64(v_["2"]):Int64(v_["NMAX"])
            v_["I-1"] = -1+I
            ig,ig_,_ = s2mpj_ii("WNES1u"*string(I),ig_)
            arrset(gtype,ig,"<=")
            arrset(pb.cnames,ig,"WNES1u"*string(I))
            iv = ix_["W"*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["S"*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            ig,ig_,_ = s2mpj_ii("WNES2u"*string(I),ig_)
            arrset(gtype,ig,"<=")
            arrset(pb.cnames,ig,"WNES2u"*string(I))
            iv = ix_["W"*string(Int64(v_["I-1"]))]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["S"*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
        end
        for I = Int64(v_["1"]):Int64(v_["NMAX"])
            for J = Int64(v_["1"]):Int64(v_["M"])
                ig,ig_,_ = s2mpj_ii("PE1"*string(I),ig_)
                arrset(gtype,ig,">=")
                arrset(pb.cnames,ig,"PE1"*string(I))
                iv = ix_["Y"*string(I)*","*string(J)]
                pbm.A[ig,iv] += Float64(1.0)
                ig,ig_,_ = s2mpj_ii("PE2"*string(I),ig_)
                arrset(gtype,ig,">=")
                arrset(pb.cnames,ig,"PE2"*string(I))
                iv = ix_["Y"*string(I)*","*string(J)]
                pbm.A[ig,iv] += Float64(1.0)
                ig,ig_,_ = s2mpj_ii("PE3"*string(I),ig_)
                arrset(gtype,ig,">=")
                arrset(pb.cnames,ig,"PE3"*string(I))
                iv = ix_["X"*string(I)*","*string(J)]
                pbm.A[ig,iv] += Float64(1.0)
                ig,ig_,_ = s2mpj_ii("PE4"*string(I),ig_)
                arrset(gtype,ig,">=")
                arrset(pb.cnames,ig,"PE4"*string(I))
                iv = ix_["X"*string(I)*","*string(J)]
                pbm.A[ig,iv] += Float64(1.0)
            end
            ig,ig_,_ = s2mpj_ii("PE1"*string(I),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"PE1"*string(I))
            iv = ix_["Z"*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            ig,ig_,_ = s2mpj_ii("PE2"*string(I),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"PE2"*string(I))
            iv = ix_["Z"*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("PE3"*string(I),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"PE3"*string(I))
            iv = ix_["Z"*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            ig,ig_,_ = s2mpj_ii("PE4"*string(I),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"PE4"*string(I))
            iv = ix_["Z"*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            for J = Int64(v_["1"]):Int64(v_["M"])
                ig,ig_,_ = s2mpj_ii("XNOT"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"<=")
                arrset(pb.cnames,ig,"XNOT"*string(I)*","*string(J))
                iv = ix_["X"*string(I)*","*string(J)]
                pbm.A[ig,iv] += Float64(1.0)
                iv = ix_["Z"*string(I)]
                pbm.A[ig,iv] += Float64(-1.0)
                ig,ig_,_ = s2mpj_ii("YNOT"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"<=")
                arrset(pb.cnames,ig,"YNOT"*string(I)*","*string(J))
                iv = ix_["Y"*string(I)*","*string(J)]
                pbm.A[ig,iv] += Float64(1.0)
                iv = ix_["Z"*string(I)]
                pbm.A[ig,iv] += Float64(-1.0)
            end
        end
        for I = Int64(v_["1"]):Int64(v_["NMAX"])
            v_["TEMP"] = -1.0*v_["AL1"]
            ig,ig_,_ = s2mpj_ii("PHEE"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"PHEE"*string(I))
            iv = ix_["X"*string(I)*","*string(Int64(v_["1"]))]
            pbm.A[ig,iv] += Float64(v_["TEMP"])
        end
        ig,ig_,_ = s2mpj_ii("DEFL",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"DEFL")
        iv = ix_["L"]
        pbm.A[ig,iv] += Float64(1.0)
        for J = Int64(v_["1"]):Int64(v_["M"])
            v_["TEMP"] = -1.0*v_["F"]
            ig,ig_,_ = s2mpj_ii("CMB1u"*string(J),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"CMB1u"*string(J))
            iv = ix_["X"*string(Int64(v_["2"]))*","*string(J)]
            pbm.A[ig,iv] += Float64(v_["TEMP"])
        end
        for I = Int64(v_["2"]):Int64(v_["NMAX"])
            for J = Int64(v_["1"]):Int64(v_["M"])
                v_["TEMP"] = -1.0*v_["BIGM"]
                ig,ig_,_ = s2mpj_ii("CMBN1"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"<=")
                arrset(pb.cnames,ig,"CMBN1"*string(I)*","*string(J))
                iv = ix_["S"*string(I)]
                pbm.A[ig,iv] += Float64(v_["BIGM"])
                ig,ig_,_ = s2mpj_ii("CMBN2"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,">=")
                arrset(pb.cnames,ig,"CMBN2"*string(I)*","*string(J))
                iv = ix_["S"*string(I)]
                pbm.A[ig,iv] += Float64(v_["TEMP"])
            end
        end
        for I = Int64(v_["2"]):Int64(v_["NMAX-1"])
            v_["TEMP1"] = v_["F"]*v_["XF"*string(Int64(v_["M"]))]
            v_["TEMP1"] = -1.0*v_["TEMP1"]
            ig,ig_,_ = s2mpj_ii("CMB1"*string(I),ig_)
            arrset(gtype,ig,"<=")
            arrset(pb.cnames,ig,"CMB1"*string(I))
            iv = ix_["S"*string(I)]
            pbm.A[ig,iv] += Float64(v_["TEMP"])
            iv = ix_["Z"*string(I)]
            pbm.A[ig,iv] += Float64(v_["BIGM"])
            iv = ix_["W"*string(I)]
            pbm.A[ig,iv] += Float64(v_["TEMP1"])
            ig,ig_,_ = s2mpj_ii("CMB2"*string(I),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"CMB2"*string(I))
            iv = ix_["S"*string(I)]
            pbm.A[ig,iv] += Float64(v_["BIGM"])
            iv = ix_["Z"*string(I)]
            pbm.A[ig,iv] += Float64(v_["TEMP"])
            iv = ix_["W"*string(I)]
            pbm.A[ig,iv] += Float64(v_["TEMP1"])
        end
        for I = Int64(v_["3"]):Int64(v_["NMAX"])
            ig,ig_,_ = s2mpj_ii("RECR"*string(I),ig_)
            arrset(gtype,ig,"<=")
            arrset(pb.cnames,ig,"RECR"*string(I))
            iv = ix_["S"*string(I)]
            pbm.A[ig,iv] += Float64(v_["BIGM"])
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
        pbm.gconst[ig_["FENTR"]] = Float64(1.0)
        pbm.gconst[ig_["NTRAY"]] = Float64(1.0)
        for I = Int64(v_["2"]):Int64(v_["NMAX"])
            pbm.gconst[ig_["WNES1u"*string(I)]] = Float64(1.0)
            pbm.gconst[ig_["WNES2u"*string(I)]] = Float64(1.0)
        end
        v_["TEMP"] = -1.0*v_["BIGM"]
        for I = Int64(v_["2"]):Int64(v_["NMAX"])
            for J = Int64(v_["1"]):Int64(v_["M"])
                pbm.gconst[ig_["CMBN1"*string(I)*","*string(J)]] = Float64(v_["BIGM"])
                pbm.gconst[ig_["CMBN2"*string(I)*","*string(J)]] = Float64(v_["TEMP"])
            end
        end
        for I = Int64(v_["2"]):Int64(v_["NMAX-1"])
            pbm.gconst[ig_["CMB1"*string(I)]] = Float64(v_["BIGM"])
            v_["TEMP"] = -1.0*v_["BIGM"]
            pbm.gconst[ig_["CMB2"*string(I)]] = Float64(v_["TEMP"])
        end
        v_["TEMP"] = v_["XF"*string(Int64(v_["1"]))]*v_["SPEC"]
        v_["TEMP1"] = v_["TEMP"]*v_["F"]
        v_["RHS"] = v_["TEMP1"]+v_["BIGM"]
        for I = Int64(v_["3"]):Int64(v_["NMAX"])
            pbm.gconst[ig_["RECR"*string(I)]] = Float64(v_["RHS"])
        end
        #%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%
        grange = Vector{Float64}(undef,ngrp)
        grange[legrps,1] = fill(Inf,pb.nle)
        grange[gegrps,1] = fill(Inf,pb.nge)
        for I = Int64(v_["1"]):Int64(v_["NMAX"])
            arrset(grange,ig_["PE1"*string(I)],Float64(2.0))
            arrset(grange,ig_["PE2"*string(I)],Float64(2.0))
            arrset(grange,ig_["PE3"*string(I)],Float64(2.0))
            arrset(grange,ig_["PE4"*string(I)],Float64(2.0))
        end
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        for I = Int64(v_["1"]):Int64(v_["NMAX"])
            pb.xupper[ix_["Z"*string(I)]] = 1.0
            pb.xupper[ix_["W"*string(I)]] = 1.0
            pb.xupper[ix_["S"*string(I)]] = 1.0
            for J = Int64(v_["1"]):Int64(v_["M"])
                pb.xupper[ix_["X"*string(I)*","*string(J)]] = 1.0
                pb.xupper[ix_["Y"*string(I)*","*string(J)]] = 1.0
            end
        end
        pb.xlower[ix_["N"]] = 3.0
        v_["TEMP"] = Float64(v_["NMAX"])
        pb.xupper[ix_["N"]] = v_["TEMP"]
        pb.xlower[ix_["P2"]] = 80.0
        pb.xupper[ix_["P2"]] = 80.0
        pb.xupper[ix_["L"]] = v_["F"]
        pb.xupper[ix_["V"]] = v_["F"]
        pb.xupper[ix_["P1"]] = v_["F"]
        pb.xupper[ix_["R"]] = 5.0
        pb.xlower[ix_["W"*string(Int64(v_["1"]))]] = 0.0
        pb.xupper[ix_["W"*string(Int64(v_["1"]))]] = 0.0
        pb.xlower[ix_["W"*string(Int64(v_["2"]))]] = 0.0
        pb.xupper[ix_["W"*string(Int64(v_["2"]))]] = 0.0
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(0.5),pb.n)
        pb.y0 = fill(Float64(0.5),pb.m)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eA2PROD", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"A")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["NMAX"])
            for K = Int64(v_["1"]):Int64(v_["M"])
                ename = "PHE"*string(I)*","*string(K)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eA2PROD")
                arrset(ielftype,ie,iet_["eA2PROD"])
                vname = "X"*string(I)*","*string(K)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
                posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "Y"*string(I)*","*string(Int64(v_["1"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
                posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                posep = findfirst(x->x=="A",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["AL"*string(K)]))
            end
        end
        ename = "DEFLE"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eA2PROD")
        arrset(ielftype,ie,iet_["eA2PROD"])
        vname = "R"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "P1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="A",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-1.0))
        for J = Int64(v_["1"]):Int64(v_["M"])
            ename = "CMB11u"*string(J)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eA2PROD")
            arrset(ielftype,ie,iet_["eA2PROD"])
            vname = "P2"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["1"]))*","*string(J)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="A",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(1.0))
            ename = "CMB12u"*string(J)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eA2PROD")
            arrset(ielftype,ie,iet_["eA2PROD"])
            vname = "V"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "Y"*string(Int64(v_["1"]))*","*string(J)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="A",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(1.0))
            ename = "CMB13u"*string(J)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eA2PROD")
            arrset(ielftype,ie,iet_["eA2PROD"])
            vname = "L"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["2"]))*","*string(J)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="A",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(-1.0))
        end
        for I = Int64(v_["2"]):Int64(v_["NMAX"])
            v_["I-1"] = -1+I
            for J = Int64(v_["1"]):Int64(v_["M"])
                ename = "CM11"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eA2PROD")
                arrset(ielftype,ie,iet_["eA2PROD"])
                vname = "L"
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
                posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
                posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                posep = findfirst(x->x=="A",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(1.0))
                ename = "CM12"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eA2PROD")
                arrset(ielftype,ie,iet_["eA2PROD"])
                vname = "P1"
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
                posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "Y"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
                posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                posep = findfirst(x->x=="A",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(1.0))
                ename = "CM13"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eA2PROD")
                arrset(ielftype,ie,iet_["eA2PROD"])
                vname = "V"
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
                posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "Y"*string(Int64(v_["I-1"]))*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
                posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                posep = findfirst(x->x=="A",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(-1.0))
                ename = "CM21"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eA2PROD")
                arrset(ielftype,ie,iet_["eA2PROD"])
                vname = "L"
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
                posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
                posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                posep = findfirst(x->x=="A",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(1.0))
                ename = "CM22"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eA2PROD")
                arrset(ielftype,ie,iet_["eA2PROD"])
                vname = "P1"
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
                posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "Y"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
                posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                posep = findfirst(x->x=="A",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(1.0))
                ename = "CM23"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eA2PROD")
                arrset(ielftype,ie,iet_["eA2PROD"])
                vname = "V"
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
                posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "Y"*string(Int64(v_["I-1"]))*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
                posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                posep = findfirst(x->x=="A",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(-1.0))
            end
        end
        for I = Int64(v_["2"]):Int64(v_["NMAX-1"])
            v_["I-1"] = -1+I
            v_["I+1"] = 1+I
            ename = "C11"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eA2PROD")
            arrset(ielftype,ie,iet_["eA2PROD"])
            vname = "L"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(I)*","*string(Int64(v_["1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="A",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(1.0))
            for K = Int64(I):Int64(v_["NMAX"])
                ename = "C12"*string(I)*","*string(K)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eA2PROD")
                arrset(ielftype,ie,iet_["eA2PROD"])
                vname = "X"*string(I)*","*string(Int64(v_["1"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
                posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "W"*string(K)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
                posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                posep = findfirst(x->x=="A",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["F"]))
            end
            ename = "C13"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eA2PROD")
            arrset(ielftype,ie,iet_["eA2PROD"])
            vname = "V"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "Y"*string(I)*","*string(Int64(v_["1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="A",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(1.0))
            ename = "C14"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eA2PROD")
            arrset(ielftype,ie,iet_["eA2PROD"])
            vname = "L"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["I+1"]))*","*string(Int64(v_["1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="A",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(-1.0))
            for K = Int64(v_["I+1"]):Int64(v_["NMAX"])
                v_["TEMP"] = -1.0*v_["F"]
                ename = "C15"*string(I)*","*string(K)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eA2PROD")
                arrset(ielftype,ie,iet_["eA2PROD"])
                vname = "X"*string(Int64(v_["I+1"]))*","*string(Int64(v_["1"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
                posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "W"*string(K)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
                posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                posep = findfirst(x->x=="A",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["TEMP"]))
            end
            ename = "C16"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eA2PROD")
            arrset(ielftype,ie,iet_["eA2PROD"])
            vname = "V"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "Y"*string(Int64(v_["I-1"]))*","*string(Int64(v_["1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="A",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(-1.0))
            ename = "C21"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eA2PROD")
            arrset(ielftype,ie,iet_["eA2PROD"])
            vname = "L"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(I)*","*string(Int64(v_["1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="A",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(1.0))
            for K = Int64(v_["1"]):Int64(v_["NMAX"])
                ename = "C22"*string(I)*","*string(K)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eA2PROD")
                arrset(ielftype,ie,iet_["eA2PROD"])
                vname = "X"*string(I)*","*string(Int64(v_["1"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
                posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "W"*string(K)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
                posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                posep = findfirst(x->x=="A",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["F"]))
            end
            ename = "C23"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eA2PROD")
            arrset(ielftype,ie,iet_["eA2PROD"])
            vname = "V"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "Y"*string(I)*","*string(Int64(v_["1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="A",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(1.0))
            ename = "C24"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eA2PROD")
            arrset(ielftype,ie,iet_["eA2PROD"])
            vname = "L"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["I+1"]))*","*string(Int64(v_["1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="A",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(-1.0))
            for K = Int64(v_["I+1"]):Int64(v_["NMAX"])
                v_["TEMP"] = -1.0*v_["F"]
                ename = "C25"*string(I)*","*string(K)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eA2PROD")
                arrset(ielftype,ie,iet_["eA2PROD"])
                vname = "X"*string(Int64(v_["I+1"]))*","*string(Int64(v_["1"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
                posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "W"*string(K)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
                posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                posep = findfirst(x->x=="A",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["TEMP"]))
            end
            ename = "C26"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eA2PROD")
            arrset(ielftype,ie,iet_["eA2PROD"])
            vname = "V"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "Y"*string(Int64(v_["I-1"]))*","*string(Int64(v_["1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="A",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(-1.0))
        end
        for I = Int64(v_["3"]):Int64(v_["NMAX"])
            ename = "REC"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eA2PROD")
            arrset(ielftype,ie,iet_["eA2PROD"])
            vname = "P1"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "Y"*string(I)*","*string(Int64(v_["1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="A",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(1.0))
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["NMAX"])
            for K = Int64(v_["1"]):Int64(v_["M"])
                ig = ig_["PHEE"*string(I)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["PHE"*string(I)*","*string(K)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,1.)
            end
        end
        ig = ig_["DEFL"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["DEFLE"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        for J = Int64(v_["1"]):Int64(v_["M"])
            ig = ig_["CMB1u"*string(J)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["CMB11u"*string(J)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["CMB12u"*string(J)])
            loaset(pbm.grelw,ig,posel, 1.)
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["CMB13u"*string(J)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        for I = Int64(v_["2"]):Int64(v_["NMAX"])
            for J = Int64(v_["1"]):Int64(v_["M"])
                ig = ig_["CMBN1"*string(I)*","*string(J)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["CM11"*string(I)*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,1.)
                posel = posel+1
                loaset(pbm.grelt,ig,posel,ie_["CM12"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel, 1.)
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["CM13"*string(I)*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,1.)
                ig = ig_["CMBN2"*string(I)*","*string(J)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["CM21"*string(I)*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,1.)
                posel = posel+1
                loaset(pbm.grelt,ig,posel,ie_["CM22"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel, 1.)
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["CM23"*string(I)*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,1.)
            end
        end
        for I = Int64(v_["2"]):Int64(v_["NMAX-1"])
            ig = ig_["CMB1"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["C11"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["C13"*string(I)])
            loaset(pbm.grelw,ig,posel, 1.)
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["C14"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["C16"*string(I)])
            loaset(pbm.grelw,ig,posel, 1.)
            ig = ig_["CMB2"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["C21"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["C23"*string(I)])
            loaset(pbm.grelw,ig,posel, 1.)
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["C24"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["C26"*string(I)])
            loaset(pbm.grelw,ig,posel, 1.)
            for K = Int64(I):Int64(v_["NMAX"])
                ig = ig_["CMB1"*string(I)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["C12"*string(I)*","*string(K)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,1.)
                ig = ig_["CMB2"*string(I)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["C22"*string(I)*","*string(K)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,1.)
            end
            v_["I+1"] = 1+I
            for K = Int64(v_["I+1"]):Int64(v_["NMAX"])
                ig = ig_["CMB1"*string(I)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["C15"*string(I)*","*string(K)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,1.)
                ig = ig_["CMB2"*string(I)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["C25"*string(I)*","*string(K)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,1.)
            end
        end
        for I = Int64(v_["3"]):Int64(v_["NMAX"])
            ig = ig_["RECR"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["REC"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[1:pb.nle] = grange[legrps]
        pb.cupper[1:pb.nle] = zeros(Float64,pb.nle)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.clower[pb.nle+pb.neq+1:pb.m] = zeros(Float64,pb.nge)
        pb.cupper[1:pb.nge] = grange[gegrps]
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-LOR2-AN-90-259"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2]
        pb.conderlvl  = pbm.conderlvl;
        return pb, pbm


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "eA2PROD"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = pbm.elpar[iel_][1]*EV_[1]*EV_[2]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = pbm.elpar[iel_][1]*EV_[2]
            g_[2] = pbm.elpar[iel_][1]*EV_[1]
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = pbm.elpar[iel_][1]
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

