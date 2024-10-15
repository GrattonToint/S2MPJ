function HAIFAS(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HAIFAS
#    *********
# 
#   Truss Topology Design (t6-9)
# 
#   Source: M. Tsibulevsky, Optimization Laboratory,
#           Faculty of Industrial Engineering, Technion,
#           Haifa, 32000, Israel.
# 
#   SIF input: Conn, Gould and Toint, May, 1992
#              minor correction by Ph. Shott, Jan 1995.
# 
#    classification = "C-LQR2-AN-13-9"
# 
#   2 * Number of nodes
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "HAIFAS"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["N"] = 12
        v_["M"] = 9
        v_["0"] = 0
        v_["1"] = 1
        v_["2"] = 2
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        iv,ix_,_ = s2mpj_ii("Z",ix_)
        arrset(pb.xnames,iv,"Z")
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["Z"]
        pbm.A[ig,iv] += Float64(1.0)
        for I = Int64(v_["1"]):Int64(v_["M"])
            ig,ig_,_ = s2mpj_ii("G"*string(I),ig_)
            arrset(gtype,ig,"<=")
            arrset(pb.cnames,ig,"G"*string(I))
            iv = ix_["Z"]
            pbm.A[ig,iv] += Float64(-1.0)
            v_["J"] = 10
            iv = ix_["X"*string(Int64(v_["J"]))]
            pbm.A[ig,iv] += Float64(-1.00000)
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
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "en2PR", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        v_["I1"] = 4
        v_["I2"] = 4
        v_["L"] = 1
        ename = "E"*string(Int64(v_["L"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "X"*string(Int64(v_["I1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["L"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "X"*string(Int64(v_["I2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        v_["I1"] = 5
        v_["I2"] = 5
        v_["L"] = 2
        ename = "E"*string(Int64(v_["L"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "X"*string(Int64(v_["I1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["L"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "X"*string(Int64(v_["I2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        v_["I1"] = 5
        v_["I2"] = 11
        v_["L"] = 3
        ename = "E"*string(Int64(v_["L"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "X"*string(Int64(v_["I1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["L"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "X"*string(Int64(v_["I2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        v_["I1"] = 11
        v_["I2"] = 11
        v_["L"] = 4
        ename = "E"*string(Int64(v_["L"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "X"*string(Int64(v_["I1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["L"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "X"*string(Int64(v_["I2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        v_["I1"] = 10
        v_["I2"] = 10
        v_["L"] = 5
        ename = "E"*string(Int64(v_["L"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "X"*string(Int64(v_["I1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["L"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "X"*string(Int64(v_["I2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        v_["I1"] = 10
        v_["I2"] = 11
        v_["L"] = 6
        ename = "E"*string(Int64(v_["L"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "X"*string(Int64(v_["I1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["L"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "X"*string(Int64(v_["I2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        v_["I1"] = 11
        v_["I2"] = 11
        v_["L"] = 7
        ename = "E"*string(Int64(v_["L"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "X"*string(Int64(v_["I1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["L"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "X"*string(Int64(v_["I2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        v_["I1"] = 4
        v_["I2"] = 4
        v_["L"] = 8
        ename = "E"*string(Int64(v_["L"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "X"*string(Int64(v_["I1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["L"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "X"*string(Int64(v_["I2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        v_["I1"] = 4
        v_["I2"] = 10
        v_["L"] = 9
        ename = "E"*string(Int64(v_["L"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "X"*string(Int64(v_["I1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["L"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "X"*string(Int64(v_["I2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        v_["I1"] = 10
        v_["I2"] = 10
        v_["L"] = 10
        ename = "E"*string(Int64(v_["L"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "X"*string(Int64(v_["I1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["L"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "X"*string(Int64(v_["I2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        v_["I1"] = 5
        v_["I2"] = 5
        v_["L"] = 11
        ename = "E"*string(Int64(v_["L"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "X"*string(Int64(v_["I1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["L"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "X"*string(Int64(v_["I2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        v_["I1"] = 6
        v_["I2"] = 6
        v_["L"] = 12
        ename = "E"*string(Int64(v_["L"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "X"*string(Int64(v_["I1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["L"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "X"*string(Int64(v_["I2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        v_["I1"] = 6
        v_["I2"] = 12
        v_["L"] = 13
        ename = "E"*string(Int64(v_["L"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "X"*string(Int64(v_["I1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["L"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "X"*string(Int64(v_["I2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        v_["I1"] = 12
        v_["I2"] = 12
        v_["L"] = 14
        ename = "E"*string(Int64(v_["L"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "X"*string(Int64(v_["I1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["L"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "X"*string(Int64(v_["I2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        v_["I1"] = 11
        v_["I2"] = 11
        v_["L"] = 15
        ename = "E"*string(Int64(v_["L"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "X"*string(Int64(v_["I1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["L"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "X"*string(Int64(v_["I2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        v_["I1"] = 11
        v_["I2"] = 12
        v_["L"] = 16
        ename = "E"*string(Int64(v_["L"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "X"*string(Int64(v_["I1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["L"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "X"*string(Int64(v_["I2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        v_["I1"] = 12
        v_["I2"] = 12
        v_["L"] = 17
        ename = "E"*string(Int64(v_["L"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "X"*string(Int64(v_["I1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["L"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "X"*string(Int64(v_["I2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        v_["I1"] = 5
        v_["I2"] = 5
        v_["L"] = 18
        ename = "E"*string(Int64(v_["L"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "X"*string(Int64(v_["I1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["L"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "X"*string(Int64(v_["I2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        v_["I1"] = 5
        v_["I2"] = 11
        v_["L"] = 19
        ename = "E"*string(Int64(v_["L"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "X"*string(Int64(v_["I1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["L"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "X"*string(Int64(v_["I2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        v_["I1"] = 11
        v_["I2"] = 11
        v_["L"] = 20
        ename = "E"*string(Int64(v_["L"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "X"*string(Int64(v_["I1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["L"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "X"*string(Int64(v_["I2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        v_["I1"] = 6
        v_["I2"] = 6
        v_["L"] = 21
        ename = "E"*string(Int64(v_["L"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "X"*string(Int64(v_["I1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["L"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "X"*string(Int64(v_["I2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        v_["J"] = 1
        v_["L"] = 1
        ig = ig_["G"*string(Int64(v_["J"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["L"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(1.00000e+01))
        v_["J"] = 2
        v_["L"] = 2
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["L"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(6.40000))
        v_["J"] = 2
        v_["L"] = 3
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["L"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(6.40000))
        v_["J"] = 2
        v_["L"] = 4
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["L"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(1.60000))
        v_["J"] = 3
        v_["L"] = 5
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["L"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(4.00000e+01))
        v_["J"] = 3
        v_["L"] = 6
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["L"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-8.00000e+01))
        v_["J"] = 3
        v_["L"] = 7
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["L"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(4.00000e+01))
        v_["J"] = 4
        v_["L"] = 8
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["L"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(6.40000))
        v_["J"] = 4
        v_["L"] = 9
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["L"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-6.40000))
        v_["J"] = 4
        v_["L"] = 10
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["L"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(1.60000))
        v_["J"] = 5
        v_["L"] = 11
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["L"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(1.00000e+01))
        v_["J"] = 6
        v_["L"] = 12
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["L"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(6.40000))
        v_["J"] = 6
        v_["L"] = 13
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["L"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(6.40000))
        v_["J"] = 6
        v_["L"] = 14
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["L"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(1.60000))
        v_["J"] = 7
        v_["L"] = 15
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["L"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(4.00000e+01))
        v_["J"] = 7
        v_["L"] = 16
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["L"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-8.00000e+01))
        v_["J"] = 7
        v_["L"] = 17
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["L"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(4.00000e+01))
        v_["J"] = 8
        v_["L"] = 18
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["L"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(6.40000))
        v_["J"] = 8
        v_["L"] = 19
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["L"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-6.40000))
        v_["J"] = 8
        v_["L"] = 20
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["L"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(1.60000))
        v_["J"] = 9
        v_["L"] = 21
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["L"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(1.00000e+01))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.cupper[1:pb.nle] = zeros(Float64,pb.nle)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-LQR2-AN-13-9"
        pb.x0          = zeros(Float64,pb.n)
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

    elseif action == "en2PR"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = 0.5*EV_[1]*EV_[2]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 0.5*EV_[2]
            g_[2] = 0.5*EV_[1]
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[2,1] = 0.5
                H_[1,2] = H_[2,1]
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

