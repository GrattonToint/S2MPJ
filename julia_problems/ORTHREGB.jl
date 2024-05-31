function ORTHREGB(action,args...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : ORTHREGB
#    *********
# 
#    An orthogonal regression problem.
# 
#    The problem is to fit (orthogonally) an ellipse to a set of 6 points
#    in the 3D space. These points are compatible with this constraint.
# 
#    Source:
#    M. Gulliksson,
#    "Algorithms for nonlinear Least-squares with Applications to
#    Orthogonal Regression",
#    UMINF-178.90, University of Umea, Sweden, 1990.
# 
#    SIF input: Ph. Toint, June 1990.
#               correction by Ph. Shott, Jan 1995.
# 
#    classification = "QQR2-AN-27-6"
# 
#    Parameters for the generation of the data points
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "ORTHREGB"

    if action == "setup"
        pbm          = PBM(name)
        pb           = PB(name)
        pb.sifpbname = "ORTHREGB"
        nargin       = length(args)
        pbm.call     = eval( Meta.parse( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["A"] = 9.0
        v_["B"] = 6.0
        v_["C"] = 7.0
        v_["CX"] = 0.5
        v_["CY"] = 0.5
        v_["CZ"] = 0.5
        v_["1"] = 1
        v_["-A"] = -1.0*v_["A"]
        v_["-B"] = -1.0*v_["B"]
        v_["-C"] = -1.0*v_["C"]
        v_["NPTS"] = 1
        v_["XZ"] = v_["CX"]
        v_["YZ"] = v_["CY"]
        v_["ZZ"] = v_["CZ"]
        v_["NPTS"] = 1
        v_["XD"*string(Int64(v_["NPTS"]))] = v_["XZ"]+v_["A"]
        v_["YD"*string(Int64(v_["NPTS"]))] = v_["YZ"]+v_["A"]
        v_["ZD"*string(Int64(v_["NPTS"]))] = v_["ZZ"]
        v_["NPTS"] = 1+v_["NPTS"]
        v_["XD"*string(Int64(v_["NPTS"]))] = v_["XZ"]+v_["B"]
        v_["YD"*string(Int64(v_["NPTS"]))] = v_["YZ"]+v_["-B"]
        v_["ZD"*string(Int64(v_["NPTS"]))] = v_["ZZ"]
        v_["NPTS"] = 1+v_["NPTS"]
        v_["XD"*string(Int64(v_["NPTS"]))] = v_["XZ"]+v_["-A"]
        v_["YD"*string(Int64(v_["NPTS"]))] = v_["YZ"]+v_["-A"]
        v_["ZD"*string(Int64(v_["NPTS"]))] = v_["ZZ"]
        v_["NPTS"] = 1+v_["NPTS"]
        v_["XD"*string(Int64(v_["NPTS"]))] = v_["XZ"]+v_["-B"]
        v_["YD"*string(Int64(v_["NPTS"]))] = v_["YZ"]+v_["B"]
        v_["ZD"*string(Int64(v_["NPTS"]))] = v_["ZZ"]
        v_["NPTS"] = 1+v_["NPTS"]
        v_["XD"*string(Int64(v_["NPTS"]))] = v_["XZ"]
        v_["YD"*string(Int64(v_["NPTS"]))] = v_["YZ"]
        v_["ZD"*string(Int64(v_["NPTS"]))] = v_["ZZ"]+v_["C"]
        v_["NPTS"] = 1+v_["NPTS"]
        v_["XD"*string(Int64(v_["NPTS"]))] = v_["XZ"]
        v_["YD"*string(Int64(v_["NPTS"]))] = v_["YZ"]
        v_["ZD"*string(Int64(v_["NPTS"]))] = v_["ZZ"]+v_["-C"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        xscale  = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        iv,ix_,_ = s2x_ii("H11",ix_)
        arrset(pb.xnames,iv,"H11")
        iv,ix_,_ = s2x_ii("H12",ix_)
        arrset(pb.xnames,iv,"H12")
        iv,ix_,_ = s2x_ii("H13",ix_)
        arrset(pb.xnames,iv,"H13")
        iv,ix_,_ = s2x_ii("H22",ix_)
        arrset(pb.xnames,iv,"H22")
        iv,ix_,_ = s2x_ii("H23",ix_)
        arrset(pb.xnames,iv,"H23")
        iv,ix_,_ = s2x_ii("H33",ix_)
        arrset(pb.xnames,iv,"H33")
        iv,ix_,_ = s2x_ii("G1",ix_)
        arrset(pb.xnames,iv,"G1")
        iv,ix_,_ = s2x_ii("G2",ix_)
        arrset(pb.xnames,iv,"G2")
        iv,ix_,_ = s2x_ii("G3",ix_)
        arrset(pb.xnames,iv,"G3")
        for I = Int64(v_["1"]):Int64(v_["NPTS"])
            iv,ix_,_ = s2x_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
            iv,ix_,_ = s2x_ii("Y"*string(I),ix_)
            arrset(pb.xnames,iv,"Y"*string(I))
            iv,ix_,_ = s2x_ii("Z"*string(I),ix_)
            arrset(pb.xnames,iv,"Z"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["1"]):Int64(v_["NPTS"])
            ig,ig_,_ = s2x_ii("OX"*string(I),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            ig,ig_,_ = s2x_ii("OY"*string(I),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["Y"*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            ig,ig_,_ = s2x_ii("OZ"*string(I),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["Z"*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            ig,ig_,_ = s2x_ii("E"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"E"*string(I))
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
        pbm.congrps = findall(x->x!="<>",gtype)
        pb.nob = ngrp-pb.m
        pbm.objgrps = findall(x->x=="<>",gtype)
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        for I = Int64(v_["1"]):Int64(v_["NPTS"])
            pbm.gconst[ig_["OX"*string(I)]] = Float64(v_["XD"*string(I)])
            pbm.gconst[ig_["OY"*string(I)]] = Float64(v_["YD"*string(I)])
            pbm.gconst[ig_["OZ"*string(I)]] = Float64(v_["ZD"*string(I)])
            pbm.gconst[ig_["E"*string(I)]] = Float64(1.0)
        end
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        if haskey(ix_,"H11")
            pb.x0[ix_["H11"]] = Float64(1.0)
        else
            pb.y0[findfirst(x->x==ig_["H11"],pbm.congrps)] = Float64(1.0)
        end
        if haskey(ix_,"H12")
            pb.x0[ix_["H12"]] = Float64(0.0)
        else
            pb.y0[findfirst(x->x==ig_["H12"],pbm.congrps)] = Float64(0.0)
        end
        if haskey(ix_,"H13")
            pb.x0[ix_["H13"]] = Float64(0.0)
        else
            pb.y0[findfirst(x->x==ig_["H13"],pbm.congrps)] = Float64(0.0)
        end
        if haskey(ix_,"H22")
            pb.x0[ix_["H22"]] = Float64(1.0)
        else
            pb.y0[findfirst(x->x==ig_["H22"],pbm.congrps)] = Float64(1.0)
        end
        if haskey(ix_,"H23")
            pb.x0[ix_["H23"]] = Float64(0.0)
        else
            pb.y0[findfirst(x->x==ig_["H23"],pbm.congrps)] = Float64(0.0)
        end
        if haskey(ix_,"H33")
            pb.x0[ix_["H33"]] = Float64(1.0)
        else
            pb.y0[findfirst(x->x==ig_["H33"],pbm.congrps)] = Float64(1.0)
        end
        if haskey(ix_,"G1")
            pb.x0[ix_["G1"]] = Float64(0.0)
        else
            pb.y0[findfirst(x->x==ig_["G1"],pbm.congrps)] = Float64(0.0)
        end
        if haskey(ix_,"G2")
            pb.x0[ix_["G2"]] = Float64(0.0)
        else
            pb.y0[findfirst(x->x==ig_["G2"],pbm.congrps)] = Float64(0.0)
        end
        if haskey(ix_,"G3")
            pb.x0[ix_["G3"]] = Float64(0.0)
        else
            pb.y0[findfirst(x->x==ig_["G3"],pbm.congrps)] = Float64(0.0)
        end
        for I = Int64(v_["1"]):Int64(v_["NPTS"])
            if haskey(ix_,"X"*string(I))
                pb.x0[ix_["X"*string(I)]] = Float64(v_["XD"*string(I)])
            else
                pb.y0[findfirst(x->x==ig_["X"*string(I)],pbm.congrps)]  = (
                      Float64(v_["XD"*string(I)]))
            end
            if haskey(ix_,"Y"*string(I))
                pb.x0[ix_["Y"*string(I)]] = Float64(v_["YD"*string(I)])
            else
                pb.y0[findfirst(x->x==ig_["Y"*string(I)],pbm.congrps)]  = (
                      Float64(v_["YD"*string(I)]))
            end
            if haskey(ix_,"Z"*string(I))
                pb.x0[ix_["Z"*string(I)]] = Float64(v_["ZD"*string(I)])
            else
                pb.y0[findfirst(x->x==ig_["Z"*string(I)],pbm.congrps)]  = (
                      Float64(v_["ZD"*string(I)]))
            end
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2x_ii( "eHXX", iet_)
        loaset(elftv,it,1,"H")
        loaset(elftv,it,2,"X")
        it,iet_,_ = s2x_ii( "eHXY", iet_)
        loaset(elftv,it,1,"H")
        loaset(elftv,it,2,"X")
        loaset(elftv,it,3,"Y")
        it,iet_,_ = s2x_ii( "eGX", iet_)
        loaset(elftv,it,1,"G")
        loaset(elftv,it,2,"X")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["NPTS"])
            ename = "EA"*string(I)
            ie,ie_,_  = s2x_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eHXX")
            arrset(ielftype, ie, iet_["eHXX"])
            vname = "H11"
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="H",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(I)
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "EB"*string(I)
            ie,ie_,_  = s2x_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eHXY")
            arrset(ielftype, ie, iet_["eHXY"])
            vname = "H12"
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="H",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(I)
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "Y"*string(I)
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "EC"*string(I)
            ie,ie_,_  = s2x_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eHXX")
            arrset(ielftype, ie, iet_["eHXX"])
            vname = "H22"
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="H",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "Y"*string(I)
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "ED"*string(I)
            ie,ie_,_  = s2x_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eGX")
            arrset(ielftype, ie, iet_["eGX"])
            vname = "G1"
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="G",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(I)
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "EE"*string(I)
            ie,ie_,_  = s2x_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eGX")
            arrset(ielftype, ie, iet_["eGX"])
            vname = "G2"
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="G",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "Y"*string(I)
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "EF"*string(I)
            ie,ie_,_  = s2x_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eHXY")
            arrset(ielftype, ie, iet_["eHXY"])
            vname = "H13"
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="H",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(I)
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "Z"*string(I)
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "EG"*string(I)
            ie,ie_,_  = s2x_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eHXY")
            arrset(ielftype, ie, iet_["eHXY"])
            vname = "H23"
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="H",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "Y"*string(I)
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "Z"*string(I)
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "EH"*string(I)
            ie,ie_,_  = s2x_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eHXX")
            arrset(ielftype, ie, iet_["eHXX"])
            vname = "H33"
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="H",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "Z"*string(I)
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "EI"*string(I)
            ie,ie_,_  = s2x_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eGX")
            arrset(ielftype, ie, iet_["eGX"])
            vname = "G3"
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="G",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "Z"*string(I)
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2x_ii("gL2",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["NPTS"])
            ig = ig_["OX"*string(I)]
            arrset(pbm.grftype,ig,"gL2")
            ig = ig_["OY"*string(I)]
            arrset(pbm.grftype,ig,"gL2")
            ig = ig_["OZ"*string(I)]
            arrset(pbm.grftype,ig,"gL2")
            ig = ig_["E"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["EA"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["EB"*string(I)])
            loaset(pbm.grelw,ig,posel,Float64(2.0))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["EC"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["ED"*string(I)])
            loaset(pbm.grelw,ig,posel,Float64(-2.0))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["EE"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(-2.0))
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["EF"*string(I)])
            loaset(pbm.grelw,ig,posel,Float64(2.0))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["EG"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(2.0))
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["EH"*string(I)])
            loaset(pbm.grelw,ig,posel, 1.)
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["EI"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(-2.0))
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
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
        lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "QQR2-AN-27-6"
        return pb, pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "eHXX"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]*EV_[2]*EV_[2]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[2]*EV_[2]
            g_[2] = 2.0*EV_[1]*EV_[2]
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = EV_[2]+EV_[2]
                H_[2,1] = H_[1,2]
                H_[2,2] = EV_[1]+EV_[1]
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eHXY"

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

    elseif action == "eGX"

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

