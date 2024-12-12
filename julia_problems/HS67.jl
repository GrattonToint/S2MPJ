function HS67(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS67
#    *********
# 
#    Source: problem 67 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    Original Source: problem 8 in
#    A.R. Colville
#    "A comparative study on nonlinear programming"
#    IBM Scientific Center Report 320-2949, New York, 1968.
# 
#    SIF input: A.R. Conn & Nick Gould, April 1991.
%#   S2MPJ adaptation: Ph. Toint, November 2024.
# 
#    classification = "C-COOI2-AN-3-14"
# 
#    Set useful parameters
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 25 XI 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "HS67"

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
        v_["4"] = 4
        v_["5"] = 5
        v_["6"] = 6
        v_["7"] = 7
        v_["8"] = 8
        v_["9"] = 9
        v_["10"] = 10
        v_["11"] = 11
        v_["12"] = 12
        v_["13"] = 13
        v_["14"] = 14
        v_["A"*string(Int64(v_["1"]))] = 0.0
        v_["A"*string(Int64(v_["2"]))] = 0.0
        v_["A"*string(Int64(v_["3"]))] = 85.0
        v_["A"*string(Int64(v_["4"]))] = 90.0
        v_["A"*string(Int64(v_["5"]))] = 3.0
        v_["A"*string(Int64(v_["6"]))] = 0.01
        v_["A"*string(Int64(v_["7"]))] = 145.0
        v_["A"*string(Int64(v_["8"]))] = 5000.0
        v_["A"*string(Int64(v_["9"]))] = 2000.0
        v_["A"*string(Int64(v_["10"]))] = 93.0
        v_["A"*string(Int64(v_["11"]))] = 95.0
        v_["A"*string(Int64(v_["12"]))] = 12.0
        v_["A"*string(Int64(v_["13"]))] = 4.0
        v_["A"*string(Int64(v_["14"]))] = 162.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        irA   = Int64[]
        icA   = Int64[]
        valA  = Float64[]
        for I = Int64(v_["1"]):Int64(v_["3"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X"*string(Int64(v_["1"]))])
        push!(valA,Float64(5.04))
        push!(irA,ig)
        push!(icA,ix_["X"*string(Int64(v_["2"]))])
        push!(valA,Float64(0.035))
        push!(irA,ig)
        push!(icA,ix_["X"*string(Int64(v_["3"]))])
        push!(valA,Float64(10.0))
        for I = Int64(v_["1"]):Int64(v_["7"])
            ig,ig_,_ = s2mpj_ii("AG"*string(I),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"AG"*string(I))
            ig,ig_,_ = s2mpj_ii("AL"*string(I),ig_)
            arrset(gtype,ig,"<=")
            arrset(pb.cnames,ig,"AL"*string(I))
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
        for I = Int64(v_["1"]):Int64(v_["7"])
            v_["I+7"] = 7+I
            pbm.gconst[ig_["AG"*string(I)]] = Float64(v_["A"*string(I)])
            pbm.gconst[ig_["AL"*string(I)]] = Float64(v_["A"*string(Int64(v_["I+7"]))])
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = fill(0.00001,pb.n)
        pb.xupper = fill(Inf,pb.n)
        pb.xupper[ix_["X"*string(Int64(v_["1"]))]] = 2000.0
        pb.xupper[ix_["X"*string(Int64(v_["2"]))]] = 16000.0
        pb.xupper[ix_["X"*string(Int64(v_["3"]))]] = 120.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        pb.x0[ix_["X"*string(Int64(v_["1"]))]] = Float64(1745.0)
        pb.x0[ix_["X"*string(Int64(v_["2"]))]] = Float64(12000.0)
        pb.x0[ix_["X"*string(Int64(v_["3"]))]] = Float64(110.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eY2Y5", iet_)
        loaset(elftv,it,1,"U1")
        loaset(elftv,it,2,"U2")
        loaset(elftv,it,3,"U3")
        it,iet_,_ = s2mpj_ii( "eY2", iet_)
        loaset(elftv,it,1,"U1")
        loaset(elftv,it,2,"U2")
        loaset(elftv,it,3,"U3")
        it,iet_,_ = s2mpj_ii( "eY3", iet_)
        loaset(elftv,it,1,"U1")
        loaset(elftv,it,2,"U2")
        loaset(elftv,it,3,"U3")
        it,iet_,_ = s2mpj_ii( "eY4", iet_)
        loaset(elftv,it,1,"U1")
        loaset(elftv,it,2,"U2")
        loaset(elftv,it,3,"U3")
        it,iet_,_ = s2mpj_ii( "eY5", iet_)
        loaset(elftv,it,1,"U1")
        loaset(elftv,it,2,"U2")
        loaset(elftv,it,3,"U3")
        it,iet_,_ = s2mpj_ii( "eY6", iet_)
        loaset(elftv,it,1,"U1")
        loaset(elftv,it,2,"U2")
        loaset(elftv,it,3,"U3")
        it,iet_,_ = s2mpj_ii( "eY7", iet_)
        loaset(elftv,it,1,"U1")
        loaset(elftv,it,2,"U2")
        loaset(elftv,it,3,"U3")
        it,iet_,_ = s2mpj_ii( "eY8", iet_)
        loaset(elftv,it,1,"U1")
        loaset(elftv,it,2,"U2")
        loaset(elftv,it,3,"U3")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "E25"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eY2Y5")
        arrset(ielftype,ie,iet_["eY2Y5"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(0.00001),nothing,nothing)
        posev = findfirst(x->x=="U1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(0.00001),nothing,nothing)
        posev = findfirst(x->x=="U2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(0.00001),nothing,nothing)
        posev = findfirst(x->x=="U3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E2"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eY2")
        arrset(ielftype,ie,iet_["eY2"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(0.00001),nothing,nothing)
        posev = findfirst(x->x=="U1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(0.00001),nothing,nothing)
        posev = findfirst(x->x=="U2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(0.00001),nothing,nothing)
        posev = findfirst(x->x=="U3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E3"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eY3")
        arrset(ielftype,ie,iet_["eY3"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(0.00001),nothing,nothing)
        posev = findfirst(x->x=="U1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(0.00001),nothing,nothing)
        posev = findfirst(x->x=="U2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(0.00001),nothing,nothing)
        posev = findfirst(x->x=="U3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E4"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eY4")
        arrset(ielftype,ie,iet_["eY4"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(0.00001),nothing,nothing)
        posev = findfirst(x->x=="U1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(0.00001),nothing,nothing)
        posev = findfirst(x->x=="U2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(0.00001),nothing,nothing)
        posev = findfirst(x->x=="U3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E5"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eY5")
        arrset(ielftype,ie,iet_["eY5"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(0.00001),nothing,nothing)
        posev = findfirst(x->x=="U1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(0.00001),nothing,nothing)
        posev = findfirst(x->x=="U2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(0.00001),nothing,nothing)
        posev = findfirst(x->x=="U3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E6"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eY6")
        arrset(ielftype,ie,iet_["eY6"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(0.00001),nothing,nothing)
        posev = findfirst(x->x=="U1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(0.00001),nothing,nothing)
        posev = findfirst(x->x=="U2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(0.00001),nothing,nothing)
        posev = findfirst(x->x=="U3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E7"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eY7")
        arrset(ielftype,ie,iet_["eY7"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(0.00001),nothing,nothing)
        posev = findfirst(x->x=="U1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(0.00001),nothing,nothing)
        posev = findfirst(x->x=="U2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(0.00001),nothing,nothing)
        posev = findfirst(x->x=="U3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E8"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eY8")
        arrset(ielftype,ie,iet_["eY8"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(0.00001),nothing,nothing)
        posev = findfirst(x->x=="U1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(0.00001),nothing,nothing)
        posev = findfirst(x->x=="U2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(0.00001),nothing,nothing)
        posev = findfirst(x->x=="U3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["OBJ"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E25"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-0.063))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E3"])
        loaset(pbm.grelw,ig,posel,Float64(3.36))
        ig = ig_["AG"*string(Int64(v_["1"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E2"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["AL"*string(Int64(v_["1"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E2"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["AG"*string(Int64(v_["2"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E3"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["AL"*string(Int64(v_["2"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E3"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["AG"*string(Int64(v_["3"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E4"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["AL"*string(Int64(v_["3"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E4"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["AG"*string(Int64(v_["4"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E5"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["AL"*string(Int64(v_["4"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E5"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["AG"*string(Int64(v_["5"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E6"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["AL"*string(Int64(v_["5"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E6"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["AG"*string(Int64(v_["6"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E7"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["AL"*string(Int64(v_["6"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E7"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["AG"*string(Int64(v_["7"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E8"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["AL"*string(Int64(v_["7"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E8"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.cupper[1:pb.nle] = zeros(Float64,pb.nle)
        pb.clower[pb.nle+pb.neq+1:pb.m] = zeros(Float64,pb.nge)
        pb.cupper[1:pb.nge] = fill(Inf,pb.nge)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-COOI2-AN-3-14"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2]
        pb.conderlvl  = pbm.conderlvl;
        return pb, pbm


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "eY2Y5"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        Y,G,H   = pbm.call("extfunc",EV_[1],EV_[2],EV_[3])
        f_      = Y[2]*Y[5]
        if nargout>1
            dim   = length(EV_)
            g_    = zeros(Float64,dim)
            g_[1] = Y[2]*G[5,1]+Y[5]*G[2,1]
            g_[2] = Y[2]*G[5,2]+Y[5]*G[2,2]
            g_[3] = Y[2]*G[5,3]+Y[5]*G[2,3]
            if nargout>2
                H_      = zeros(Float64,3,3)
                H_[1,1] = Y[2]*H[5,1]+2.0*G[5,1]*G[2,1]+Y[5]*H[2,1]
                H_[1,2] = Y[2]*H[5,2]+G[5,1]*G[2,2]+G[5,2]*G[2,1]+Y[5]*H[2,2]
                H_[2,1] = H_[1,2]
                H_[1,3] = Y[2]*H[5,3]+G[5,1]*G[2,3]+G[5,3]*G[2,1]+Y[5]*H[2,3]
                H_[3,1] = H_[1,3]
                H_[2,2] = Y[2]*H[5,4]+2.0*G[5,2]*G[2,2]+Y[5]*H[2,4]
                H_[2,3] = Y[2]*H[5,5]+G[5,2]*G[2,3]+G[5,3]*G[2,2]+Y[5]*H[2,5]
                H_[3,2] = H_[2,3]
                H_[3,3] = Y[2]*H[5,6]+2.0*G[5,3]*G[2,3]+Y[5]*H[2,6]
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eY2"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        Y,G,H   = pbm.call("extfunc",EV_[1],EV_[2],EV_[3])
        f_      = Y[2]
        if nargout>1
            dim   = length(EV_)
            g_    = zeros(Float64,dim)
            g_[1] = G[2,1]
            g_[2] = G[2,2]
            g_[3] = 0.0
            if nargout>2
                H_      = zeros(Float64,3,3)
                H_[1,1] = H[2,1]
                H_[2,1] = H[2,2]
                H_[1,2] = H_[2,1]
                H_[2,2] = H[2,4]
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eY3"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        Y,G,H   = pbm.call("extfunc",EV_[1],EV_[2],EV_[3])
        f_      = Y[3]
        if nargout>1
            dim   = length(EV_)
            g_    = zeros(Float64,dim)
            g_[1] = G[3,1]
            g_[2] = G[3,2]
            g_[3] = 0.0
            if nargout>2
                H_      = zeros(Float64,3,3)
                H_[1,1] = H[3,1]
                H_[2,1] = H[3,2]
                H_[1,2] = H_[2,1]
                H_[2,2] = H[3,4]
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eY4"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        Y,G,H   = pbm.call("extfunc",EV_[1],EV_[2],EV_[3])
        f_      = Y[4]
        if nargout>1
            dim   = length(EV_)
            g_    = zeros(Float64,dim)
            g_[1] = G[4,1]
            g_[2] = G[4,2]
            g_[3] = G[4,3]
            if nargout>2
                H_      = zeros(Float64,3,3)
                H_[1,1] = H[4,1]
                H_[2,1] = H[4,2]
                H_[1,2] = H_[2,1]
                H_[3,1] = H[4,3]
                H_[1,3] = H_[3,1]
                H_[2,2] = H[4,4]
                H_[3,2] = H[4,5]
                H_[2,3] = H_[3,2]
                H_[3,3] = H[4,6]
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eY5"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        Y,G,H   = pbm.call("extfunc",EV_[1],EV_[2],EV_[3])
        f_      = Y[5]
        if nargout>1
            dim =   length(EV_)
            g_    = zeros(Float64,dim)
            g_[1] = G[5,1]
            g_[2] = G[5,2]
            g_[3] = G[5,3]
            if nargout>2
                H_      = zeros(Float64,3,3)
                H_[1,1] = H[5,1]
                H_[2,1] = H[5,2]
                H_[1,2] = H_[2,1]
                H_[3,1] = H[5,3]
                H_[1,3] = H_[3,1]
                H_[2,2] = H[5,4]
                H_[3,2] = H[5,5]
                H_[2,3] = H_[3,2]
                H_[3,3] = H[5,6]
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eY6"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        Y,G,H   = pbm.call("extfunc",EV_[1],EV_[2],EV_[3])
        f_      = Y[6]
        if nargout>1
            dim   = length(EV_)
            g_    = zeros(Float64,dim)
            g_[1] = G[6,1]
            g_[2] = G[6,2]
            g_[3] = 0.0
            if nargout>2
                H_      = zeros(Float64,3,3)
                H_[1,1] = H[6,1]
                H_[2,1] = H[6,2]
                H_[1,2] = H_[2,1]
                H_[2,2] = H[6,4]
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eY7"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        Y,G,H   = pbm.call("extfunc",EV_[1],EV_[2],EV_[3])
        f_      = Y[7]
        if nargout>1
            dim   = length(EV_)
            g_    = zeros(Float64,dim)
            g_[1] = G[7,1]
            g_[2] = G[7,2]
            g_[3] = G[7,3]
            if nargout>2
                H_      = zeros(Float64,3,3)
                H_[1,1] = H[7,1]
                H_[2,1] = H[7,2]
                H_[1,2] = H_[2,1]
                H_[3,1] = H[7,3]
                H_[1,3] = H_[3,1]
                H_[2,2] = H[7,4]
                H_[3,2] = H[7,5]
                H_[2,3] = H_[3,2]
                H_[3,3] = H[7,6]
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eY8"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        Y,G,H   = pbm.call("extfunc",EV_[1],EV_[2],EV_[3])
        f_      = Y[8]
        if nargout>1
            dim   = length(EV_)
            g_    = zeros(Float64,dim)
            g_[1] = G[8,1]
            g_[2] = G[8,2]
            g_[3] = G[8,3]
            if nargout>2
                H_      = zeros(Float64,3,3)
                H_[1,1] = H[8,1]
                H_[2,1] = H[8,2]
                H_[1,2] = H_[2,1]
                H_[3,1] = H[8,3]
                H_[1,3] = H_[3,1]
                H_[2,2] = H[8,4]
                H_[3,2] = H[8,5]
                H_[2,3] = H_[3,2]
                H_[3,3] = H[8,6]
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "extfunc"
    
    #   PAGE 129 Hock and Schittkowski.

        X1 = args[1]
        X2 = args[2]
        X3 = args[3]
        
        Y = zeros(Float64,8,1);
        G = zeros(Float64,8,3);
        H = zeros(Float64,8,6);

        #  First approximation to Y2.

        Y[2]   = 1.6 * X1

        #  Loop until Y2 converges.

        for iloop = 1:1000

        #  Y3.

           Y[3]   = 1.22 * Y[2] - X1
           G[3,1] = 1.22 * G[2,1] - 1.0
           G[3,2] = 1.22 * G[2,2]
           G[3,3] = 1.22 * G[2,3]
           H[3,1] = 1.22 * H[2,1]
           H[3,2] = 1.22 * H[2,2]
           H[3,3] = 1.22 * H[2,3]
           H[3,4] = 1.22 * H[2,4]
           H[3,5] = 1.22 * H[2,5]
           H[3,6] = 1.22 * H[2,6]

        #  Y6.

           Y[6]   = ( X2 + Y[3] ) / X1
           G[6,1] = - Y[6] / X1 + G[3,1] / X1
           G[6,2] =   1.0  / X1 + G[3,2] / X1
           G[6,3] =  G[3,3] / X1
           H[6,1] = - G[6,1] / X1 + Y[6]/ X1^2 - G[3,1] / X1^2 + H[3,1] / X1
           H[6,2] = - G[6,2] / X1 + H[3,2] / X1
           H[6,3] = - G[6,3] / X1 + H[3,3] / X1
           H[6,4] = H[3,4] / X1
           H[6,5] = H[3,5] / X1
           H[6,6] = H[3,6] / X1
           Y2C    = 0.01 * X1 * ( 112.0 + 13.167 * Y[6] - 0.6667 * Y[6]^2 )
           
        #  Y2.

           if abs( Y2C - Y[2] ) > 0.001 
           
              Y[2]   = Y2C
              G[2,1] = ( 0.01 * ( 112.0 + 13.167 * Y[6] - 0.6667 * Y[6]^2 )  
                           +  X1 * 0.13167 * G[6,1] -  X1 * 0.013334 * Y[6] * G[6,1] )
              G[2,2] = X1 * ( 0.13167 * G[6,2] - 0.013334 * Y[6] * G[6,2] )
              G[2,3] = X1 * ( 0.13167 * G[6,3] - 0.013334 * Y[6] * G[6,3] )
              H[2,1] = (  0.13167 * G[6,1] - 0.013334 * Y[6] * G[6,1] + 0.13167 * G[6,1]                
                       - 0.013334 * Y[6] * G[6,1] + X1 * 0.13167 * H[6,1]    - X1 * 0.013334 * G[6,1]^2 
                       - X1 * 0.013334 * Y[6] * H[6,1] )
              H[2,2] = (  0.13167 * G[6,2] - 0.013334 * Y[6] * G[6,2]         + X1 * 0.13167 * H[6,2] 
                       - X1 * 0.013334 * G[6,2] * G[6,1]                     - X1 * 0.013334 * Y[6] * H[6,2] )
              H[2,3] = ( 0.13167 * G[6,3] - 0.013334 * Y[6] *  G[6,3]          + X1 * 0.13167 * H[6,3] 
                       - X1 * 0.013334 * G[6,3] * G[6,1]                     - X1 * 0.013334 * Y[6] * H[6,3] )
              H[2,4] = X1 * ( 0.13167 * H[6,4] - 0.013334 * G[6,2]^2         - 0.013334 * Y[6] * H[6,4] )
              H[2,5] = X1 * ( 0.13167 * H[6,5] - 0.013334 * G[6,3] * G[6,2]  - 0.013334 * Y[6] * H[6,5] )
              H[2,6] = X1 * ( 0.13167 * H[6,6] - 0.013334 * G[6,3]^2         - 0.013334 * Y[6] * H[6,6] )
           else
              break
           end
        end

        #Y#D
        #G#D
        #H#D
        #return

        #  First approximation to Y4.

        Y[4]    = 93.0

        #  Loop until Y4 converges.

        for iloop = 1:1000

        #  Y5.

           Y[5]   = 86.35 + 1.098 * Y[6] - 0.038 * Y[6]^2 + 0.325 * ( Y[4] - 89.0 )
           G[5,1] = 1.098 * G[6,1] - 0.076 * Y[6] * G[6,1] + 0.325 * G[4,1]
           G[5,2] = 1.098 * G[6,2] - 0.076 * Y[6] * G[6,2] + 0.325 * G[4,2]
           G[5,3] = 1.098 * G[6,3] - 0.076 * Y[6] * G[6,3] + 0.325 * G[4,3]
           H[5,1] = 1.098 * H[6,1] - 0.076 * G[6,1] * G[6,1] - 0.076 * Y[6] * H[6,1] + 0.325 * H[4,1]
           H[5,2] = 1.098 * H[6,2] - 0.076 * G[6,1] * G[6,2] - 0.076 * Y[6] * H[6,2] + 0.325 * H[4,2]
           H[5,3] = 1.098 * H[6,3] - 0.076 * G[6,1] * G[6,3] - 0.076 * Y[6] * H[6,3] + 0.325 * H[4,3]
           H[5,4] = 1.098 * H[6,4] - 0.076 * G[6,2] * G[6,2] - 0.076 * Y[6] * H[6,4] + 0.325 * H[4,4]
           H[5,5] = 1.098 * H[6,5] - 0.076 * G[6,2] * G[6,3] - 0.076 * Y[6] * H[6,5] + 0.325 * H[4,5]
           H[5,6] = 1.098 * H[6,6] - 0.076 * G[6,3] * G[6,3] - 0.076 * Y[6] * H[6,6] + 0.325 * H[4,6]

        #  Y8.

           Y[8]   = 3.0 * Y[5] - 133.0
           G[8,1] = 3.0 * G[5,1]
           G[8,2] = 3.0 * G[5,2]
           G[8,3] = 3.0 * G[5,3]
           H[8,1] = 3.0 * H[5,1]
           H[8,2] = 3.0 * H[5,2]
           H[8,3] = 3.0 * H[5,3]
           H[8,4] = 3.0 * H[5,4]
           H[8,5] = 3.0 * H[5,5]
           H[8,6] = 3.0 * H[5,6]

        #  Y7.

           Y[7]   = 35.82 - 0.222 * Y[8]
           G[7,1] = - 0.222 * G[8,1]
           G[7,2] = - 0.222 * G[8,2]
           G[7,3] = - 0.222 * G[8,3]
           H[7,1] = - 0.222 * H[8,1]
           H[7,2] = - 0.222 * H[8,2]
           H[7,3] = - 0.222 * H[8,3]
           H[7,4] = - 0.222 * H[8,4]
           H[7,5] = - 0.222 * H[8,5]
           H[7,6] = - 0.222 * H[8,6]
           Y2Y7X3 = Y[2] * Y[7] + 1000.0 * X3
           Y4C    = 98000.0 * X3 / Y2Y7X3


        #H#D
        #'===1==='
        #return

           
        #  Y4.

           if abs( Y4C - Y[4] ) > 0.001
              Y[4]    = Y4C
              G[4,1] =   -  98000.0 * X3 * ( G[2,1] * Y[7] + Y[2] * G[7,1] ) / Y2Y7X3^2
              G[4,2] =   -  98000.0 * X3 * ( G[2,2] * Y[7] + Y[2] * G[7,2] ) / Y2Y7X3^2
              G[4,3] =      98000.0 / Y2Y7X3  - 98000.0 * X3 * ( G[2,3] * Y[7] + Y[2] * G[7,3] + 1000.0 ) / Y2Y7X3^2
              H[4,1] = ( -  98000.0 * X3 * ( H[2,1] * Y[7] + 2.0 * G[2,1] * G[7,1] + Y[2] * H[7,1] ) / Y2Y7X3^2
                         + 196000.0 * X3 * ( G[2,1] * Y[7] + Y[2] * G[7,1] )^2 / Y2Y7X3^3 )
              H[4,2] = ( -  98000.0 * X3 * ( H[2,2] * Y[7] + G[2,2] * G[7,1] + G[2,1] * G[7,2]  + Y[2] * H[7,2] ) / Y2Y7X3^2 
                       + 196000.0 * X3 * ( G[2,2] * Y[7] + Y[2] * G[7,2] ) * ( G[2,1] * Y[7]  + Y[2] * G[7,1] ) / Y2Y7X3^3 )
              H[4,3] = ( -  98000.0 * ( Y[2] * G[7,1] + Y[7] *  G[2,1] ) / Y2Y7X3^2 
                         -  98000.0 * X3 * ( G[2,3] * G[7,1] + Y[2] * H[7,3] + G[2,1] * G[7,3] + H[2,3] * Y[7] ) / Y2Y7X3^2 
                         + 196000.0 * X3 * ( Y[2] * G[7,1] + Y[7] * G[2,1] ) *
                            ( G[2,3 ] * Y[7] + Y[2] * G[7,3 ] + 1000.0 ) / Y2Y7X3^3 )
              H[4,4] = ( -  98000.0 * X3 * ( H[2,4] * Y[7] + 2.0 * G[2,2] * G[7,2] + Y[2] * H[7,4] ) / Y2Y7X3^2
                         + 196000.0 * X3 * ( G[2,2] * Y[7] + Y[2] * G[7,2] )^2 / Y2Y7X3^3 )
              H[4,5] = ( -  98000.0 * ( Y[2] * G[7,2] + Y[7] *  G[2,2] ) / Y2Y7X3^2 
                         -  98000.0 * X3 * ( G[2,3] * G[7,2] + Y[2] * H[ 7, 5 ]       
                         + G[2,2] * G[7,3] + H[2,5] * Y[7] ) / Y2Y7X3^2  + 196000.0 * X3 * ( Y[2] * G[7,2]
                         + Y[7] * G[2,2] ) * ( G[2,3] * Y[7] + Y[2] *  G[7,3] + 1000.0 ) / Y2Y7X3^3 )
              H[4,6] = ( - 196000.0 * ( Y[2] * G[7,3] + Y[7] * G[2,3] + 1000.0 ) / Y2Y7X3^2  
                         -  98000.0 * X3 * ( H[2,6] * Y[7] + 2.0 * G[2,3] * G[7,3] + Y[2] * H[7,6] ) / Y2Y7X3^2
                         + 196000.0 * X3 * ( G[2,3] * Y[7] + Y[2] * G[7,3] + 1000.0 )^2 / Y2Y7X3^3 )
           else
              break
           end
        end


#Y#D
#G#D
#H#D
#'============='

        return Y, G, H

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

