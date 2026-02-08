function HS88(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS88
#    *********
# 
#    A time-optimal heat conduction problem.
# 
#    Source: problem 88 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: Nick Gould, September 1991.
#        Julia coding: Cunxin Huang, 2025.
# 
#    classification = "C-CQOR2-MN-6-1"
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "HS88"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["N"] = 6
        v_["EPS"] = 0.01
        v_["EPSSQR"] = v_["EPS"]*v_["EPS"]
        v_["-EPSSQR"] = -1.0*v_["EPSSQR"]
        v_["1"] = 1
        v_["2"] = 2
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        irA   = Int64[]
        icA   = Int64[]
        valA  = Float64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        ig,ig_,_ = s2mpj_ii("CON",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"CON")
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
        pbm.gconst[ig_["CON"]] = Float64(v_["-EPSSQR"])
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        for I = Int64(v_["1"]):Int64(v_["2"]):Int64(v_["N"])
            if haskey(ix_,"X"*string(I))
                pb.x0[ix_["X"*string(I)]] = Float64(0.5)
            else
                pb.y0[findfirst(x->x==ig_["X"*string(I)],pbm.congrps)] = Float64(0.5)
            end
        end
        for I = Int64(v_["2"]):Int64(v_["2"]):Int64(v_["N"])
            if haskey(ix_,"X"*string(I))
                pb.x0[ix_["X"*string(I)]] = Float64(-0.5)
            else
                pb.y0[findfirst(x->x==ig_["X"*string(I)],pbm.congrps)] = Float64(-0.5)
            end
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eSQR", iet_)
        loaset(elftv,it,1,"X")
        it,iet_,_ = s2mpj_ii( "eH", iet_)
        loaset(elftv,it,1,"X1")
        loaset(elftv,it,2,"X2")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["N"])
            ename = "O"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
            vname = "X"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        ename = "H"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eH")
        arrset(ielftype,ie,iet_["eH"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            ig = ig_["OBJ"]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["O"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        ig = ig_["CON"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["H"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+pb.neq+1:pb.m] = zeros(Float64,pb.nge)
        pb.cupper[1:pb.nge] = fill(Inf,pb.nge)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-CQOR2-MN-2-1"
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

    elseif action == "eSQR"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]*EV_[1]
        if nargout>1
            dim = length(EV_)
            g_  = zeros(Float64,dim)
            g_[1] = EV_[1]+EV_[1]
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 2.0e+0
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eH"

        nargout  = args[3]
        pbm      = args[4]
        f_,g_,H_ = pbm.call("extfunc", args[1] )
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "extfunc"

        x     = args[1]
        n     = length(x);
        g     = zeros(Float64,n,1);
        H     = zeros(Float64,n,n);
        A     = zeros(Float64,30,1);
        R     = zeros(Float64,30,30);
        S     = zeros(Float64,30,1);
        RHO   = zeros(Float64,30,1);
        DRHO  = zeros(Float64,30,n);
        D2RHO = zeros(Float64,30,n,n);
        P     = zeros(Float64,n+1,1);
        
        mu = [ 8.6033358901938017e-01,  3.4256184594817283e+00, 6.4372981791719468e+00,  9.5293344053619631e+00,
               1.2645287223856643e+01,  1.5771284874815882e+01, 1.8902409956860023e+01,  2.2036496727938566e+01,
               2.5172446326646664e+01,  2.8309642854452012e+01, 3.1447714637546234e+01,  3.4586424215288922e+01,
               3.7725612827776501e+01,  4.0865170330488070e+01, 4.4005017920830845e+01,  4.7145097736761031e+01,
               5.0285366337773652e+01,  5.3425790477394663e+01, 5.6566344279821521e+01,  5.9707007305335459e+01,
               6.2847763194454451e+01,  6.5988598698490392e+01, 6.9129502973895256e+01,  7.2270467060308960e+01,
               7.5411483488848148e+01,  7.8552545984242926e+01, 8.1693649235601683e+01,  8.4834788718042290e+01,
               8.7975960552493220e+01,  9.1117161394464745e+01 ];

        T = 2.0 / 15.0;
        for i = 1:30
            MUI  = mu[i]
            SMUI = sin(MUI)
            CMUI = cos(MUI)
            AI   = 2.0*SMUI/(MUI+SMUI*CMUI)
            A[i] = AI
            S[i] = 2.0*AI*(CMUI-SMUI/MUI)
            AIMUI2 = AI * MUI^2
            for j = 1:i
                if ( i == j ) 
                    R[i,i] = 0.5*(1.0+0.5*sin(MUI+MUI)/MUI)*AIMUI2^2;
                else
                    MUJ    = mu[j]
                    R[i,j] = 0.5*(sin(MUI+MUJ )/(MUI+MUJ)+sin(MUI-MUJ )/(MUI-MUJ))*AIMUI2*A[j]*MUJ^2
                    R[j,i] = R[i,j];
                end
            end
        end
        
    #                                  n   2
    #  Calculate the functions p(x) = SUM x .
    #                           j     i=j  i

        for k = n:-1:1
            P[k] = P[k+1]+x[k]^2;
        end

    #  Calculate the functions rho.

        for j = 1:30
            MUJ2 = mu[j]*mu[j]
            U    = exp(-MUJ2*P[1])
            for k =1:n
                DRHO[j,k] = 2.0*U*x[k]
                for l = k:n
                    D2RHO[j,k,l] = -4.0*MUJ2*U*x[k]*x[l]
                    if ( l == k )
                        D2RHO[j,k,l] = D2RHO[j,k,l]+2.0*U
                    end
                end
            end
            ALPHA = -2.0
            for i = 2:n
                EU = ALPHA*exp(-MUJ2*P[i])
                U  = U+EU
                for k = i:n
                    DRHO[j,k] = DRHO[j,k]+2.0*EU*x[k]
                    for l = k:n
                        D2RHO[j,k,l] = D2RHO[j,k,l]-4.0*MUJ2*EU*x[k]*x[l]
                        if ( l == k )
                             D2RHO[j,k,l] = D2RHO[j,k,l]+2.0*EU
                        end
                    end
                end
                ALPHA = - ALPHA
            end
            U      = U+0.5*ALPHA
            RHO[j] = -U/MUJ2;
        end

    #  Evaluate the function and derivatives.

        f = T;
        for i = 1:30
            SI   = S[i]
            RHOI = RHO[i]
            f    = f+SI*RHOI;
            for k = 1:n
                g[k] = g[k]+SI*DRHO[i,k]
                for l = k:n
                   H[k,l] = H[k,l]+SI*D2RHO[i,k,l]
                end
            end
            for j = 1:30
                RIJ  = R[i,j]
                RHOJ = RHO[j]
                f    = f+RIJ*RHOI*RHOJ;
                for k = 1:n
                    g[k] = g[k]+RIJ*(RHOI*DRHO[j,k]+RHOJ*DRHO[i,k])
                    for l = k:n
                        H[k,l] = H[k,l]+RIJ*(RHOI*D2RHO[j,k,l]+RHOJ*D2RHO[i,k,l]+DRHO[i,k]*DRHO[j,l]+DRHO[j,k]*DRHO[i,l])
                    end
                end
            end
        end

    #   Symmetrize the Hessian.

        for k = 1:n
            for l = k+1:n
                H[l,k] = H[k,l]
            end
        end
        
        return f, g, H

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

