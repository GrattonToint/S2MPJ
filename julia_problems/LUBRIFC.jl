function LUBRIFC(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : LUBRIFC
#    *********
# 
#    Corrected version of LUBRIF which contained an error
#    in the definition of the Reynold's equation (ELEMENT USES)
#    mixing H & P, see line 298ff below (or search for ***).
#    Fix by: Sven Leyffer, U. Dundee, September 2000
# 
#    The elastodynamic lubrification problem by Kostreva.
# 
#    Source:
#    M.M. Kostreva,
#    "Elasto-hydrodynamic lubrification: a non-linear
#    complementarity problem",
#    International Journal for Numerical Methods in Fluids,
#    4: 377-397, 1984.
# 
#    This problem is problem #5 in More's test set.
# 
#    SIF input: Ph. Toint, June 1990.
# 
#    classification = "C-QOR2-MN-V-V"
# 
#    Number of discretized points per unit length
# 
#       Alternative values for the SIF file parameters:
# IE NN                  10             $-PARAMETER n = 151    original value
# IE NN                  50             $-PARAMETER n = 751
# IE NN                  250            $-PARAMETER n = 3751
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "LUBRIFC"

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
            v_["NN"] = Int64(5);  #  SIF file default value
        else
            v_["NN"] = Int64(args[1]);
        end
        v_["ALPHA"] = 1.838
        v_["LAMBDA"] = 1.642
        v_["XA"] = -3.0
        v_["XF"] = 2.0
        v_["N"] = 5*v_["NN"]
        v_["0"] = 0
        v_["1"] = 1
        v_["2"] = 2
        v_["3"] = 3
        v_["4"] = 4
        v_["PI"] = 3.1415926535
        v_["2N"] = 2*v_["N"]
        v_["2N-2"] = -2+v_["2N"]
        v_["2N-1"] = -1+v_["2N"]
        v_["2N+2"] = 2+v_["2N"]
        v_["-XA"] = -1.0*v_["XA"]
        v_["LEN"] = v_["XF"]+v_["-XA"]
        v_["1/PI"] = 1.0/v_["PI"]
        v_["1/2PI"] = 0.5*v_["1/PI"]
        v_["RN"] = Float64(v_["N"])
        v_["1/N"] = 1.0/v_["RN"]
        v_["DX"] = v_["LEN"]*v_["1/N"]
        v_["1/DX"] = 1.0/v_["DX"]
        v_["L/DX"] = v_["LAMBDA"]*v_["1/DX"]
        v_["-L/DX"] = -1.0*v_["L/DX"]
        v_["1/DX2"] = v_["1/DX"]*v_["1/DX"]
        v_["-1/DX2"] = -1.0*v_["1/DX2"]
        v_["DX/PI"] = v_["DX"]*v_["1/PI"]
        v_["2DX/PI"] = 2.0*v_["DX/PI"]
        v_["DX/2"] = 0.5*v_["DX"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        iv,ix_,_ = s2mpj_ii("K",ix_)
        arrset(pb.xnames,iv,"K")
        for I = Int64(v_["0"]):Int64(v_["2"]):Int64(v_["2N"])
            iv,ix_,_ = s2mpj_ii("P"*string(I),ix_)
            arrset(pb.xnames,iv,"P"*string(I))
        end
        for J = Int64(v_["1"]):Int64(v_["2"]):Int64(v_["2N-1"])
            iv,ix_,_ = s2mpj_ii("H"*string(J),ix_)
            arrset(pb.xnames,iv,"H"*string(J))
        end
        for I = Int64(v_["2"]):Int64(v_["2"]):Int64(v_["2N-2"])
            iv,ix_,_ = s2mpj_ii("R"*string(I),ix_)
            arrset(pb.xnames,iv,"R"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["2"]):Int64(v_["2"]):Int64(v_["2N-2"])
            ig,ig_,_ = s2mpj_ii("R"*string(Int64(v_["0"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"R"*string(Int64(v_["0"])))
            iv = ix_["P"*string(I)]
            pbm.A[ig,iv] += Float64(v_["2DX/PI"])
        end
        ig,ig_,_ = s2mpj_ii("R"*string(Int64(v_["0"])),ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R"*string(Int64(v_["0"])))
        iv = ix_["P"*string(Int64(v_["2N"]))]
        pbm.A[ig,iv] += Float64(v_["DX/PI"])
        ig,ig_,_ = s2mpj_ii("COMPL",ig_)
        arrset(gtype,ig,"<>")
        for I = Int64(v_["2"]):Int64(v_["2"]):Int64(v_["2N-2"])
            v_["I-1"] = -1+I
            v_["I+1"] = 1+I
            ig,ig_,_ = s2mpj_ii("DR"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"DR"*string(I))
            iv = ix_["H"*string(Int64(v_["I+1"]))]
            pbm.A[ig,iv] += Float64(v_["L/DX"])
            iv = ix_["H"*string(Int64(v_["I-1"]))]
            pbm.A[ig,iv] += Float64(v_["-L/DX"])
            iv = ix_["R"*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
        end
        for J = Int64(v_["1"]):Int64(v_["2"]):Int64(v_["2N-1"])
            v_["-J"] = -1*J
            ig,ig_,_ = s2mpj_ii("DH"*string(J),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"DH"*string(J))
            iv = ix_["K"]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["H"*string(J)]
            pbm.A[ig,iv] += Float64(-1.0)
            for I = Int64(v_["2"]):Int64(v_["2N"])
                v_["C"*string(I)] = 0.0
            end
            v_["RI-J"] = Float64(v_["-J"])
            v_["I-JDX"] = v_["RI-J"]*v_["DX/2"]
            v_["ALN"] = abs(v_["I-JDX"])
            v_["LN"] = log(v_["ALN"])
            v_["T1"] = v_["I-JDX"]*v_["LN"]
            v_["COEFF"] = v_["T1"]*v_["1/2PI"]
            v_["C"*string(Int64(v_["2"]))] = v_["C"*string(Int64(v_["2"]))]+v_["COEFF"]
            v_["I-J"] = 2+v_["-J"]
            v_["RI-J"] = Float64(v_["I-J"])
            v_["I-JDX"] = v_["RI-J"]*v_["DX/2"]
            v_["ALN"] = abs(v_["I-JDX"])
            v_["LN"] = log(v_["ALN"])
            v_["T1"] = v_["I-JDX"]*v_["LN"]
            v_["COEFF"] = v_["T1"]*v_["1/PI"]
            v_["C"*string(Int64(v_["4"]))] = v_["C"*string(Int64(v_["4"]))]+v_["COEFF"]
            for I = Int64(v_["4"]):Int64(v_["2"]):Int64(v_["2N-2"])
                v_["I-2"] = -2+I
                v_["I+2"] = 2+I
                v_["I-J"] = I+v_["-J"]
                v_["RI-J"] = Float64(v_["I-J"])
                v_["I-JDX"] = v_["RI-J"]*v_["DX/2"]
                v_["ALN"] = abs(v_["I-JDX"])
                v_["LN"] = log(v_["ALN"])
                v_["T1"] = v_["I-JDX"]*v_["LN"]
                v_["COEFF"] = v_["T1"]*v_["1/PI"]
                v_["C"*string(Int64(v_["I+2"]))] = (v_["C"*string(Int64(v_["I+2"]))]+
                     v_["COEFF"])
                v_["-COEFF"] = -1.0*v_["COEFF"]
                v_["C"*string(Int64(v_["I-2"]))] = (v_["C"*string(Int64(v_["I-2"]))]+
                     v_["-COEFF"])
            end
            v_["I-J"] = v_["2N"]+v_["-J"]
            v_["RI-J"] = Float64(v_["I-J"])
            v_["I-JDX"] = v_["RI-J"]*v_["DX/2"]
            v_["ALN"] = abs(v_["I-JDX"])
            v_["LN"] = log(v_["ALN"])
            v_["T1"] = v_["I-JDX"]*v_["LN"]
            v_["COEFF"] = v_["T1"]*v_["1/2PI"]
            v_["-COEFF"] = -1.0*v_["COEFF"]
            v_["C"*string(Int64(v_["2N-2"]))] = (v_["C"*string(Int64(v_["2N-2"]))]+
                 v_["-COEFF"])
            for I = Int64(v_["2"]):Int64(v_["2"]):Int64(v_["2N-2"])
                ig,ig_,_ = s2mpj_ii("DH"*string(J),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"DH"*string(J))
                iv = ix_["P"*string(I)]
                pbm.A[ig,iv] += Float64(v_["C"*string(I)])
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
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        pbm.gconst[ig_["R"*string(Int64(v_["0"]))]] = Float64(1.0)
        for J = Int64(v_["1"]):Int64(v_["2"]):Int64(v_["2N-1"])
            v_["RJ"] = Float64(J)
            v_["JDX"] = v_["RJ"]*v_["DX/2"]
            v_["XJ"] = v_["XA"]+v_["JDX"]
            v_["XJSQ"] = v_["XJ"]*v_["XJ"]
            v_["XJSQ+1"] = 1.0+v_["XJSQ"]
            v_["RHS"] = -1.0*v_["XJSQ+1"]
            pbm.gconst[ig_["DH"*string(J)]] = Float64(v_["RHS"])
        end
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower[ix_["K"]] = -Inf
        pb.xupper[ix_["K"]] = +Inf
        for I = Int64(v_["2"]):Int64(v_["2"]):Int64(v_["2N-2"])
            pb.xupper[ix_["P"*string(I)]] = 3.0
            pb.xlower[ix_["P"*string(I)]] = 0.0
        end
        for I = Int64(v_["1"]):Int64(v_["2"]):Int64(v_["2N-1"])
            pb.xlower[ix_["H"*string(I)]] = -Inf
            pb.xupper[ix_["H"*string(I)]] = +Inf
        end
        pb.xlower[ix_["P"*string(Int64(v_["0"]))]] = 0.0
        pb.xupper[ix_["P"*string(Int64(v_["0"]))]] = 0.0
        pb.xlower[ix_["P"*string(Int64(v_["2N"]))]] = 0.0
        pb.xupper[ix_["P"*string(Int64(v_["2N"]))]] = 0.0
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(0.0),pb.n)
        v_["2NN"] = v_["NN"]+v_["NN"]
        v_["4NN"] = 4*v_["NN"]
        for I = Int64(v_["2"]):Int64(v_["2"]):Int64(v_["4NN"])
            v_["RI"] = Float64(I)
            v_["IDX"] = v_["RI"]*v_["DX/2"]
            v_["XI"] = v_["XA"]+v_["IDX"]
            v_["LIN"] = 0.02*v_["XI"]
            v_["PI0"] = 0.06+v_["LIN"]
            if haskey(ix_,"P"*string(I))
                pb.x0[ix_["P"*string(I)]] = Float64(v_["PI0"])
            else
                pb.y0[findfirst(x->x==ig_["P"*string(I)],pbm.congrps)] = Float64(v_["PI0"])
            end
        end
        v_["4NN+2"] = 2+v_["4NN"]
        v_["8NN"] = 8*v_["NN"]
        for I = Int64(v_["4NN+2"]):Int64(v_["2"]):Int64(v_["8NN"])
            v_["RI"] = Float64(I)
            v_["IDX"] = v_["RI"]*v_["DX/2"]
            v_["XI"] = v_["XA"]+v_["IDX"]
            v_["XISQ"] = v_["XI"]*v_["XI"]
            v_["-XISQ"] = -1.0*v_["XISQ"]
            v_["1-XISQ"] = 1.0+v_["-XISQ"]
            v_["PI0"] = sqrt(v_["1-XISQ"])
            if haskey(ix_,"P"*string(I))
                pb.x0[ix_["P"*string(I)]] = Float64(v_["PI0"])
            else
                pb.y0[findfirst(x->x==ig_["P"*string(I)],pbm.congrps)] = Float64(v_["PI0"])
            end
        end
        v_["8NN+2"] = 2+v_["8NN"]
        for I = Int64(v_["8NN+2"]):Int64(v_["2"]):Int64(v_["2N"])
            if haskey(ix_,"P"*string(I))
                pb.x0[ix_["P"*string(I)]] = Float64(0.0)
            else
                pb.y0[findfirst(x->x==ig_["P"*string(I)],pbm.congrps)] = Float64(0.0)
            end
        end
        for J = Int64(v_["1"]):Int64(v_["2"]):Int64(v_["2N-1"])
            v_["RJ"] = Float64(J)
            v_["JDX"] = v_["RJ"]*v_["DX/2"]
            v_["XJ"] = v_["XA"]+v_["JDX"]
            v_["XJSQ"] = v_["XJ"]*v_["XJ"]
            if haskey(ix_,"H"*string(J))
                pb.x0[ix_["H"*string(J)]] = Float64(v_["XJSQ"])
            else
                pb.y0[findfirst(x->x==ig_["H"*string(J)],pbm.congrps)] = Float64(v_["XJSQ"])
            end
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eREY", iet_)
        loaset(elftv,it,1,"PA")
        loaset(elftv,it,2,"PB")
        loaset(elftv,it,3,"H")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"A")
        it,iet_,_ = s2mpj_ii( "en2PR", iet_)
        loaset(elftv,it,1,"P")
        loaset(elftv,it,2,"R")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for J = Int64(v_["1"]):Int64(v_["2"]):Int64(v_["2N-1"])
            v_["I+"] = 1+J
            v_["I-"] = -1+J
            ename = "ER"*string(J)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eREY")
            arrset(ielftype,ie,iet_["eREY"])
            vname = "P"*string(Int64(v_["I-"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
            posev = findfirst(x->x=="PA",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "H"*string(J)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
            posev = findfirst(x->x=="H",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "P"*string(Int64(v_["I+"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
            posev = findfirst(x->x=="PB",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="A",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["ALPHA"]))
        end
        for I = Int64(v_["2"]):Int64(v_["2"]):Int64(v_["2N-2"])
            ename = "EC"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
            vname = "P"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
            posev = findfirst(x->x=="P",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "R"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
            posev = findfirst(x->x=="R",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["2"]):Int64(v_["2"]):Int64(v_["2N-2"])
            ig = ig_["COMPL"]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["EC"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        for I = Int64(v_["2"]):Int64(v_["2"]):Int64(v_["2N-2"])
            v_["I-1"] = -1+I
            v_["I+1"] = 1+I
            ig = ig_["DR"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["ER"*string(Int64(v_["I-1"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["1/DX2"]))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["ER"*string(Int64(v_["I+1"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["-1/DX2"]))
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN                0.0
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
        pb.pbclass = "C-QOR2-MN-V-V"
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

    elseif action == "eREY"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        HA = -0.5*pbm.elpar[iel_][1]
        EARG = HA*(EV_[1]+EV_[2])
        E = exp(EARG)
        PAMPB = EV_[1]-EV_[2]
        T1 = PAMPB*HA+1.0
        T2 = PAMPB*HA-1.0
        HSQ = EV_[3]*EV_[3]
        HCB = HSQ*EV_[3]
        f_   = PAMPB*HCB*E
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = T1*HCB*E
            g_[2] = T2*HCB*E
            g_[3] = 3.0*PAMPB*HSQ*E
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,1] = HCB*E*HA*(T1+1.0)
                H_[1,2] = HCB*E*HA*(T1-1.0)
                H_[2,1] = H_[1,2]
                H_[1,3] = 3.0*T1*HSQ*E
                H_[3,1] = H_[1,3]
                H_[2,2] = HCB*E*HA*(T2-1.0)
                H_[2,3] = 3.0*T2*HSQ*E
                H_[3,2] = H_[2,3]
                H_[3,3] = 6.0*EV_[3]*PAMPB*E
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

