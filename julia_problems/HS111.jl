function HS111(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS111
#    *********
# 
#    This problem is a chemical equilibrium problem involving 3 linear
#    equality constraints.
# 
#    Source: problem 111 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: Nick Gould, August 1991.
# 
#    classification = "C-OOR2-AN-10-3"
# 
#    N is the number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "HS111"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["N"] = 10
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
        v_["C1"] = -6.089
        v_["C2"] = -17.164
        v_["C3"] = -34.054
        v_["C4"] = -5.914
        v_["C5"] = -24.721
        v_["C6"] = -14.986
        v_["C7"] = -24.100
        v_["C8"] = -10.708
        v_["C9"] = -26.662
        v_["C10"] = -22.179
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        ig,ig_,_ = s2mpj_ii("CON1",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"CON1")
        ig,ig_,_ = s2mpj_ii("CON2",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"CON2")
        ig,ig_,_ = s2mpj_ii("CON3",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"CON3")
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
        pbm.gconst[ig_["CON1"]] = Float64(2.0)
        pbm.gconst[ig_["CON2"]] = Float64(1.0)
        pbm.gconst[ig_["CON3"]] = Float64(1.0)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = fill(-100.0,pb.n)
        pb.xupper = fill(100.0,pb.n)
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(-2.3),pb.n)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eOBJ", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        loaset(elftv,it,3,"V3")
        loaset(elftv,it,4,"V4")
        loaset(elftv,it,5,"V5")
        loaset(elftv,it,6,"V6")
        loaset(elftv,it,7,"V7")
        loaset(elftv,it,8,"V8")
        loaset(elftv,it,9,"V9")
        loaset(elftv,it,10,"V10")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"C")
        it,iet_,_ = s2mpj_ii( "eEXP", iet_)
        loaset(elftv,it,1,"X")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["N"])
            v_["RI"*string(I)] = Float64(I)
            v_["RI"*string(I)] = 0.1+v_["RI"*string(I)]
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            ename = "O"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eOBJ")
            arrset(ielftype,ie,iet_["eOBJ"])
            v_["TEMP"] = v_["RI"*string(Int64(v_["1"]))]
            v_["RI"*string(Int64(v_["1"]))] = v_["RI"*string(I)]
            v_["RI"*string(I)] = v_["TEMP"]
            v_["R"] = v_["RI"*string(Int64(v_["1"]))]
            v_["J"] = trunc(Int,v_["R"])
            vname = "X"*string(Int64(v_["J"]))
            iv,ix_,pb  = (
                  s2mpj_nlx(vname,ix_,pb,1,Float64(-100.0),Float64(100.0),Float64(-2.3)))
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            v_["R"] = v_["RI"*string(Int64(v_["2"]))]
            v_["J"] = trunc(Int,v_["R"])
            vname = "X"*string(Int64(v_["J"]))
            iv,ix_,pb  = (
                  s2mpj_nlx(vname,ix_,pb,1,Float64(-100.0),Float64(100.0),Float64(-2.3)))
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            v_["R"] = v_["RI"*string(Int64(v_["3"]))]
            v_["J"] = trunc(Int,v_["R"])
            vname = "X"*string(Int64(v_["J"]))
            iv,ix_,pb  = (
                  s2mpj_nlx(vname,ix_,pb,1,Float64(-100.0),Float64(100.0),Float64(-2.3)))
            posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            v_["R"] = v_["RI"*string(Int64(v_["4"]))]
            v_["J"] = trunc(Int,v_["R"])
            vname = "X"*string(Int64(v_["J"]))
            iv,ix_,pb  = (
                  s2mpj_nlx(vname,ix_,pb,1,Float64(-100.0),Float64(100.0),Float64(-2.3)))
            posev = findfirst(x->x=="V4",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            v_["R"] = v_["RI"*string(Int64(v_["5"]))]
            v_["J"] = trunc(Int,v_["R"])
            vname = "X"*string(Int64(v_["J"]))
            iv,ix_,pb  = (
                  s2mpj_nlx(vname,ix_,pb,1,Float64(-100.0),Float64(100.0),Float64(-2.3)))
            posev = findfirst(x->x=="V5",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            v_["R"] = v_["RI"*string(Int64(v_["6"]))]
            v_["J"] = trunc(Int,v_["R"])
            vname = "X"*string(Int64(v_["J"]))
            iv,ix_,pb  = (
                  s2mpj_nlx(vname,ix_,pb,1,Float64(-100.0),Float64(100.0),Float64(-2.3)))
            posev = findfirst(x->x=="V6",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            v_["R"] = v_["RI"*string(Int64(v_["7"]))]
            v_["J"] = trunc(Int,v_["R"])
            vname = "X"*string(Int64(v_["J"]))
            iv,ix_,pb  = (
                  s2mpj_nlx(vname,ix_,pb,1,Float64(-100.0),Float64(100.0),Float64(-2.3)))
            posev = findfirst(x->x=="V7",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            v_["R"] = v_["RI"*string(Int64(v_["8"]))]
            v_["J"] = trunc(Int,v_["R"])
            vname = "X"*string(Int64(v_["J"]))
            iv,ix_,pb  = (
                  s2mpj_nlx(vname,ix_,pb,1,Float64(-100.0),Float64(100.0),Float64(-2.3)))
            posev = findfirst(x->x=="V8",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            v_["R"] = v_["RI"*string(Int64(v_["9"]))]
            v_["J"] = trunc(Int,v_["R"])
            vname = "X"*string(Int64(v_["J"]))
            iv,ix_,pb  = (
                  s2mpj_nlx(vname,ix_,pb,1,Float64(-100.0),Float64(100.0),Float64(-2.3)))
            posev = findfirst(x->x=="V9",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            v_["R"] = v_["RI"*string(Int64(v_["10"]))]
            v_["J"] = trunc(Int,v_["R"])
            vname = "X"*string(Int64(v_["J"]))
            iv,ix_,pb  = (
                  s2mpj_nlx(vname,ix_,pb,1,Float64(-100.0),Float64(100.0),Float64(-2.3)))
            posev = findfirst(x->x=="V10",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="C",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["C"*string(I)]))
            ename = "E"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eEXP")
            arrset(ielftype,ie,iet_["eEXP"])
            vname = "X"*string(I)
            iv,ix_,pb  = (
                  s2mpj_nlx(vname,ix_,pb,1,Float64(-100.0),Float64(100.0),Float64(-2.3)))
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
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
        ig = ig_["CON1"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(1.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E2"])
        loaset(pbm.grelw,ig,posel,Float64(2.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E3"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(2.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E6"])
        loaset(pbm.grelw,ig,posel,Float64(1.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E10"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(1.0))
        ig = ig_["CON2"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E4"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(1.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E5"])
        loaset(pbm.grelw,ig,posel,Float64(2.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E6"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(1.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E7"])
        loaset(pbm.grelw,ig,posel,Float64(1.0))
        ig = ig_["CON3"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E3"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(1.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E7"])
        loaset(pbm.grelw,ig,posel,Float64(1.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E8"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(1.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E9"])
        loaset(pbm.grelw,ig,posel,Float64(2.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E10"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               -47.707579
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
        pb.pbclass = "C-OOR2-AN-10-3"
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

    elseif action == "eEXP"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        EX = exp(EV_[1])
        f_   = EX
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EX
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = EX
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eOBJ"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        E1 = exp(EV_[1])
        E2 = exp(EV_[2])
        E3 = exp(EV_[3])
        E4 = exp(EV_[4])
        E5 = exp(EV_[5])
        E6 = exp(EV_[6])
        E7 = exp(EV_[7])
        E8 = exp(EV_[8])
        E9 = exp(EV_[9])
        E10 = exp(EV_[10])
        SUM = E1+E2+E3+E4+E5+E6+E7+E8+E9+E10
        f_   = E1*(pbm.elpar[iel_][1]+EV_[1]-log(SUM))
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = E1*(pbm.elpar[iel_][1]+EV_[1]-log(SUM))+E1*(1.0e+0-E1/SUM)
            g_[2] = -E1*E2/SUM
            g_[3] = -E1*E3/SUM
            g_[4] = -E1*E4/SUM
            g_[5] = -E1*E5/SUM
            g_[6] = -E1*E6/SUM
            g_[7] = -E1*E7/SUM
            g_[8] = -E1*E8/SUM
            g_[9] = -E1*E9/SUM
            g_[10] = -E1*E10/SUM
            if nargout>2
                H_ = zeros(Float64,10,10)
                H_[1,1] = (E1*(pbm.elpar[iel_][1]+EV_[1]-log(SUM))+E1*(1.0e+0-E1/SUM)+
                     E1*(1.0e+0-E1/SUM)+E1*(-E1/SUM)+E1*(E1^2/SUM^2))
                H_[1,2] = (-1.0e+0+E1/SUM)*E1*E2/SUM
                H_[2,1] = H_[1,2]
                H_[2,2] = (-1.0e+0+E2/SUM)*E1*E2/SUM
                H_[1,3] = (-1.0e+0+E1/SUM)*E1*E3/SUM
                H_[3,1] = H_[1,3]
                H_[2,3] = E1*E2*E3/SUM^2
                H_[3,2] = H_[2,3]
                H_[3,3] = (-1.0e+0+E3/SUM)*E1*E3/SUM
                H_[1,4] = (-1.0e+0+E1/SUM)*E1*E4/SUM
                H_[4,1] = H_[1,4]
                H_[2,4] = E1*E2*E4/SUM^2
                H_[4,2] = H_[2,4]
                H_[3,4] = E1*E3*E4/SUM^2
                H_[4,3] = H_[3,4]
                H_[4,4] = (-1.0e+0+E4/SUM)*E1*E4/SUM
                H_[1,5] = (-1.0e+0+E1/SUM)*E1*E5/SUM
                H_[5,1] = H_[1,5]
                H_[2,5] = E1*E2*E5/SUM^2
                H_[5,2] = H_[2,5]
                H_[3,5] = E1*E3*E5/SUM^2
                H_[5,3] = H_[3,5]
                H_[4,5] = E1*E4*E5/SUM^2
                H_[5,4] = H_[4,5]
                H_[5,5] = (-1.0e+0+E5/SUM)*E1*E5/SUM
                H_[1,6] = (-1.0e+0+E1/SUM)*E1*E6/SUM
                H_[6,1] = H_[1,6]
                H_[2,6] = E1*E2*E6/SUM^2
                H_[6,2] = H_[2,6]
                H_[3,6] = E1*E3*E6/SUM^2
                H_[6,3] = H_[3,6]
                H_[4,6] = E1*E4*E6/SUM^2
                H_[6,4] = H_[4,6]
                H_[5,6] = E1*E5*E6/SUM^2
                H_[6,5] = H_[5,6]
                H_[6,6] = (-1.0e+0+E6/SUM)*E1*E6/SUM
                H_[1,7] = (-1.0e+0+E1/SUM)*E1*E7/SUM
                H_[7,1] = H_[1,7]
                H_[2,7] = E1*E2*E7/SUM^2
                H_[7,2] = H_[2,7]
                H_[3,7] = E1*E3*E7/SUM^2
                H_[7,3] = H_[3,7]
                H_[4,7] = E1*E4*E7/SUM^2
                H_[7,4] = H_[4,7]
                H_[5,7] = E1*E5*E7/SUM^2
                H_[7,5] = H_[5,7]
                H_[6,7] = E1*E6*E7/SUM^2
                H_[7,6] = H_[6,7]
                H_[7,7] = (-1.0e+0+E7/SUM)*E1*E7/SUM
                H_[1,8] = (-1.0e+0+E1/SUM)*E1*E8/SUM
                H_[8,1] = H_[1,8]
                H_[2,8] = E1*E2*E8/SUM^2
                H_[8,2] = H_[2,8]
                H_[3,8] = E1*E3*E8/SUM^2
                H_[8,3] = H_[3,8]
                H_[4,8] = E1*E4*E8/SUM^2
                H_[8,4] = H_[4,8]
                H_[5,8] = E1*E5*E8/SUM^2
                H_[8,5] = H_[5,8]
                H_[6,8] = E1*E6*E8/SUM^2
                H_[8,6] = H_[6,8]
                H_[7,8] = E1*E7*E8/SUM^2
                H_[8,7] = H_[7,8]
                H_[8,8] = (-1.0e+0+E8/SUM)*E1*E8/SUM
                H_[1,9] = (-1.0e+0+E1/SUM)*E1*E9/SUM
                H_[9,1] = H_[1,9]
                H_[2,9] = E1*E2*E9/SUM^2
                H_[9,2] = H_[2,9]
                H_[3,9] = E1*E3*E9/SUM^2
                H_[9,3] = H_[3,9]
                H_[4,9] = E1*E4*E9/SUM^2
                H_[9,4] = H_[4,9]
                H_[5,9] = E1*E5*E9/SUM^2
                H_[9,5] = H_[5,9]
                H_[6,9] = E1*E6*E9/SUM^2
                H_[9,6] = H_[6,9]
                H_[7,9] = E1*E7*E9/SUM^2
                H_[9,7] = H_[7,9]
                H_[8,9] = E1*E8*E9/SUM^2
                H_[9,8] = H_[8,9]
                H_[9,9] = (-1.0e+0+E9/SUM)*E1*E9/SUM
                H_[1,10] = (-1.0e+0+E1/SUM)*E1*E10/SUM
                H_[10,1] = H_[1,10]
                H_[2,10] = E1*E2*E10/SUM^2
                H_[10,2] = H_[2,10]
                H_[3,10] = E1*E3*E10/SUM^2
                H_[10,3] = H_[3,10]
                H_[4,10] = E1*E4*E10/SUM^2
                H_[10,4] = H_[4,10]
                H_[5,10] = E1*E5*E10/SUM^2
                H_[10,5] = H_[5,10]
                H_[6,10] = E1*E6*E10/SUM^2
                H_[10,6] = H_[6,10]
                H_[7,10] = E1*E7*E10/SUM^2
                H_[10,7] = H_[7,10]
                H_[8,10] = E1*E8*E10/SUM^2
                H_[10,8] = H_[8,10]
                H_[9,10] = E1*E9*E10/SUM^2
                H_[10,9] = H_[9,10]
                H_[10,10] = (-1.0e+0+E10/SUM)*E1*E10/SUM
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

