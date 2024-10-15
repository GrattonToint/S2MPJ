function DITTERT(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DITTERT
#    *********
# 
#    Minimize the Dittert function.
#    Source: See Minc, Linear and Multilinear Algebra 21, 1987
# 
#    SIF input: N. Gould, March 1992.
#               minor correction by Ph. Shott, Jan 1995.
# 
#    classification = "C-OQR2-AN-V-V"
# 
#    Size of matrix
# 
#       Alternative values for the SIF file parameters:
# IE N                   2              $-PARAMETER
# IE N                   3              $-PARAMETER
# IE N                   4              $-PARAMETER
# IE N                   5              $-PARAMETER
# IE N                   6              $-PARAMETER
# IE N                   7              $-PARAMETER
# IE N                   8              $-PARAMETER     original value
# IE N                   9              $-PARAMETER
# IE N                   10             $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "DITTERT"

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
            v_["N"] = Int64(5);  #  SIF file default value
        else
            v_["N"] = Int64(args[1]);
        end
        v_["0"] = 0
        v_["1"] = 1
        v_["2"] = 2
        v_["RN"] = Float64(v_["N"])
        v_["N+1"] = 1+v_["N"]
        v_["2**N"] = 1
        for I = Int64(v_["1"]):Int64(v_["N"])
            v_["N-I+1"] = v_["N+1"]-I
            v_["I-1"] = - 1+I
            v_["R2**N"] = Float64(v_["2**N"])
            v_["S"*string(Int64(v_["N-I+1"]))] = 0.1+v_["R2**N"]
            v_["T"*string(Int64(v_["I-1"]))] = 0.1+v_["R2**N"]
            v_["2**N"] = v_["2**N"]*v_["2"]
        end
        v_["N-1"] = - 1+v_["N"]
        v_["2**N-1"] = - 1+v_["2**N"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for M = Int64(v_["1"]):Int64(v_["N-1"])
            v_["RK1"] = v_["T"*string(M)]
            v_["K1"] = trunc(Int,v_["RK1"])
            v_["K2"] = 2*v_["K1"]
            v_["K1"] = 1+v_["K1"]
            v_["K2"] = - 1+v_["K2"]
            for K = Int64(v_["K1"]):Int64(v_["K2"])
                iv,ix_,_ = s2mpj_ii("P"*string(K),ix_)
                arrset(pb.xnames,iv,"P"*string(K))
            end
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            for J = Int64(v_["1"]):Int64(v_["N"])
                iv,ix_,_ = s2mpj_ii("A"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"A"*string(I)*","*string(J))
            end
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("R"*string(I),ix_)
            arrset(pb.xnames,iv,"R"*string(I))
            iv,ix_,_ = s2mpj_ii("C"*string(I),ix_)
            arrset(pb.xnames,iv,"C"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("ROWPROD",ig_)
        arrset(gtype,ig,"<>")
        ig,ig_,_ = s2mpj_ii("COLPROD",ig_)
        arrset(gtype,ig,"<>")
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["P"*string(Int64(v_["2**N-1"]))]
        pbm.A[ig,iv] += Float64(1.0)
        for M = Int64(v_["1"]):Int64(v_["N-1"])
            v_["RK1"] = v_["T"*string(M)]
            v_["K1"] = trunc(Int,v_["RK1"])
            v_["K2"] = 2*v_["K1"]
            v_["K1"] = 1+v_["K1"]
            v_["K2"] = - 1+v_["K2"]
            for K = Int64(v_["K1"]):Int64(v_["K2"])
                ig,ig_,_ = s2mpj_ii("PE"*string(K),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"PE"*string(K))
                iv = ix_["P"*string(K)]
                pbm.A[ig,iv] += Float64(- 1.0)
            end
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            for J = Int64(v_["1"]):Int64(v_["N"])
                ig,ig_,_ = s2mpj_ii("R"*string(I),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"R"*string(I))
                iv = ix_["A"*string(I)*","*string(J)]
                pbm.A[ig,iv] += Float64(1.0)
                ig,ig_,_ = s2mpj_ii("C"*string(J),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"C"*string(J))
                iv = ix_["A"*string(I)*","*string(J)]
                pbm.A[ig,iv] += Float64(1.0)
            end
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            ig,ig_,_ = s2mpj_ii("R"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"R"*string(I))
            iv = ix_["R"*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("C"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"C"*string(I))
            iv = ix_["C"*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            ig,ig_,_ = s2mpj_ii("SUM",ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"SUM")
            iv = ix_["R"*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
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
        pbm.gconst[ig_["SUM"]] = Float64(v_["RN"])
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        for I = Int64(v_["1"]):Int64(v_["N"])
            pb.xlower[ix_["R"*string(I)]] = 1.0e-6
            pb.xlower[ix_["C"*string(I)]] = 1.0e-6
            for J = Int64(v_["1"]):Int64(v_["N"])
                pb.xupper[ix_["A"*string(I)*","*string(J)]] = 1.0
            end
        end
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "en2PR", iet_)
        loaset(elftv,it,1,"A")
        loaset(elftv,it,2,"P")
        it,iet_,_ = s2mpj_ii( "eLOG", iet_)
        loaset(elftv,it,1,"Y")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for M = Int64(v_["1"]):Int64(v_["N-1"])
            v_["RK1"] = v_["T"*string(M)]
            v_["K1"] = trunc(Int,v_["RK1"])
            v_["K2"] = 2*v_["K1"]
            v_["K1"] = 1+v_["K1"]
            v_["K2"] = - 1+v_["K2"]
            for K = Int64(v_["K1"]):Int64(v_["K2"])
                v_["ID"] = 0
                v_["PT"] = 1
                v_["KK"] = K
                for I = Int64(v_["1"]):Int64(v_["N"])
                    v_["SI"] = v_["S"*string(I)]
                    v_["ISI"] = trunc(Int,v_["SI"])
                    v_["BI"] = trunc(Int,(v_["KK"]/v_["ISI"]))
                    v_["ID"] = v_["ID"]+v_["BI"]
                    v_["BISI"] = v_["BI"]*v_["ISI"]
                    v_["KK"] = v_["KK"]-v_["BISI"]
                    v_["RI"] = Float64(I)
                    v_["RNZ"*string(Int64(v_["PT"]))] = 0.1+v_["RI"]
                    v_["PT"] = v_["PT"]+v_["BI"]
                end
                v_["I1"] = v_["0"]
                v_["I2"] = v_["1"]
                v_["ID-2"] = - 2+v_["ID"]
                for I = Int64(v_["1"]):Int64(v_["ID-2"])
                    v_["I1"] = v_["ID"]
                    v_["I2"] = v_["0"]
                end
                for I = Int64(v_["1"]):Int64(v_["I1"])
                    v_["RJ"] = v_["RNZ"*string(I)]
                    v_["J"] = trunc(Int,v_["RJ"])
                    v_["SI"] = v_["S"*string(Int64(v_["J"]))]
                    v_["ISI"] = trunc(Int,v_["SI"])
                    v_["IPP"] = K-v_["ISI"]
                    ename = "E"*string(K)*","*string(I)
                    ie,ie_,_  = s2mpj_ii(ename,ie_)
                    arrset(pbm.elftype,ie,"en2PR")
                    arrset(ielftype,ie,iet_["en2PR"])
                    vname = "A"*string(Int64(v_["ID"]))*","*string(Int64(v_["J"]))
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="A",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                    vname = "P"*string(Int64(v_["IPP"]))
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="P",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                end
                for I = Int64(v_["1"]):Int64(v_["I2"])
                    v_["RJ"] = v_["RNZ"*string(Int64(v_["1"]))]
                    v_["J"] = trunc(Int,v_["RJ"])
                    v_["RJJ"] = v_["RNZ"*string(Int64(v_["2"]))]
                    v_["JJ"] = trunc(Int,v_["RJJ"])
                    ename = "E"*string(K)*","*string(Int64(v_["1"]))
                    ie,ie_,_  = s2mpj_ii(ename,ie_)
                    arrset(pbm.elftype,ie,"en2PR")
                    arrset(ielftype,ie,iet_["en2PR"])
                    ename = "E"*string(K)*","*string(Int64(v_["1"]))
                    ie,ie_,_  = s2mpj_ii(ename,ie_)
                    vname = "A"*string(Int64(v_["2"]))*","*string(Int64(v_["J"]))
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="A",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                    ename = "E"*string(K)*","*string(Int64(v_["1"]))
                    ie,ie_,_  = s2mpj_ii(ename,ie_)
                    vname = "A"*string(Int64(v_["1"]))*","*string(Int64(v_["JJ"]))
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="P",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                    ename = "E"*string(K)*","*string(Int64(v_["2"]))
                    ie,ie_,_  = s2mpj_ii(ename,ie_)
                    arrset(pbm.elftype,ie,"en2PR")
                    arrset(ielftype,ie,iet_["en2PR"])
                    ename = "E"*string(K)*","*string(Int64(v_["2"]))
                    ie,ie_,_  = s2mpj_ii(ename,ie_)
                    vname = "A"*string(Int64(v_["2"]))*","*string(Int64(v_["JJ"]))
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="A",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                    ename = "E"*string(K)*","*string(Int64(v_["2"]))
                    ie,ie_,_  = s2mpj_ii(ename,ie_)
                    vname = "A"*string(Int64(v_["1"]))*","*string(Int64(v_["J"]))
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="P",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                end
                v_["RD"] = Float64(v_["ID"])
                v_["D"*string(K)] = 0.1+v_["RD"]
            end
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            ename = "LOGC"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eLOG")
            arrset(ielftype,ie,iet_["eLOG"])
            vname = "C"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "LOGR"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eLOG")
            arrset(ielftype,ie,iet_["eLOG"])
            vname = "R"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gEXP",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for M = Int64(v_["1"]):Int64(v_["N-1"])
            v_["RK1"] = v_["T"*string(M)]
            v_["K1"] = trunc(Int,v_["RK1"])
            v_["K2"] = 2*v_["K1"]
            v_["K1"] = 1+v_["K1"]
            v_["K2"] = - 1+v_["K2"]
            for K = Int64(v_["K1"]):Int64(v_["K2"])
                v_["RD"] = v_["D"*string(K)]
                v_["ID"] = trunc(Int,v_["RD"])
                for I = Int64(v_["1"]):Int64(v_["ID"])
                    ig = ig_["PE"*string(K)]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["E"*string(K)*","*string(I)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,1.)
                end
            end
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            ig = ig_["ROWPROD"]
            arrset(pbm.grftype,ig,"gEXP")
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["LOGR"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
            ig = ig_["COLPROD"]
            arrset(pbm.grftype,ig,"gEXP")
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["LOGC"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN(2)           5.0D-1
# LO SOLTN(3)           2.22222222D-1
# LO SOLTN(4)           9.375-2
# LO SOLTN(5)           3.84D-2
# LO SOLTN(6)           1.54321098D-2
# LO SOLTN(7)           6.11989902D-3
# LO SOLTN(8)           2.40325927D-3
# LO SOLTN(9)           9.36656708D-4
# LO SOLTN(10)          3.6288D-4
# LO SOLTN(11)          1.39905948D-4
# LO SOLTN(12)          5.37232170D-5
# LO SOLTN(13)          2.05596982D-5
# LO SOLTN(14)          7.84541375D-6
# LO SOLTN(15)          2.98628137D-6
# LO SOLTN(16)          1.13422671D-6
# LO SOLTN(17)          4.29968709D-7
# LO SOLTN(18)          1.62718123D-7
# LO SOLTN(19)          6.14859946D-8
# LO SOLTN(20)          2.32019615D-8
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
        pb.pbclass = "C-OQR2-AN-V-V"
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

    elseif action == "eLOG"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = log(EV_[1])
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 1.0/EV_[1]
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = -1.0/(EV_[1]*EV_[1])
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

    elseif action == "gEXP"

        GVAR_   = args[1]
        igr_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        EXZ = -exp(GVAR_)
        f_= EXZ
        if nargout>1
            g_ = EXZ
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_ = EXZ
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

