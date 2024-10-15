function BLOWEYC(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : BLOWEYC
#    *********
# 
#    A nonconvex quadratic program proposed by 
#    James Blowey (University of Durham)
# 
#    Given function v(s) and f(s) = v(s) + A(inv) v(s), s in [0,1],
#    minimize 
# 
#         (u(s) - v(s))(trans) ( A + A(inv) ) u(s) - (u(s) - v(s))(trans)f(s)
# 
#    where 
# 
#       u(s) in [-1,1] and int[0,1] u(s) ds = int[0,1] v(s) ds
# 
#    and A is the 
# 
#       "- Laplacian with Neumann boundary conditions on a uniform mesh"
# 
#    The troublesome term A(inv) u(s) is replaced by the additional 
#    variable w(s) and the constraint A w(s) = u(s)
# 
#    The function v(s) is chosen to be 
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 6 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "BLOWEYC"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
#    classification = "C-QLR2-MN-V-V"
#       Alternative values for the SIF file parameters:
# IE N                   10             $-PARAMETER  n = 22, m = 12
# IE N                   100            $-PARAMETER  n = 202, m = 102
# IE N                   1000           $-PARAMETER  n = 2002, m = 1002
# IE N                   2000           $-PARAMETER  n = 4002, m = 2002
        if nargin<1
            v_["N"] = Int64(10);  #  SIF file default value
        else
            v_["N"] = Int64(args[1]);
        end
# IE N                   4000           $-PARAMETER  n = 8002, m = 4002
# IE N                   8000           $-PARAMETER  n = 16002, m = 8002
        v_["0"] = 0
        v_["1"] = 1
        v_["2"] = 2
        v_["3"] = 3
        v_["4"] = 4
        v_["5"] = 5
        v_["ONE"] = 1.0
        v_["-ONE"] = -1.0
        v_["TWO"] = 2.0
        v_["-TWO"] = -2.0
        v_["RN"] = Float64(v_["N"])
        v_["N**2"] = v_["RN"]*v_["RN"]
        v_["N-1"] = -1+v_["N"]
        v_["1/N**2"] = v_["ONE"]/v_["N**2"]
        v_["-1/N**2"] = v_["-ONE"]/v_["N**2"]
        v_["-2/N**2"] = 2.0*v_["-1/N**2"]
        v_["2N**2"] = 2.0*v_["N**2"]
        v_["-2N**2"] = -2.0*v_["N**2"]
        v_["N/2"] = trunc(Int,(v_["N"]/v_["2"]))
        v_["N/5"] = trunc(Int,(v_["N"]/v_["5"]))
        v_["NA"] = v_["N/5"]
        v_["A"] = Float64(v_["NA"])
        v_["A"] = v_["A"]/v_["RN"]
        v_["NA+1"] = 1+v_["NA"]
        v_["NB"] = v_["N/2"]
        v_["B"] = Float64(v_["NB"])
        v_["B"] = v_["B"]/v_["RN"]
        v_["NB+1"] = 1+v_["NB"]
        v_["NC"] = v_["N/2"]
        v_["C"] = Float64(v_["NC"])
        v_["C"] = v_["C"]/v_["RN"]
        v_["NC+1"] = 1+v_["NC"]
        v_["ND"] = v_["N/5"]*v_["4"]
        v_["D"] = Float64(v_["ND"])
        v_["D"] = v_["D"]/v_["RN"]
        v_["ND+1"] = 1+v_["ND"]
        v_["INT"] = v_["ONE"]
        v_["INT"] = v_["INT"]+v_["A"]
        v_["INT"] = v_["INT"]+v_["B"]
        v_["INT"] = v_["INT"]-v_["C"]
        v_["INT"] = v_["INT"]-v_["D"]
        v_["INT"] = v_["INT"]*v_["RN"]
        for I = Int64(v_["0"]):Int64(v_["NA"])
            v_["V"*string(I)] = 1.0
        end
        v_["STEP"] = v_["B"]-v_["A"]
        v_["STEP"] = v_["STEP"]*v_["RN"]
        v_["STEP"] = v_["TWO"]/v_["STEP"]
        for I = Int64(v_["NA+1"]):Int64(v_["NB"])
            v_["J"] = I-v_["NA"]
            v_["RJ"] = Float64(v_["J"])
            v_["VAL"] = v_["RJ"]*v_["STEP"]
            v_["VAL"] = v_["ONE"]-v_["VAL"]
            v_["V"*string(I)] = v_["VAL"]
        end
        for I = Int64(v_["NB+1"]):Int64(v_["NC"])
            v_["V"*string(I)] = -1.0
        end
        v_["STEP"] = v_["D"]-v_["C"]
        v_["STEP"] = v_["STEP"]*v_["RN"]
        v_["STEP"] = v_["TWO"]/v_["STEP"]
        for I = Int64(v_["NC+1"]):Int64(v_["ND"])
            v_["J"] = I-v_["NC"]
            v_["RJ"] = Float64(v_["J"])
            v_["VAL"] = v_["RJ"]*v_["STEP"]
            v_["VAL"] = v_["-ONE"]+v_["VAL"]
            v_["V"*string(I)] = v_["VAL"]
        end
        for I = Int64(v_["ND"]):Int64(v_["N"])
            v_["V"*string(I)] = 1.0
        end
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["0"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("U"*string(I),ix_)
            arrset(pb.xnames,iv,"U"*string(I))
            iv,ix_,_ = s2mpj_ii("W"*string(I),ix_)
            arrset(pb.xnames,iv,"W"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["0"]):Int64(v_["N"])
            v_["VAL"] = v_["V"*string(I)]*v_["-1/N**2"]
            ig,ig_,_ = s2mpj_ii("OBJ",ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["U"*string(I)]
            pbm.A[ig,iv] += Float64(v_["VAL"])
            v_["VAL"] = v_["V"*string(I)]*v_["-2/N**2"]
            iv = ix_["W"*string(I)]
            pbm.A[ig,iv] += Float64(v_["VAL"])
        end
        v_["VAL"] = v_["V"*string(Int64(v_["1"]))]-v_["V"*string(Int64(v_["0"]))]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["U"*string(Int64(v_["0"]))]
        pbm.A[ig,iv] += Float64(v_["VAL"])
        for I = Int64(v_["1"]):Int64(v_["N-1"])
            v_["I+1"] = 1+I
            v_["I-1"] = -1+I
            v_["VAL"] = -2.0*v_["V"*string(I)]
            v_["VAL"] = v_["VAL"]+v_["V"*string(Int64(v_["I-1"]))]
            v_["VAL"] = v_["VAL"]+v_["V"*string(Int64(v_["I+1"]))]
            ig,ig_,_ = s2mpj_ii("OBJ",ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["U"*string(I)]
            pbm.A[ig,iv] += Float64(v_["VAL"])
        end
        v_["VAL"] = v_["V"*string(Int64(v_["N-1"]))]-v_["V"*string(Int64(v_["N"]))]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["U"*string(Int64(v_["N"]))]
        pbm.A[ig,iv] += Float64(v_["VAL"])
        ig,ig_,_ = s2mpj_ii("INT",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"INT")
        iv = ix_["U"*string(Int64(v_["0"]))]
        pbm.A[ig,iv] += Float64(0.5)
        ig,ig_,_ = s2mpj_ii("CON"*string(Int64(v_["0"])),ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"CON"*string(Int64(v_["0"])))
        iv = ix_["U"*string(Int64(v_["0"]))]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["U"*string(Int64(v_["1"]))]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("CON"*string(Int64(v_["0"])),ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"CON"*string(Int64(v_["0"])))
        iv = ix_["W"*string(Int64(v_["0"]))]
        pbm.A[ig,iv] += Float64(v_["-1/N**2"])
        for I = Int64(v_["1"]):Int64(v_["N-1"])
            v_["I+1"] = 1+I
            v_["I-1"] = -1+I
            ig,ig_,_ = s2mpj_ii("CON"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"CON"*string(I))
            iv = ix_["U"*string(I)]
            pbm.A[ig,iv] += Float64(2.0)
            iv = ix_["U"*string(Int64(v_["I+1"]))]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["U"*string(Int64(v_["I-1"]))]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["W"*string(I)]
            pbm.A[ig,iv] += Float64(v_["-1/N**2"])
            ig,ig_,_ = s2mpj_ii("INT",ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"INT")
            iv = ix_["U"*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
        end
        ig,ig_,_ = s2mpj_ii("CON"*string(Int64(v_["N"])),ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"CON"*string(Int64(v_["N"])))
        iv = ix_["U"*string(Int64(v_["N"]))]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["U"*string(Int64(v_["N-1"]))]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("CON"*string(Int64(v_["N"])),ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"CON"*string(Int64(v_["N"])))
        iv = ix_["W"*string(Int64(v_["N"]))]
        pbm.A[ig,iv] += Float64(v_["-1/N**2"])
        ig,ig_,_ = s2mpj_ii("INT",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"INT")
        iv = ix_["U"*string(Int64(v_["N"]))]
        pbm.A[ig,iv] += Float64(0.5)
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
        pbm.gconst[ig_["INT"]] = Float64(v_["INT"])
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        for I = Int64(v_["0"]):Int64(v_["N"])
            pb.xlower[ix_["U"*string(I)]] = -1.0
            pb.xupper[ix_["U"*string(I)]] = 1.0
        end
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(0.0),pb.n)
        for I = Int64(v_["0"]):Int64(v_["N"])
            pb.x0[ix_["U"*string(I)]] = Float64(v_["V"*string(I)])
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eSQ", iet_)
        loaset(elftv,it,1,"Z")
        it,iet_,_ = s2mpj_ii( "ePROD", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["0"]):Int64(v_["N"])
            ename = "C"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"ePROD")
            arrset(ielftype,ie,iet_["ePROD"])
            vname = "U"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "W"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
            posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        for I = Int64(v_["0"]):Int64(v_["N-1"])
            v_["I+1"] = 1+I
            ename = "D"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSQ")
            arrset(ielftype,ie,iet_["eSQ"])
            vname = "U"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
            posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "O"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"ePROD")
            arrset(ielftype,ie,iet_["ePROD"])
            vname = "U"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "U"*string(Int64(v_["I+1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
            posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        ename = "D"*string(Int64(v_["N"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        ename = "D"*string(Int64(v_["N"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "U"*string(Int64(v_["N"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["OBJ"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["O"*string(Int64(v_["0"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["-TWO"]))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["D"*string(Int64(v_["0"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["ONE"]))
        for I = Int64(v_["1"]):Int64(v_["N-1"])
            ig = ig_["OBJ"]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["O"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["-TWO"]))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["D"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["TWO"]))
        end
        ig = ig_["OBJ"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["D"*string(Int64(v_["N"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["ONE"]))
        for I = Int64(v_["0"]):Int64(v_["N"])
            ig = ig_["OBJ"]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["C"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["1/N**2"]))
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# XL SOLUTION            -2.38816D+02   $ N = 10 
# XL SOLUTION            -2.65340D+03   $ N = 100
# XL SOLUTION            -2.67211D+04   $ N = 1000
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
        pb.pbclass = "C-QLR2-MN-V-V"
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

    elseif action == "eSQ"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]*EV_[1]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[1]+EV_[1]
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 2.0
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "ePROD"

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

