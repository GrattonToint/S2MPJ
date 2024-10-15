function HS84(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    *********
#    Problem : HS84
#    *********
# 
#    Source: problem 84 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: A.R. Conn, March 1991.
# 
#    classification = "C-QQR2-AN-5-3"
# 
#    Set useful parameters
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "HS84"

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
        v_["15"] = 15
        v_["16"] = 16
        v_["17"] = 17
        v_["18"] = 18
        v_["19"] = 19
        v_["20"] = 20
        v_["21"] = 21
        v_["A"*string(Int64(v_["1"]))] = -24345.0
        v_["A"*string(Int64(v_["2"]))] = -8720288.849
        v_["MA"*string(Int64(v_["2"]))] = -1.0*v_["A"*string(Int64(v_["2"]))]
        v_["A"*string(Int64(v_["3"]))] = 150512.5253
        v_["MA"*string(Int64(v_["3"]))] = -1.0*v_["A"*string(Int64(v_["3"]))]
        v_["A"*string(Int64(v_["4"]))] = -156.6950325
        v_["MA"*string(Int64(v_["4"]))] = -1.0*v_["A"*string(Int64(v_["4"]))]
        v_["A"*string(Int64(v_["5"]))] = 476470.3222
        v_["MA"*string(Int64(v_["5"]))] = -1.0*v_["A"*string(Int64(v_["5"]))]
        v_["A"*string(Int64(v_["6"]))] = 729482.8271
        v_["MA"*string(Int64(v_["6"]))] = -1.0*v_["A"*string(Int64(v_["6"]))]
        v_["A"*string(Int64(v_["7"]))] = -145421.402
        v_["A"*string(Int64(v_["8"]))] = 2931.1506
        v_["A"*string(Int64(v_["9"]))] = -40.427932
        v_["A"*string(Int64(v_["10"]))] = 5106.192
        v_["A"*string(Int64(v_["11"]))] = 15711.36
        v_["A"*string(Int64(v_["12"]))] = -155011.1084
        v_["A"*string(Int64(v_["13"]))] = 4360.53352
        v_["A"*string(Int64(v_["14"]))] = 12.9492344
        v_["A"*string(Int64(v_["15"]))] = 10236.884
        v_["A"*string(Int64(v_["16"]))] = 13176.786
        v_["A"*string(Int64(v_["17"]))] = -326669.5104
        v_["A"*string(Int64(v_["18"]))] = 7390.68412
        v_["A"*string(Int64(v_["19"]))] = -27.8986976
        v_["A"*string(Int64(v_["20"]))] = 16643.076
        v_["A"*string(Int64(v_["21"]))] = 30988.146
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["5"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X"*string(Int64(v_["1"]))]
        pbm.A[ig,iv] += Float64(v_["MA"*string(Int64(v_["2"]))])
        ig,ig_,_ = s2mpj_ii("CON1",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"CON1")
        iv = ix_["X"*string(Int64(v_["1"]))]
        pbm.A[ig,iv] += Float64(v_["A"*string(Int64(v_["7"]))])
        ig,ig_,_ = s2mpj_ii("CON2",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"CON2")
        iv = ix_["X"*string(Int64(v_["1"]))]
        pbm.A[ig,iv] += Float64(v_["A"*string(Int64(v_["12"]))])
        ig,ig_,_ = s2mpj_ii("CON3",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"CON3")
        iv = ix_["X"*string(Int64(v_["1"]))]
        pbm.A[ig,iv] += Float64(v_["A"*string(Int64(v_["17"]))])
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
        pbm.gconst[ig_["OBJ"]] = Float64(v_["A"*string(Int64(v_["1"]))])
        #%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%
        grange = Vector{Float64}(undef,ngrp)
        grange[gegrps,1] = fill(Inf,pb.nge)
        arrset(grange,ig_["CON1"],Float64(294000.0))
        arrset(grange,ig_["CON2"],Float64(294000.0))
        arrset(grange,ig_["CON3"],Float64(277200.0))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower[ix_["X"*string(Int64(v_["1"]))]] = 0.0
        pb.xlower[ix_["X"*string(Int64(v_["2"]))]] = 1.2
        pb.xlower[ix_["X"*string(Int64(v_["3"]))]] = 20.0
        pb.xlower[ix_["X"*string(Int64(v_["4"]))]] = 9.0
        pb.xlower[ix_["X"*string(Int64(v_["5"]))]] = 6.5
        pb.xupper[ix_["X"*string(Int64(v_["1"]))]] = 1000.0
        pb.xupper[ix_["X"*string(Int64(v_["2"]))]] = 2.4
        pb.xupper[ix_["X"*string(Int64(v_["3"]))]] = 60.0
        pb.xupper[ix_["X"*string(Int64(v_["4"]))]] = 9.3
        pb.xupper[ix_["X"*string(Int64(v_["5"]))]] = 7.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        pb.x0[ix_["X"*string(Int64(v_["1"]))]] = Float64(2.52)
        pb.x0[ix_["X"*string(Int64(v_["2"]))]] = Float64(2.0)
        pb.x0[ix_["X"*string(Int64(v_["3"]))]] = Float64(37.5)
        pb.x0[ix_["X"*string(Int64(v_["4"]))]] = Float64(9.25)
        pb.x0[ix_["X"*string(Int64(v_["5"]))]] = Float64(6.8)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "ePROD", iet_)
        loaset(elftv,it,1,"U1")
        loaset(elftv,it,2,"U2")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["4"])
            ename = "E"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"ePROD")
            arrset(ielftype,ie,iet_["ePROD"])
            v_["IP1"] = 1+I
            vname = "X"*string(Int64(v_["1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="U1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["IP1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="U2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["OBJ"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["1"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["MA"*string(Int64(v_["3"]))]))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["2"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["MA"*string(Int64(v_["4"]))]))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["3"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["MA"*string(Int64(v_["5"]))]))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["4"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["MA"*string(Int64(v_["6"]))]))
        ig = ig_["CON1"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["1"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["A"*string(Int64(v_["8"]))]))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["2"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["A"*string(Int64(v_["9"]))]))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["3"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["A"*string(Int64(v_["10"]))]))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["4"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["A"*string(Int64(v_["11"]))]))
        ig = ig_["CON2"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["1"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["A"*string(Int64(v_["13"]))]))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["2"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["A"*string(Int64(v_["14"]))]))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["3"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["A"*string(Int64(v_["15"]))]))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["4"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["A"*string(Int64(v_["16"]))]))
        ig = ig_["CON3"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["1"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["A"*string(Int64(v_["18"]))]))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["2"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["A"*string(Int64(v_["19"]))]))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["3"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["A"*string(Int64(v_["20"]))]))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["4"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["A"*string(Int64(v_["21"]))]))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+pb.neq+1:pb.m] = zeros(Float64,pb.nge)
        pb.cupper[1:pb.nge] = grange[gegrps]
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-QQR2-AN-5-3"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2]
        pb.conderlvl  = pbm.conderlvl;
        return pb, pbm


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

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

