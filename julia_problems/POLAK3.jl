function POLAK3(action,args...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : POLAK3
#    *********
# 
#    A nonlinear minmax problem in eleven variables.
# 
#    Note: the original statement of the problem contains an inconsistant
#          index i.  This has been replaced by 1, assuming a very common typo.
#          But the optimal solution of the resulting problem differs from that
#          quoted in the source. 
# 
#    Source: 
#    E. Polak, D.H. Mayne and J.E. Higgins,
#    "Superlinearly convergent algorithm for min-max problems"
#    JOTA 69, pp. 407-439, 1991.
# 
#    SIF input: Ph. Toint, Nov 1993.
# 
#    classification = "LOR2-AN-12-10"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "POLAK3"

    if action == "setup"
        pbm          = PBM(name)
        pb           = PB(name)
        pb.sifpbname = "POLAK3"
        nargin       = length(args)
        pbm.call     = eval( Meta.parse( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["1"] = 1
        v_["10"] = 10
        v_["11"] = 11
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        xscale  = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["11"])
            iv,ix_,_ = s2x_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        iv,ix_,_ = s2x_ii("U",ix_)
        arrset(pb.xnames,iv,"U")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2x_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["U"]
        pbm.A[ig,iv] += Float64(1.0)
        for I = Int64(v_["1"]):Int64(v_["10"])
            ig,ig_,_ = s2x_ii("F"*string(I),ig_)
            arrset(gtype,ig,"<=")
            arrset(pb.cnames,ig,"F"*string(I))
            iv = ix_["U"]
            pbm.A[ig,iv] += Float64(-1.0)
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
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(1.0),pb.n)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2x_ii( "eEL", iet_)
        loaset(elftv,it,1,"XX")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"A")
        loaset(elftp,it,2,"B")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["10"])
            v_["I-1"] = -1+I
            v_["RI-1"] = Float64(v_["I-1"])
            for J = Int64(v_["1"]):Int64(v_["11"])
                v_["J-1"] = -1+J
                v_["RJ-1"] = Float64(J)
                ename = "E"*string(I)*","*string(J)
                ie,ie_,_  = s2x_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eEL")
                arrset(ielftype, ie, iet_["eEL"])
                vname = "X"*string(J)
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,1.0)
                posev = findfirst(x->x=="XX",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                posep = findfirst(x->x=="A",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["RI-1"]))
                posep = findfirst(x->x=="B",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["RJ-1"]))
            end
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["10"])
            for J = Int64(v_["1"]):Int64(v_["11"])
                v_["RJ"] = Float64(J)
                v_["1/J"] = 1.0/v_["RJ"]
                ig = ig_["F"*string(I)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["E"*string(I)*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["1/J"]))
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.cupper[1:pb.nle] = zeros(Float64,pb.nle)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "LOR2-AN-12-10"
        return pb, pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "eEL"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        V = EV_[1]-sin(pbm.elpar[iel_][1]+pbm.elpar[iel_][2]+pbm.elpar[iel_][2])
        EV = exp(V*V)
        f_   = EV
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 2.0*V*EV
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = (4.0*V*V+2.0)*EV
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
