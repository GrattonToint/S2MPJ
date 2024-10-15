function CHARDIS0(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : CHARDIS0
#    *********
# 
#    Distribution of (equal)charges on [-R,R]x[-R,R] (2D)
# 
#    SIF input: R. Felkel, Jun 1999.
#               incorrectly decoded version (see CHARDIS0 for correction)
# 
#    classification = "C-OBR2-AY-V-V"
# 
#    Number of positive (or negative) charges -> Number of variables 2*NP1
# 
#       Alternative values for the SIF file parameters:
# IE NP1                 5              $-PARAMETER
# IE NP1                 9              $-PARAMETER
# IE NP1                 20             $-PARAMETER
# IE NP1                 30             $-PARAMETER
# IE NP1                 50             $-PARAMETER     original value
# IE NP1                 100            $-PARAMETER
# IE NP1                 200            $-PARAMETER
# IE NP1                 500            $-PARAMETER
# IE NP1                 1000           $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "CHARDIS0"

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
            v_["NP1"] = Int64(20);  #  SIF file default value
        else
            v_["NP1"] = Int64(args[1]);
        end
# IE NP1                 2000           $-PARAMETER
# IE NP1                 5000           $-PARAMETER
        v_["R"] = 10.0
        v_["R-"] = -10.0
        v_["N"] = -1+v_["NP1"]
        v_["NReal"] = Float64(v_["N"])
        v_["NP1Real"] = Float64(v_["NP1"])
        v_["halfPI"] = asin(1.0)
        v_["PI"] = 2.0*v_["halfPI"]
        v_["2PI"] = 4.0*v_["halfPI"]
        v_["4PI"] = 8.0*v_["halfPI"]
        v_["4PIqN"] = v_["4PI"]/v_["NReal"]
        v_["2PIqN"] = v_["2PI"]/v_["NReal"]
        v_["PIqN"] = v_["PI"]/v_["NReal"]
        v_["RqN"] = v_["R"]/v_["NReal"]
        v_["1"] = 1
        v_["2"] = 2
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["NP1"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
            iv,ix_,_ = s2mpj_ii("Y"*string(I),ix_)
            arrset(pb.xnames,iv,"Y"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["1"]):Int64(v_["NP1"])
            v_["I+"] = 1+I
            for J = Int64(v_["I+"]):Int64(v_["NP1"])
                ig,ig_,_ = s2mpj_ii("O"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"<>")
                arrset(pbm.gscale,ig,Float64(0.01))
            end
        end
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        for I = Int64(v_["1"]):Int64(v_["NP1"])
            pb.xlower[ix_["X"*string(I)]] = v_["R-"]
            pb.xupper[ix_["X"*string(I)]] = v_["R"]
            pb.xlower[ix_["Y"*string(I)]] = v_["R-"]
            pb.xupper[ix_["Y"*string(I)]] = v_["R"]
        end
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        for I = Int64(v_["1"]):Int64(v_["NP1"])
            v_["RealI-"] = Float64(I)
            v_["RealNP1-I"] = v_["NP1Real"]-v_["RealI-"]
            v_["PHII-"] = v_["2PIqN"]*v_["RealI-"]
            v_["RI-"] = v_["RqN"]*v_["RealNP1-I"]
            v_["XSTT"] = cos(v_["PHII-"])
            v_["YSTT"] = sin(v_["PHII-"])
            v_["XST"] = v_["XSTT"]*v_["RI-"]
            v_["YST"] = v_["YSTT"]*v_["RI-"]
            v_["XS"] = 0.5*v_["XST"]
            v_["YS"] = 0.5*v_["YST"]
            pb.x0[ix_["X"*string(I)]] = Float64(v_["XS"])
            pb.x0[ix_["Y"*string(I)]] = Float64(v_["YS"])
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eDIFSQR", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["NP1"])
            v_["I+"] = 1+I
            for J = Int64(v_["I+"]):Int64(v_["NP1"])
                ename = "X"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eDIFSQR")
                arrset(ielftype,ie,iet_["eDIFSQR"])
                vname = "X"*string(I)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                ename = "Y"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eDIFSQR")
                arrset(ielftype,ie,iet_["eDIFSQR"])
                vname = "Y"*string(I)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "Y"*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
        end
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gREZIP",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["NP1"])
            v_["I+"] = 1+I
            for J = Int64(v_["I+"]):Int64(v_["NP1"])
                ig = ig_["O"*string(I)*","*string(J)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["X"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel,1.)
                posel = posel+1
                loaset(pbm.grelt,ig,posel,ie_["Y"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel, 1.)
            end
        end
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.A = spzeros(Float64,0,0)
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-OBR2-AY-V-V"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        return pb, pbm


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "eDIFSQR"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = (EV_[1]-EV_[2])*(EV_[1]-EV_[2])
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 2.0*(EV_[1]-EV_[2])
            g_[2] = -2.0*(EV_[1]-EV_[2])
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = 2.0
                H_[1,2] = -2.0
                H_[2,1] = H_[1,2]
                H_[2,2] = 2.0
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

    elseif action == "gREZIP"

        GVAR_   = args[1]
        igr_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_= 1.0/GVAR_
        if nargout>1
            g_ = -1.0/(GVAR_*GVAR_)
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_ = 2.0/(GVAR_*GVAR_*GVAR_)
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

