function MINSURFO(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : MINSURFO
#    *********
# 
#    Find the surface with minimal area, given boundary conditions, 
#    and above an obstacle.
# 
#    This is problem 17 in the COPS (Version 2) collection of 
#    E. Dolan and J. More'
#    see "Benchmarking Optimization Software with COPS"
#    Argonne National Labs Technical Report ANL/MCS-246 (2000)
# 
#    SIF input: Nick Gould, December 2000
# 
#    classification = "C-OBR2-AN-V-V"
# 
#  grid points in x direction (fixed at 50 in COPS)
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "MINSURFO"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["NX"] = 5
        v_["NY"] = 10
        v_["0"] = 0
        v_["1"] = 1
        v_["ONE"] = 1.0
        v_["NX+1"] = 1+v_["NX"]
        v_["NY+1"] = 1+v_["NY"]
        v_["RNX+1"] = Float64(v_["NX+1"])
        v_["RNY+1"] = Float64(v_["NY+1"])
        v_["HX"] = 1.0/v_["RNX+1"]
        v_["HY"] = 1.0/v_["RNY+1"]
        v_["AREA"] = v_["HX"]*v_["HY"]
        v_["AREA"] = 0.5*v_["AREA"]
        v_["1/AREA"] = 1.0/v_["AREA"]
        v_["1/HX"] = 1.0/v_["HX"]
        v_["1/HX2"] = v_["1/HX"]*v_["1/HX"]
        v_["1/HY"] = 1.0/v_["HY"]
        v_["1/HY2"] = v_["1/HY"]*v_["1/HY"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["0"]):Int64(v_["NX+1"])
            for J = Int64(v_["0"]):Int64(v_["NY+1"])
                iv,ix_,_ = s2mpj_ii("V"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"V"*string(I)*","*string(J))
            end
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["0"]):Int64(v_["NX"])
            for J = Int64(v_["0"]):Int64(v_["NY"])
                ig,ig_,_ = s2mpj_ii("A"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"<>")
                arrset(pbm.gscale,ig,Float64(v_["1/AREA"]))
            end
        end
        for I = Int64(v_["1"]):Int64(v_["NX+1"])
            for J = Int64(v_["1"]):Int64(v_["NY+1"])
                ig,ig_,_ = s2mpj_ii("B"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"<>")
                arrset(pbm.gscale,ig,Float64(v_["1/AREA"]))
            end
        end
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        for I = Int64(v_["0"]):Int64(v_["NX"])
            for J = Int64(v_["0"]):Int64(v_["NY"])
                pbm.gconst[ig_["A"*string(I)*","*string(J)]] = Float64(-1.0)
            end
        end
        for I = Int64(v_["1"]):Int64(v_["NX+1"])
            for J = Int64(v_["1"]):Int64(v_["NY+1"])
                pbm.gconst[ig_["B"*string(I)*","*string(J)]] = Float64(-1.0)
            end
        end
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        v_["1/4HX"] = 0.25/v_["HX"]
        v_["3/4HX"] = 0.75/v_["HX"]
        v_["1/4HY"] = 0.25/v_["HY"]
        v_["3/4HY"] = 0.75/v_["HY"]
        v_["3/4HX"] = 0.9999999999+v_["3/4HX"]
        v_["3/4HY"] = 0.9999999999+v_["3/4HY"]
        v_["1/4HX"] = trunc(Int,v_["1/4HX"])
        v_["1/4HY"] = trunc(Int,v_["1/4HY"])
        v_["3/4HX"] = trunc(Int,v_["3/4HX"])
        v_["3/4HY"] = trunc(Int,v_["3/4HY"])
        for I = Int64(v_["1/4HX"]):Int64(v_["3/4HX"])
            for J = Int64(v_["1/4HY"]):Int64(v_["3/4HY"])
                pb.xlower[ix_["V"*string(I)*","*string(J)]] = 1.0
            end
        end
        for J = Int64(v_["0"]):Int64(v_["NY+1"])
            pb.xlower[ix_["V"*string(Int64(v_["0"]))*","*string(J)]] = 0.0
            pb.xupper[ix_["V"*string(Int64(v_["0"]))*","*string(J)]] = 0.0
            pb.xlower[ix_["V"*string(Int64(v_["NX+1"]))*","*string(J)]] = 0.0
            pb.xupper[ix_["V"*string(Int64(v_["NX+1"]))*","*string(J)]] = 0.0
        end
        for I = Int64(v_["0"]):Int64(v_["NX+1"])
            v_["I"] = Float64(I)
            v_["VIJ"] = 2.0*I
            v_["VIJ"] = v_["VIJ"]*v_["HX"]
            v_["VIJ"] = -1.0+v_["VIJ"]
            v_["VIJ"] = v_["VIJ"]*v_["VIJ"]
            v_["VIJ"] = v_["ONE"]-v_["VIJ"]
            pb.xlower[ix_["V"*string(I)*","*string(Int64(v_["0"]))]] = v_["VIJ"]
            pb.xupper[ix_["V"*string(I)*","*string(Int64(v_["0"]))]] = v_["VIJ"]
            pb.xlower[ix_["V"*string(I)*","*string(Int64(v_["NY+1"]))]] = v_["VIJ"]
            pb.xupper[ix_["V"*string(I)*","*string(Int64(v_["NY+1"]))]] = v_["VIJ"]
        end
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        for I = Int64(v_["0"]):Int64(v_["NX+1"])
            v_["I"] = Float64(I)
            v_["VIJ"] = 2.0*I
            v_["VIJ"] = v_["VIJ"]*v_["HX"]
            v_["VIJ"] = -1.0+v_["VIJ"]
            v_["VIJ"] = v_["VIJ"]*v_["VIJ"]
            v_["VIJ"] = v_["ONE"]-v_["VIJ"]
            for J = Int64(v_["0"]):Int64(v_["NY+1"])
                pb.x0[ix_["V"*string(I)*","*string(J)]] = Float64(v_["VIJ"])
            end
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eISQ", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["0"]):Int64(v_["NX"])
            v_["I+1"] = 1+I
            for J = Int64(v_["0"]):Int64(v_["NY"])
                v_["J+1"] = 1+J
                ename = "I"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eISQ")
                arrset(ielftype,ie,iet_["eISQ"])
                vname = "V"*string(Int64(v_["I+1"]))*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "V"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                ename = "J"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eISQ")
                arrset(ielftype,ie,iet_["eISQ"])
                vname = "V"*string(I)*","*string(Int64(v_["J+1"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "V"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
        end
        for J = Int64(v_["0"]):Int64(v_["NY+1"])
            v_["J1"] = 1+J
            ename = "J"*string(Int64(v_["NX+1"]))*","*string(J)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eISQ")
            arrset(ielftype,ie,iet_["eISQ"])
            ename = "J"*string(Int64(v_["NX+1"]))*","*string(J)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "V"*string(Int64(v_["NX+1"]))*","*string(Int64(v_["J1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "J"*string(Int64(v_["NX+1"]))*","*string(J)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "V"*string(Int64(v_["NX+1"]))*","*string(J)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        for I = Int64(v_["0"]):Int64(v_["NX+1"])
            v_["I1"] = 1+I
            ename = "I"*string(I)*","*string(Int64(v_["NY+1"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eISQ")
            arrset(ielftype,ie,iet_["eISQ"])
            ename = "I"*string(I)*","*string(Int64(v_["NY+1"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "V"*string(Int64(v_["I1"]))*","*string(Int64(v_["NY+1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "I"*string(I)*","*string(Int64(v_["NY+1"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "V"*string(I)*","*string(Int64(v_["NY+1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gROOT",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["0"]):Int64(v_["NX"])
            for J = Int64(v_["0"]):Int64(v_["NY"])
                ig = ig_["A"*string(I)*","*string(J)]
                arrset(pbm.grftype,ig,"gROOT")
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["I"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel,Float64(v_["1/HX2"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["J"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel,Float64(v_["1/HY2"]))
            end
        end
        for I = Int64(v_["1"]):Int64(v_["NX+1"])
            v_["I-1"] = -1+I
            for J = Int64(v_["1"]):Int64(v_["NY+1"])
                v_["J-1"] = -1+J
                ig = ig_["B"*string(I)*","*string(J)]
                arrset(pbm.grftype,ig,"gROOT")
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["I"*string(Int64(v_["I-1"]))*","*string(J)])
                loaset(pbm.grelw,ig,posel,Float64(v_["1/HX2"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["J"*string(I)*","*string(Int64(v_["J-1"]))])
                loaset(pbm.grelw,ig,posel,Float64(v_["1/HY2"]))
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLUTION            2.51948D+00    $ (NX=50,NY=25)
# LO SOLUTION            2.51488D+00    $ (NX=50,NY=50)
# LO SOLUTION            2.50568D+00    $ (NX=50,NY=75)
# LO SOLUTION            2.50694D+00    $ (NX=50,NY=100)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.A = spzeros(Float64,0,0)
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-OBR2-AN-V-V"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        return pb, pbm

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "eISQ"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,1,2)
        IV_ =  zeros(Float64,1)
        U_[1,1] = U_[1,1]+1
        U_[1,2] = U_[1,2]-1
        IV_[1] = dot(U_[1,:],EV_)
        f_   = IV_[1]*IV_[1]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = IV_[1]+IV_[1]
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 2.0
                H_ = U_'*H_*U_
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

    elseif action == "gROOT"

        GVAR_   = args[1]
        igr_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        ROOTAL = sqrt(GVAR_)
        f_= ROOTAL
        if nargout>1
            g_ = 0.5/ROOTAL
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_ = -0.25/ROOTAL^3
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

