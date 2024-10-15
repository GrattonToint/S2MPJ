function VIBRBEAM(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    A nonlinear least-squares problem arising from laser-Doppler
#    measurements of a vibrating beam.  The data correspond to a simulated
#    experiment where two laser-Doppler velocimeters take measurements
#    at random points along the centreline of the beam.  These measurements
#    consist of a position (x), an incident angle (p) and the magnitude
#    of the velocity along the line of sight (v).
#    The problem is then to fit
# 
#                          2      3                    2     3
#        v = (c + c x + c x  + c x ) cos[ d + d x + d x + d x  - p ]
#              0   1     2      3          0   1     2     3
#            <---- magnitude ----->       <------ phase ----->
# 
#    in the least-squares sense.
# 
#    Source: 
#    a modification of an exercize for L. Watson course on LANCELOT in
#    the Spring 1993. Compared to the original proposal, the unnecessary
#    elements were removed as well as an unnecessary constraint on the phase.
# 
#    SIF input: Ph. L. Toint, May 1993, based on a proposal by
#               D. E. Montgomery, Virginia Tech., April 1993.
# 
#    classification = "C-SUR2-MN-8-0"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "VIBRBEAM"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["0"] = 0
        v_["1"] = 1
        v_["3"] = 3
        v_["m"] = 30
        v_["x1"] = 39.1722
        v_["x2"] = 53.9707
        v_["x3"] = 47.9829
        v_["x4"] = 12.5925
        v_["x5"] = 16.5414
        v_["x6"] = 18.9548
        v_["x7"] = 27.7168
        v_["x8"] = 31.9201
        v_["x9"] = 45.6830
        v_["x10"] = 22.2524
        v_["x11"] = 33.9805
        v_["x12"] = 6.8425
        v_["x13"] = 35.1677
        v_["x14"] = 33.5682
        v_["x15"] = 43.3659
        v_["x16"] = 13.3835
        v_["x17"] = 25.7273
        v_["x18"] = 21.0230
        v_["x19"] = 10.9755
        v_["x20"] = 1.5323
        v_["x21"] = 45.4416
        v_["x22"] = 14.5431
        v_["x23"] = 22.4313
        v_["x24"] = 29.0144
        v_["x25"] = 25.2675
        v_["x26"] = 15.5095
        v_["x27"] = 9.6297
        v_["x28"] = 8.3009
        v_["x29"] = 30.8694
        v_["x30"] = 43.3299
        v_["v1"] = -1.2026
        v_["v2"] = 1.7053
        v_["v3"] = 0.5410
        v_["v4"] = 1.1477
        v_["v5"] = 1.2447
        v_["v6"] = 0.9428
        v_["v7"] = -0.1360
        v_["v8"] = -0.7542
        v_["v9"] = -0.3396
        v_["v10"] = 0.7057
        v_["v11"] = -0.8509
        v_["v12"] = -0.1201
        v_["v13"] = -1.2193
        v_["v14"] = -1.0448
        v_["v15"] = -0.7723
        v_["v16"] = 0.4342
        v_["v17"] = 0.1154
        v_["v18"] = 0.2868
        v_["v19"] = 0.3558
        v_["v20"] = -0.5090
        v_["v21"] = -0.0842
        v_["v22"] = 0.6021
        v_["v23"] = 0.1197
        v_["v24"] = -0.1827
        v_["v25"] = 0.1806
        v_["v26"] = 0.5395
        v_["v27"] = 0.2072
        v_["v28"] = 0.1466
        v_["v29"] = -0.2672
        v_["v30"] = -0.3038
        v_["p1"] = 2.5736
        v_["p2"] = 2.7078
        v_["p3"] = 2.6613
        v_["p4"] = 2.0374
        v_["p5"] = 2.1553
        v_["p6"] = 2.2195
        v_["p7"] = 2.4077
        v_["p8"] = 2.4772
        v_["p9"] = 2.6409
        v_["p10"] = 2.2981
        v_["p11"] = 2.5073
        v_["p12"] = 1.8380
        v_["p13"] = 2.5236
        v_["p14"] = 2.5015
        v_["p15"] = 2.6186
        v_["p16"] = 0.4947
        v_["p17"] = 0.6062
        v_["p18"] = 0.5588
        v_["p19"] = 0.4772
        v_["p20"] = 0.4184
        v_["p21"] = 0.9051
        v_["p22"] = 0.5035
        v_["p23"] = 0.5723
        v_["p24"] = 0.6437
        v_["p25"] = 0.6013
        v_["p26"] = 0.5111
        v_["p27"] = 0.4679
        v_["p28"] = 0.4590
        v_["p29"] = 0.6666
        v_["p30"] = 0.8630
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        iv,ix_,_ = s2mpj_ii("c0",ix_)
        arrset(pb.xnames,iv,"c0")
        iv,ix_,_ = s2mpj_ii("c1",ix_)
        arrset(pb.xnames,iv,"c1")
        iv,ix_,_ = s2mpj_ii("c2",ix_)
        arrset(pb.xnames,iv,"c2")
        iv,ix_,_ = s2mpj_ii("c3",ix_)
        arrset(pb.xnames,iv,"c3")
        iv,ix_,_ = s2mpj_ii("d0",ix_)
        arrset(pb.xnames,iv,"d0")
        iv,ix_,_ = s2mpj_ii("d1",ix_)
        arrset(pb.xnames,iv,"d1")
        iv,ix_,_ = s2mpj_ii("d2",ix_)
        arrset(pb.xnames,iv,"d2")
        iv,ix_,_ = s2mpj_ii("d3",ix_)
        arrset(pb.xnames,iv,"d3")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for i = Int64(v_["1"]):Int64(v_["m"])
            ig,ig_,_ = s2mpj_ii("f"*string(i),ig_)
            arrset(gtype,ig,"<>")
        end
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        for i = Int64(v_["1"]):Int64(v_["m"])
            pbm.gconst[ig_["f"*string(i)]] = Float64(v_["v"*string(i)])
        end
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower[ix_["c0"]] = -Inf
        pb.xupper[ix_["c0"]] = +Inf
        pb.xlower[ix_["c1"]] = -Inf
        pb.xupper[ix_["c1"]] = +Inf
        pb.xlower[ix_["c2"]] = -Inf
        pb.xupper[ix_["c2"]] = +Inf
        pb.xlower[ix_["c3"]] = -Inf
        pb.xupper[ix_["c3"]] = +Inf
        pb.xlower[ix_["d0"]] = -Inf
        pb.xupper[ix_["d0"]] = +Inf
        pb.xlower[ix_["d1"]] = -Inf
        pb.xupper[ix_["d1"]] = +Inf
        pb.xlower[ix_["d2"]] = -Inf
        pb.xupper[ix_["d2"]] = +Inf
        pb.xlower[ix_["d3"]] = -Inf
        pb.xupper[ix_["d3"]] = +Inf
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.x0[ix_["c0"]] = Float64(-3.5)
        pb.x0[ix_["c1"]] = Float64(1.0)
        pb.x0[ix_["d0"]] = Float64(1.7)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "efun", iet_)
        loaset(elftv,it,1,"a0")
        loaset(elftv,it,2,"a1")
        loaset(elftv,it,3,"a2")
        loaset(elftv,it,4,"a3")
        loaset(elftv,it,5,"b")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"y")
        loaset(elftp,it,2,"q")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for i = Int64(v_["1"]):Int64(v_["m"])
            for j = Int64(v_["0"]):Int64(v_["3"])
                ename = "fu"*string(i)*","*string(j)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"efun")
                arrset(ielftype,ie,iet_["efun"])
                vname = "d0"
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="a0",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "d1"
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="a1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "d2"
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="a2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "d3"
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="a3",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "c"*string(j)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="b",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                posep = findfirst(x->x=="y",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["x"*string(i)]))
                posep = findfirst(x->x=="q",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["p"*string(i)]))
            end
        end
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gsquare",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for i = Int64(v_["1"]):Int64(v_["m"])
            ig = ig_["f"*string(i)]
            arrset(pbm.grftype,ig,"gsquare")
            v_["y"] = 1.0
            for j = Int64(v_["0"]):Int64(v_["3"])
                ig = ig_["f"*string(i)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["fu"*string(i)*","*string(j)])
                loaset(pbm.grelw,ig,posel,Float64(v_["y"]))
                v_["y"] = v_["y"]*v_["x"*string(i)]
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
# LO SOLUTION             0.15644607137
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.A = spzeros(Float64,0,0)
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-SUR2-MN-8-0"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        return pb, pbm

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "efun"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        y2 = pbm.elpar[iel_][1]*pbm.elpar[iel_][1]
        y3 = pbm.elpar[iel_][1]*y2
        y4 = y2*y2
        y5 = y2*y3
        y6 = y3*y3
        phi  = (
              EV_[1]+pbm.elpar[iel_][1]*(EV_[2]+pbm.elpar[iel_][1]*(EV_[3]+pbm.elpar[iel_][1]*EV_[4]))-pbm.elpar[iel_][2])
        cosphi = cos(phi)
        sinphi = sin(phi)
        bcos = EV_[5]*cosphi
        bsin = EV_[5]*sinphi
        f_   = bcos
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = -bsin
            g_[2] = -bsin*pbm.elpar[iel_][1]
            g_[3] = -bsin*y2
            g_[4] = -bsin*y3
            g_[5] = cosphi
            if nargout>2
                H_ = zeros(Float64,5,5)
                H_[1,1] = -bcos
                H_[1,2] = -bcos*pbm.elpar[iel_][1]
                H_[2,1] = H_[1,2]
                H_[1,3] = -bcos*y2
                H_[3,1] = H_[1,3]
                H_[1,4] = -bcos*y3
                H_[4,1] = H_[1,4]
                H_[1,5] = -sinphi
                H_[5,1] = H_[1,5]
                H_[2,2] = -bcos*y2
                H_[2,3] = -bcos*y3
                H_[3,2] = H_[2,3]
                H_[2,4] = -bcos*y4
                H_[4,2] = H_[2,4]
                H_[2,5] = -sinphi*pbm.elpar[iel_][1]
                H_[5,2] = H_[2,5]
                H_[3,3] = -bcos*y4
                H_[3,4] = -bcos*y5
                H_[4,3] = H_[3,4]
                H_[3,5] = -sinphi*y2
                H_[5,3] = H_[3,5]
                H_[4,4] = -bcos*y6
                H_[4,5] = -sinphi*y3
                H_[5,4] = H_[4,5]
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

    elseif action == "gsquare"

        GVAR_   = args[1]
        igr_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_= GVAR_*GVAR_
        if nargout>1
            g_ = GVAR_+GVAR_
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_ = 2.0
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

