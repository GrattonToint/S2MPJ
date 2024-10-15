function MAXLIKA(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : MAXLIKA
#    *********
# 
#    A variant of Hock and Schittkowski problem 105, where the
#    (inactive) inequality constraint is dropped.
# 
#    Source:
#    Ph. Toint and A. Griewank.
# 
#    SIF input: Ph. Toint, June 1990.
# 
#    classification = "C-OBR2-AY-8-0"
# 
#    Constants
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "MAXLIKA"

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
        v_["235"] = 235
        v_["Y1"] = 95.0
        v_["Y2"] = 105.0
        v_["Y3"] = 110.0
        v_["Y4"] = 110.0
        v_["Y5"] = 110.0
        v_["Y6"] = 110.0
        v_["Y7"] = 115.0
        v_["Y8"] = 115.0
        v_["Y9"] = 115.0
        v_["Y10"] = 115.0
        v_["Y11"] = 120.0
        v_["Y12"] = 120.0
        v_["Y13"] = 120.0
        v_["Y14"] = 120.0
        v_["Y15"] = 120.0
        v_["Y16"] = 120.0
        v_["Y17"] = 120.0
        v_["Y18"] = 120.0
        v_["Y19"] = 120.0
        v_["Y20"] = 120.0
        v_["Y21"] = 120.0
        v_["Y22"] = 120.0
        v_["Y23"] = 120.0
        v_["Y24"] = 120.0
        v_["Y25"] = 120.0
        v_["Y26"] = 125.0
        v_["Y27"] = 125.0
        v_["Y28"] = 125.0
        v_["Y29"] = 125.0
        v_["Y30"] = 125.0
        v_["Y31"] = 125.0
        v_["Y32"] = 125.0
        v_["Y33"] = 125.0
        v_["Y34"] = 125.0
        v_["Y35"] = 125.0
        v_["Y36"] = 125.0
        v_["Y37"] = 125.0
        v_["Y38"] = 125.0
        v_["Y39"] = 125.0
        v_["Y40"] = 125.0
        v_["Y41"] = 130.0
        v_["Y42"] = 130.0
        v_["Y43"] = 130.0
        v_["Y44"] = 130.0
        v_["Y45"] = 130.0
        v_["Y46"] = 130.0
        v_["Y47"] = 130.0
        v_["Y48"] = 130.0
        v_["Y49"] = 130.0
        v_["Y50"] = 130.0
        v_["Y51"] = 130.0
        v_["Y52"] = 130.0
        v_["Y53"] = 130.0
        v_["Y54"] = 130.0
        v_["Y55"] = 130.0
        v_["Y56"] = 135.0
        v_["Y57"] = 135.0
        v_["Y58"] = 135.0
        v_["Y59"] = 135.0
        v_["Y60"] = 135.0
        v_["Y61"] = 135.0
        v_["Y62"] = 135.0
        v_["Y63"] = 135.0
        v_["Y64"] = 135.0
        v_["Y65"] = 135.0
        v_["Y66"] = 135.0
        v_["Y67"] = 135.0
        v_["Y68"] = 135.0
        v_["Y69"] = 140.0
        v_["Y70"] = 140.0
        v_["Y71"] = 140.0
        v_["Y72"] = 140.0
        v_["Y73"] = 140.0
        v_["Y74"] = 140.0
        v_["Y75"] = 140.0
        v_["Y76"] = 140.0
        v_["Y77"] = 140.0
        v_["Y78"] = 140.0
        v_["Y79"] = 140.0
        v_["Y80"] = 140.0
        v_["Y81"] = 140.0
        v_["Y82"] = 140.0
        v_["Y83"] = 140.0
        v_["Y84"] = 140.0
        v_["Y85"] = 140.0
        v_["Y86"] = 140.0
        v_["Y87"] = 140.0
        v_["Y88"] = 140.0
        v_["Y89"] = 140.0
        v_["Y90"] = 145.0
        v_["Y91"] = 145.0
        v_["Y92"] = 145.0
        v_["Y93"] = 145.0
        v_["Y94"] = 145.0
        v_["Y95"] = 145.0
        v_["Y96"] = 145.0
        v_["Y97"] = 145.0
        v_["Y98"] = 145.0
        v_["Y99"] = 145.0
        v_["Y100"] = 145.0
        v_["Y101"] = 145.0
        v_["Y102"] = 150.0
        v_["Y103"] = 150.0
        v_["Y104"] = 150.0
        v_["Y105"] = 150.0
        v_["Y106"] = 150.0
        v_["Y107"] = 150.0
        v_["Y108"] = 150.0
        v_["Y109"] = 150.0
        v_["Y110"] = 150.0
        v_["Y111"] = 150.0
        v_["Y112"] = 150.0
        v_["Y113"] = 150.0
        v_["Y114"] = 150.0
        v_["Y115"] = 150.0
        v_["Y116"] = 150.0
        v_["Y117"] = 150.0
        v_["Y118"] = 150.0
        v_["Y119"] = 155.0
        v_["Y120"] = 155.0
        v_["Y121"] = 155.0
        v_["Y122"] = 155.0
        v_["Y123"] = 160.0
        v_["Y124"] = 160.0
        v_["Y125"] = 160.0
        v_["Y126"] = 160.0
        v_["Y127"] = 160.0
        v_["Y128"] = 160.0
        v_["Y129"] = 160.0
        v_["Y130"] = 160.0
        v_["Y131"] = 160.0
        v_["Y132"] = 160.0
        v_["Y133"] = 160.0
        v_["Y134"] = 160.0
        v_["Y135"] = 160.0
        v_["Y136"] = 160.0
        v_["Y137"] = 160.0
        v_["Y138"] = 160.0
        v_["Y139"] = 160.0
        v_["Y140"] = 160.0
        v_["Y141"] = 160.0
        v_["Y142"] = 160.0
        v_["Y143"] = 165.0
        v_["Y144"] = 165.0
        v_["Y145"] = 165.0
        v_["Y146"] = 165.0
        v_["Y147"] = 165.0
        v_["Y148"] = 165.0
        v_["Y149"] = 165.0
        v_["Y150"] = 165.0
        v_["Y151"] = 170.0
        v_["Y152"] = 170.0
        v_["Y153"] = 170.0
        v_["Y154"] = 170.0
        v_["Y155"] = 170.0
        v_["Y156"] = 170.0
        v_["Y157"] = 170.0
        v_["Y158"] = 170.0
        v_["Y159"] = 170.0
        v_["Y160"] = 170.0
        v_["Y161"] = 170.0
        v_["Y162"] = 170.0
        v_["Y163"] = 170.0
        v_["Y164"] = 170.0
        v_["Y165"] = 170.0
        v_["Y166"] = 170.0
        v_["Y167"] = 170.0
        v_["Y168"] = 175.0
        v_["Y169"] = 175.0
        v_["Y170"] = 175.0
        v_["Y171"] = 175.0
        v_["Y172"] = 175.0
        v_["Y173"] = 175.0
        v_["Y174"] = 175.0
        v_["Y175"] = 175.0
        v_["Y176"] = 180.0
        v_["Y177"] = 180.0
        v_["Y178"] = 180.0
        v_["Y179"] = 180.0
        v_["Y180"] = 180.0
        v_["Y181"] = 180.0
        v_["Y182"] = 185.0
        v_["Y183"] = 185.0
        v_["Y184"] = 185.0
        v_["Y185"] = 185.0
        v_["Y186"] = 185.0
        v_["Y187"] = 185.0
        v_["Y188"] = 190.0
        v_["Y189"] = 190.0
        v_["Y190"] = 190.0
        v_["Y191"] = 190.0
        v_["Y192"] = 190.0
        v_["Y193"] = 190.0
        v_["Y194"] = 190.0
        v_["Y195"] = 195.0
        v_["Y196"] = 195.0
        v_["Y197"] = 195.0
        v_["Y198"] = 195.0
        v_["Y199"] = 200.0
        v_["Y200"] = 200.0
        v_["Y201"] = 200.0
        v_["Y202"] = 205.0
        v_["Y203"] = 205.0
        v_["Y204"] = 205.0
        v_["Y205"] = 210.0
        v_["Y206"] = 210.0
        v_["Y207"] = 210.0
        v_["Y208"] = 210.0
        v_["Y209"] = 210.0
        v_["Y210"] = 210.0
        v_["Y211"] = 210.0
        v_["Y212"] = 210.0
        v_["Y213"] = 215.0
        v_["Y214"] = 220.0
        v_["Y215"] = 220.0
        v_["Y216"] = 220.0
        v_["Y217"] = 220.0
        v_["Y218"] = 220.0
        v_["Y219"] = 220.0
        v_["Y220"] = 230.0
        v_["Y221"] = 230.0
        v_["Y222"] = 230.0
        v_["Y223"] = 230.0
        v_["Y224"] = 230.0
        v_["Y225"] = 235.0
        v_["Y226"] = 240.0
        v_["Y227"] = 240.0
        v_["Y228"] = 240.0
        v_["Y229"] = 240.0
        v_["Y230"] = 240.0
        v_["Y231"] = 240.0
        v_["Y232"] = 240.0
        v_["Y233"] = 245.0
        v_["Y234"] = 250.0
        v_["Y235"] = 250.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        iv,ix_,_ = s2mpj_ii("X1",ix_)
        arrset(pb.xnames,iv,"X1")
        iv,ix_,_ = s2mpj_ii("X2",ix_)
        arrset(pb.xnames,iv,"X2")
        iv,ix_,_ = s2mpj_ii("X3",ix_)
        arrset(pb.xnames,iv,"X3")
        iv,ix_,_ = s2mpj_ii("X4",ix_)
        arrset(pb.xnames,iv,"X4")
        iv,ix_,_ = s2mpj_ii("X5",ix_)
        arrset(pb.xnames,iv,"X5")
        iv,ix_,_ = s2mpj_ii("X6",ix_)
        arrset(pb.xnames,iv,"X6")
        iv,ix_,_ = s2mpj_ii("X7",ix_)
        arrset(pb.xnames,iv,"X7")
        iv,ix_,_ = s2mpj_ii("X8",ix_)
        arrset(pb.xnames,iv,"X8")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["1"]):Int64(v_["235"])
            ig,ig_,_ = s2mpj_ii("L"*string(I),ig_)
            arrset(gtype,ig,"<>")
            arrset(pbm.gscale,ig,Float64(-1.0))
        end
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower[ix_["X1"]] = 0.001
        pb.xupper[ix_["X1"]] = 0.499
        pb.xlower[ix_["X2"]] = 0.001
        pb.xupper[ix_["X2"]] = 0.499
        pb.xlower[ix_["X3"]] = 100.0
        pb.xupper[ix_["X3"]] = 180.0
        pb.xlower[ix_["X4"]] = 130.0
        pb.xupper[ix_["X4"]] = 210.0
        pb.xlower[ix_["X5"]] = 170.0
        pb.xupper[ix_["X5"]] = 240.0
        pb.xlower[ix_["X6"]] = 5.0
        pb.xupper[ix_["X6"]] = 25.0
        pb.xlower[ix_["X7"]] = 5.0
        pb.xupper[ix_["X7"]] = 25.0
        pb.xlower[ix_["X8"]] = 5.0
        pb.xupper[ix_["X8"]] = 25.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.x0[ix_["X1"]] = Float64(0.1)
        pb.x0[ix_["X2"]] = Float64(0.2)
        pb.x0[ix_["X3"]] = Float64(100.0)
        pb.x0[ix_["X4"]] = Float64(125.0)
        pb.x0[ix_["X5"]] = Float64(175.0)
        pb.x0[ix_["X6"]] = Float64(11.2)
        pb.x0[ix_["X7"]] = Float64(13.2)
        pb.x0[ix_["X8"]] = Float64(15.8)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eAB", iet_)
        loaset(elftv,it,1,"U")
        loaset(elftv,it,2,"V")
        loaset(elftv,it,3,"W")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"Y")
        it,iet_,_ = s2mpj_ii( "eC", iet_)
        loaset(elftv,it,1,"U")
        loaset(elftv,it,2,"X")
        loaset(elftv,it,3,"S")
        loaset(elftv,it,4,"T")
        loaset(elftp,it,1,"Y")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["235"])
            ename = "A"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eAB")
            arrset(ielftype,ie,iet_["eAB"])
            vname = "X1"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="U",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X6"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X3"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="W",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="Y",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["Y"*string(I)]))
            ename = "B"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eAB")
            arrset(ielftype,ie,iet_["eAB"])
            vname = "X2"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="U",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X7"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X4"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="W",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="Y",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["Y"*string(I)]))
            ename = "C"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eC")
            arrset(ielftype,ie,iet_["eC"])
            vname = "X2"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="U",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X1"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X8"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="S",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X5"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="T",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="Y",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["Y"*string(I)]))
        end
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gLN",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["235"])
            ig = ig_["L"*string(I)]
            arrset(pbm.grftype,ig,"gLN")
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["A"*string(I)])
            loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["B"*string(I)])
            loaset(pbm.grelw,ig,posel, 1.)
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["C"*string(I)])
            loaset(pbm.grelw,ig,posel,1.)
        end
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.A = spzeros(Float64,0,0)
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-OBR2-AY-8-0"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        return pb, pbm


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "eAB"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        YMW = pbm.elpar[iel_][1]-EV_[3]
        YMWSQ = YMW*YMW
        VSQ = EV_[2]*EV_[2]
        VCB = VSQ*EV_[2]
        A = -YMWSQ/(2.0*VSQ)
        DADV = YMWSQ/VCB
        DADW = YMW/VSQ
        D2ADV2 = -3.0*YMWSQ/(VSQ*VSQ)
        D2ADVW = -2.0*YMW/VCB
        D2ADW2 = -1.0/VSQ
        E = exp(A)
        DEDV = E*DADV
        DEDW = E*DADW
        B = EV_[1]*E
        DBDV = B*DADV
        DBDW = B*DADW
        D2BDV2 = DBDV*DADV+B*D2ADV2
        D2BDVW = DBDW*DADV+B*D2ADVW
        D2BDW2 = DBDW*DADW+B*D2ADW2
        f_   = B/EV_[2]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = E/EV_[2]
            g_[2] = (DBDV-B/EV_[2])/EV_[2]
            g_[3] = DBDW/EV_[2]
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,2] = (DEDV-E/EV_[2])/EV_[2]
                H_[2,1] = H_[1,2]
                H_[1,3] = DEDW/EV_[2]
                H_[3,1] = H_[1,3]
                H_[2,2] = (D2BDV2-DBDV/EV_[2]+B/VSQ)/EV_[2]-(DBDV-B/EV_[2])/VSQ
                H_[2,3] = (D2BDVW-DBDW/EV_[2])/EV_[2]
                H_[3,2] = H_[2,3]
                H_[3,3] = D2BDW2/EV_[2]
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eC"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,3,4)
        IV_ =  zeros(Float64,3)
        U_[1,1] = U_[1,1]+1
        U_[1,2] = U_[1,2]+1
        U_[2,3] = U_[2,3]+1
        U_[3,4] = U_[3,4]+1
        IV_[1] = dot(U_[1,:],EV_)
        IV_[2] = dot(U_[2,:],EV_)
        IV_[3] = dot(U_[3,:],EV_)
        YMW = pbm.elpar[iel_][1]-IV_[3]
        YMWSQ = YMW*YMW
        VSQ = IV_[2]*IV_[2]
        VCB = VSQ*IV_[2]
        A = -YMWSQ/(2.0*VSQ)
        DADV = YMWSQ/VCB
        DADW = YMW/VSQ
        D2ADV2 = -3.0*YMWSQ/(VSQ*VSQ)
        D2ADVW = -2.0*YMW/VCB
        D2ADW2 = -1.0/VSQ
        E = exp(A)
        DEDV = E*DADV
        DEDW = E*DADW
        B = (1.0-IV_[1])*E
        DBDZ = -E
        DBDV = B*DADV
        DBDW = B*DADW
        D2BDVZ = -DEDV
        D2BDWZ = -DEDW
        D2BDV2 = DBDV*DADV+B*D2ADV2
        D2BDVW = DBDW*DADV+B*D2ADVW
        D2BDW2 = DBDW*DADW+B*D2ADW2
        f_   = B/IV_[2]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = DBDZ/IV_[2]
            g_[2] = (DBDV-B/IV_[2])/IV_[2]
            g_[3] = DBDW/IV_[2]
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,2] = (DEDV-E/IV_[2])/IV_[2]
                H_[2,1] = H_[1,2]
                H_[1,3] = D2BDWZ/IV_[2]
                H_[3,1] = H_[1,3]
                H_[2,2] = (D2BDV2-DBDV/IV_[2]+B/VSQ)/IV_[2]-(DBDV-B/IV_[2])/VSQ
                H_[2,3] = (D2BDVW-DBDW/IV_[2])/IV_[2]
                H_[3,2] = H_[2,3]
                H_[3,3] = D2BDW2/IV_[2]
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

    elseif action == "gLN"

        GVAR_   = args[1]
        igr_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_= log(GVAR_*0.39894228)
        if nargout>1
            g_ = 1.0/GVAR_
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_ = -1.0/GVAR_^2
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

