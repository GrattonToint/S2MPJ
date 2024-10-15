function ROBOT(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    This program solves the displacement optimization problem in 
#    REDUNDANT ROBOTS. A redundant robot is one which has more links than 
#    the dimensions it moves in.  Because this redundancy allows almost
#    infinite combinations of joint angles for a particular orientation of
#    end-effector of a robot, choosing an optimum combination has always been a
#    problem of research. 
#    The ROBOT considered here is a 7 link robot moving in 2 dimensional space.
# 
#    Source: an exercize for L. Watson course on LANCELOT in the Spring 1993.
#    B.Benhabib, R.G.Fenton and A.A.Goldberg, 
#    "Analytical trajectory optimization of seven degrees of freedom redundant
#    robot",  
#    Transactions of the Canadian Society for Mechanical Engineering,
#    vol.11(4), 1987, pp 197-200.
# 
#    SIF input: Manish Sabu at Virginia Tech., Spring 1993.
#               Minor modifications by Ph. L. Toint, April 1993.
# 
#    classification = "C-QOR2-MY-14-2"
# 
#  This segment describes the initial values of angles (by THnIN)
#   and final position of the end effector (by XPOS and YPOS)
#   these values can be changed here according to the needs of the user.
#  The segment also defines the upper and lower bounds of the various joint 
#   angles (by HIGH and DOWN)
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "ROBOT"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["TH1IN"] = 0.0
        v_["TH2IN"] = 0.0
        v_["TH3IN"] = 0.0
        v_["TH4IN"] = 0.0
        v_["TH5IN"] = 0.0
        v_["TH6IN"] = 0.0
        v_["TH7IN"] = 0.0
        v_["XPOS"] = 4.0
        v_["YPOS"] = 4.0
        v_["HIGH"] = 2.356194
        v_["DOWN"] = -2.356194
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        iv,ix_,_ = s2mpj_ii("TH1",ix_)
        arrset(pb.xnames,iv,"TH1")
        iv,ix_,_ = s2mpj_ii("TH2",ix_)
        arrset(pb.xnames,iv,"TH2")
        iv,ix_,_ = s2mpj_ii("TH3",ix_)
        arrset(pb.xnames,iv,"TH3")
        iv,ix_,_ = s2mpj_ii("TH4",ix_)
        arrset(pb.xnames,iv,"TH4")
        iv,ix_,_ = s2mpj_ii("TH5",ix_)
        arrset(pb.xnames,iv,"TH5")
        iv,ix_,_ = s2mpj_ii("TH6",ix_)
        arrset(pb.xnames,iv,"TH6")
        iv,ix_,_ = s2mpj_ii("TH7",ix_)
        arrset(pb.xnames,iv,"TH7")
        iv,ix_,_ = s2mpj_ii("TH1I",ix_)
        arrset(pb.xnames,iv,"TH1I")
        iv,ix_,_ = s2mpj_ii("TH2I",ix_)
        arrset(pb.xnames,iv,"TH2I")
        iv,ix_,_ = s2mpj_ii("TH3I",ix_)
        arrset(pb.xnames,iv,"TH3I")
        iv,ix_,_ = s2mpj_ii("TH4I",ix_)
        arrset(pb.xnames,iv,"TH4I")
        iv,ix_,_ = s2mpj_ii("TH5I",ix_)
        arrset(pb.xnames,iv,"TH5I")
        iv,ix_,_ = s2mpj_ii("TH6I",ix_)
        arrset(pb.xnames,iv,"TH6I")
        iv,ix_,_ = s2mpj_ii("TH7I",ix_)
        arrset(pb.xnames,iv,"TH7I")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        ig,ig_,_ = s2mpj_ii("CONSTR1",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"CONSTR1")
        ig,ig_,_ = s2mpj_ii("CONSTR2",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"CONSTR2")
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
        pbm.gconst[ig_["CONSTR1"]] = Float64(v_["XPOS"])
        pbm.gconst[ig_["CONSTR2"]] = Float64(v_["YPOS"])
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower[ix_["TH1"]] = v_["DOWN"]
        pb.xupper[ix_["TH1"]] = v_["HIGH"]
        pb.xlower[ix_["TH2"]] = v_["DOWN"]
        pb.xupper[ix_["TH2"]] = v_["HIGH"]
        pb.xlower[ix_["TH3"]] = v_["DOWN"]
        pb.xupper[ix_["TH3"]] = v_["HIGH"]
        pb.xlower[ix_["TH4"]] = v_["DOWN"]
        pb.xupper[ix_["TH4"]] = v_["HIGH"]
        pb.xlower[ix_["TH5"]] = v_["DOWN"]
        pb.xupper[ix_["TH5"]] = v_["HIGH"]
        pb.xlower[ix_["TH6"]] = v_["DOWN"]
        pb.xupper[ix_["TH6"]] = v_["HIGH"]
        pb.xlower[ix_["TH7"]] = v_["DOWN"]
        pb.xupper[ix_["TH7"]] = v_["HIGH"]
        pb.xlower[ix_["TH1I"]] = v_["TH1IN"]
        pb.xupper[ix_["TH1I"]] = v_["TH1IN"]
        pb.xlower[ix_["TH2I"]] = v_["TH2IN"]
        pb.xupper[ix_["TH2I"]] = v_["TH2IN"]
        pb.xlower[ix_["TH3I"]] = v_["TH3IN"]
        pb.xupper[ix_["TH3I"]] = v_["TH3IN"]
        pb.xlower[ix_["TH4I"]] = v_["TH4IN"]
        pb.xupper[ix_["TH4I"]] = v_["TH4IN"]
        pb.xlower[ix_["TH5I"]] = v_["TH5IN"]
        pb.xupper[ix_["TH5I"]] = v_["TH5IN"]
        pb.xlower[ix_["TH6I"]] = v_["TH6IN"]
        pb.xupper[ix_["TH6I"]] = v_["TH6IN"]
        pb.xlower[ix_["TH7I"]] = v_["TH7IN"]
        pb.xupper[ix_["TH7I"]] = v_["TH7IN"]
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        pb.x0[ix_["TH1"]] = Float64(0.0)
        pb.x0[ix_["TH2"]] = Float64(0.0)
        pb.x0[ix_["TH3"]] = Float64(0.0)
        pb.x0[ix_["TH4"]] = Float64(0.0)
        pb.x0[ix_["TH5"]] = Float64(0.0)
        pb.x0[ix_["TH6"]] = Float64(0.0)
        pb.x0[ix_["TH7"]] = Float64(0.0)
        pb.x0[ix_["TH1I"]] = Float64(0.0)
        pb.x0[ix_["TH2I"]] = Float64(0.0)
        pb.x0[ix_["TH3I"]] = Float64(0.0)
        pb.x0[ix_["TH4I"]] = Float64(0.0)
        pb.x0[ix_["TH5I"]] = Float64(0.0)
        pb.x0[ix_["TH6I"]] = Float64(0.0)
        pb.x0[ix_["TH7I"]] = Float64(0.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eISQ", iet_)
        loaset(elftv,it,1,"V")
        loaset(elftv,it,2,"W")
        it,iet_,_ = s2mpj_ii( "eCOSTH", iet_)
        loaset(elftv,it,1,"THETAC")
        it,iet_,_ = s2mpj_ii( "eSINTH", iet_)
        loaset(elftv,it,1,"THETAS")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "TH1SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eISQ")
        arrset(ielftype,ie,iet_["eISQ"])
        vname = "TH1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TH1I"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="W",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "TH2SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eISQ")
        arrset(ielftype,ie,iet_["eISQ"])
        vname = "TH2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TH2I"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="W",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "TH3SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eISQ")
        arrset(ielftype,ie,iet_["eISQ"])
        vname = "TH3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TH3I"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="W",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "TH4SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eISQ")
        arrset(ielftype,ie,iet_["eISQ"])
        vname = "TH4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TH4I"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="W",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "TH5SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eISQ")
        arrset(ielftype,ie,iet_["eISQ"])
        vname = "TH5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TH5I"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="W",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "TH6SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eISQ")
        arrset(ielftype,ie,iet_["eISQ"])
        vname = "TH6"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TH6I"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="W",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "TH7SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eISQ")
        arrset(ielftype,ie,iet_["eISQ"])
        vname = "TH7"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TH7I"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="W",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "C1TH"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eCOSTH")
        arrset(ielftype,ie,iet_["eCOSTH"])
        vname = "TH1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="THETAC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "C2TH"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eCOSTH")
        arrset(ielftype,ie,iet_["eCOSTH"])
        vname = "TH2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="THETAC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "C3TH"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eCOSTH")
        arrset(ielftype,ie,iet_["eCOSTH"])
        vname = "TH3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="THETAC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "C4TH"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eCOSTH")
        arrset(ielftype,ie,iet_["eCOSTH"])
        vname = "TH4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="THETAC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "C5TH"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eCOSTH")
        arrset(ielftype,ie,iet_["eCOSTH"])
        vname = "TH5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="THETAC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "C6TH"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eCOSTH")
        arrset(ielftype,ie,iet_["eCOSTH"])
        vname = "TH6"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="THETAC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "C7TH"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eCOSTH")
        arrset(ielftype,ie,iet_["eCOSTH"])
        vname = "TH7"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="THETAC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "S1TH"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSINTH")
        arrset(ielftype,ie,iet_["eSINTH"])
        vname = "TH1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="THETAS",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "S2TH"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSINTH")
        arrset(ielftype,ie,iet_["eSINTH"])
        vname = "TH2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="THETAS",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "S3TH"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSINTH")
        arrset(ielftype,ie,iet_["eSINTH"])
        vname = "TH3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="THETAS",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "S4TH"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSINTH")
        arrset(ielftype,ie,iet_["eSINTH"])
        vname = "TH4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="THETAS",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "S5TH"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSINTH")
        arrset(ielftype,ie,iet_["eSINTH"])
        vname = "TH5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="THETAS",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "S6TH"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSINTH")
        arrset(ielftype,ie,iet_["eSINTH"])
        vname = "TH6"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="THETAS",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "S7TH"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSINTH")
        arrset(ielftype,ie,iet_["eSINTH"])
        vname = "TH7"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="THETAS",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["OBJ"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["TH1SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["TH2SQ"])
        loaset(pbm.grelw,ig,posel, 1.)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["TH3SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["TH4SQ"])
        loaset(pbm.grelw,ig,posel, 1.)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["TH5SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["TH6SQ"])
        loaset(pbm.grelw,ig,posel, 1.)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["TH7SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["CONSTR1"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["C1TH"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["C2TH"])
        loaset(pbm.grelw,ig,posel, 1.)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["C3TH"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["C4TH"])
        loaset(pbm.grelw,ig,posel, 1.)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["C5TH"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["C6TH"])
        loaset(pbm.grelw,ig,posel, 1.)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["C7TH"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(0.5))
        ig = ig_["CONSTR2"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["S1TH"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["S2TH"])
        loaset(pbm.grelw,ig,posel, 1.)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["S3TH"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["S4TH"])
        loaset(pbm.grelw,ig,posel, 1.)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["S5TH"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["S6TH"])
        loaset(pbm.grelw,ig,posel, 1.)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["S7TH"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(0.5))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
# LO SOLUTION            5.46283877
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
        pb.pbclass = "C-QOR2-MY-14-2"
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

    elseif action == "eCOSTH"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        TMP = cos(EV_[1])
        f_   = TMP
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = -sin(EV_[1])
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = -TMP
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eSINTH"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        TMP = sin(EV_[1])
        f_   = TMP
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = cos(EV_[1])
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = -TMP
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

