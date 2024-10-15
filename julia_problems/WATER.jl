function WATER(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    A small nonlinear network problem.
#    The problem is to compute the flows in a water distribution network
#    with 7 nodes and 8 links, subject to known supply/demand at the nodes 
#    and a unique reservoir at node 1.
# 
#    The problem is convex.
# 
#    Source:
#    an exercize for L. Watson course on LANCELOT in the Spring 1993.
# 
#    SIF input: E. P. Smith, Virginia Tech., Spring 1993.
# 
#    classification = "C-ONR2-MN-31-10"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "WATER"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("obj0102",ig_)
        arrset(gtype,ig,"<>")
        arrset(pbm.gscale,ig,Float64(105665.6))
        ig,ig_,_ = s2mpj_ii("obj0203",ig_)
        arrset(gtype,ig,"<>")
        arrset(pbm.gscale,ig,Float64(3613.412))
        ig,ig_,_ = s2mpj_ii("obj0204",ig_)
        arrset(gtype,ig,"<>")
        arrset(pbm.gscale,ig,Float64(105665.6))
        ig,ig_,_ = s2mpj_ii("obj0305",ig_)
        arrset(gtype,ig,"<>")
        arrset(pbm.gscale,ig,Float64(890.1553))
        ig,ig_,_ = s2mpj_ii("obj0405",ig_)
        arrset(gtype,ig,"<>")
        arrset(pbm.gscale,ig,Float64(76.66088))
        ig,ig_,_ = s2mpj_ii("obj0406",ig_)
        arrset(gtype,ig,"<>")
        arrset(pbm.gscale,ig,Float64(55145.82))
        ig,ig_,_ = s2mpj_ii("obj0607",ig_)
        arrset(gtype,ig,"<>")
        arrset(pbm.gscale,ig,Float64(26030.46))
        ig,ig_,_ = s2mpj_ii("obj0705",ig_)
        arrset(gtype,ig,"<>")
        arrset(pbm.gscale,ig,Float64(890.1553))
        ig,ig_,_ = s2mpj_ii("obj",ig_)
        arrset(gtype,ig,"<>")
        ig,ig_,_ = s2mpj_ii("c1",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"c1")
        ig,ig_,_ = s2mpj_ii("c2",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"c2")
        ig,ig_,_ = s2mpj_ii("c3",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"c3")
        ig,ig_,_ = s2mpj_ii("c4",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"c4")
        ig,ig_,_ = s2mpj_ii("c5",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"c5")
        ig,ig_,_ = s2mpj_ii("c6",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"c6")
        ig,ig_,_ = s2mpj_ii("c7",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"c7")
        ig,ig_,_ = s2mpj_ii("c8",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"c8")
        ig,ig_,_ = s2mpj_ii("c9",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"c9")
        ig,ig_,_ = s2mpj_ii("c10",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"c10")
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        ngrp   = length(ig_)
        iv,ix_,_ = s2mpj_ii("Q0102",ix_)
        arrset(pb.xnames,iv,"Q0102")
        ig = ig_["obj0102"]
        pbm.A[ig,iv] += Float64(1)
        ig = ig_["c1"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("Q0102",ix_)
        arrset(pb.xnames,iv,"Q0102")
        ig = ig_["c2"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("Q0203",ix_)
        arrset(pb.xnames,iv,"Q0203")
        ig = ig_["obj0203"]
        pbm.A[ig,iv] += Float64(1)
        ig = ig_["c2"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("Q0203",ix_)
        arrset(pb.xnames,iv,"Q0203")
        ig = ig_["c3"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("Q0204",ix_)
        arrset(pb.xnames,iv,"Q0204")
        ig = ig_["obj0204"]
        pbm.A[ig,iv] += Float64(1)
        ig = ig_["c2"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("Q0204",ix_)
        arrset(pb.xnames,iv,"Q0204")
        ig = ig_["c4"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("Q0305",ix_)
        arrset(pb.xnames,iv,"Q0305")
        ig = ig_["obj0305"]
        pbm.A[ig,iv] += Float64(1)
        ig = ig_["c3"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("Q0305",ix_)
        arrset(pb.xnames,iv,"Q0305")
        ig = ig_["c5"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("Q0405",ix_)
        arrset(pb.xnames,iv,"Q0405")
        ig = ig_["obj0405"]
        pbm.A[ig,iv] += Float64(1)
        ig = ig_["c4"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("Q0405",ix_)
        arrset(pb.xnames,iv,"Q0405")
        ig = ig_["c5"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("Q0406",ix_)
        arrset(pb.xnames,iv,"Q0406")
        ig = ig_["obj0406"]
        pbm.A[ig,iv] += Float64(1)
        ig = ig_["c4"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("Q0406",ix_)
        arrset(pb.xnames,iv,"Q0406")
        ig = ig_["c6"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("Q0607",ix_)
        arrset(pb.xnames,iv,"Q0607")
        ig = ig_["obj0607"]
        pbm.A[ig,iv] += Float64(1)
        ig = ig_["c6"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("Q0607",ix_)
        arrset(pb.xnames,iv,"Q0607")
        ig = ig_["c7"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("Q0705",ix_)
        arrset(pb.xnames,iv,"Q0705")
        ig = ig_["obj0705"]
        pbm.A[ig,iv] += Float64(1)
        ig = ig_["c5"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("Q0705",ix_)
        arrset(pb.xnames,iv,"Q0705")
        ig = ig_["c7"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("Q01u0",ix_)
        arrset(pb.xnames,iv,"Q01u0")
        ig = ig_["obj"]
        pbm.A[ig,iv] += Float64(210)
        ig = ig_["c1"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("Q01u0",ix_)
        arrset(pb.xnames,iv,"Q01u0")
        ig = ig_["c8"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("y02up",ix_)
        arrset(pb.xnames,iv,"y02up")
        ig = ig_["obj"]
        pbm.A[ig,iv] += Float64(210)
        ig = ig_["c2"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("y02up",ix_)
        arrset(pb.xnames,iv,"y02up")
        ig = ig_["c9"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("y03up",ix_)
        arrset(pb.xnames,iv,"y03up")
        ig = ig_["obj"]
        pbm.A[ig,iv] += Float64(210)
        ig = ig_["c3"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("y03up",ix_)
        arrset(pb.xnames,iv,"y03up")
        ig = ig_["c9"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("y04up",ix_)
        arrset(pb.xnames,iv,"y04up")
        ig = ig_["obj"]
        pbm.A[ig,iv] += Float64(210)
        ig = ig_["c4"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("y04up",ix_)
        arrset(pb.xnames,iv,"y04up")
        ig = ig_["c9"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("y05up",ix_)
        arrset(pb.xnames,iv,"y05up")
        ig = ig_["obj"]
        pbm.A[ig,iv] += Float64(210)
        ig = ig_["c5"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("y05up",ix_)
        arrset(pb.xnames,iv,"y05up")
        ig = ig_["c9"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("y06up",ix_)
        arrset(pb.xnames,iv,"y06up")
        ig = ig_["obj"]
        pbm.A[ig,iv] += Float64(210)
        ig = ig_["c6"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("y06up",ix_)
        arrset(pb.xnames,iv,"y06up")
        ig = ig_["c9"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("y07up",ix_)
        arrset(pb.xnames,iv,"y07up")
        ig = ig_["obj"]
        pbm.A[ig,iv] += Float64(210)
        ig = ig_["c7"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("y07up",ix_)
        arrset(pb.xnames,iv,"y07up")
        ig = ig_["c9"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("yqu02",ix_)
        arrset(pb.xnames,iv,"yqu02")
        ig = ig_["obj"]
        pbm.A[ig,iv] += Float64(-175)
        ig = ig_["c2"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("yqu02",ix_)
        arrset(pb.xnames,iv,"yqu02")
        ig = ig_["c10"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("yqu03",ix_)
        arrset(pb.xnames,iv,"yqu03")
        ig = ig_["obj"]
        pbm.A[ig,iv] += Float64(-190)
        ig = ig_["c3"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("yqu03",ix_)
        arrset(pb.xnames,iv,"yqu03")
        ig = ig_["c10"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("yqu04",ix_)
        arrset(pb.xnames,iv,"yqu04")
        ig = ig_["obj"]
        pbm.A[ig,iv] += Float64(-185)
        ig = ig_["c4"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("yqu04",ix_)
        arrset(pb.xnames,iv,"yqu04")
        ig = ig_["c10"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("yqu05",ix_)
        arrset(pb.xnames,iv,"yqu05")
        ig = ig_["obj"]
        pbm.A[ig,iv] += Float64(-180)
        ig = ig_["c5"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("yqu05",ix_)
        arrset(pb.xnames,iv,"yqu05")
        ig = ig_["c10"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("yqu06",ix_)
        arrset(pb.xnames,iv,"yqu06")
        ig = ig_["obj"]
        pbm.A[ig,iv] += Float64(-195)
        ig = ig_["c6"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("yqu06",ix_)
        arrset(pb.xnames,iv,"yqu06")
        ig = ig_["c10"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("yqu07",ix_)
        arrset(pb.xnames,iv,"yqu07")
        ig = ig_["obj"]
        pbm.A[ig,iv] += Float64(-190)
        ig = ig_["c7"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("yqu07",ix_)
        arrset(pb.xnames,iv,"yqu07")
        ig = ig_["c10"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("Q0201",ix_)
        arrset(pb.xnames,iv,"Q0201")
        ig = ig_["c1"]
        pbm.A[ig,iv] += Float64(-1)
        ig = ig_["c2"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("Q0302",ix_)
        arrset(pb.xnames,iv,"Q0302")
        ig = ig_["c2"]
        pbm.A[ig,iv] += Float64(-1)
        ig = ig_["c3"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("Q0402",ix_)
        arrset(pb.xnames,iv,"Q0402")
        ig = ig_["c2"]
        pbm.A[ig,iv] += Float64(-1)
        ig = ig_["c4"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("Q0503",ix_)
        arrset(pb.xnames,iv,"Q0503")
        ig = ig_["c3"]
        pbm.A[ig,iv] += Float64(-1)
        ig = ig_["c5"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("Q0504",ix_)
        arrset(pb.xnames,iv,"Q0504")
        ig = ig_["c4"]
        pbm.A[ig,iv] += Float64(-1)
        ig = ig_["c5"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("Q0604",ix_)
        arrset(pb.xnames,iv,"Q0604")
        ig = ig_["c4"]
        pbm.A[ig,iv] += Float64(-1)
        ig = ig_["c6"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("Q0507",ix_)
        arrset(pb.xnames,iv,"Q0507")
        ig = ig_["c5"]
        pbm.A[ig,iv] += Float64(1)
        ig = ig_["c7"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("Q0706",ix_)
        arrset(pb.xnames,iv,"Q0706")
        ig = ig_["c6"]
        pbm.A[ig,iv] += Float64(-1)
        ig = ig_["c7"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("yupu0",ix_)
        arrset(pb.xnames,iv,"yupu0")
        ig = ig_["c8"]
        pbm.A[ig,iv] += Float64(-1)
        ig = ig_["c9"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("yu0uq",ix_)
        arrset(pb.xnames,iv,"yu0uq")
        ig = ig_["c8"]
        pbm.A[ig,iv] += Float64(1)
        ig = ig_["c10"]
        pbm.A[ig,iv] += Float64(-1)
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
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
        pbm.gconst[ig_["c1"]] = Float64(1120)
        pbm.gconst[ig_["c2"]] = Float64(-100)
        pbm.gconst[ig_["c3"]] = Float64(-100)
        pbm.gconst[ig_["c4"]] = Float64(-120)
        pbm.gconst[ig_["c5"]] = Float64(-270)
        pbm.gconst[ig_["c6"]] = Float64(-330)
        pbm.gconst[ig_["c7"]] = Float64(-200)
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xupper[ix_["Q0102"]] = 1200
        pb.xupper[ix_["Q0203"]] = 1200
        pb.xupper[ix_["Q0204"]] = 1200
        pb.xupper[ix_["Q0305"]] = 1200
        pb.xupper[ix_["Q0405"]] = 1200
        pb.xupper[ix_["Q0406"]] = 1200
        pb.xupper[ix_["Q0607"]] = 1200
        pb.xupper[ix_["Q0705"]] = 1200
        pb.xupper[ix_["Q01u0"]] = 1200
        pb.xupper[ix_["y02up"]] = 1200
        pb.xupper[ix_["y03up"]] = 1200
        pb.xupper[ix_["y04up"]] = 1200
        pb.xupper[ix_["y05up"]] = 1200
        pb.xupper[ix_["y06up"]] = 1200
        pb.xupper[ix_["y07up"]] = 1200
        pb.xupper[ix_["yqu02"]] = 1200
        pb.xupper[ix_["yqu03"]] = 1200
        pb.xupper[ix_["yqu04"]] = 1200
        pb.xupper[ix_["yqu05"]] = 1200
        pb.xupper[ix_["yqu06"]] = 1200
        pb.xupper[ix_["yqu07"]] = 1200
        pb.xupper[ix_["Q0201"]] = 1200
        pb.xupper[ix_["Q0302"]] = 1200
        pb.xupper[ix_["Q0402"]] = 1200
        pb.xupper[ix_["Q0503"]] = 1200
        pb.xupper[ix_["Q0504"]] = 1200
        pb.xupper[ix_["Q0604"]] = 1200
        pb.xupper[ix_["Q0507"]] = 1200
        pb.xupper[ix_["Q0706"]] = 1200
        pb.xupper[ix_["yupu0"]] = 1200
        pb.xupper[ix_["yu0uq"]] = 1200
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gPOWER",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["obj0102"]
        arrset(pbm.grftype,ig,"gPOWER")
        ig = ig_["obj0203"]
        arrset(pbm.grftype,ig,"gPOWER")
        ig = ig_["obj0204"]
        arrset(pbm.grftype,ig,"gPOWER")
        ig = ig_["obj0305"]
        arrset(pbm.grftype,ig,"gPOWER")
        ig = ig_["obj0405"]
        arrset(pbm.grftype,ig,"gPOWER")
        ig = ig_["obj0406"]
        arrset(pbm.grftype,ig,"gPOWER")
        ig = ig_["obj0607"]
        arrset(pbm.grftype,ig,"gPOWER")
        ig = ig_["obj0705"]
        arrset(pbm.grftype,ig,"gPOWER")
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
# LO SOLUTION           1.054938D+04
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
        pb.lincons   = collect(1:length(pbm.congrps))
        pb.pbclass = "C-ONR2-MN-31-10"
        pb.x0          = zeros(Float64,pb.n)
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2]
        pb.conderlvl  = pbm.conderlvl;
        return pb, pbm

# ********************
#  SET UP THE GROUPS *
#  ROUTINE           *
# ********************

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    elseif action == "gPOWER"

        GVAR_   = args[1]
        igr_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_= GVAR_^2.852
        if nargout>1
            g_ = 2.852*GVAR_^1.852
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_ = 5.282*GVAR_^.852
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

