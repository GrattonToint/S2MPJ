function NASH(action::String,args::Union{Any}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : NASH
#    *********
# 
#    A quadratic programming reformulation of a linear
#    complementarity problem arising from Nash equilibrium
#    provided by Michael Ferris
# 
#    classification = "C-CQLR2-AN-72-24"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 21 VI 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "NASH"
    if ( !isdefined(@__MODULE__, :s2mpj_ii) )
        error( "Please include(\"s2mpjlib.jl\") using \"s2mpjlib.jl\" from the S2MPJ distribution before calling NASH.")
    end

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
        v_["N"] = 72
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        irA   = Int64[]
        icA   = Int64[]
        valA  = Float64[]
        for J = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("X"*string(J),ix_)
            arrset(pb.xnames,iv,"X"*string(J))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X56"])
        push!(valA,Float64(1000.0))
        push!(irA,ig)
        push!(icA,ix_["X57"])
        push!(valA,Float64(500.0))
        push!(irA,ig)
        push!(icA,ix_["X58"])
        push!(valA,Float64(1000.0))
        ig,ig_,_ = s2mpj_ii("C1",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C1")
        push!(irA,ig)
        push!(icA,ix_["X25"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X49"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X1"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X2"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X3"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X4"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X5"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X6"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X7"])
        push!(valA,Float64(-1.0))
        ig,ig_,_ = s2mpj_ii("C2",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C2")
        push!(irA,ig)
        push!(icA,ix_["X26"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X50"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X1"])
        push!(valA,Float64(0.02309))
        push!(irA,ig)
        push!(icA,ix_["X2"])
        push!(valA,Float64(0.02309))
        push!(irA,ig)
        push!(icA,ix_["X4"])
        push!(valA,Float64(0.02309))
        push!(irA,ig)
        push!(icA,ix_["X6"])
        push!(valA,Float64(0.02309))
        push!(irA,ig)
        push!(icA,ix_["X11"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X22"])
        push!(valA,Float64(0.288626))
        push!(irA,ig)
        push!(icA,ix_["X23"])
        push!(valA,Float64(0.263887))
        push!(irA,ig)
        push!(icA,ix_["X24"])
        push!(valA,Float64(0.447486))
        ig,ig_,_ = s2mpj_ii("C3",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C3")
        push!(irA,ig)
        push!(icA,ix_["X27"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X51"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X1"])
        push!(valA,Float64(0.02309))
        push!(irA,ig)
        push!(icA,ix_["X3"])
        push!(valA,Float64(0.02309))
        push!(irA,ig)
        push!(icA,ix_["X5"])
        push!(valA,Float64(0.02309))
        push!(irA,ig)
        push!(icA,ix_["X7"])
        push!(valA,Float64(0.02309))
        push!(irA,ig)
        push!(icA,ix_["X12"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X22"])
        push!(valA,Float64(0.288626))
        push!(irA,ig)
        push!(icA,ix_["X23"])
        push!(valA,Float64(0.263887))
        push!(irA,ig)
        push!(icA,ix_["X24"])
        push!(valA,Float64(0.447486))
        ig,ig_,_ = s2mpj_ii("C4",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C4")
        push!(irA,ig)
        push!(icA,ix_["X28"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X52"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X1"])
        push!(valA,Float64(0.02309))
        push!(irA,ig)
        push!(icA,ix_["X2"])
        push!(valA,Float64(0.02309))
        push!(irA,ig)
        push!(icA,ix_["X4"])
        push!(valA,Float64(0.02309))
        push!(irA,ig)
        push!(icA,ix_["X6"])
        push!(valA,Float64(0.02309))
        push!(irA,ig)
        push!(icA,ix_["X11"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X22"])
        push!(valA,Float64(0.288626))
        push!(irA,ig)
        push!(icA,ix_["X23"])
        push!(valA,Float64(0.263887))
        push!(irA,ig)
        push!(icA,ix_["X24"])
        push!(valA,Float64(0.447486))
        ig,ig_,_ = s2mpj_ii("C5",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C5")
        push!(irA,ig)
        push!(icA,ix_["X29"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X53"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X1"])
        push!(valA,Float64(0.02309))
        push!(irA,ig)
        push!(icA,ix_["X3"])
        push!(valA,Float64(0.02309))
        push!(irA,ig)
        push!(icA,ix_["X5"])
        push!(valA,Float64(0.02309))
        push!(irA,ig)
        push!(icA,ix_["X7"])
        push!(valA,Float64(0.02309))
        push!(irA,ig)
        push!(icA,ix_["X12"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X22"])
        push!(valA,Float64(0.288626))
        push!(irA,ig)
        push!(icA,ix_["X23"])
        push!(valA,Float64(0.263887))
        push!(irA,ig)
        push!(icA,ix_["X24"])
        push!(valA,Float64(0.447486))
        ig,ig_,_ = s2mpj_ii("C6",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C6")
        push!(irA,ig)
        push!(icA,ix_["X30"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X54"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X1"])
        push!(valA,Float64(0.02309))
        push!(irA,ig)
        push!(icA,ix_["X2"])
        push!(valA,Float64(0.02309))
        push!(irA,ig)
        push!(icA,ix_["X4"])
        push!(valA,Float64(0.02309))
        push!(irA,ig)
        push!(icA,ix_["X6"])
        push!(valA,Float64(0.02309))
        push!(irA,ig)
        push!(icA,ix_["X11"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X22"])
        push!(valA,Float64(0.288626))
        push!(irA,ig)
        push!(icA,ix_["X23"])
        push!(valA,Float64(0.263887))
        push!(irA,ig)
        push!(icA,ix_["X24"])
        push!(valA,Float64(0.447486))
        ig,ig_,_ = s2mpj_ii("C7",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C7")
        push!(irA,ig)
        push!(icA,ix_["X31"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X55"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X1"])
        push!(valA,Float64(0.02309))
        push!(irA,ig)
        push!(icA,ix_["X3"])
        push!(valA,Float64(0.02309))
        push!(irA,ig)
        push!(icA,ix_["X5"])
        push!(valA,Float64(0.02309))
        push!(irA,ig)
        push!(icA,ix_["X7"])
        push!(valA,Float64(0.02309))
        push!(irA,ig)
        push!(icA,ix_["X12"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X22"])
        push!(valA,Float64(0.288626))
        push!(irA,ig)
        push!(icA,ix_["X23"])
        push!(valA,Float64(0.263887))
        push!(irA,ig)
        push!(icA,ix_["X24"])
        push!(valA,Float64(0.447486))
        ig,ig_,_ = s2mpj_ii("C8",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C8")
        push!(irA,ig)
        push!(icA,ix_["X32"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X56"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X11"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X22"])
        push!(valA,Float64(-1.0))
        ig,ig_,_ = s2mpj_ii("C9",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C9")
        push!(irA,ig)
        push!(icA,ix_["X33"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X57"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X11"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X23"])
        push!(valA,Float64(-1.0))
        ig,ig_,_ = s2mpj_ii("C10",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C10")
        push!(irA,ig)
        push!(icA,ix_["X34"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X58"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X12"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X23"])
        push!(valA,Float64(-1.0))
        ig,ig_,_ = s2mpj_ii("C11",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C11")
        push!(irA,ig)
        push!(icA,ix_["X35"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X59"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X2"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X4"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X6"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X8"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X9"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("C12",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C12")
        push!(irA,ig)
        push!(icA,ix_["X36"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X60"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X3"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X5"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X7"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X10"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("C13",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C13")
        push!(irA,ig)
        push!(icA,ix_["X37"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X61"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X16"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X19"])
        push!(valA,Float64(-.33))
        push!(irA,ig)
        push!(icA,ix_["X20"])
        push!(valA,Float64(0.33))
        ig,ig_,_ = s2mpj_ii("C14",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C14")
        push!(irA,ig)
        push!(icA,ix_["X38"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X62"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X17"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X19"])
        push!(valA,Float64(-.67))
        push!(irA,ig)
        push!(icA,ix_["X20"])
        push!(valA,Float64(-.33))
        ig,ig_,_ = s2mpj_ii("C15",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C15")
        push!(irA,ig)
        push!(icA,ix_["X39"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X63"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X18"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X19"])
        push!(valA,Float64(-.33))
        push!(irA,ig)
        push!(icA,ix_["X20"])
        push!(valA,Float64(-.67))
        ig,ig_,_ = s2mpj_ii("C16",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C16")
        push!(irA,ig)
        push!(icA,ix_["X40"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X64"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X13"])
        push!(valA,Float64(-1.0))
        ig,ig_,_ = s2mpj_ii("C17",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C17")
        push!(irA,ig)
        push!(icA,ix_["X41"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X65"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X14"])
        push!(valA,Float64(-1.0))
        ig,ig_,_ = s2mpj_ii("C18",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C18")
        push!(irA,ig)
        push!(icA,ix_["X42"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X66"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X15"])
        push!(valA,Float64(-1.0))
        ig,ig_,_ = s2mpj_ii("C19",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C19")
        push!(irA,ig)
        push!(icA,ix_["X43"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X67"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X13"])
        push!(valA,Float64(0.33))
        push!(irA,ig)
        push!(icA,ix_["X14"])
        push!(valA,Float64(0.67))
        push!(irA,ig)
        push!(icA,ix_["X15"])
        push!(valA,Float64(0.33))
        push!(irA,ig)
        push!(icA,ix_["X22"])
        push!(valA,Float64(-1.0))
        ig,ig_,_ = s2mpj_ii("C20",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C20")
        push!(irA,ig)
        push!(icA,ix_["X44"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X68"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X13"])
        push!(valA,Float64(-.33))
        push!(irA,ig)
        push!(icA,ix_["X14"])
        push!(valA,Float64(0.33))
        push!(irA,ig)
        push!(icA,ix_["X15"])
        push!(valA,Float64(0.67))
        push!(irA,ig)
        push!(icA,ix_["X23"])
        push!(valA,Float64(-1.0))
        ig,ig_,_ = s2mpj_ii("C21",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C21")
        push!(irA,ig)
        push!(icA,ix_["X45"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X69"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X24"])
        push!(valA,Float64(-1.0))
        ig,ig_,_ = s2mpj_ii("C22",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C22")
        push!(irA,ig)
        push!(icA,ix_["X46"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X70"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X1"])
        push!(valA,Float64(-.288626))
        push!(irA,ig)
        push!(icA,ix_["X8"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X19"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X22"])
        push!(valA,Float64(8.892169))
        push!(irA,ig)
        push!(icA,ix_["X23"])
        push!(valA,Float64(-3.298588))
        push!(irA,ig)
        push!(icA,ix_["X24"])
        push!(valA,Float64(-5.593581))
        ig,ig_,_ = s2mpj_ii("C23",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C23")
        push!(irA,ig)
        push!(icA,ix_["X47"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X71"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X1"])
        push!(valA,Float64(-.263887))
        push!(irA,ig)
        push!(icA,ix_["X9"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X10"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X20"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X22"])
        push!(valA,Float64(-3.298588))
        push!(irA,ig)
        push!(icA,ix_["X23"])
        push!(valA,Float64(8.412719))
        push!(irA,ig)
        push!(icA,ix_["X24"])
        push!(valA,Float64(-5.114131))
        ig,ig_,_ = s2mpj_ii("C24",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C24")
        push!(irA,ig)
        push!(icA,ix_["X48"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X72"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X1"])
        push!(valA,Float64(-.447486))
        push!(irA,ig)
        push!(icA,ix_["X21"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X22"])
        push!(valA,Float64(-5.593581))
        push!(irA,ig)
        push!(icA,ix_["X23"])
        push!(valA,Float64(-5.114131))
        push!(irA,ig)
        push!(icA,ix_["X24"])
        push!(valA,Float64(10.707712))
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
        pbm.gconst[ig_["C2"]] = Float64(35.100673)
        pbm.gconst[ig_["C3"]] = Float64(35.100673)
        pbm.gconst[ig_["C4"]] = Float64(35.100673)
        pbm.gconst[ig_["C5"]] = Float64(35.100673)
        pbm.gconst[ig_["C6"]] = Float64(35.100673)
        pbm.gconst[ig_["C7"]] = Float64(35.100673)
        pbm.gconst[ig_["C8"]] = Float64(-15.0)
        pbm.gconst[ig_["C9"]] = Float64(-15.0)
        pbm.gconst[ig_["C10"]] = Float64(-20.0)
        pbm.gconst[ig_["C22"]] = Float64(61.241589)
        pbm.gconst[ig_["C23"]] = Float64(-1.150548)
        pbm.gconst[ig_["C24"]] = Float64(-60.091041)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = fill(0.E+00,pb.n)
        pb.xupper = fill(0.E+00,pb.n)
        pb.xlower[ix_["X1"]] = -Inf
        pb.xupper[ix_["X1"]] = +Inf
        pb.xupper[ix_["X8"]] = 1000.0
        pb.xupper[ix_["X9"]] = 500.0
        pb.xupper[ix_["X10"]] = 1000.0
        pb.xlower[ix_["X11"]] = -Inf
        pb.xupper[ix_["X11"]] = +Inf
        pb.xlower[ix_["X12"]] = -Inf
        pb.xupper[ix_["X12"]] = +Inf
        pb.xlower[ix_["X13"]] = -Inf
        pb.xupper[ix_["X13"]] = +Inf
        pb.xlower[ix_["X14"]] = -Inf
        pb.xupper[ix_["X14"]] = +Inf
        pb.xlower[ix_["X15"]] = -Inf
        pb.xupper[ix_["X15"]] = +Inf
        pb.xlower[ix_["X16"]] = -Inf
        pb.xupper[ix_["X16"]] = +Inf
        pb.xlower[ix_["X17"]] = -Inf
        pb.xupper[ix_["X17"]] = +Inf
        pb.xlower[ix_["X18"]] = -Inf
        pb.xupper[ix_["X18"]] = +Inf
        pb.xlower[ix_["X19"]] = -Inf
        pb.xupper[ix_["X19"]] = +Inf
        pb.xlower[ix_["X20"]] = -Inf
        pb.xupper[ix_["X20"]] = +Inf
        pb.xlower[ix_["X21"]] = -Inf
        pb.xupper[ix_["X21"]] = +Inf
        pb.xlower[ix_["X22"]] = -Inf
        pb.xupper[ix_["X22"]] = +Inf
        pb.xlower[ix_["X23"]] = -Inf
        pb.xupper[ix_["X23"]] = +Inf
        pb.xlower[ix_["X24"]] = -Inf
        pb.xupper[ix_["X24"]] = +Inf
        #%%%%%%%%%%%%%%%%%%%% QUADRATIC %%%%%%%%%%%%%%%%%%%
        irH  = Int64[]
        icH  = Int64[]
        valH = Float64[]
        push!(irH,ix_["X25"])
        push!(icH,ix_["X1"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X1"])
        push!(icH,ix_["X25"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X26"])
        push!(icH,ix_["X2"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X2"])
        push!(icH,ix_["X26"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X27"])
        push!(icH,ix_["X3"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X3"])
        push!(icH,ix_["X27"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X28"])
        push!(icH,ix_["X4"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X4"])
        push!(icH,ix_["X28"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X29"])
        push!(icH,ix_["X5"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X5"])
        push!(icH,ix_["X29"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X30"])
        push!(icH,ix_["X6"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X6"])
        push!(icH,ix_["X30"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X31"])
        push!(icH,ix_["X7"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X7"])
        push!(icH,ix_["X31"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X32"])
        push!(icH,ix_["X8"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X8"])
        push!(icH,ix_["X32"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X33"])
        push!(icH,ix_["X9"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X9"])
        push!(icH,ix_["X33"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X34"])
        push!(icH,ix_["X10"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X10"])
        push!(icH,ix_["X34"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X35"])
        push!(icH,ix_["X11"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X11"])
        push!(icH,ix_["X35"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X36"])
        push!(icH,ix_["X12"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X12"])
        push!(icH,ix_["X36"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X37"])
        push!(icH,ix_["X13"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X13"])
        push!(icH,ix_["X37"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X38"])
        push!(icH,ix_["X14"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X14"])
        push!(icH,ix_["X38"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X39"])
        push!(icH,ix_["X15"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X15"])
        push!(icH,ix_["X39"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X40"])
        push!(icH,ix_["X16"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X16"])
        push!(icH,ix_["X40"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X41"])
        push!(icH,ix_["X17"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X17"])
        push!(icH,ix_["X41"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X42"])
        push!(icH,ix_["X18"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X18"])
        push!(icH,ix_["X42"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X43"])
        push!(icH,ix_["X19"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X19"])
        push!(icH,ix_["X43"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X44"])
        push!(icH,ix_["X20"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X20"])
        push!(icH,ix_["X44"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X45"])
        push!(icH,ix_["X21"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X21"])
        push!(icH,ix_["X45"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X46"])
        push!(icH,ix_["X22"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X22"])
        push!(icH,ix_["X46"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X47"])
        push!(icH,ix_["X23"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X23"])
        push!(icH,ix_["X47"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X48"])
        push!(icH,ix_["X24"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X24"])
        push!(icH,ix_["X48"])
        push!(valH,Float64(1.0))
        push!(irH,ix_["X49"])
        push!(icH,ix_["X1"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X1"])
        push!(icH,ix_["X49"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X50"])
        push!(icH,ix_["X2"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X2"])
        push!(icH,ix_["X50"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X51"])
        push!(icH,ix_["X3"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X3"])
        push!(icH,ix_["X51"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X52"])
        push!(icH,ix_["X4"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X4"])
        push!(icH,ix_["X52"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X53"])
        push!(icH,ix_["X5"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X5"])
        push!(icH,ix_["X53"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X54"])
        push!(icH,ix_["X6"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X6"])
        push!(icH,ix_["X54"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X55"])
        push!(icH,ix_["X7"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X7"])
        push!(icH,ix_["X55"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X56"])
        push!(icH,ix_["X8"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X8"])
        push!(icH,ix_["X56"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X57"])
        push!(icH,ix_["X9"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X9"])
        push!(icH,ix_["X57"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X58"])
        push!(icH,ix_["X10"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X10"])
        push!(icH,ix_["X58"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X59"])
        push!(icH,ix_["X11"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X11"])
        push!(icH,ix_["X59"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X60"])
        push!(icH,ix_["X12"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X12"])
        push!(icH,ix_["X60"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X61"])
        push!(icH,ix_["X13"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X13"])
        push!(icH,ix_["X61"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X62"])
        push!(icH,ix_["X14"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X14"])
        push!(icH,ix_["X62"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X63"])
        push!(icH,ix_["X15"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X15"])
        push!(icH,ix_["X63"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X64"])
        push!(icH,ix_["X16"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X16"])
        push!(icH,ix_["X64"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X65"])
        push!(icH,ix_["X17"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X17"])
        push!(icH,ix_["X65"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X66"])
        push!(icH,ix_["X18"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X18"])
        push!(icH,ix_["X66"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X67"])
        push!(icH,ix_["X19"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X19"])
        push!(icH,ix_["X67"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X68"])
        push!(icH,ix_["X20"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X20"])
        push!(icH,ix_["X68"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X69"])
        push!(icH,ix_["X21"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X21"])
        push!(icH,ix_["X69"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X70"])
        push!(icH,ix_["X22"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X22"])
        push!(icH,ix_["X70"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X71"])
        push!(icH,ix_["X23"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X23"])
        push!(icH,ix_["X71"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X72"])
        push!(icH,ix_["X24"])
        push!(valH,Float64(-1.0))
        push!(irH,ix_["X24"])
        push!(icH,ix_["X72"])
        push!(valH,Float64(-1.0))
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n)
        pbm.H = sparse(irH,icH,valH,pb.n,pb.n)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons   = collect(1:length(pbm.congrps))
        pb.pbclass = "C-CQLR2-AN-72-24"
        pb.x0          = zeros(Float64,pb.n)
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2]
        pb.conderlvl  = pbm.conderlvl;
        return pb, pbm


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

