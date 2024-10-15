function FCCU(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : FCCU
#    *********
# 
#    A simple data reconciliation for a fluid catalytic cracker.
# 
#                      +--------------------------+
#                      | FCCU data reconciliation |
#     +-----------+    +--------------------------+
#  1) | Flowsheet |                      Off_gas                      |------->
#     +-----------+                     |----------->                 | Propane
#                              MF_ohd   |                  |-------->F7
#                              |------>F4<-------|         |          |
#                              |        |        |         |DC3_feed  | Butane
#  Feed      Effluent          ^        |        |         |          |------->
#  ------->F1-------->F2......>|        |---------------->F5                  
#           ^                  v       DC4_feed  |         |DC4_btms  |------->
#           |                  |                 |         |          |  LCN
#           |                  |                 |<-------F6-------->F8
#           |                  | HCN              Lean_oil   C8spl_fd |  MCN
#           |                  |-------->                             |------->
#           |                  | LCO          
#           |                  |-------->              
#           |                  | HCO          
#           |                  |-------->              
#           |                  | MF_btms     
#           |                  v 
#           |<----------------F3-------->
#              Dec_recy         Decant
#     +--------------------+
#  2) | Objective function |
#     +--------------------+
#    Obj = sum 1->i [W_i*(C_flow_i - M_flow_i)**2]
#                           
#     Where: W_i       Weight on term i of objective function
#            C_flow_i  Computed flow i (a variable for this problem)
#            M_flow_i  Measrued flow i (a constant for this problem)
#     +-------------+
#  3) | Constraints |
#     +-------------+
#     These represent the linear mass balances around each
#     node, where a node (Fx) represents a single unit operation
#     in a fluid catalytics cracker.
#     +---------------+
#  4) | Initial point |
#     +---------------+
#     Feed       1.0
#     Effluent   1.0
#     MF_ohd     1.0
#     HCN        1.0
#     LCO        1.0
#     HCO        1.0
#     MF_btms    1.0
#     Decant     1.0
#     Dec_recy   1.0
#     Off_gas    1.0
#     DC4_feed   1.0
#     DC3_feed   1.0
#     DC4_btms   1.0
#     Lean_oil   1.0
#     Propane    1.0
#     Butane     1.0
#     C8spl_fd   1.0
#     LCN        1.0
#     MCN        1.0
#     Obj        7.36259000271320D+03
#     +------------------+
#  5) | Optimal solution |
#     +------------------+
#     Feed       3.11639D+01
#     Effluent   3.53528D+01
#     MF_ohd     1.94669D+01
#     HCN        2.94255D+00
#     LCO        4.94255D+00
#     HCO        3.44255D+00
#     MF_btms    4.55828D+00
#     Decant     3.69371D-01
#     Dec_recy   4.18891D+00
#     Off_gas    2.56075D+00
#     DC4_feed   2.41207D+01
#     DC3_feed   5.15601D+00
#     DC4_btms   1.89647D+01
#     Lean_oil   7.21458D+00
#     Propane    2.42801D+00
#     Butane     2.72801D+00
#     C8spl_fd   1.17501D+01
#     LCN        5.87506D+00
#     MCN        5.87506D+00
#     Obj        1.11491D+01
#     +-----------------------------------------------+
#  6) | SPEC.SPC (remove 1st * of every line to use). |
#     +-----------------------------------------------+
# BEGIN
# * maximizer-sought
# *  check-derivatives
#   ignore-derivative-bugs
# * use-scalings
# * print-scalings
#   finite-difference-gradients
# *  exact-second-derivatives-used
# * bfgs-approximate-second-derivatives-used
# * sr1-approximate-second-derivatives-used
#   bandsolver-preconditioned-cg-solver-used   5
# * diagonal-preconditioned-cg-solver-used
# * gill-murray-ponceleon-saunders-preconditioned-cg-solver-used
# * schnabel-eskow-preconditioned-cg-solver-used
# * munksgaards-preconditioned-cg-solver-used
#   exact-cauchy-point-required
# * inexact-cauchy-point-required
# * solve-bqp-accurately
# * two-norm-trust-region
# * gradient-tolerance    1.0D-5
# * constraint-tolerance  1.0D-5
#   trust-region-radius   1.0D+0
#   maximum-number-of-iterations   1000
#   print-level                    1
#   start-printing-at-iteration    0
#   stop-printing-at-iteration     1000
# END
# 
#    Source:
#    W. J. Korchinski, Profimatics, Inc,
#    325 Rolling Oaks Drive, Thousand Oaks, California, USA 91361-1200
#    Telephone: 1-805 496 6661, Fax: 1-805 373 5108
# 
#    SIF input: W. Korchinski, Spring 1993.
# 
#    classification = "C-SLR2-MN-19-8"
# 
# ***************************************************************
#  PROBLEM SPECIFICATION BEGINS HERE.
#  **********************************
#  **********************************           
# *************************************
#  Define objective function weights. *
# *************************************
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "FCCU"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["W1"] = 0.2
        v_["W2"] = 1.0
        v_["W3"] = 1.0
        v_["W4"] = 0.33333333
        v_["W5"] = 0.33333333
        v_["W6"] = 0.33333333
        v_["W7"] = 1.0
        v_["W8"] = 1.0
        v_["W9"] = 1.0
        v_["W10"] = 1.0
        v_["W11"] = 1.0
        v_["W12"] = 1.0
        v_["W13"] = 1.0
        v_["W14"] = 1.0
        v_["W15"] = 0.33333333
        v_["W16"] = 0.33333333
        v_["W17"] = 1.0
        v_["W18"] = 0.33333333
        v_["W19"] = 0.33333333
        v_["M1"] = 31.0
        v_["M2"] = 36.0
        v_["M3"] = 20.0
        v_["M4"] = 3.0
        v_["M5"] = 5.0
        v_["M6"] = 3.5
        v_["M7"] = 4.2
        v_["M8"] = 0.9
        v_["M9"] = 3.9
        v_["M10"] = 2.2
        v_["M11"] = 22.8
        v_["M12"] = 6.8
        v_["M13"] = 19.0
        v_["M14"] = 8.5
        v_["M15"] = 2.2
        v_["M16"] = 2.5
        v_["M17"] = 10.8
        v_["M18"] = 6.5
        v_["M19"] = 6.5
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        iv,ix_,_ = s2mpj_ii("Feed",ix_)
        arrset(pb.xnames,iv,"Feed")
        iv,ix_,_ = s2mpj_ii("Effluent",ix_)
        arrset(pb.xnames,iv,"Effluent")
        iv,ix_,_ = s2mpj_ii("MFuohd",ix_)
        arrset(pb.xnames,iv,"MFuohd")
        iv,ix_,_ = s2mpj_ii("HCN",ix_)
        arrset(pb.xnames,iv,"HCN")
        iv,ix_,_ = s2mpj_ii("LCO",ix_)
        arrset(pb.xnames,iv,"LCO")
        iv,ix_,_ = s2mpj_ii("HCO",ix_)
        arrset(pb.xnames,iv,"HCO")
        iv,ix_,_ = s2mpj_ii("MFubtms",ix_)
        arrset(pb.xnames,iv,"MFubtms")
        iv,ix_,_ = s2mpj_ii("Decant",ix_)
        arrset(pb.xnames,iv,"Decant")
        iv,ix_,_ = s2mpj_ii("Decurecy",ix_)
        arrset(pb.xnames,iv,"Decurecy")
        iv,ix_,_ = s2mpj_ii("Offugas",ix_)
        arrset(pb.xnames,iv,"Offugas")
        iv,ix_,_ = s2mpj_ii("DC4ufeed",ix_)
        arrset(pb.xnames,iv,"DC4ufeed")
        iv,ix_,_ = s2mpj_ii("DC3ufeed",ix_)
        arrset(pb.xnames,iv,"DC3ufeed")
        iv,ix_,_ = s2mpj_ii("DC4ubtms",ix_)
        arrset(pb.xnames,iv,"DC4ubtms")
        iv,ix_,_ = s2mpj_ii("Leanuoil",ix_)
        arrset(pb.xnames,iv,"Leanuoil")
        iv,ix_,_ = s2mpj_ii("Propane",ix_)
        arrset(pb.xnames,iv,"Propane")
        iv,ix_,_ = s2mpj_ii("Butane",ix_)
        arrset(pb.xnames,iv,"Butane")
        iv,ix_,_ = s2mpj_ii("C8splufd",ix_)
        arrset(pb.xnames,iv,"C8splufd")
        iv,ix_,_ = s2mpj_ii("LCN",ix_)
        arrset(pb.xnames,iv,"LCN")
        iv,ix_,_ = s2mpj_ii("MCN",ix_)
        arrset(pb.xnames,iv,"MCN")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("F1",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"F1")
        iv = ix_["Feed"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["Decurecy"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["Effluent"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("F2",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"F2")
        iv = ix_["Effluent"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["MFuohd"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["HCN"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["LCO"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["HCO"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["MFubtms"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("F3",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"F3")
        iv = ix_["MFubtms"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["Decant"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["Decurecy"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("F4",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"F4")
        iv = ix_["MFuohd"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["Leanuoil"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["Offugas"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["DC4ufeed"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("F5",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"F5")
        iv = ix_["DC4ufeed"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["DC3ufeed"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["DC4ubtms"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("F6",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"F6")
        iv = ix_["DC4ubtms"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["Leanuoil"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["C8splufd"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("F7",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"F7")
        iv = ix_["DC3ufeed"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["Propane"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["Butane"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("F8",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"F8")
        iv = ix_["C8splufd"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["LCN"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["MCN"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("Obj1",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["Feed"]
        pbm.A[ig,iv] += Float64(1.0)
        arrset(pbm.gscale,ig,Float64(v_["W1"]))
        ig,ig_,_ = s2mpj_ii("Obj2",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["Effluent"]
        pbm.A[ig,iv] += Float64(1.0)
        arrset(pbm.gscale,ig,Float64(v_["W2"]))
        ig,ig_,_ = s2mpj_ii("Obj3",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["MFuohd"]
        pbm.A[ig,iv] += Float64(1.0)
        arrset(pbm.gscale,ig,Float64(v_["W3"]))
        ig,ig_,_ = s2mpj_ii("Obj4",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["HCN"]
        pbm.A[ig,iv] += Float64(1.0)
        arrset(pbm.gscale,ig,Float64(v_["W4"]))
        ig,ig_,_ = s2mpj_ii("Obj5",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["LCO"]
        pbm.A[ig,iv] += Float64(1.0)
        arrset(pbm.gscale,ig,Float64(v_["W5"]))
        ig,ig_,_ = s2mpj_ii("Obj6",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["HCO"]
        pbm.A[ig,iv] += Float64(1.0)
        arrset(pbm.gscale,ig,Float64(v_["W6"]))
        ig,ig_,_ = s2mpj_ii("Obj7",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["MFubtms"]
        pbm.A[ig,iv] += Float64(1.0)
        arrset(pbm.gscale,ig,Float64(v_["W7"]))
        ig,ig_,_ = s2mpj_ii("Obj8",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["Decant"]
        pbm.A[ig,iv] += Float64(1.0)
        arrset(pbm.gscale,ig,Float64(v_["W8"]))
        ig,ig_,_ = s2mpj_ii("Obj9",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["Decurecy"]
        pbm.A[ig,iv] += Float64(1.0)
        arrset(pbm.gscale,ig,Float64(v_["W9"]))
        ig,ig_,_ = s2mpj_ii("Obj10",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["Offugas"]
        pbm.A[ig,iv] += Float64(1.0)
        arrset(pbm.gscale,ig,Float64(v_["W10"]))
        ig,ig_,_ = s2mpj_ii("Obj11",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["DC4ufeed"]
        pbm.A[ig,iv] += Float64(1.0)
        arrset(pbm.gscale,ig,Float64(v_["W11"]))
        ig,ig_,_ = s2mpj_ii("Obj12",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["DC3ufeed"]
        pbm.A[ig,iv] += Float64(1.0)
        arrset(pbm.gscale,ig,Float64(v_["W12"]))
        ig,ig_,_ = s2mpj_ii("Obj13",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["DC4ubtms"]
        pbm.A[ig,iv] += Float64(1.0)
        arrset(pbm.gscale,ig,Float64(v_["W13"]))
        ig,ig_,_ = s2mpj_ii("Obj14",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["Leanuoil"]
        pbm.A[ig,iv] += Float64(1.0)
        arrset(pbm.gscale,ig,Float64(v_["W14"]))
        ig,ig_,_ = s2mpj_ii("Obj15",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["Propane"]
        pbm.A[ig,iv] += Float64(1.0)
        arrset(pbm.gscale,ig,Float64(v_["W15"]))
        ig,ig_,_ = s2mpj_ii("Obj16",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["Butane"]
        pbm.A[ig,iv] += Float64(1.0)
        arrset(pbm.gscale,ig,Float64(v_["W16"]))
        ig,ig_,_ = s2mpj_ii("Obj17",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["C8splufd"]
        pbm.A[ig,iv] += Float64(1.0)
        arrset(pbm.gscale,ig,Float64(v_["W17"]))
        ig,ig_,_ = s2mpj_ii("Obj18",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["LCN"]
        pbm.A[ig,iv] += Float64(1.0)
        arrset(pbm.gscale,ig,Float64(v_["W18"]))
        ig,ig_,_ = s2mpj_ii("Obj19",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["MCN"]
        pbm.A[ig,iv] += Float64(1.0)
        arrset(pbm.gscale,ig,Float64(v_["W19"]))
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
        pbm.gconst[ig_["Obj1"]] = Float64(v_["M1"])
        pbm.gconst[ig_["Obj2"]] = Float64(v_["M2"])
        pbm.gconst[ig_["Obj3"]] = Float64(v_["M3"])
        pbm.gconst[ig_["Obj4"]] = Float64(v_["M4"])
        pbm.gconst[ig_["Obj5"]] = Float64(v_["M5"])
        pbm.gconst[ig_["Obj6"]] = Float64(v_["M6"])
        pbm.gconst[ig_["Obj7"]] = Float64(v_["M7"])
        pbm.gconst[ig_["Obj8"]] = Float64(v_["M8"])
        pbm.gconst[ig_["Obj9"]] = Float64(v_["M9"])
        pbm.gconst[ig_["Obj10"]] = Float64(v_["M10"])
        pbm.gconst[ig_["Obj11"]] = Float64(v_["M11"])
        pbm.gconst[ig_["Obj12"]] = Float64(v_["M12"])
        pbm.gconst[ig_["Obj13"]] = Float64(v_["M13"])
        pbm.gconst[ig_["Obj14"]] = Float64(v_["M14"])
        pbm.gconst[ig_["Obj15"]] = Float64(v_["M15"])
        pbm.gconst[ig_["Obj16"]] = Float64(v_["M16"])
        pbm.gconst[ig_["Obj17"]] = Float64(v_["M17"])
        pbm.gconst[ig_["Obj18"]] = Float64(v_["M18"])
        pbm.gconst[ig_["Obj19"]] = Float64(v_["M19"])
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        pb.x0[ix_["Feed"]] = Float64(1.0)
        pb.x0[ix_["Effluent"]] = Float64(1.0)
        pb.x0[ix_["MFuohd"]] = Float64(1.0)
        pb.x0[ix_["HCN"]] = Float64(1.0)
        pb.x0[ix_["LCO"]] = Float64(1.0)
        pb.x0[ix_["HCO"]] = Float64(1.0)
        pb.x0[ix_["MFubtms"]] = Float64(1.0)
        pb.x0[ix_["Decant"]] = Float64(1.0)
        pb.x0[ix_["Decurecy"]] = Float64(1.0)
        pb.x0[ix_["Offugas"]] = Float64(1.0)
        pb.x0[ix_["DC4ufeed"]] = Float64(1.0)
        pb.x0[ix_["DC3ufeed"]] = Float64(1.0)
        pb.x0[ix_["DC4ubtms"]] = Float64(1.0)
        pb.x0[ix_["Leanuoil"]] = Float64(1.0)
        pb.x0[ix_["Propane"]] = Float64(1.0)
        pb.x0[ix_["Butane"]] = Float64(1.0)
        pb.x0[ix_["C8splufd"]] = Float64(1.0)
        pb.x0[ix_["LCN"]] = Float64(1.0)
        pb.x0[ix_["MCN"]] = Float64(1.0)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gSQUARE",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["Obj1"]
        arrset(pbm.grftype,ig,"gSQUARE")
        ig = ig_["Obj2"]
        arrset(pbm.grftype,ig,"gSQUARE")
        ig = ig_["Obj3"]
        arrset(pbm.grftype,ig,"gSQUARE")
        ig = ig_["Obj4"]
        arrset(pbm.grftype,ig,"gSQUARE")
        ig = ig_["Obj5"]
        arrset(pbm.grftype,ig,"gSQUARE")
        ig = ig_["Obj6"]
        arrset(pbm.grftype,ig,"gSQUARE")
        ig = ig_["Obj7"]
        arrset(pbm.grftype,ig,"gSQUARE")
        ig = ig_["Obj8"]
        arrset(pbm.grftype,ig,"gSQUARE")
        ig = ig_["Obj9"]
        arrset(pbm.grftype,ig,"gSQUARE")
        ig = ig_["Obj10"]
        arrset(pbm.grftype,ig,"gSQUARE")
        ig = ig_["Obj11"]
        arrset(pbm.grftype,ig,"gSQUARE")
        ig = ig_["Obj12"]
        arrset(pbm.grftype,ig,"gSQUARE")
        ig = ig_["Obj13"]
        arrset(pbm.grftype,ig,"gSQUARE")
        ig = ig_["Obj14"]
        arrset(pbm.grftype,ig,"gSQUARE")
        ig = ig_["Obj15"]
        arrset(pbm.grftype,ig,"gSQUARE")
        ig = ig_["Obj16"]
        arrset(pbm.grftype,ig,"gSQUARE")
        ig = ig_["Obj17"]
        arrset(pbm.grftype,ig,"gSQUARE")
        ig = ig_["Obj18"]
        arrset(pbm.grftype,ig,"gSQUARE")
        ig = ig_["Obj19"]
        arrset(pbm.grftype,ig,"gSQUARE")
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
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
        pb.pbclass = "C-SLR2-MN-19-8"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2]
        pb.conderlvl  = pbm.conderlvl;
        return pb, pbm


    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    elseif action == "gSQUARE"

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

