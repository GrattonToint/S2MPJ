function CORE1(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : CORE1
#    *********
# 
#    A problem from the exploitation of a gas transmission network
#    with consideration of the head loses in the pipes. The aim is
#    to satisfy the demand at several points in the network at a
#    minimal pressure, pumping the gas from a number of different
#    entry points.
# 
#    Sources:
#    D. De Wolf, "Optimisation de reseaux de transport de gas avec
#                 consideration des pertes de charge dans les gazoducs",
#                Ph. D. dissertation, CORE, Belgium, 1992, and
#    D. De Wolf, O. Janssens de Bisthoven and Y. Smeers,
#                "The simplex algorithm extended to piecewise linearly
#                 constrained problems II; an application to the gas
#                 transmission problem", CORE discussion paper 9103, 1991.
# 
# 
#    SDIF input: E. Loute and D. De Wolf, September 1992.
# 
#    classification = "C-LQI2-RN-65-59"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "CORE1"

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
        ig,ig_,_ = s2mpj_ii("COST",ig_)
        arrset(gtype,ig,"<>")
        ig,ig_,_ = s2mpj_ii("NODE0001",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0001")
        ig,ig_,_ = s2mpj_ii("NODE0002",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0002")
        ig,ig_,_ = s2mpj_ii("NODE0003",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0003")
        ig,ig_,_ = s2mpj_ii("NODE0004",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0004")
        ig,ig_,_ = s2mpj_ii("NODE0005",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0005")
        ig,ig_,_ = s2mpj_ii("NODE0006",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0006")
        ig,ig_,_ = s2mpj_ii("NODE0007",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0007")
        ig,ig_,_ = s2mpj_ii("NODE0008",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0008")
        ig,ig_,_ = s2mpj_ii("NODE0009",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0009")
        ig,ig_,_ = s2mpj_ii("NODE0010",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0010")
        ig,ig_,_ = s2mpj_ii("NODE0011",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0011")
        ig,ig_,_ = s2mpj_ii("NODE0012",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0012")
        ig,ig_,_ = s2mpj_ii("NODE0013",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0013")
        ig,ig_,_ = s2mpj_ii("NODE0014",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0014")
        ig,ig_,_ = s2mpj_ii("NODE0015",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0015")
        ig,ig_,_ = s2mpj_ii("NODE0016",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0016")
        ig,ig_,_ = s2mpj_ii("NODE0017",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0017")
        ig,ig_,_ = s2mpj_ii("NODE0018",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0018")
        ig,ig_,_ = s2mpj_ii("NODE0019",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0019")
        ig,ig_,_ = s2mpj_ii("NODE0020",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0020")
        ig,ig_,_ = s2mpj_ii("ARC00001",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00001")
        ig,ig_,_ = s2mpj_ii("ARC00002",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00002")
        ig,ig_,_ = s2mpj_ii("ARC00003",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00003")
        ig,ig_,_ = s2mpj_ii("ARC00004",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00004")
        ig,ig_,_ = s2mpj_ii("ARC00005",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00005")
        ig,ig_,_ = s2mpj_ii("ARC00006",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00006")
        ig,ig_,_ = s2mpj_ii("ARC00007",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00007")
        ig,ig_,_ = s2mpj_ii("ARC00008",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00008")
        ig,ig_,_ = s2mpj_ii("ARC00009",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00009")
        ig,ig_,_ = s2mpj_ii("ARC00010",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"ARC00010")
        ig,ig_,_ = s2mpj_ii("ARC00011",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"ARC00011")
        ig,ig_,_ = s2mpj_ii("ARC00012",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00012")
        ig,ig_,_ = s2mpj_ii("ARC00013",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00013")
        ig,ig_,_ = s2mpj_ii("ARC00014",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00014")
        ig,ig_,_ = s2mpj_ii("ARC00015",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00015")
        ig,ig_,_ = s2mpj_ii("ARC00016",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00016")
        ig,ig_,_ = s2mpj_ii("ARC00017",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00017")
        ig,ig_,_ = s2mpj_ii("ARC00018",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00018")
        ig,ig_,_ = s2mpj_ii("ARC00019",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00019")
        ig,ig_,_ = s2mpj_ii("ARC00020",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00020")
        ig,ig_,_ = s2mpj_ii("ARC00021",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00021")
        ig,ig_,_ = s2mpj_ii("ARC00022",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"ARC00022")
        ig,ig_,_ = s2mpj_ii("ARC00023",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00023")
        ig,ig_,_ = s2mpj_ii("ARC00024",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00024")
        ig,ig_,_ = s2mpj_ii("REGIO001",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"REGIO001")
        ig,ig_,_ = s2mpj_ii("REGIO002",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"REGIO002")
        ig,ig_,_ = s2mpj_ii("REGIO003",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"REGIO003")
        ig,ig_,_ = s2mpj_ii("REGIO004",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"REGIO004")
        ig,ig_,_ = s2mpj_ii("REGIO005",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"REGIO005")
        ig,ig_,_ = s2mpj_ii("REGIO006",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"REGIO006")
        ig,ig_,_ = s2mpj_ii("REGIO007",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"REGIO007")
        ig,ig_,_ = s2mpj_ii("REGIO008",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"REGIO008")
        ig,ig_,_ = s2mpj_ii("REGIO009",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"REGIO009")
        ig,ig_,_ = s2mpj_ii("PROD0001",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"PROD0001")
        ig,ig_,_ = s2mpj_ii("PROD0002",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"PROD0002")
        ig,ig_,_ = s2mpj_ii("PROD0003",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"PROD0003")
        ig,ig_,_ = s2mpj_ii("PROD0004",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"PROD0004")
        ig,ig_,_ = s2mpj_ii("PROD0005",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"PROD0005")
        ig,ig_,_ = s2mpj_ii("PROD0006",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"PROD0006")
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        ngrp   = length(ig_)
        iv,ix_,_ = s2mpj_ii("FLOW0001",ix_)
        arrset(pb.xnames,iv,"FLOW0001")
        ig = ig_["NODE0001"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0002"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0002",ix_)
        arrset(pb.xnames,iv,"FLOW0002")
        ig = ig_["NODE0001"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0002"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0003",ix_)
        arrset(pb.xnames,iv,"FLOW0003")
        ig = ig_["NODE0002"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0003"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0004",ix_)
        arrset(pb.xnames,iv,"FLOW0004")
        ig = ig_["NODE0002"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0003"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0005",ix_)
        arrset(pb.xnames,iv,"FLOW0005")
        ig = ig_["NODE0003"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0004"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0006",ix_)
        arrset(pb.xnames,iv,"FLOW0006")
        ig = ig_["NODE0005"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0006"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0007",ix_)
        arrset(pb.xnames,iv,"FLOW0007")
        ig = ig_["NODE0006"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0007"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0008",ix_)
        arrset(pb.xnames,iv,"FLOW0008")
        ig = ig_["NODE0007"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0004"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0009",ix_)
        arrset(pb.xnames,iv,"FLOW0009")
        ig = ig_["NODE0004"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0014"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0010",ix_)
        arrset(pb.xnames,iv,"FLOW0010")
        ig = ig_["NODE0008"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0009"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0011",ix_)
        arrset(pb.xnames,iv,"FLOW0011")
        ig = ig_["NODE0008"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0009"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0012",ix_)
        arrset(pb.xnames,iv,"FLOW0012")
        ig = ig_["NODE0009"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0010"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0013",ix_)
        arrset(pb.xnames,iv,"FLOW0013")
        ig = ig_["NODE0009"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0010"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0014",ix_)
        arrset(pb.xnames,iv,"FLOW0014")
        ig = ig_["NODE0010"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0011"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0015",ix_)
        arrset(pb.xnames,iv,"FLOW0015")
        ig = ig_["NODE0010"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0011"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0016",ix_)
        arrset(pb.xnames,iv,"FLOW0016")
        ig = ig_["NODE0011"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0012"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0017",ix_)
        arrset(pb.xnames,iv,"FLOW0017")
        ig = ig_["NODE0012"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0013"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0018",ix_)
        arrset(pb.xnames,iv,"FLOW0018")
        ig = ig_["NODE0013"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0014"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0019",ix_)
        arrset(pb.xnames,iv,"FLOW0019")
        ig = ig_["NODE0014"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0015"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0020",ix_)
        arrset(pb.xnames,iv,"FLOW0020")
        ig = ig_["NODE0015"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0016"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0021",ix_)
        arrset(pb.xnames,iv,"FLOW0021")
        ig = ig_["NODE0011"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0017"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0022",ix_)
        arrset(pb.xnames,iv,"FLOW0022")
        ig = ig_["NODE0017"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0018"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0023",ix_)
        arrset(pb.xnames,iv,"FLOW0023")
        ig = ig_["NODE0018"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0019"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0024",ix_)
        arrset(pb.xnames,iv,"FLOW0024")
        ig = ig_["NODE0019"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0020"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00001",ix_)
        arrset(pb.xnames,iv,"DEM00001")
        ig = ig_["NODE0003"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO001"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00002",ix_)
        arrset(pb.xnames,iv,"DEM00002")
        ig = ig_["NODE0006"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO002"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00003",ix_)
        arrset(pb.xnames,iv,"DEM00003")
        ig = ig_["NODE0007"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO003"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00004",ix_)
        arrset(pb.xnames,iv,"DEM00004")
        ig = ig_["NODE0010"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO004"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00005",ix_)
        arrset(pb.xnames,iv,"DEM00005")
        ig = ig_["NODE0012"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO005"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00006",ix_)
        arrset(pb.xnames,iv,"DEM00006")
        ig = ig_["NODE0015"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO006"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00007",ix_)
        arrset(pb.xnames,iv,"DEM00007")
        ig = ig_["NODE0016"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO007"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00008",ix_)
        arrset(pb.xnames,iv,"DEM00008")
        ig = ig_["NODE0019"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO008"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00009",ix_)
        arrset(pb.xnames,iv,"DEM00009")
        ig = ig_["NODE0020"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO009"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("SUPP0001",ix_)
        arrset(pb.xnames,iv,"SUPP0001")
        ig = ig_["PROD0001"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0001"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("SUPP0002",ix_)
        arrset(pb.xnames,iv,"SUPP0002")
        ig = ig_["PROD0002"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0002"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("SUPP0003",ix_)
        arrset(pb.xnames,iv,"SUPP0003")
        ig = ig_["PROD0003"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0005"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("SUPP0004",ix_)
        arrset(pb.xnames,iv,"SUPP0004")
        ig = ig_["PROD0004"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0008"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("SUPP0005",ix_)
        arrset(pb.xnames,iv,"SUPP0005")
        ig = ig_["PROD0005"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0013"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("SUPP0006",ix_)
        arrset(pb.xnames,iv,"SUPP0006")
        ig = ig_["PROD0006"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0014"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("PROD0001",ix_)
        arrset(pb.xnames,iv,"PROD0001")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(2.28000E+00)
        ig = ig_["PROD0001"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("PROD0002",ix_)
        arrset(pb.xnames,iv,"PROD0002")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(2.28000E+00)
        ig = ig_["PROD0002"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("PROD0003",ix_)
        arrset(pb.xnames,iv,"PROD0003")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(2.28000E+00)
        ig = ig_["PROD0003"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("PROD0004",ix_)
        arrset(pb.xnames,iv,"PROD0004")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(1.68000E+00)
        ig = ig_["PROD0004"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("PROD0005",ix_)
        arrset(pb.xnames,iv,"PROD0005")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(1.68000E+00)
        ig = ig_["PROD0005"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("PROD0006",ix_)
        arrset(pb.xnames,iv,"PROD0006")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(1.68000E+00)
        ig = ig_["PROD0006"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("PI000001",ix_)
        arrset(pb.xnames,iv,"PI000001")
        ig = ig_["ARC00001"]
        pbm.A[ig,iv] += Float64(-9.07027E+00)
        iv,ix_,_ = s2mpj_ii("PI000001",ix_)
        arrset(pb.xnames,iv,"PI000001")
        ig = ig_["ARC00002"]
        pbm.A[ig,iv] += Float64(-9.07027E+00)
        iv,ix_,_ = s2mpj_ii("PI000002",ix_)
        arrset(pb.xnames,iv,"PI000002")
        ig = ig_["ARC00001"]
        pbm.A[ig,iv] += Float64(9.07027E+00)
        iv,ix_,_ = s2mpj_ii("PI000002",ix_)
        arrset(pb.xnames,iv,"PI000002")
        ig = ig_["ARC00002"]
        pbm.A[ig,iv] += Float64(9.07027E+00)
        iv,ix_,_ = s2mpj_ii("PI000002",ix_)
        arrset(pb.xnames,iv,"PI000002")
        ig = ig_["ARC00003"]
        pbm.A[ig,iv] += Float64(-6.04685E+00)
        iv,ix_,_ = s2mpj_ii("PI000002",ix_)
        arrset(pb.xnames,iv,"PI000002")
        ig = ig_["ARC00004"]
        pbm.A[ig,iv] += Float64(-6.04685E+00)
        iv,ix_,_ = s2mpj_ii("PI000003",ix_)
        arrset(pb.xnames,iv,"PI000003")
        ig = ig_["ARC00003"]
        pbm.A[ig,iv] += Float64(6.04685E+00)
        iv,ix_,_ = s2mpj_ii("PI000003",ix_)
        arrset(pb.xnames,iv,"PI000003")
        ig = ig_["ARC00004"]
        pbm.A[ig,iv] += Float64(6.04685E+00)
        iv,ix_,_ = s2mpj_ii("PI000003",ix_)
        arrset(pb.xnames,iv,"PI000003")
        ig = ig_["ARC00005"]
        pbm.A[ig,iv] += Float64(-1.39543E+00)
        iv,ix_,_ = s2mpj_ii("PI000004",ix_)
        arrset(pb.xnames,iv,"PI000004")
        ig = ig_["ARC00005"]
        pbm.A[ig,iv] += Float64(1.39543E+00)
        iv,ix_,_ = s2mpj_ii("PI000004",ix_)
        arrset(pb.xnames,iv,"PI000004")
        ig = ig_["ARC00008"]
        pbm.A[ig,iv] += Float64(2.26895E-01)
        iv,ix_,_ = s2mpj_ii("PI000004",ix_)
        arrset(pb.xnames,iv,"PI000004")
        ig = ig_["ARC00009"]
        pbm.A[ig,iv] += Float64(-6.59656E-01)
        iv,ix_,_ = s2mpj_ii("PI000005",ix_)
        arrset(pb.xnames,iv,"PI000005")
        ig = ig_["ARC00006"]
        pbm.A[ig,iv] += Float64(-1.00256E-01)
        iv,ix_,_ = s2mpj_ii("PI000006",ix_)
        arrset(pb.xnames,iv,"PI000006")
        ig = ig_["ARC00006"]
        pbm.A[ig,iv] += Float64(1.00256E-01)
        iv,ix_,_ = s2mpj_ii("PI000006",ix_)
        arrset(pb.xnames,iv,"PI000006")
        ig = ig_["ARC00007"]
        pbm.A[ig,iv] += Float64(-1.48655E-01)
        iv,ix_,_ = s2mpj_ii("PI000007",ix_)
        arrset(pb.xnames,iv,"PI000007")
        ig = ig_["ARC00007"]
        pbm.A[ig,iv] += Float64(1.48655E-01)
        iv,ix_,_ = s2mpj_ii("PI000007",ix_)
        arrset(pb.xnames,iv,"PI000007")
        ig = ig_["ARC00008"]
        pbm.A[ig,iv] += Float64(-2.26895E-01)
        iv,ix_,_ = s2mpj_ii("PI000008",ix_)
        arrset(pb.xnames,iv,"PI000008")
        ig = ig_["ARC00010"]
        pbm.A[ig,iv] += Float64(-7.25622E+00)
        iv,ix_,_ = s2mpj_ii("PI000008",ix_)
        arrset(pb.xnames,iv,"PI000008")
        ig = ig_["ARC00011"]
        pbm.A[ig,iv] += Float64(-1.08033E-01)
        iv,ix_,_ = s2mpj_ii("PI000009",ix_)
        arrset(pb.xnames,iv,"PI000009")
        ig = ig_["ARC00010"]
        pbm.A[ig,iv] += Float64(7.25622E+00)
        iv,ix_,_ = s2mpj_ii("PI000009",ix_)
        arrset(pb.xnames,iv,"PI000009")
        ig = ig_["ARC00011"]
        pbm.A[ig,iv] += Float64(1.08033E-01)
        iv,ix_,_ = s2mpj_ii("PI000009",ix_)
        arrset(pb.xnames,iv,"PI000009")
        ig = ig_["ARC00012"]
        pbm.A[ig,iv] += Float64(-1.81405E+00)
        iv,ix_,_ = s2mpj_ii("PI000009",ix_)
        arrset(pb.xnames,iv,"PI000009")
        ig = ig_["ARC00013"]
        pbm.A[ig,iv] += Float64(-2.70084E-02)
        iv,ix_,_ = s2mpj_ii("PI000010",ix_)
        arrset(pb.xnames,iv,"PI000010")
        ig = ig_["ARC00012"]
        pbm.A[ig,iv] += Float64(1.81405E+00)
        iv,ix_,_ = s2mpj_ii("PI000010",ix_)
        arrset(pb.xnames,iv,"PI000010")
        ig = ig_["ARC00013"]
        pbm.A[ig,iv] += Float64(2.70084E-02)
        iv,ix_,_ = s2mpj_ii("PI000010",ix_)
        arrset(pb.xnames,iv,"PI000010")
        ig = ig_["ARC00014"]
        pbm.A[ig,iv] += Float64(-1.45124E+00)
        iv,ix_,_ = s2mpj_ii("PI000010",ix_)
        arrset(pb.xnames,iv,"PI000010")
        ig = ig_["ARC00015"]
        pbm.A[ig,iv] += Float64(-2.16067E-02)
        iv,ix_,_ = s2mpj_ii("PI000011",ix_)
        arrset(pb.xnames,iv,"PI000011")
        ig = ig_["ARC00014"]
        pbm.A[ig,iv] += Float64(1.45124E+00)
        iv,ix_,_ = s2mpj_ii("PI000011",ix_)
        arrset(pb.xnames,iv,"PI000011")
        ig = ig_["ARC00015"]
        pbm.A[ig,iv] += Float64(2.16067E-02)
        iv,ix_,_ = s2mpj_ii("PI000011",ix_)
        arrset(pb.xnames,iv,"PI000011")
        ig = ig_["ARC00016"]
        pbm.A[ig,iv] += Float64(-8.63836E-01)
        iv,ix_,_ = s2mpj_ii("PI000011",ix_)
        arrset(pb.xnames,iv,"PI000011")
        ig = ig_["ARC00021"]
        pbm.A[ig,iv] += Float64(-5.14445E-02)
        iv,ix_,_ = s2mpj_ii("PI000012",ix_)
        arrset(pb.xnames,iv,"PI000012")
        ig = ig_["ARC00016"]
        pbm.A[ig,iv] += Float64(8.63836E-01)
        iv,ix_,_ = s2mpj_ii("PI000012",ix_)
        arrset(pb.xnames,iv,"PI000012")
        ig = ig_["ARC00017"]
        pbm.A[ig,iv] += Float64(-9.07027E-01)
        iv,ix_,_ = s2mpj_ii("PI000013",ix_)
        arrset(pb.xnames,iv,"PI000013")
        ig = ig_["ARC00017"]
        pbm.A[ig,iv] += Float64(9.07027E-01)
        iv,ix_,_ = s2mpj_ii("PI000013",ix_)
        arrset(pb.xnames,iv,"PI000013")
        ig = ig_["ARC00018"]
        pbm.A[ig,iv] += Float64(-7.25622E+00)
        iv,ix_,_ = s2mpj_ii("PI000014",ix_)
        arrset(pb.xnames,iv,"PI000014")
        ig = ig_["ARC00009"]
        pbm.A[ig,iv] += Float64(6.59656E-01)
        iv,ix_,_ = s2mpj_ii("PI000014",ix_)
        arrset(pb.xnames,iv,"PI000014")
        ig = ig_["ARC00018"]
        pbm.A[ig,iv] += Float64(7.25622E+00)
        iv,ix_,_ = s2mpj_ii("PI000014",ix_)
        arrset(pb.xnames,iv,"PI000014")
        ig = ig_["ARC00019"]
        pbm.A[ig,iv] += Float64(-3.62811E+00)
        iv,ix_,_ = s2mpj_ii("PI000015",ix_)
        arrset(pb.xnames,iv,"PI000015")
        ig = ig_["ARC00019"]
        pbm.A[ig,iv] += Float64(3.62811E+00)
        iv,ix_,_ = s2mpj_ii("PI000015",ix_)
        arrset(pb.xnames,iv,"PI000015")
        ig = ig_["ARC00020"]
        pbm.A[ig,iv] += Float64(-1.45124E+00)
        iv,ix_,_ = s2mpj_ii("PI000016",ix_)
        arrset(pb.xnames,iv,"PI000016")
        ig = ig_["ARC00020"]
        pbm.A[ig,iv] += Float64(1.45124E+00)
        iv,ix_,_ = s2mpj_ii("PI000017",ix_)
        arrset(pb.xnames,iv,"PI000017")
        ig = ig_["ARC00021"]
        pbm.A[ig,iv] += Float64(5.14445E-02)
        iv,ix_,_ = s2mpj_ii("PI000017",ix_)
        arrset(pb.xnames,iv,"PI000017")
        ig = ig_["ARC00022"]
        pbm.A[ig,iv] += Float64(-6.41977E-03)
        iv,ix_,_ = s2mpj_ii("PI000018",ix_)
        arrset(pb.xnames,iv,"PI000018")
        ig = ig_["ARC00022"]
        pbm.A[ig,iv] += Float64(6.41977E-03)
        iv,ix_,_ = s2mpj_ii("PI000018",ix_)
        arrset(pb.xnames,iv,"PI000018")
        ig = ig_["ARC00023"]
        pbm.A[ig,iv] += Float64(-1.70320E-03)
        iv,ix_,_ = s2mpj_ii("PI000019",ix_)
        arrset(pb.xnames,iv,"PI000019")
        ig = ig_["ARC00023"]
        pbm.A[ig,iv] += Float64(1.70320E-03)
        iv,ix_,_ = s2mpj_ii("PI000019",ix_)
        arrset(pb.xnames,iv,"PI000019")
        ig = ig_["ARC00024"]
        pbm.A[ig,iv] += Float64(-2.78190E-02)
        iv,ix_,_ = s2mpj_ii("PI000020",ix_)
        arrset(pb.xnames,iv,"PI000020")
        ig = ig_["ARC00024"]
        pbm.A[ig,iv] += Float64(2.78190E-02)
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
        pbm.gconst[ig_["REGIO001"]] = Float64(-3.91800E+00)
        pbm.gconst[ig_["REGIO002"]] = Float64(-4.03400E+00)
        pbm.gconst[ig_["REGIO003"]] = Float64(-5.25600E+00)
        pbm.gconst[ig_["REGIO004"]] = Float64(-6.36500E+00)
        pbm.gconst[ig_["REGIO005"]] = Float64(-2.12000E+00)
        pbm.gconst[ig_["REGIO006"]] = Float64(-6.84800E+00)
        pbm.gconst[ig_["REGIO007"]] = Float64(-1.56160E+01)
        pbm.gconst[ig_["REGIO008"]] = Float64(-2.22000E-01)
        pbm.gconst[ig_["REGIO009"]] = Float64(-1.91900E+00)
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower[ix_["FLOW0001"]] = -2.20120E+02
        pb.xupper[ix_["FLOW0001"]] = 2.20120E+02
        pb.xlower[ix_["FLOW0002"]] = -2.20120E+02
        pb.xupper[ix_["FLOW0002"]] = 2.20120E+02
        pb.xlower[ix_["FLOW0003"]] = -2.20120E+02
        pb.xupper[ix_["FLOW0003"]] = 2.20120E+02
        pb.xlower[ix_["FLOW0004"]] = -2.20120E+02
        pb.xupper[ix_["FLOW0004"]] = 2.20120E+02
        pb.xlower[ix_["FLOW0005"]] = -2.20120E+02
        pb.xupper[ix_["FLOW0005"]] = 2.20120E+02
        pb.xlower[ix_["FLOW0006"]] = -2.20120E+02
        pb.xupper[ix_["FLOW0006"]] = 2.20120E+02
        pb.xlower[ix_["FLOW0007"]] = -2.20120E+02
        pb.xupper[ix_["FLOW0007"]] = 2.20120E+02
        pb.xlower[ix_["FLOW0008"]] = -2.20120E+02
        pb.xupper[ix_["FLOW0008"]] = 2.20120E+02
        pb.xlower[ix_["FLOW0009"]] = -2.20120E+02
        pb.xupper[ix_["FLOW0009"]] = 2.20120E+02
        pb.xlower[ix_["FLOW0010"]] = 0.00000E+00
        pb.xupper[ix_["FLOW0010"]] = 2.20120E+02
        pb.xlower[ix_["FLOW0011"]] = 0.00000E+00
        pb.xupper[ix_["FLOW0011"]] = 2.20120E+02
        pb.xlower[ix_["FLOW0012"]] = -2.20120E+02
        pb.xupper[ix_["FLOW0012"]] = 2.20120E+02
        pb.xlower[ix_["FLOW0013"]] = -2.20120E+02
        pb.xupper[ix_["FLOW0013"]] = 2.20120E+02
        pb.xlower[ix_["FLOW0014"]] = -2.20120E+02
        pb.xupper[ix_["FLOW0014"]] = 2.20120E+02
        pb.xlower[ix_["FLOW0015"]] = -2.20120E+02
        pb.xupper[ix_["FLOW0015"]] = 2.20120E+02
        pb.xlower[ix_["FLOW0016"]] = -2.20120E+02
        pb.xupper[ix_["FLOW0016"]] = 2.20120E+02
        pb.xlower[ix_["FLOW0017"]] = -2.20120E+02
        pb.xupper[ix_["FLOW0017"]] = 2.20120E+02
        pb.xlower[ix_["FLOW0018"]] = -2.20120E+02
        pb.xupper[ix_["FLOW0018"]] = 2.20120E+02
        pb.xlower[ix_["FLOW0019"]] = -2.20120E+02
        pb.xupper[ix_["FLOW0019"]] = 2.20120E+02
        pb.xlower[ix_["FLOW0020"]] = -2.20120E+02
        pb.xupper[ix_["FLOW0020"]] = 2.20120E+02
        pb.xlower[ix_["FLOW0021"]] = -2.20120E+02
        pb.xupper[ix_["FLOW0021"]] = 2.20120E+02
        pb.xlower[ix_["FLOW0022"]] = 0.00000E+00
        pb.xupper[ix_["FLOW0022"]] = 2.20120E+02
        pb.xlower[ix_["FLOW0023"]] = -2.20120E+02
        pb.xupper[ix_["FLOW0023"]] = 2.20120E+02
        pb.xlower[ix_["FLOW0024"]] = -2.20120E+02
        pb.xupper[ix_["FLOW0024"]] = 2.20120E+02
        pb.xupper[ix_["SUPP0001"]] = 1.15940E+01
        pb.xupper[ix_["SUPP0002"]] = 8.40000E+00
        pb.xupper[ix_["SUPP0003"]] = 4.80000E+00
        pb.xupper[ix_["SUPP0004"]] = 2.20120E+01
        pb.xupper[ix_["SUPP0005"]] = 1.20000E+00
        pb.xupper[ix_["SUPP0006"]] = 9.60000E-01
        pb.xlower[ix_["PROD0001"]] = 8.87000E+00
        pb.xupper[ix_["PROD0001"]] = 1.15940E+01
        pb.xlower[ix_["PROD0002"]] = 0.00000E+00
        pb.xupper[ix_["PROD0002"]] = 8.40000E+00
        pb.xlower[ix_["PROD0003"]] = 0.00000E+00
        pb.xupper[ix_["PROD0003"]] = 4.80000E+00
        pb.xlower[ix_["PROD0004"]] = 2.03440E+01
        pb.xupper[ix_["PROD0004"]] = 2.20120E+01
        pb.xlower[ix_["PROD0005"]] = 0.00000E+00
        pb.xupper[ix_["PROD0005"]] = 1.20000E+00
        pb.xlower[ix_["PROD0006"]] = 0.00000E+00
        pb.xupper[ix_["PROD0006"]] = 9.60000E-01
        pb.xlower[ix_["PI000001"]] = 0.00000E+00
        pb.xupper[ix_["PI000001"]] = 5.92900E+03
        pb.xlower[ix_["PI000002"]] = 0.00000E+00
        pb.xupper[ix_["PI000002"]] = 5.92900E+03
        pb.xlower[ix_["PI000003"]] = 9.00000E+02
        pb.xupper[ix_["PI000003"]] = 6.40000E+03
        pb.xlower[ix_["PI000004"]] = 0.00000E+00
        pb.xupper[ix_["PI000004"]] = 6.40000E+03
        pb.xlower[ix_["PI000005"]] = 0.00000E+00
        pb.xupper[ix_["PI000005"]] = 5.92900E+03
        pb.xlower[ix_["PI000006"]] = 9.00000E+02
        pb.xupper[ix_["PI000006"]] = 6.40000E+03
        pb.xlower[ix_["PI000007"]] = 9.00000E+02
        pb.xupper[ix_["PI000007"]] = 6.40000E+03
        pb.xlower[ix_["PI000008"]] = 2.50000E+03
        pb.xupper[ix_["PI000008"]] = 4.38244E+03
        pb.xlower[ix_["PI000009"]] = 0.00000E+00
        pb.xupper[ix_["PI000009"]] = 4.38244E+03
        pb.xlower[ix_["PI000010"]] = 9.00000E+02
        pb.xupper[ix_["PI000010"]] = 4.38244E+03
        pb.xlower[ix_["PI000011"]] = 0.00000E+00
        pb.xupper[ix_["PI000011"]] = 4.38244E+03
        pb.xlower[ix_["PI000012"]] = 9.00000E+02
        pb.xupper[ix_["PI000012"]] = 4.38244E+03
        pb.xlower[ix_["PI000013"]] = 0.00000E+00
        pb.xupper[ix_["PI000013"]] = 4.38244E+03
        pb.xlower[ix_["PI000014"]] = 0.00000E+00
        pb.xupper[ix_["PI000014"]] = 4.38244E+03
        pb.xlower[ix_["PI000015"]] = 9.00000E+02
        pb.xupper[ix_["PI000015"]] = 4.38244E+03
        pb.xlower[ix_["PI000016"]] = 2.50000E+03
        pb.xupper[ix_["PI000016"]] = 4.38244E+03
        pb.xlower[ix_["PI000017"]] = 0.00000E+00
        pb.xupper[ix_["PI000017"]] = 4.38244E+03
        pb.xlower[ix_["PI000018"]] = 0.00000E+00
        pb.xupper[ix_["PI000018"]] = 3.96900E+03
        pb.xlower[ix_["PI000019"]] = 6.25000E+02
        pb.xupper[ix_["PI000019"]] = 4.38244E+03
        pb.xlower[ix_["PI000020"]] = 6.25000E+02
        pb.xupper[ix_["PI000020"]] = 4.38244E+03
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        pb.x0[ix_["FLOW0001"]] = Float64(5.79700E+00)
        pb.x0[ix_["FLOW0002"]] = Float64(5.79700E+00)
        pb.x0[ix_["FLOW0003"]] = Float64(9.99700E+00)
        pb.x0[ix_["FLOW0004"]] = Float64(9.99700E+00)
        pb.x0[ix_["FLOW0005"]] = Float64(1.60760E+01)
        pb.x0[ix_["FLOW0006"]] = Float64(4.80000E+00)
        pb.x0[ix_["FLOW0007"]] = Float64(7.66000E-01)
        pb.x0[ix_["FLOW0008"]] = Float64(-4.49000E+00)
        pb.x0[ix_["FLOW0009"]] = Float64(1.15860E+01)
        pb.x0[ix_["FLOW0010"]] = Float64(1.72404E+01)
        pb.x0[ix_["FLOW0011"]] = Float64(2.10363E+00)
        pb.x0[ix_["FLOW0012"]] = Float64(1.72404E+01)
        pb.x0[ix_["FLOW0013"]] = Float64(2.10363E+00)
        pb.x0[ix_["FLOW0014"]] = Float64(1.15676E+01)
        pb.x0[ix_["FLOW0015"]] = Float64(1.41145E+00)
        pb.x0[ix_["FLOW0016"]] = Float64(1.08380E+01)
        pb.x0[ix_["FLOW0017"]] = Float64(8.71800E+00)
        pb.x0[ix_["FLOW0018"]] = Float64(9.91800E+00)
        pb.x0[ix_["FLOW0019"]] = Float64(2.24640E+01)
        pb.x0[ix_["FLOW0020"]] = Float64(1.56160E+01)
        pb.x0[ix_["FLOW0021"]] = Float64(2.14100E+00)
        pb.x0[ix_["FLOW0022"]] = Float64(2.14100E+00)
        pb.x0[ix_["FLOW0023"]] = Float64(2.14100E+00)
        pb.x0[ix_["FLOW0024"]] = Float64(1.91900E+00)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eSQR", iet_)
        loaset(elftv,it,1,"X")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "F00001SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0001"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00002SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0002"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00003SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0003"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00004SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0004"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00005SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0005"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00006SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0006"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00007SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0007"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00008SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0008"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00009SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0009"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00010SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0010"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00011SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0011"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00012SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0012"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00013SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0013"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00014SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0014"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00015SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0015"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00016SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0016"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00017SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0017"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00018SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0018"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00019SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0019"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00020SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0020"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00021SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0021"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00022SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0022"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00023SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0023"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00024SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0024"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["ARC00001"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00001SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00002"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00002SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00003"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00003SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00004"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00004SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00005"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00005SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00006"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00006SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00007"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00007SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00008"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00008SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00009"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00009SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00010"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00010SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00011"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00011SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00012"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00012SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00013"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00013SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00014"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00014SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00015"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00015SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00016"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00016SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00017"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00017SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00018"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00018SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00019"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00019SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00020"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00020SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00021"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00021SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00022"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00022SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00023"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00023SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00024"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00024SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.cupper[1:pb.nle] = zeros(Float64,pb.nle)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.clower[pb.nle+pb.neq+1:pb.m] = zeros(Float64,pb.nge)
        pb.cupper[1:pb.nge] = fill(Inf,pb.nge)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-LQI2-RN-65-59"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2]
        pb.conderlvl  = pbm.conderlvl;
        return pb, pbm


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "eSQR"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        XGE0 = EV_[1]>=0.0e+0
        if XGE0
            G = 2.0e+0*EV_[1]
        end
        if !XGE0
            G = -2.0e+0*EV_[1]
        end
        if XGE0
            H = 2.0e+0
        end
        if !XGE0
            H = -2.0e+0
        end
        f_   = EV_[1]*abs(EV_[1])
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = G
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = H
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

