function CORE2(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : CORE2
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
#    classification = "C-CLQI2-RN-157-134"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 9 XI 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "CORE2"

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
        ig,ig_,_ = s2mpj_ii("NODE0021",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0021")
        ig,ig_,_ = s2mpj_ii("NODE0022",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0022")
        ig,ig_,_ = s2mpj_ii("NODE0023",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0023")
        ig,ig_,_ = s2mpj_ii("NODE0024",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0024")
        ig,ig_,_ = s2mpj_ii("NODE0025",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0025")
        ig,ig_,_ = s2mpj_ii("NODE0026",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0026")
        ig,ig_,_ = s2mpj_ii("NODE0027",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0027")
        ig,ig_,_ = s2mpj_ii("NODE0028",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0028")
        ig,ig_,_ = s2mpj_ii("NODE0029",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0029")
        ig,ig_,_ = s2mpj_ii("NODE0030",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0030")
        ig,ig_,_ = s2mpj_ii("NODE0031",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0031")
        ig,ig_,_ = s2mpj_ii("NODE0032",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0032")
        ig,ig_,_ = s2mpj_ii("NODE0033",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0033")
        ig,ig_,_ = s2mpj_ii("NODE0034",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0034")
        ig,ig_,_ = s2mpj_ii("NODE0035",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0035")
        ig,ig_,_ = s2mpj_ii("NODE0036",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0036")
        ig,ig_,_ = s2mpj_ii("NODE0037",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0037")
        ig,ig_,_ = s2mpj_ii("NODE0038",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0038")
        ig,ig_,_ = s2mpj_ii("NODE0039",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0039")
        ig,ig_,_ = s2mpj_ii("NODE0040",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0040")
        ig,ig_,_ = s2mpj_ii("NODE0041",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0041")
        ig,ig_,_ = s2mpj_ii("NODE0042",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0042")
        ig,ig_,_ = s2mpj_ii("NODE0043",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0043")
        ig,ig_,_ = s2mpj_ii("NODE0044",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0044")
        ig,ig_,_ = s2mpj_ii("NODE0045",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0045")
        ig,ig_,_ = s2mpj_ii("NODE0046",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0046")
        ig,ig_,_ = s2mpj_ii("NODE0047",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0047")
        ig,ig_,_ = s2mpj_ii("NODE0048",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"NODE0048")
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
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00010")
        ig,ig_,_ = s2mpj_ii("ARC00011",ig_)
        arrset(gtype,ig,"==")
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
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00022")
        ig,ig_,_ = s2mpj_ii("ARC00023",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00023")
        ig,ig_,_ = s2mpj_ii("ARC00024",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00024")
        ig,ig_,_ = s2mpj_ii("ARC00025",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00025")
        ig,ig_,_ = s2mpj_ii("ARC00026",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00026")
        ig,ig_,_ = s2mpj_ii("ARC00027",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00027")
        ig,ig_,_ = s2mpj_ii("ARC00028",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00028")
        ig,ig_,_ = s2mpj_ii("ARC00029",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00029")
        ig,ig_,_ = s2mpj_ii("ARC00030",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00030")
        ig,ig_,_ = s2mpj_ii("ARC00031",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00031")
        ig,ig_,_ = s2mpj_ii("ARC00032",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00032")
        ig,ig_,_ = s2mpj_ii("ARC00033",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00033")
        ig,ig_,_ = s2mpj_ii("ARC00034",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00034")
        ig,ig_,_ = s2mpj_ii("ARC00035",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00035")
        ig,ig_,_ = s2mpj_ii("ARC00036",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00036")
        ig,ig_,_ = s2mpj_ii("ARC00037",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00037")
        ig,ig_,_ = s2mpj_ii("ARC00038",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00038")
        ig,ig_,_ = s2mpj_ii("ARC00039",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00039")
        ig,ig_,_ = s2mpj_ii("ARC00040",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00040")
        ig,ig_,_ = s2mpj_ii("ARC00041",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00041")
        ig,ig_,_ = s2mpj_ii("ARC00042",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00042")
        ig,ig_,_ = s2mpj_ii("ARC00043",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00043")
        ig,ig_,_ = s2mpj_ii("ARC00044",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00044")
        ig,ig_,_ = s2mpj_ii("ARC00045",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00045")
        ig,ig_,_ = s2mpj_ii("ARC00046",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00046")
        ig,ig_,_ = s2mpj_ii("ARC00047",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00047")
        ig,ig_,_ = s2mpj_ii("ARC00048",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00048")
        ig,ig_,_ = s2mpj_ii("ARC00049",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00049")
        ig,ig_,_ = s2mpj_ii("ARC00050",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00050")
        ig,ig_,_ = s2mpj_ii("ARC00051",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00051")
        ig,ig_,_ = s2mpj_ii("ARC00052",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00052")
        ig,ig_,_ = s2mpj_ii("ARC00053",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00053")
        ig,ig_,_ = s2mpj_ii("ARC00054",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00054")
        ig,ig_,_ = s2mpj_ii("ARC00055",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00055")
        ig,ig_,_ = s2mpj_ii("ARC00056",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00056")
        ig,ig_,_ = s2mpj_ii("ARC00057",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00057")
        ig,ig_,_ = s2mpj_ii("ARC00058",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00058")
        ig,ig_,_ = s2mpj_ii("ARC00059",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00059")
        ig,ig_,_ = s2mpj_ii("ARC00060",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"ARC00060")
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
        ig,ig_,_ = s2mpj_ii("REGIO010",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"REGIO010")
        ig,ig_,_ = s2mpj_ii("REGIO011",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"REGIO011")
        ig,ig_,_ = s2mpj_ii("REGIO012",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"REGIO012")
        ig,ig_,_ = s2mpj_ii("REGIO013",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"REGIO013")
        ig,ig_,_ = s2mpj_ii("REGIO014",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"REGIO014")
        ig,ig_,_ = s2mpj_ii("REGIO015",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"REGIO015")
        ig,ig_,_ = s2mpj_ii("REGIO016",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"REGIO016")
        ig,ig_,_ = s2mpj_ii("REGIO017",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"REGIO017")
        ig,ig_,_ = s2mpj_ii("REGIO018",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"REGIO018")
        ig,ig_,_ = s2mpj_ii("REGIO019",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"REGIO019")
        ig,ig_,_ = s2mpj_ii("REGIO020",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"REGIO020")
        ig,ig_,_ = s2mpj_ii("REGIO021",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"REGIO021")
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
        ig = ig_["NODE0003"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0004"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0003",ix_)
        arrset(pb.xnames,iv,"FLOW0003")
        ig = ig_["NODE0005"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0006"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0004",ix_)
        arrset(pb.xnames,iv,"FLOW0004")
        ig = ig_["NODE0006"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0007"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0005",ix_)
        arrset(pb.xnames,iv,"FLOW0005")
        ig = ig_["NODE0008"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0009"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0006",ix_)
        arrset(pb.xnames,iv,"FLOW0006")
        ig = ig_["NODE0009"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0010"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0007",ix_)
        arrset(pb.xnames,iv,"FLOW0007")
        ig = ig_["NODE0010"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0011"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0008",ix_)
        arrset(pb.xnames,iv,"FLOW0008")
        ig = ig_["NODE0012"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0010"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0009",ix_)
        arrset(pb.xnames,iv,"FLOW0009")
        ig = ig_["NODE0008"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0013"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0010",ix_)
        arrset(pb.xnames,iv,"FLOW0010")
        ig = ig_["NODE0013"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0012"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0011",ix_)
        arrset(pb.xnames,iv,"FLOW0011")
        ig = ig_["NODE0012"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0014"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0012",ix_)
        arrset(pb.xnames,iv,"FLOW0012")
        ig = ig_["NODE0014"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0015"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0013",ix_)
        arrset(pb.xnames,iv,"FLOW0013")
        ig = ig_["NODE0016"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0017"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0014",ix_)
        arrset(pb.xnames,iv,"FLOW0014")
        ig = ig_["NODE0017"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0018"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0015",ix_)
        arrset(pb.xnames,iv,"FLOW0015")
        ig = ig_["NODE0018"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0015"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0016",ix_)
        arrset(pb.xnames,iv,"FLOW0016")
        ig = ig_["NODE0015"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0001"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0017",ix_)
        arrset(pb.xnames,iv,"FLOW0017")
        ig = ig_["NODE0019"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0020"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0018",ix_)
        arrset(pb.xnames,iv,"FLOW0018")
        ig = ig_["NODE0020"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0021"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0019",ix_)
        arrset(pb.xnames,iv,"FLOW0019")
        ig = ig_["NODE0021"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0022"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0020",ix_)
        arrset(pb.xnames,iv,"FLOW0020")
        ig = ig_["NODE0023"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0024"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0021",ix_)
        arrset(pb.xnames,iv,"FLOW0021")
        ig = ig_["NODE0024"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0025"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0022",ix_)
        arrset(pb.xnames,iv,"FLOW0022")
        ig = ig_["NODE0025"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0026"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0023",ix_)
        arrset(pb.xnames,iv,"FLOW0023")
        ig = ig_["NODE0026"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0027"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0024",ix_)
        arrset(pb.xnames,iv,"FLOW0024")
        ig = ig_["NODE0028"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0027"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0025",ix_)
        arrset(pb.xnames,iv,"FLOW0025")
        ig = ig_["NODE0029"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0028"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0026",ix_)
        arrset(pb.xnames,iv,"FLOW0026")
        ig = ig_["NODE0028"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0030"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0027",ix_)
        arrset(pb.xnames,iv,"FLOW0027")
        ig = ig_["NODE0019"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0031"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0028",ix_)
        arrset(pb.xnames,iv,"FLOW0028")
        ig = ig_["NODE0031"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0032"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0029",ix_)
        arrset(pb.xnames,iv,"FLOW0029")
        ig = ig_["NODE0032"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0033"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0030",ix_)
        arrset(pb.xnames,iv,"FLOW0030")
        ig = ig_["NODE0033"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0030"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0031",ix_)
        arrset(pb.xnames,iv,"FLOW0031")
        ig = ig_["NODE0030"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0034"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0032",ix_)
        arrset(pb.xnames,iv,"FLOW0032")
        ig = ig_["NODE0034"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0035"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0033",ix_)
        arrset(pb.xnames,iv,"FLOW0033")
        ig = ig_["NODE0035"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0024"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0034",ix_)
        arrset(pb.xnames,iv,"FLOW0034")
        ig = ig_["NODE0027"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0034"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0035",ix_)
        arrset(pb.xnames,iv,"FLOW0035")
        ig = ig_["NODE0036"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0007"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0036",ix_)
        arrset(pb.xnames,iv,"FLOW0036")
        ig = ig_["NODE0007"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0016"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0037",ix_)
        arrset(pb.xnames,iv,"FLOW0037")
        ig = ig_["NODE0037"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0038"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0038",ix_)
        arrset(pb.xnames,iv,"FLOW0038")
        ig = ig_["NODE0038"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0036"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0039",ix_)
        arrset(pb.xnames,iv,"FLOW0039")
        ig = ig_["NODE0028"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0030"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0040",ix_)
        arrset(pb.xnames,iv,"FLOW0040")
        ig = ig_["NODE0019"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0031"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0041",ix_)
        arrset(pb.xnames,iv,"FLOW0041")
        ig = ig_["NODE0031"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0032"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0042",ix_)
        arrset(pb.xnames,iv,"FLOW0042")
        ig = ig_["NODE0032"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0033"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0043",ix_)
        arrset(pb.xnames,iv,"FLOW0043")
        ig = ig_["NODE0033"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0030"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0044",ix_)
        arrset(pb.xnames,iv,"FLOW0044")
        ig = ig_["NODE0030"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0034"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0045",ix_)
        arrset(pb.xnames,iv,"FLOW0045")
        ig = ig_["NODE0034"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0035"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0046",ix_)
        arrset(pb.xnames,iv,"FLOW0046")
        ig = ig_["NODE0035"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0024"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0047",ix_)
        arrset(pb.xnames,iv,"FLOW0047")
        ig = ig_["NODE0036"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0007"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0048",ix_)
        arrset(pb.xnames,iv,"FLOW0048")
        ig = ig_["NODE0007"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0016"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0049",ix_)
        arrset(pb.xnames,iv,"FLOW0049")
        ig = ig_["NODE0037"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0039"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0050",ix_)
        arrset(pb.xnames,iv,"FLOW0050")
        ig = ig_["NODE0039"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0040"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0051",ix_)
        arrset(pb.xnames,iv,"FLOW0051")
        ig = ig_["NODE0041"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0040"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0052",ix_)
        arrset(pb.xnames,iv,"FLOW0052")
        ig = ig_["NODE0040"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0029"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0053",ix_)
        arrset(pb.xnames,iv,"FLOW0053")
        ig = ig_["NODE0029"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0042"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0054",ix_)
        arrset(pb.xnames,iv,"FLOW0054")
        ig = ig_["NODE0042"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0016"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0055",ix_)
        arrset(pb.xnames,iv,"FLOW0055")
        ig = ig_["NODE0003"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0043"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0056",ix_)
        arrset(pb.xnames,iv,"FLOW0056")
        ig = ig_["NODE0043"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0001"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0057",ix_)
        arrset(pb.xnames,iv,"FLOW0057")
        ig = ig_["NODE0044"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0045"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0058",ix_)
        arrset(pb.xnames,iv,"FLOW0058")
        ig = ig_["NODE0045"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0046"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0059",ix_)
        arrset(pb.xnames,iv,"FLOW0059")
        ig = ig_["NODE0044"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0047"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("FLOW0060",ix_)
        arrset(pb.xnames,iv,"FLOW0060")
        ig = ig_["NODE0047"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0048"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00001",ix_)
        arrset(pb.xnames,iv,"DEM00001")
        ig = ig_["NODE0044"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO001"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00002",ix_)
        arrset(pb.xnames,iv,"DEM00002")
        ig = ig_["NODE0038"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO002"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00003",ix_)
        arrset(pb.xnames,iv,"DEM00003")
        ig = ig_["NODE0045"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO002"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00004",ix_)
        arrset(pb.xnames,iv,"DEM00004")
        ig = ig_["NODE0046"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO002"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00005",ix_)
        arrset(pb.xnames,iv,"DEM00005")
        ig = ig_["NODE0004"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO003"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00006",ix_)
        arrset(pb.xnames,iv,"DEM00006")
        ig = ig_["NODE0005"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO003"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00007",ix_)
        arrset(pb.xnames,iv,"DEM00007")
        ig = ig_["NODE0002"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO004"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00008",ix_)
        arrset(pb.xnames,iv,"DEM00008")
        ig = ig_["NODE0006"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO004"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00009",ix_)
        arrset(pb.xnames,iv,"DEM00009")
        ig = ig_["NODE0011"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO004"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00010",ix_)
        arrset(pb.xnames,iv,"DEM00010")
        ig = ig_["NODE0016"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO004"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00011",ix_)
        arrset(pb.xnames,iv,"DEM00011")
        ig = ig_["NODE0036"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO004"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00012",ix_)
        arrset(pb.xnames,iv,"DEM00012")
        ig = ig_["NODE0042"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO005"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00013",ix_)
        arrset(pb.xnames,iv,"DEM00013")
        ig = ig_["NODE0047"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO005"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00014",ix_)
        arrset(pb.xnames,iv,"DEM00014")
        ig = ig_["NODE0039"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO006"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00015",ix_)
        arrset(pb.xnames,iv,"DEM00015")
        ig = ig_["NODE0040"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO006"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00016",ix_)
        arrset(pb.xnames,iv,"DEM00016")
        ig = ig_["NODE0048"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO006"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00017",ix_)
        arrset(pb.xnames,iv,"DEM00017")
        ig = ig_["NODE0041"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO007"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00018",ix_)
        arrset(pb.xnames,iv,"DEM00018")
        ig = ig_["NODE0030"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO008"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00019",ix_)
        arrset(pb.xnames,iv,"DEM00019")
        ig = ig_["NODE0028"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO009"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00020",ix_)
        arrset(pb.xnames,iv,"DEM00020")
        ig = ig_["NODE0043"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO010"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00021",ix_)
        arrset(pb.xnames,iv,"DEM00021")
        ig = ig_["NODE0008"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO011"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00022",ix_)
        arrset(pb.xnames,iv,"DEM00022")
        ig = ig_["NODE0009"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO012"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00023",ix_)
        arrset(pb.xnames,iv,"DEM00023")
        ig = ig_["NODE0013"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO012"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00024",ix_)
        arrset(pb.xnames,iv,"DEM00024")
        ig = ig_["NODE0014"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO013"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00025",ix_)
        arrset(pb.xnames,iv,"DEM00025")
        ig = ig_["NODE0017"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO013"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00026",ix_)
        arrset(pb.xnames,iv,"DEM00026")
        ig = ig_["NODE0018"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO013"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00027",ix_)
        arrset(pb.xnames,iv,"DEM00027")
        ig = ig_["NODE0020"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO014"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00028",ix_)
        arrset(pb.xnames,iv,"DEM00028")
        ig = ig_["NODE0021"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO014"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00029",ix_)
        arrset(pb.xnames,iv,"DEM00029")
        ig = ig_["NODE0031"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO015"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00030",ix_)
        arrset(pb.xnames,iv,"DEM00030")
        ig = ig_["NODE0022"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO016"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00031",ix_)
        arrset(pb.xnames,iv,"DEM00031")
        ig = ig_["NODE0022"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO017"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00032",ix_)
        arrset(pb.xnames,iv,"DEM00032")
        ig = ig_["NODE0033"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO018"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00033",ix_)
        arrset(pb.xnames,iv,"DEM00033")
        ig = ig_["NODE0025"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO019"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00034",ix_)
        arrset(pb.xnames,iv,"DEM00034")
        ig = ig_["NODE0026"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO019"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00035",ix_)
        arrset(pb.xnames,iv,"DEM00035")
        ig = ig_["NODE0035"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO019"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00036",ix_)
        arrset(pb.xnames,iv,"DEM00036")
        ig = ig_["NODE0023"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO020"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("DEM00037",ix_)
        arrset(pb.xnames,iv,"DEM00037")
        ig = ig_["NODE0023"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["REGIO021"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("SUPP0001",ix_)
        arrset(pb.xnames,iv,"SUPP0001")
        ig = ig_["PROD0001"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0021"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("SUPP0002",ix_)
        arrset(pb.xnames,iv,"SUPP0002")
        ig = ig_["PROD0002"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0044"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("SUPP0003",ix_)
        arrset(pb.xnames,iv,"SUPP0003")
        ig = ig_["PROD0003"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0041"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("SUPP0004",ix_)
        arrset(pb.xnames,iv,"SUPP0004")
        ig = ig_["PROD0004"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0003"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("SUPP0005",ix_)
        arrset(pb.xnames,iv,"SUPP0005")
        ig = ig_["PROD0004"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0008"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("SUPP0006",ix_)
        arrset(pb.xnames,iv,"SUPP0006")
        ig = ig_["PROD0004"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0023"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("SUPP0007",ix_)
        arrset(pb.xnames,iv,"SUPP0007")
        ig = ig_["PROD0005"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["NODE0037"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("PROD0001",ix_)
        arrset(pb.xnames,iv,"PROD0001")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["PROD0001"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("PROD0002",ix_)
        arrset(pb.xnames,iv,"PROD0002")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["PROD0002"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("PROD0003",ix_)
        arrset(pb.xnames,iv,"PROD0003")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["PROD0003"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("PROD0004",ix_)
        arrset(pb.xnames,iv,"PROD0004")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["PROD0004"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("PROD0005",ix_)
        arrset(pb.xnames,iv,"PROD0005")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(1.00000E+00)
        ig = ig_["PROD0005"]
        pbm.A[ig,iv] += Float64(-1.00000E+00)
        iv,ix_,_ = s2mpj_ii("PI000001",ix_)
        arrset(pb.xnames,iv,"PI000001")
        ig = ig_["ARC00001"]
        pbm.A[ig,iv] += Float64(-6.97133E-01)
        iv,ix_,_ = s2mpj_ii("PI000001",ix_)
        arrset(pb.xnames,iv,"PI000001")
        ig = ig_["ARC00016"]
        pbm.A[ig,iv] += Float64(1.81479E-02)
        iv,ix_,_ = s2mpj_ii("PI000001",ix_)
        arrset(pb.xnames,iv,"PI000001")
        ig = ig_["ARC00056"]
        pbm.A[ig,iv] += Float64(2.85181E-02)
        iv,ix_,_ = s2mpj_ii("PI000002",ix_)
        arrset(pb.xnames,iv,"PI000002")
        ig = ig_["ARC00001"]
        pbm.A[ig,iv] += Float64(6.97133E-01)
        iv,ix_,_ = s2mpj_ii("PI000003",ix_)
        arrset(pb.xnames,iv,"PI000003")
        ig = ig_["ARC00002"]
        pbm.A[ig,iv] += Float64(-2.65927E-03)
        iv,ix_,_ = s2mpj_ii("PI000003",ix_)
        arrset(pb.xnames,iv,"PI000003")
        ig = ig_["ARC00055"]
        pbm.A[ig,iv] += Float64(-3.99253E-02)
        iv,ix_,_ = s2mpj_ii("PI000004",ix_)
        arrset(pb.xnames,iv,"PI000004")
        ig = ig_["ARC00002"]
        pbm.A[ig,iv] += Float64(2.65927E-03)
        iv,ix_,_ = s2mpj_ii("PI000005",ix_)
        arrset(pb.xnames,iv,"PI000005")
        ig = ig_["ARC00003"]
        pbm.A[ig,iv] += Float64(-1.72371E-01)
        iv,ix_,_ = s2mpj_ii("PI000006",ix_)
        arrset(pb.xnames,iv,"PI000006")
        ig = ig_["ARC00003"]
        pbm.A[ig,iv] += Float64(1.72371E-01)
        iv,ix_,_ = s2mpj_ii("PI000006",ix_)
        arrset(pb.xnames,iv,"PI000006")
        ig = ig_["ARC00004"]
        pbm.A[ig,iv] += Float64(-2.58556E-01)
        iv,ix_,_ = s2mpj_ii("PI000007",ix_)
        arrset(pb.xnames,iv,"PI000007")
        ig = ig_["ARC00004"]
        pbm.A[ig,iv] += Float64(2.58556E-01)
        iv,ix_,_ = s2mpj_ii("PI000007",ix_)
        arrset(pb.xnames,iv,"PI000007")
        ig = ig_["ARC00035"]
        pbm.A[ig,iv] += Float64(9.98133E-02)
        iv,ix_,_ = s2mpj_ii("PI000007",ix_)
        arrset(pb.xnames,iv,"PI000007")
        ig = ig_["ARC00036"]
        pbm.A[ig,iv] += Float64(-4.99067E-02)
        iv,ix_,_ = s2mpj_ii("PI000007",ix_)
        arrset(pb.xnames,iv,"PI000007")
        ig = ig_["ARC00047"]
        pbm.A[ig,iv] += Float64(2.09022E+00)
        iv,ix_,_ = s2mpj_ii("PI000007",ix_)
        arrset(pb.xnames,iv,"PI000007")
        ig = ig_["ARC00048"]
        pbm.A[ig,iv] += Float64(-6.96739E-01)
        iv,ix_,_ = s2mpj_ii("PI000008",ix_)
        arrset(pb.xnames,iv,"PI000008")
        ig = ig_["ARC00005"]
        pbm.A[ig,iv] += Float64(-1.90020E-01)
        iv,ix_,_ = s2mpj_ii("PI000008",ix_)
        arrset(pb.xnames,iv,"PI000008")
        ig = ig_["ARC00009"]
        pbm.A[ig,iv] += Float64(-6.11442E-03)
        iv,ix_,_ = s2mpj_ii("PI000009",ix_)
        arrset(pb.xnames,iv,"PI000009")
        ig = ig_["ARC00005"]
        pbm.A[ig,iv] += Float64(1.90020E-01)
        iv,ix_,_ = s2mpj_ii("PI000009",ix_)
        arrset(pb.xnames,iv,"PI000009")
        ig = ig_["ARC00006"]
        pbm.A[ig,iv] += Float64(-4.18044E+00)
        iv,ix_,_ = s2mpj_ii("PI000010",ix_)
        arrset(pb.xnames,iv,"PI000010")
        ig = ig_["ARC00006"]
        pbm.A[ig,iv] += Float64(4.18044E+00)
        iv,ix_,_ = s2mpj_ii("PI000010",ix_)
        arrset(pb.xnames,iv,"PI000010")
        ig = ig_["ARC00007"]
        pbm.A[ig,iv] += Float64(-3.21572E-01)
        iv,ix_,_ = s2mpj_ii("PI000010",ix_)
        arrset(pb.xnames,iv,"PI000010")
        ig = ig_["ARC00008"]
        pbm.A[ig,iv] += Float64(4.32263E-02)
        iv,ix_,_ = s2mpj_ii("PI000011",ix_)
        arrset(pb.xnames,iv,"PI000011")
        ig = ig_["ARC00007"]
        pbm.A[ig,iv] += Float64(3.21572E-01)
        iv,ix_,_ = s2mpj_ii("PI000012",ix_)
        arrset(pb.xnames,iv,"PI000012")
        ig = ig_["ARC00008"]
        pbm.A[ig,iv] += Float64(-4.32263E-02)
        iv,ix_,_ = s2mpj_ii("PI000012",ix_)
        arrset(pb.xnames,iv,"PI000012")
        ig = ig_["ARC00010"]
        pbm.A[ig,iv] += Float64(1.16174E-02)
        iv,ix_,_ = s2mpj_ii("PI000012",ix_)
        arrset(pb.xnames,iv,"PI000012")
        ig = ig_["ARC00011"]
        pbm.A[ig,iv] += Float64(-2.59358E-01)
        iv,ix_,_ = s2mpj_ii("PI000013",ix_)
        arrset(pb.xnames,iv,"PI000013")
        ig = ig_["ARC00009"]
        pbm.A[ig,iv] += Float64(6.11442E-03)
        iv,ix_,_ = s2mpj_ii("PI000013",ix_)
        arrset(pb.xnames,iv,"PI000013")
        ig = ig_["ARC00010"]
        pbm.A[ig,iv] += Float64(-1.16174E-02)
        iv,ix_,_ = s2mpj_ii("PI000014",ix_)
        arrset(pb.xnames,iv,"PI000014")
        ig = ig_["ARC00011"]
        pbm.A[ig,iv] += Float64(2.59358E-01)
        iv,ix_,_ = s2mpj_ii("PI000014",ix_)
        arrset(pb.xnames,iv,"PI000014")
        ig = ig_["ARC00012"]
        pbm.A[ig,iv] += Float64(-2.59358E-01)
        iv,ix_,_ = s2mpj_ii("PI000015",ix_)
        arrset(pb.xnames,iv,"PI000015")
        ig = ig_["ARC00012"]
        pbm.A[ig,iv] += Float64(2.59358E-01)
        iv,ix_,_ = s2mpj_ii("PI000015",ix_)
        arrset(pb.xnames,iv,"PI000015")
        ig = ig_["ARC00015"]
        pbm.A[ig,iv] += Float64(2.85181E-02)
        iv,ix_,_ = s2mpj_ii("PI000015",ix_)
        arrset(pb.xnames,iv,"PI000015")
        ig = ig_["ARC00016"]
        pbm.A[ig,iv] += Float64(-1.81479E-02)
        iv,ix_,_ = s2mpj_ii("PI000016",ix_)
        arrset(pb.xnames,iv,"PI000016")
        ig = ig_["ARC00013"]
        pbm.A[ig,iv] += Float64(-2.85181E-02)
        iv,ix_,_ = s2mpj_ii("PI000016",ix_)
        arrset(pb.xnames,iv,"PI000016")
        ig = ig_["ARC00036"]
        pbm.A[ig,iv] += Float64(4.99067E-02)
        iv,ix_,_ = s2mpj_ii("PI000016",ix_)
        arrset(pb.xnames,iv,"PI000016")
        ig = ig_["ARC00048"]
        pbm.A[ig,iv] += Float64(6.96739E-01)
        iv,ix_,_ = s2mpj_ii("PI000016",ix_)
        arrset(pb.xnames,iv,"PI000016")
        ig = ig_["ARC00054"]
        pbm.A[ig,iv] += Float64(1.36747E-01)
        iv,ix_,_ = s2mpj_ii("PI000017",ix_)
        arrset(pb.xnames,iv,"PI000017")
        ig = ig_["ARC00013"]
        pbm.A[ig,iv] += Float64(2.85181E-02)
        iv,ix_,_ = s2mpj_ii("PI000017",ix_)
        arrset(pb.xnames,iv,"PI000017")
        ig = ig_["ARC00014"]
        pbm.A[ig,iv] += Float64(-3.99253E-02)
        iv,ix_,_ = s2mpj_ii("PI000018",ix_)
        arrset(pb.xnames,iv,"PI000018")
        ig = ig_["ARC00014"]
        pbm.A[ig,iv] += Float64(3.99253E-02)
        iv,ix_,_ = s2mpj_ii("PI000018",ix_)
        arrset(pb.xnames,iv,"PI000018")
        ig = ig_["ARC00015"]
        pbm.A[ig,iv] += Float64(-2.85181E-02)
        iv,ix_,_ = s2mpj_ii("PI000019",ix_)
        arrset(pb.xnames,iv,"PI000019")
        ig = ig_["ARC00017"]
        pbm.A[ig,iv] += Float64(-8.64525E-02)
        iv,ix_,_ = s2mpj_ii("PI000019",ix_)
        arrset(pb.xnames,iv,"PI000019")
        ig = ig_["ARC00027"]
        pbm.A[ig,iv] += Float64(-3.32711E-02)
        iv,ix_,_ = s2mpj_ii("PI000019",ix_)
        arrset(pb.xnames,iv,"PI000019")
        ig = ig_["ARC00040"]
        pbm.A[ig,iv] += Float64(-1.93623E-02)
        iv,ix_,_ = s2mpj_ii("PI000020",ix_)
        arrset(pb.xnames,iv,"PI000020")
        ig = ig_["ARC00017"]
        pbm.A[ig,iv] += Float64(8.64525E-02)
        iv,ix_,_ = s2mpj_ii("PI000020",ix_)
        arrset(pb.xnames,iv,"PI000020")
        ig = ig_["ARC00018"]
        pbm.A[ig,iv] += Float64(-7.41022E-02)
        iv,ix_,_ = s2mpj_ii("PI000021",ix_)
        arrset(pb.xnames,iv,"PI000021")
        ig = ig_["ARC00018"]
        pbm.A[ig,iv] += Float64(7.41022E-02)
        iv,ix_,_ = s2mpj_ii("PI000021",ix_)
        arrset(pb.xnames,iv,"PI000021")
        ig = ig_["ARC00019"]
        pbm.A[ig,iv] += Float64(-1.90020E-01)
        iv,ix_,_ = s2mpj_ii("PI000022",ix_)
        arrset(pb.xnames,iv,"PI000022")
        ig = ig_["ARC00019"]
        pbm.A[ig,iv] += Float64(1.90020E-01)
        iv,ix_,_ = s2mpj_ii("PI000023",ix_)
        arrset(pb.xnames,iv,"PI000023")
        ig = ig_["ARC00020"]
        pbm.A[ig,iv] += Float64(-2.35780E-02)
        iv,ix_,_ = s2mpj_ii("PI000024",ix_)
        arrset(pb.xnames,iv,"PI000024")
        ig = ig_["ARC00020"]
        pbm.A[ig,iv] += Float64(2.35780E-02)
        iv,ix_,_ = s2mpj_ii("PI000024",ix_)
        arrset(pb.xnames,iv,"PI000024")
        ig = ig_["ARC00021"]
        pbm.A[ig,iv] += Float64(-4.18044E+00)
        iv,ix_,_ = s2mpj_ii("PI000024",ix_)
        arrset(pb.xnames,iv,"PI000024")
        ig = ig_["ARC00033"]
        pbm.A[ig,iv] += Float64(6.07700E-03)
        iv,ix_,_ = s2mpj_ii("PI000024",ix_)
        arrset(pb.xnames,iv,"PI000024")
        ig = ig_["ARC00046"]
        pbm.A[ig,iv] += Float64(6.07700E-03)
        iv,ix_,_ = s2mpj_ii("PI000025",ix_)
        arrset(pb.xnames,iv,"PI000025")
        ig = ig_["ARC00021"]
        pbm.A[ig,iv] += Float64(4.18044E+00)
        iv,ix_,_ = s2mpj_ii("PI000025",ix_)
        arrset(pb.xnames,iv,"PI000025")
        ig = ig_["ARC00022"]
        pbm.A[ig,iv] += Float64(-5.22555E-01)
        iv,ix_,_ = s2mpj_ii("PI000026",ix_)
        arrset(pb.xnames,iv,"PI000026")
        ig = ig_["ARC00022"]
        pbm.A[ig,iv] += Float64(5.22555E-01)
        iv,ix_,_ = s2mpj_ii("PI000026",ix_)
        arrset(pb.xnames,iv,"PI000026")
        ig = ig_["ARC00023"]
        pbm.A[ig,iv] += Float64(-1.39348E+00)
        iv,ix_,_ = s2mpj_ii("PI000027",ix_)
        arrset(pb.xnames,iv,"PI000027")
        ig = ig_["ARC00023"]
        pbm.A[ig,iv] += Float64(1.39348E+00)
        iv,ix_,_ = s2mpj_ii("PI000027",ix_)
        arrset(pb.xnames,iv,"PI000027")
        ig = ig_["ARC00024"]
        pbm.A[ig,iv] += Float64(6.96739E-01)
        iv,ix_,_ = s2mpj_ii("PI000027",ix_)
        arrset(pb.xnames,iv,"PI000027")
        ig = ig_["ARC00034"]
        pbm.A[ig,iv] += Float64(-4.10240E-01)
        iv,ix_,_ = s2mpj_ii("PI000028",ix_)
        arrset(pb.xnames,iv,"PI000028")
        ig = ig_["ARC00024"]
        pbm.A[ig,iv] += Float64(-6.96739E-01)
        iv,ix_,_ = s2mpj_ii("PI000028",ix_)
        arrset(pb.xnames,iv,"PI000028")
        ig = ig_["ARC00025"]
        pbm.A[ig,iv] += Float64(4.18044E-01)
        iv,ix_,_ = s2mpj_ii("PI000028",ix_)
        arrset(pb.xnames,iv,"PI000028")
        ig = ig_["ARC00026"]
        pbm.A[ig,iv] += Float64(-6.33373E-04)
        iv,ix_,_ = s2mpj_ii("PI000028",ix_)
        arrset(pb.xnames,iv,"PI000028")
        ig = ig_["ARC00039"]
        pbm.A[ig,iv] += Float64(-1.99627E-02)
        iv,ix_,_ = s2mpj_ii("PI000029",ix_)
        arrset(pb.xnames,iv,"PI000029")
        ig = ig_["ARC00025"]
        pbm.A[ig,iv] += Float64(-4.18044E-01)
        iv,ix_,_ = s2mpj_ii("PI000029",ix_)
        arrset(pb.xnames,iv,"PI000029")
        ig = ig_["ARC00052"]
        pbm.A[ig,iv] += Float64(3.21572E-01)
        iv,ix_,_ = s2mpj_ii("PI000029",ix_)
        arrset(pb.xnames,iv,"PI000029")
        ig = ig_["ARC00053"]
        pbm.A[ig,iv] += Float64(-1.64096E-01)
        iv,ix_,_ = s2mpj_ii("PI000030",ix_)
        arrset(pb.xnames,iv,"PI000030")
        ig = ig_["ARC00026"]
        pbm.A[ig,iv] += Float64(6.33373E-04)
        iv,ix_,_ = s2mpj_ii("PI000030",ix_)
        arrset(pb.xnames,iv,"PI000030")
        ig = ig_["ARC00030"]
        pbm.A[ig,iv] += Float64(3.32711E-02)
        iv,ix_,_ = s2mpj_ii("PI000030",ix_)
        arrset(pb.xnames,iv,"PI000030")
        ig = ig_["ARC00031"]
        pbm.A[ig,iv] += Float64(-9.54957E-03)
        iv,ix_,_ = s2mpj_ii("PI000030",ix_)
        arrset(pb.xnames,iv,"PI000030")
        ig = ig_["ARC00039"]
        pbm.A[ig,iv] += Float64(1.99627E-02)
        iv,ix_,_ = s2mpj_ii("PI000030",ix_)
        arrset(pb.xnames,iv,"PI000030")
        ig = ig_["ARC00043"]
        pbm.A[ig,iv] += Float64(1.93623E-02)
        iv,ix_,_ = s2mpj_ii("PI000030",ix_)
        arrset(pb.xnames,iv,"PI000030")
        ig = ig_["ARC00044"]
        pbm.A[ig,iv] += Float64(-9.54957E-03)
        iv,ix_,_ = s2mpj_ii("PI000031",ix_)
        arrset(pb.xnames,iv,"PI000031")
        ig = ig_["ARC00027"]
        pbm.A[ig,iv] += Float64(3.32711E-02)
        iv,ix_,_ = s2mpj_ii("PI000031",ix_)
        arrset(pb.xnames,iv,"PI000031")
        ig = ig_["ARC00028"]
        pbm.A[ig,iv] += Float64(-9.98133E-02)
        iv,ix_,_ = s2mpj_ii("PI000031",ix_)
        arrset(pb.xnames,iv,"PI000031")
        ig = ig_["ARC00040"]
        pbm.A[ig,iv] += Float64(1.93623E-02)
        iv,ix_,_ = s2mpj_ii("PI000031",ix_)
        arrset(pb.xnames,iv,"PI000031")
        ig = ig_["ARC00041"]
        pbm.A[ig,iv] += Float64(-5.80870E-02)
        iv,ix_,_ = s2mpj_ii("PI000032",ix_)
        arrset(pb.xnames,iv,"PI000032")
        ig = ig_["ARC00028"]
        pbm.A[ig,iv] += Float64(9.98133E-02)
        iv,ix_,_ = s2mpj_ii("PI000032",ix_)
        arrset(pb.xnames,iv,"PI000032")
        ig = ig_["ARC00029"]
        pbm.A[ig,iv] += Float64(-3.32711E-02)
        iv,ix_,_ = s2mpj_ii("PI000032",ix_)
        arrset(pb.xnames,iv,"PI000032")
        ig = ig_["ARC00041"]
        pbm.A[ig,iv] += Float64(5.80870E-02)
        iv,ix_,_ = s2mpj_ii("PI000032",ix_)
        arrset(pb.xnames,iv,"PI000032")
        ig = ig_["ARC00042"]
        pbm.A[ig,iv] += Float64(-1.93623E-02)
        iv,ix_,_ = s2mpj_ii("PI000033",ix_)
        arrset(pb.xnames,iv,"PI000033")
        ig = ig_["ARC00029"]
        pbm.A[ig,iv] += Float64(3.32711E-02)
        iv,ix_,_ = s2mpj_ii("PI000033",ix_)
        arrset(pb.xnames,iv,"PI000033")
        ig = ig_["ARC00030"]
        pbm.A[ig,iv] += Float64(-3.32711E-02)
        iv,ix_,_ = s2mpj_ii("PI000033",ix_)
        arrset(pb.xnames,iv,"PI000033")
        ig = ig_["ARC00042"]
        pbm.A[ig,iv] += Float64(1.93623E-02)
        iv,ix_,_ = s2mpj_ii("PI000033",ix_)
        arrset(pb.xnames,iv,"PI000033")
        ig = ig_["ARC00043"]
        pbm.A[ig,iv] += Float64(-1.93623E-02)
        iv,ix_,_ = s2mpj_ii("PI000034",ix_)
        arrset(pb.xnames,iv,"PI000034")
        ig = ig_["ARC00031"]
        pbm.A[ig,iv] += Float64(9.54957E-03)
        iv,ix_,_ = s2mpj_ii("PI000034",ix_)
        arrset(pb.xnames,iv,"PI000034")
        ig = ig_["ARC00032"]
        pbm.A[ig,iv] += Float64(-6.68470E-02)
        iv,ix_,_ = s2mpj_ii("PI000034",ix_)
        arrset(pb.xnames,iv,"PI000034")
        ig = ig_["ARC00034"]
        pbm.A[ig,iv] += Float64(4.10240E-01)
        iv,ix_,_ = s2mpj_ii("PI000034",ix_)
        arrset(pb.xnames,iv,"PI000034")
        ig = ig_["ARC00044"]
        pbm.A[ig,iv] += Float64(9.54957E-03)
        iv,ix_,_ = s2mpj_ii("PI000034",ix_)
        arrset(pb.xnames,iv,"PI000034")
        ig = ig_["ARC00045"]
        pbm.A[ig,iv] += Float64(-6.68470E-02)
        iv,ix_,_ = s2mpj_ii("PI000035",ix_)
        arrset(pb.xnames,iv,"PI000035")
        ig = ig_["ARC00032"]
        pbm.A[ig,iv] += Float64(6.68470E-02)
        iv,ix_,_ = s2mpj_ii("PI000035",ix_)
        arrset(pb.xnames,iv,"PI000035")
        ig = ig_["ARC00033"]
        pbm.A[ig,iv] += Float64(-6.07700E-03)
        iv,ix_,_ = s2mpj_ii("PI000035",ix_)
        arrset(pb.xnames,iv,"PI000035")
        ig = ig_["ARC00045"]
        pbm.A[ig,iv] += Float64(6.68470E-02)
        iv,ix_,_ = s2mpj_ii("PI000035",ix_)
        arrset(pb.xnames,iv,"PI000035")
        ig = ig_["ARC00046"]
        pbm.A[ig,iv] += Float64(-6.07700E-03)
        iv,ix_,_ = s2mpj_ii("PI000036",ix_)
        arrset(pb.xnames,iv,"PI000036")
        ig = ig_["ARC00035"]
        pbm.A[ig,iv] += Float64(-9.98133E-02)
        iv,ix_,_ = s2mpj_ii("PI000036",ix_)
        arrset(pb.xnames,iv,"PI000036")
        ig = ig_["ARC00038"]
        pbm.A[ig,iv] += Float64(6.96739E-01)
        iv,ix_,_ = s2mpj_ii("PI000036",ix_)
        arrset(pb.xnames,iv,"PI000036")
        ig = ig_["ARC00047"]
        pbm.A[ig,iv] += Float64(-2.09022E+00)
        iv,ix_,_ = s2mpj_ii("PI000037",ix_)
        arrset(pb.xnames,iv,"PI000037")
        ig = ig_["ARC00037"]
        pbm.A[ig,iv] += Float64(-3.48370E-01)
        iv,ix_,_ = s2mpj_ii("PI000037",ix_)
        arrset(pb.xnames,iv,"PI000037")
        ig = ig_["ARC00049"]
        pbm.A[ig,iv] += Float64(-1.43027E-02)
        iv,ix_,_ = s2mpj_ii("PI000038",ix_)
        arrset(pb.xnames,iv,"PI000038")
        ig = ig_["ARC00037"]
        pbm.A[ig,iv] += Float64(3.48370E-01)
        iv,ix_,_ = s2mpj_ii("PI000038",ix_)
        arrset(pb.xnames,iv,"PI000038")
        ig = ig_["ARC00038"]
        pbm.A[ig,iv] += Float64(-6.96739E-01)
        iv,ix_,_ = s2mpj_ii("PI000039",ix_)
        arrset(pb.xnames,iv,"PI000039")
        ig = ig_["ARC00049"]
        pbm.A[ig,iv] += Float64(1.43027E-02)
        iv,ix_,_ = s2mpj_ii("PI000039",ix_)
        arrset(pb.xnames,iv,"PI000039")
        ig = ig_["ARC00050"]
        pbm.A[ig,iv] += Float64(-1.09654E-01)
        iv,ix_,_ = s2mpj_ii("PI000040",ix_)
        arrset(pb.xnames,iv,"PI000040")
        ig = ig_["ARC00050"]
        pbm.A[ig,iv] += Float64(1.09654E-01)
        iv,ix_,_ = s2mpj_ii("PI000040",ix_)
        arrset(pb.xnames,iv,"PI000040")
        ig = ig_["ARC00051"]
        pbm.A[ig,iv] += Float64(5.22555E-01)
        iv,ix_,_ = s2mpj_ii("PI000040",ix_)
        arrset(pb.xnames,iv,"PI000040")
        ig = ig_["ARC00052"]
        pbm.A[ig,iv] += Float64(-3.21572E-01)
        iv,ix_,_ = s2mpj_ii("PI000041",ix_)
        arrset(pb.xnames,iv,"PI000041")
        ig = ig_["ARC00051"]
        pbm.A[ig,iv] += Float64(-5.22555E-01)
        iv,ix_,_ = s2mpj_ii("PI000042",ix_)
        arrset(pb.xnames,iv,"PI000042")
        ig = ig_["ARC00053"]
        pbm.A[ig,iv] += Float64(1.64096E-01)
        iv,ix_,_ = s2mpj_ii("PI000042",ix_)
        arrset(pb.xnames,iv,"PI000042")
        ig = ig_["ARC00054"]
        pbm.A[ig,iv] += Float64(-1.36747E-01)
        iv,ix_,_ = s2mpj_ii("PI000043",ix_)
        arrset(pb.xnames,iv,"PI000043")
        ig = ig_["ARC00055"]
        pbm.A[ig,iv] += Float64(3.99253E-02)
        iv,ix_,_ = s2mpj_ii("PI000043",ix_)
        arrset(pb.xnames,iv,"PI000043")
        ig = ig_["ARC00056"]
        pbm.A[ig,iv] += Float64(-2.85181E-02)
        iv,ix_,_ = s2mpj_ii("PI000044",ix_)
        arrset(pb.xnames,iv,"PI000044")
        ig = ig_["ARC00057"]
        pbm.A[ig,iv] += Float64(-1.82329E-01)
        iv,ix_,_ = s2mpj_ii("PI000044",ix_)
        arrset(pb.xnames,iv,"PI000044")
        ig = ig_["ARC00059"]
        pbm.A[ig,iv] += Float64(-3.32711E-02)
        iv,ix_,_ = s2mpj_ii("PI000045",ix_)
        arrset(pb.xnames,iv,"PI000045")
        ig = ig_["ARC00057"]
        pbm.A[ig,iv] += Float64(1.82329E-01)
        iv,ix_,_ = s2mpj_ii("PI000045",ix_)
        arrset(pb.xnames,iv,"PI000045")
        ig = ig_["ARC00058"]
        pbm.A[ig,iv] += Float64(-1.82329E-01)
        iv,ix_,_ = s2mpj_ii("PI000046",ix_)
        arrset(pb.xnames,iv,"PI000046")
        ig = ig_["ARC00058"]
        pbm.A[ig,iv] += Float64(1.82329E-01)
        iv,ix_,_ = s2mpj_ii("PI000047",ix_)
        arrset(pb.xnames,iv,"PI000047")
        ig = ig_["ARC00059"]
        pbm.A[ig,iv] += Float64(3.32711E-02)
        iv,ix_,_ = s2mpj_ii("PI000047",ix_)
        arrset(pb.xnames,iv,"PI000047")
        ig = ig_["ARC00060"]
        pbm.A[ig,iv] += Float64(-9.98133E-03)
        iv,ix_,_ = s2mpj_ii("PI000048",ix_)
        arrset(pb.xnames,iv,"PI000048")
        ig = ig_["ARC00060"]
        pbm.A[ig,iv] += Float64(9.98133E-03)
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
        pbm.gconst[ig_["REGIO001"]] = Float64(-8.90000E+00)
        pbm.gconst[ig_["REGIO002"]] = Float64(-3.30000E+00)
        pbm.gconst[ig_["REGIO003"]] = Float64(-5.40000E+00)
        pbm.gconst[ig_["REGIO004"]] = Float64(-1.45000E+01)
        pbm.gconst[ig_["REGIO005"]] = Float64(-2.20000E+00)
        pbm.gconst[ig_["REGIO006"]] = Float64(-5.80000E+00)
        pbm.gconst[ig_["REGIO007"]] = Float64(-3.50000E+00)
        pbm.gconst[ig_["REGIO008"]] = Float64(-1.70000E+00)
        pbm.gconst[ig_["REGIO009"]] = Float64(-1.50000E+00)
        pbm.gconst[ig_["REGIO010"]] = Float64(-8.00000E-01)
        pbm.gconst[ig_["REGIO011"]] = Float64(-1.30000E+00)
        pbm.gconst[ig_["REGIO012"]] = Float64(-2.30000E+00)
        pbm.gconst[ig_["REGIO013"]] = Float64(-2.20000E+00)
        pbm.gconst[ig_["REGIO014"]] = Float64(-1.10000E+00)
        pbm.gconst[ig_["REGIO015"]] = Float64(-7.00000E-01)
        pbm.gconst[ig_["REGIO016"]] = Float64(-1.50000E+00)
        pbm.gconst[ig_["REGIO017"]] = Float64(-4.30000E+00)
        pbm.gconst[ig_["REGIO018"]] = Float64(-1.50000E+00)
        pbm.gconst[ig_["REGIO019"]] = Float64(-4.30000E+00)
        pbm.gconst[ig_["REGIO020"]] = Float64(-1.10000E+00)
        pbm.gconst[ig_["REGIO021"]] = Float64(-5.00000E+00)
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower[ix_["FLOW0001"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0001"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0002"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0002"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0003"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0003"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0004"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0004"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0005"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0005"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0006"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0006"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0007"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0007"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0008"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0008"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0009"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0009"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0010"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0010"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0011"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0011"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0012"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0012"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0013"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0013"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0014"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0014"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0015"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0015"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0016"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0016"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0017"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0017"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0018"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0018"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0019"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0019"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0020"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0020"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0021"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0021"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0022"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0022"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0023"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0023"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0024"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0024"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0025"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0025"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0026"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0026"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0027"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0027"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0028"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0028"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0029"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0029"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0030"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0030"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0031"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0031"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0032"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0032"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0033"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0033"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0034"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0034"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0035"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0035"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0036"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0036"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0037"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0037"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0038"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0038"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0039"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0039"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0040"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0040"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0041"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0041"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0042"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0042"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0043"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0043"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0044"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0044"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0045"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0045"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0046"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0046"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0047"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0047"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0048"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0048"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0049"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0049"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0050"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0050"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0051"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0051"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0052"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0052"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0053"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0053"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0054"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0054"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0055"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0055"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0056"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0056"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0057"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0057"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0058"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0058"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0059"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0059"]] = 3.61167E+02
        pb.xlower[ix_["FLOW0060"]] = -3.61167E+02
        pb.xupper[ix_["FLOW0060"]] = 3.61167E+02
        pb.xupper[ix_["SUPP0001"]] = 1.00000E+30
        pb.xupper[ix_["SUPP0002"]] = 1.00000E+30
        pb.xupper[ix_["SUPP0003"]] = 1.00000E+30
        pb.xupper[ix_["SUPP0004"]] = 4.80000E+00
        pb.xupper[ix_["SUPP0005"]] = 2.45000E+01
        pb.xupper[ix_["SUPP0006"]] = 1.32000E+01
        pb.xupper[ix_["SUPP0007"]] = 1.00000E+30
        pb.xlower[ix_["PROD0001"]] = 2.33333E+00
        pb.xupper[ix_["PROD0001"]] = 6.41667E+00
        pb.xlower[ix_["PROD0002"]] = 4.92222E+00
        pb.xupper[ix_["PROD0002"]] = 1.35361E+01
        pb.xlower[ix_["PROD0003"]] = 1.31333E+01
        pb.xupper[ix_["PROD0003"]] = 3.61167E+01
        pb.xlower[ix_["PROD0004"]] = 9.55556E+00
        pb.xupper[ix_["PROD0004"]] = 2.62778E+01
        pb.xlower[ix_["PROD0005"]] = 6.06667E+00
        pb.xupper[ix_["PROD0005"]] = 1.66833E+01
        pb.xlower[ix_["PI000001"]] = 0.00000E+00
        pb.xupper[ix_["PI000001"]] = 6.40000E+03
        pb.xlower[ix_["PI000002"]] = 1.60000E+03
        pb.xupper[ix_["PI000002"]] = 6.40000E+03
        pb.xlower[ix_["PI000003"]] = 0.00000E+00
        pb.xupper[ix_["PI000003"]] = 6.40000E+03
        pb.xlower[ix_["PI000004"]] = 1.60000E+03
        pb.xupper[ix_["PI000004"]] = 6.40000E+03
        pb.xlower[ix_["PI000005"]] = 1.60000E+03
        pb.xupper[ix_["PI000005"]] = 6.40000E+03
        pb.xlower[ix_["PI000006"]] = 1.60000E+03
        pb.xupper[ix_["PI000006"]] = 6.40000E+03
        pb.xlower[ix_["PI000007"]] = 0.00000E+00
        pb.xupper[ix_["PI000007"]] = 6.40000E+03
        pb.xlower[ix_["PI000008"]] = 1.60000E+03
        pb.xupper[ix_["PI000008"]] = 6.40000E+03
        pb.xlower[ix_["PI000009"]] = 1.60000E+03
        pb.xupper[ix_["PI000009"]] = 6.40000E+03
        pb.xlower[ix_["PI000010"]] = 0.00000E+00
        pb.xupper[ix_["PI000010"]] = 6.40000E+03
        pb.xlower[ix_["PI000011"]] = 1.60000E+03
        pb.xupper[ix_["PI000011"]] = 6.40000E+03
        pb.xlower[ix_["PI000012"]] = 0.00000E+00
        pb.xupper[ix_["PI000012"]] = 6.40000E+03
        pb.xlower[ix_["PI000013"]] = 1.60000E+03
        pb.xupper[ix_["PI000013"]] = 6.40000E+03
        pb.xlower[ix_["PI000014"]] = 1.60000E+03
        pb.xupper[ix_["PI000014"]] = 6.40000E+03
        pb.xlower[ix_["PI000015"]] = 0.00000E+00
        pb.xupper[ix_["PI000015"]] = 6.40000E+03
        pb.xlower[ix_["PI000016"]] = 1.60000E+03
        pb.xupper[ix_["PI000016"]] = 6.40000E+03
        pb.xlower[ix_["PI000017"]] = 1.60000E+03
        pb.xupper[ix_["PI000017"]] = 6.40000E+03
        pb.xlower[ix_["PI000018"]] = 1.60000E+03
        pb.xupper[ix_["PI000018"]] = 6.40000E+03
        pb.xlower[ix_["PI000019"]] = 0.00000E+00
        pb.xupper[ix_["PI000019"]] = 6.40000E+03
        pb.xlower[ix_["PI000020"]] = 1.60000E+03
        pb.xupper[ix_["PI000020"]] = 6.40000E+03
        pb.xlower[ix_["PI000021"]] = 1.60000E+03
        pb.xupper[ix_["PI000021"]] = 4.48900E+03
        pb.xlower[ix_["PI000022"]] = 1.60000E+03
        pb.xupper[ix_["PI000022"]] = 6.40000E+03
        pb.xlower[ix_["PI000023"]] = 1.60000E+03
        pb.xupper[ix_["PI000023"]] = 6.40000E+03
        pb.xlower[ix_["PI000024"]] = 0.00000E+00
        pb.xupper[ix_["PI000024"]] = 6.40000E+03
        pb.xlower[ix_["PI000025"]] = 1.60000E+03
        pb.xupper[ix_["PI000025"]] = 6.40000E+03
        pb.xlower[ix_["PI000026"]] = 1.60000E+03
        pb.xupper[ix_["PI000026"]] = 6.40000E+03
        pb.xlower[ix_["PI000027"]] = 0.00000E+00
        pb.xupper[ix_["PI000027"]] = 6.40000E+03
        pb.xlower[ix_["PI000028"]] = 1.60000E+03
        pb.xupper[ix_["PI000028"]] = 6.40000E+03
        pb.xlower[ix_["PI000029"]] = 0.00000E+00
        pb.xupper[ix_["PI000029"]] = 6.40000E+03
        pb.xlower[ix_["PI000030"]] = 1.60000E+03
        pb.xupper[ix_["PI000030"]] = 6.40000E+03
        pb.xlower[ix_["PI000031"]] = 1.60000E+03
        pb.xupper[ix_["PI000031"]] = 6.40000E+03
        pb.xlower[ix_["PI000032"]] = 0.00000E+00
        pb.xupper[ix_["PI000032"]] = 6.40000E+03
        pb.xlower[ix_["PI000033"]] = 1.60000E+03
        pb.xupper[ix_["PI000033"]] = 6.40000E+03
        pb.xlower[ix_["PI000034"]] = 0.00000E+00
        pb.xupper[ix_["PI000034"]] = 6.40000E+03
        pb.xlower[ix_["PI000035"]] = 1.60000E+03
        pb.xupper[ix_["PI000035"]] = 6.40000E+03
        pb.xlower[ix_["PI000036"]] = 1.60000E+03
        pb.xupper[ix_["PI000036"]] = 6.40000E+03
        pb.xlower[ix_["PI000037"]] = 0.00000E+00
        pb.xupper[ix_["PI000037"]] = 4.48900E+03
        pb.xlower[ix_["PI000038"]] = 1.60000E+03
        pb.xupper[ix_["PI000038"]] = 6.40000E+03
        pb.xlower[ix_["PI000039"]] = 1.60000E+03
        pb.xupper[ix_["PI000039"]] = 6.40000E+03
        pb.xlower[ix_["PI000040"]] = 1.60000E+03
        pb.xupper[ix_["PI000040"]] = 6.40000E+03
        pb.xlower[ix_["PI000041"]] = 1.60000E+03
        pb.xupper[ix_["PI000041"]] = 6.40000E+03
        pb.xlower[ix_["PI000042"]] = 1.60000E+03
        pb.xupper[ix_["PI000042"]] = 6.40000E+03
        pb.xlower[ix_["PI000043"]] = 1.60000E+03
        pb.xupper[ix_["PI000043"]] = 6.40000E+03
        pb.xlower[ix_["PI000044"]] = 1.60000E+03
        pb.xupper[ix_["PI000044"]] = 4.48900E+03
        pb.xlower[ix_["PI000045"]] = 1.60000E+03
        pb.xupper[ix_["PI000045"]] = 6.40000E+03
        pb.xlower[ix_["PI000046"]] = 1.60000E+03
        pb.xupper[ix_["PI000046"]] = 6.40000E+03
        pb.xlower[ix_["PI000047"]] = 1.60000E+03
        pb.xupper[ix_["PI000047"]] = 6.40000E+03
        pb.xlower[ix_["PI000048"]] = 1.60000E+03
        pb.xupper[ix_["PI000048"]] = 6.40000E+03
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        pb.x0[ix_["FLOW0001"]] = Float64(2.28416E+00)
        pb.x0[ix_["FLOW0002"]] = Float64(1.28924E+00)
        pb.x0[ix_["FLOW0003"]] = Float64(-4.11076E+00)
        pb.x0[ix_["FLOW0004"]] = Float64(-4.11076E+00)
        pb.x0[ix_["FLOW0005"]] = Float64(8.34083E+00)
        pb.x0[ix_["FLOW0006"]] = Float64(6.67356E+00)
        pb.x0[ix_["FLOW0007"]] = Float64(5.15250E+00)
        pb.x0[ix_["FLOW0008"]] = Float64(-1.52106E+00)
        pb.x0[ix_["FLOW0009"]] = Float64(1.49619E+00)
        pb.x0[ix_["FLOW0010"]] = Float64(8.63464E-01)
        pb.x0[ix_["FLOW0011"]] = Float64(2.38452E+00)
        pb.x0[ix_["FLOW0012"]] = Float64(1.84522E-01)
        pb.x0[ix_["FLOW0013"]] = Float64(-2.80167E-01)
        pb.x0[ix_["FLOW0014"]] = Float64(-2.80167E-01)
        pb.x0[ix_["FLOW0015"]] = Float64(-2.80167E-01)
        pb.x0[ix_["FLOW0016"]] = Float64(-9.56457E-02)
        pb.x0[ix_["FLOW0017"]] = Float64(4.83333E-01)
        pb.x0[ix_["FLOW0018"]] = Float64(0.00000E+00)
        pb.x0[ix_["FLOW0019"]] = Float64(5.80000E+00)
        pb.x0[ix_["FLOW0020"]] = Float64(3.62164E+00)
        pb.x0[ix_["FLOW0021"]] = Float64(3.15920E+00)
        pb.x0[ix_["FLOW0022"]] = Float64(0.00000E+00)
        pb.x0[ix_["FLOW0023"]] = Float64(-1.14080E+00)
        pb.x0[ix_["FLOW0024"]] = Float64(2.96863E+00)
        pb.x0[ix_["FLOW0025"]] = Float64(6.56169E+00)
        pb.x0[ix_["FLOW0026"]] = Float64(3.16455E-01)
        pb.x0[ix_["FLOW0027"]] = Float64(-2.74176E-01)
        pb.x0[ix_["FLOW0028"]] = Float64(-6.71257E-01)
        pb.x0[ix_["FLOW0029"]] = Float64(-6.71257E-01)
        pb.x0[ix_["FLOW0030"]] = Float64(-1.52215E+00)
        pb.x0[ix_["FLOW0031"]] = Float64(-1.14514E+00)
        pb.x0[ix_["FLOW0032"]] = Float64(-2.31222E-01)
        pb.x0[ix_["FLOW0033"]] = Float64(-2.31221E-01)
        pb.x0[ix_["FLOW0034"]] = Float64(1.82783E+00)
        pb.x0[ix_["FLOW0035"]] = Float64(4.80828E-01)
        pb.x0[ix_["FLOW0036"]] = Float64(-3.01828E-01)
        pb.x0[ix_["FLOW0037"]] = Float64(1.06084E+01)
        pb.x0[ix_["FLOW0038"]] = Float64(9.74452E+00)
        pb.x0[ix_["FLOW0039"]] = Float64(1.77661E+00)
        pb.x0[ix_["FLOW0040"]] = Float64(-2.09158E-01)
        pb.x0[ix_["FLOW0041"]] = Float64(-5.12076E-01)
        pb.x0[ix_["FLOW0042"]] = Float64(-5.12076E-01)
        pb.x0[ix_["FLOW0043"]] = Float64(-1.16119E+00)
        pb.x0[ix_["FLOW0044"]] = Float64(-1.14514E+00)
        pb.x0[ix_["FLOW0045"]] = Float64(-2.31221E-01)
        pb.x0[ix_["FLOW0046"]] = Float64(-2.31221E-01)
        pb.x0[ix_["FLOW0047"]] = Float64(2.20035E+00)
        pb.x0[ix_["FLOW0048"]] = Float64(-1.12775E+00)
        pb.x0[ix_["FLOW0049"]] = Float64(1.91798E+00)
        pb.x0[ix_["FLOW0050"]] = Float64(0.00000E+00)
        pb.x0[ix_["FLOW0051"]] = Float64(1.15931E+01)
        pb.x0[ix_["FLOW0052"]] = Float64(7.71110E+00)
        pb.x0[ix_["FLOW0053"]] = Float64(1.14942E+00)
        pb.x0[ix_["FLOW0054"]] = Float64(1.14942E+00)
        pb.x0[ix_["FLOW0055"]] = Float64(3.17981E+00)
        pb.x0[ix_["FLOW0056"]] = Float64(2.37981E+00)
        pb.x0[ix_["FLOW0057"]] = Float64(2.43611E+00)
        pb.x0[ix_["FLOW0058"]] = Float64(0.00000E+00)
        pb.x0[ix_["FLOW0059"]] = Float64(2.20000E+00)
        pb.x0[ix_["FLOW0060"]] = Float64(0.00000E+00)
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
        ename = "F00025SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0025"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00026SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0026"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00027SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0027"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00028SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0028"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00029SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0029"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00030SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0030"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00031SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0031"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00032SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0032"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00033SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0033"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00034SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0034"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00035SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0035"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00036SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0036"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00037SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0037"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00038SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0038"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00039SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0039"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00040SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0040"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00041SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0041"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00042SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0042"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00043SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0043"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00044SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0044"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00045SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0045"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00046SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0046"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00047SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0047"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00048SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0048"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00049SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0049"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00050SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0050"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00051SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0051"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00052SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0052"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00053SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0053"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00054SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0054"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00055SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0055"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00056SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0056"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00057SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0057"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00058SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0058"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00059SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0059"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "F00060SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "FLOW0060"
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
        ig = ig_["ARC00025"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00025SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00026"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00026SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00027"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00027SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00028"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00028SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00029"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00029SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00030"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00030SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00031"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00031SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00032"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00032SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00033"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00033SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00034"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00034SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00035"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00035SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00036"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00036SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00037"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00037SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00038"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00038SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00039"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00039SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00040"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00040SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00041"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00041SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00042"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00042SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00043"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00043SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00044"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00044SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00045"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00045SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00046"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00046SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00047"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00047SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00048"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00048SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00049"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00049SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00050"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00050SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00051"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00051SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00052"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00052SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00053"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00053SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00054"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00054SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00055"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00055SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00056"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00056SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00057"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00057SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00058"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00058SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00059"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00059SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["ARC00060"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F00060SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.cupper[1:pb.nle] = zeros(Float64,pb.nle)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-CLQI2-RN-157-134"
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

