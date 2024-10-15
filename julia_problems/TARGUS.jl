function TARGUS(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : TARGUS
#    *********
# 
#    A two-norm fitted formulation of the problem of finding the 
#    smallest perturbation of data that fits a linear model
#    arising in large-scale tabular data protection.
# 
#    Source:
#    J. Castro, 
#    Minimum-distance controlled perturbation methods for 
#    large-scale tabular data protection, 
#    European Journal of Operational Research 171 (2006) pp 39-52.
# 
#    SIF input: Jordi Castro, 2006 as L2_targus.mps
#    see http://www-eio.upc.es/~jcastro/data.html
# 
#    classification = "C-QLR2-RN-162-63"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "TARGUS"

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
        ig,ig_,_ = s2mpj_ii("R0001",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0001")
        ig,ig_,_ = s2mpj_ii("R0002",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0002")
        ig,ig_,_ = s2mpj_ii("R0003",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0003")
        ig,ig_,_ = s2mpj_ii("R0004",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0004")
        ig,ig_,_ = s2mpj_ii("R0005",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0005")
        ig,ig_,_ = s2mpj_ii("R0006",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0006")
        ig,ig_,_ = s2mpj_ii("R0007",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0007")
        ig,ig_,_ = s2mpj_ii("R0008",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0008")
        ig,ig_,_ = s2mpj_ii("R0009",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0009")
        ig,ig_,_ = s2mpj_ii("R0010",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0010")
        ig,ig_,_ = s2mpj_ii("R0011",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0011")
        ig,ig_,_ = s2mpj_ii("R0012",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0012")
        ig,ig_,_ = s2mpj_ii("R0013",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0013")
        ig,ig_,_ = s2mpj_ii("R0014",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0014")
        ig,ig_,_ = s2mpj_ii("R0015",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0015")
        ig,ig_,_ = s2mpj_ii("R0016",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0016")
        ig,ig_,_ = s2mpj_ii("R0017",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0017")
        ig,ig_,_ = s2mpj_ii("R0018",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0018")
        ig,ig_,_ = s2mpj_ii("R0019",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0019")
        ig,ig_,_ = s2mpj_ii("R0020",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0020")
        ig,ig_,_ = s2mpj_ii("R0021",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0021")
        ig,ig_,_ = s2mpj_ii("R0022",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0022")
        ig,ig_,_ = s2mpj_ii("R0023",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0023")
        ig,ig_,_ = s2mpj_ii("R0024",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0024")
        ig,ig_,_ = s2mpj_ii("R0025",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0025")
        ig,ig_,_ = s2mpj_ii("R0026",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0026")
        ig,ig_,_ = s2mpj_ii("R0027",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0027")
        ig,ig_,_ = s2mpj_ii("R0028",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0028")
        ig,ig_,_ = s2mpj_ii("R0029",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0029")
        ig,ig_,_ = s2mpj_ii("R0030",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0030")
        ig,ig_,_ = s2mpj_ii("R0031",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0031")
        ig,ig_,_ = s2mpj_ii("R0032",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0032")
        ig,ig_,_ = s2mpj_ii("R0033",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0033")
        ig,ig_,_ = s2mpj_ii("R0034",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0034")
        ig,ig_,_ = s2mpj_ii("R0035",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0035")
        ig,ig_,_ = s2mpj_ii("R0036",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0036")
        ig,ig_,_ = s2mpj_ii("R0037",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0037")
        ig,ig_,_ = s2mpj_ii("R0038",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0038")
        ig,ig_,_ = s2mpj_ii("R0039",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0039")
        ig,ig_,_ = s2mpj_ii("R0040",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0040")
        ig,ig_,_ = s2mpj_ii("R0041",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0041")
        ig,ig_,_ = s2mpj_ii("R0042",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0042")
        ig,ig_,_ = s2mpj_ii("R0043",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0043")
        ig,ig_,_ = s2mpj_ii("R0044",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0044")
        ig,ig_,_ = s2mpj_ii("R0045",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0045")
        ig,ig_,_ = s2mpj_ii("R0046",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0046")
        ig,ig_,_ = s2mpj_ii("R0047",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0047")
        ig,ig_,_ = s2mpj_ii("R0048",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0048")
        ig,ig_,_ = s2mpj_ii("R0049",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0049")
        ig,ig_,_ = s2mpj_ii("R0050",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0050")
        ig,ig_,_ = s2mpj_ii("R0051",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0051")
        ig,ig_,_ = s2mpj_ii("R0052",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0052")
        ig,ig_,_ = s2mpj_ii("R0053",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0053")
        ig,ig_,_ = s2mpj_ii("R0054",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0054")
        ig,ig_,_ = s2mpj_ii("R0055",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0055")
        ig,ig_,_ = s2mpj_ii("R0056",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0056")
        ig,ig_,_ = s2mpj_ii("R0057",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0057")
        ig,ig_,_ = s2mpj_ii("R0058",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0058")
        ig,ig_,_ = s2mpj_ii("R0059",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0059")
        ig,ig_,_ = s2mpj_ii("R0060",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0060")
        ig,ig_,_ = s2mpj_ii("R0061",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0061")
        ig,ig_,_ = s2mpj_ii("R0062",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0062")
        ig,ig_,_ = s2mpj_ii("R0063",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0063")
        ig,ig_,_ = s2mpj_ii("R0064",ig_)
        arrset(gtype,ig,"<>")
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        ngrp   = length(ig_)
        iv,ix_,_ = s2mpj_ii("C0001",ix_)
        arrset(pb.xnames,iv,"C0001")
        ig = ig_["R0001"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0001",ix_)
        arrset(pb.xnames,iv,"C0001")
        ig = ig_["R0019"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0002",ix_)
        arrset(pb.xnames,iv,"C0002")
        ig = ig_["R0002"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0002",ix_)
        arrset(pb.xnames,iv,"C0002")
        ig = ig_["R0019"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0002",ix_)
        arrset(pb.xnames,iv,"C0002")
        ig = ig_["R0020"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0003",ix_)
        arrset(pb.xnames,iv,"C0003")
        ig = ig_["R0003"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0003",ix_)
        arrset(pb.xnames,iv,"C0003")
        ig = ig_["R0020"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0004",ix_)
        arrset(pb.xnames,iv,"C0004")
        ig = ig_["R0004"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0004",ix_)
        arrset(pb.xnames,iv,"C0004")
        ig = ig_["R0020"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0005",ix_)
        arrset(pb.xnames,iv,"C0005")
        ig = ig_["R0005"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0005",ix_)
        arrset(pb.xnames,iv,"C0005")
        ig = ig_["R0020"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0006",ix_)
        arrset(pb.xnames,iv,"C0006")
        ig = ig_["R0006"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0006",ix_)
        arrset(pb.xnames,iv,"C0006")
        ig = ig_["R0019"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0006",ix_)
        arrset(pb.xnames,iv,"C0006")
        ig = ig_["R0021"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0007",ix_)
        arrset(pb.xnames,iv,"C0007")
        ig = ig_["R0007"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0007",ix_)
        arrset(pb.xnames,iv,"C0007")
        ig = ig_["R0021"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0008",ix_)
        arrset(pb.xnames,iv,"C0008")
        ig = ig_["R0008"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0008",ix_)
        arrset(pb.xnames,iv,"C0008")
        ig = ig_["R0021"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0009",ix_)
        arrset(pb.xnames,iv,"C0009")
        ig = ig_["R0009"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0009",ix_)
        arrset(pb.xnames,iv,"C0009")
        ig = ig_["R0021"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0010",ix_)
        arrset(pb.xnames,iv,"C0010")
        ig = ig_["R0010"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0010",ix_)
        arrset(pb.xnames,iv,"C0010")
        ig = ig_["R0021"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0011",ix_)
        arrset(pb.xnames,iv,"C0011")
        ig = ig_["R0011"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0011",ix_)
        arrset(pb.xnames,iv,"C0011")
        ig = ig_["R0019"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0011",ix_)
        arrset(pb.xnames,iv,"C0011")
        ig = ig_["R0022"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0012",ix_)
        arrset(pb.xnames,iv,"C0012")
        ig = ig_["R0012"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0012",ix_)
        arrset(pb.xnames,iv,"C0012")
        ig = ig_["R0022"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0013",ix_)
        arrset(pb.xnames,iv,"C0013")
        ig = ig_["R0013"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0013",ix_)
        arrset(pb.xnames,iv,"C0013")
        ig = ig_["R0022"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0014",ix_)
        arrset(pb.xnames,iv,"C0014")
        ig = ig_["R0014"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0014",ix_)
        arrset(pb.xnames,iv,"C0014")
        ig = ig_["R0022"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0015",ix_)
        arrset(pb.xnames,iv,"C0015")
        ig = ig_["R0015"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0015",ix_)
        arrset(pb.xnames,iv,"C0015")
        ig = ig_["R0019"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0015",ix_)
        arrset(pb.xnames,iv,"C0015")
        ig = ig_["R0023"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0016",ix_)
        arrset(pb.xnames,iv,"C0016")
        ig = ig_["R0016"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0016",ix_)
        arrset(pb.xnames,iv,"C0016")
        ig = ig_["R0023"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0017",ix_)
        arrset(pb.xnames,iv,"C0017")
        ig = ig_["R0017"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0017",ix_)
        arrset(pb.xnames,iv,"C0017")
        ig = ig_["R0023"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0018",ix_)
        arrset(pb.xnames,iv,"C0018")
        ig = ig_["R0018"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0018",ix_)
        arrset(pb.xnames,iv,"C0018")
        ig = ig_["R0019"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0019",ix_)
        arrset(pb.xnames,iv,"C0019")
        ig = ig_["R0001"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0019",ix_)
        arrset(pb.xnames,iv,"C0019")
        ig = ig_["R0024"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0020",ix_)
        arrset(pb.xnames,iv,"C0020")
        ig = ig_["R0002"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0020",ix_)
        arrset(pb.xnames,iv,"C0020")
        ig = ig_["R0024"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0020",ix_)
        arrset(pb.xnames,iv,"C0020")
        ig = ig_["R0025"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0021",ix_)
        arrset(pb.xnames,iv,"C0021")
        ig = ig_["R0003"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0021",ix_)
        arrset(pb.xnames,iv,"C0021")
        ig = ig_["R0025"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0022",ix_)
        arrset(pb.xnames,iv,"C0022")
        ig = ig_["R0004"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0022",ix_)
        arrset(pb.xnames,iv,"C0022")
        ig = ig_["R0025"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0023",ix_)
        arrset(pb.xnames,iv,"C0023")
        ig = ig_["R0005"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0023",ix_)
        arrset(pb.xnames,iv,"C0023")
        ig = ig_["R0025"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0024",ix_)
        arrset(pb.xnames,iv,"C0024")
        ig = ig_["R0006"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0024",ix_)
        arrset(pb.xnames,iv,"C0024")
        ig = ig_["R0024"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0024",ix_)
        arrset(pb.xnames,iv,"C0024")
        ig = ig_["R0026"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0025",ix_)
        arrset(pb.xnames,iv,"C0025")
        ig = ig_["R0007"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0025",ix_)
        arrset(pb.xnames,iv,"C0025")
        ig = ig_["R0026"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0026",ix_)
        arrset(pb.xnames,iv,"C0026")
        ig = ig_["R0008"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0026",ix_)
        arrset(pb.xnames,iv,"C0026")
        ig = ig_["R0026"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0027",ix_)
        arrset(pb.xnames,iv,"C0027")
        ig = ig_["R0009"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0027",ix_)
        arrset(pb.xnames,iv,"C0027")
        ig = ig_["R0026"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0028",ix_)
        arrset(pb.xnames,iv,"C0028")
        ig = ig_["R0010"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0028",ix_)
        arrset(pb.xnames,iv,"C0028")
        ig = ig_["R0026"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0029",ix_)
        arrset(pb.xnames,iv,"C0029")
        ig = ig_["R0011"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0029",ix_)
        arrset(pb.xnames,iv,"C0029")
        ig = ig_["R0024"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0029",ix_)
        arrset(pb.xnames,iv,"C0029")
        ig = ig_["R0027"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0030",ix_)
        arrset(pb.xnames,iv,"C0030")
        ig = ig_["R0012"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0030",ix_)
        arrset(pb.xnames,iv,"C0030")
        ig = ig_["R0027"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0031",ix_)
        arrset(pb.xnames,iv,"C0031")
        ig = ig_["R0013"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0031",ix_)
        arrset(pb.xnames,iv,"C0031")
        ig = ig_["R0027"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0032",ix_)
        arrset(pb.xnames,iv,"C0032")
        ig = ig_["R0014"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0032",ix_)
        arrset(pb.xnames,iv,"C0032")
        ig = ig_["R0027"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0033",ix_)
        arrset(pb.xnames,iv,"C0033")
        ig = ig_["R0015"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0033",ix_)
        arrset(pb.xnames,iv,"C0033")
        ig = ig_["R0024"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0033",ix_)
        arrset(pb.xnames,iv,"C0033")
        ig = ig_["R0028"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0034",ix_)
        arrset(pb.xnames,iv,"C0034")
        ig = ig_["R0016"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0034",ix_)
        arrset(pb.xnames,iv,"C0034")
        ig = ig_["R0028"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0035",ix_)
        arrset(pb.xnames,iv,"C0035")
        ig = ig_["R0017"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0035",ix_)
        arrset(pb.xnames,iv,"C0035")
        ig = ig_["R0028"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0036",ix_)
        arrset(pb.xnames,iv,"C0036")
        ig = ig_["R0018"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0036",ix_)
        arrset(pb.xnames,iv,"C0036")
        ig = ig_["R0024"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0037",ix_)
        arrset(pb.xnames,iv,"C0037")
        ig = ig_["R0001"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0037",ix_)
        arrset(pb.xnames,iv,"C0037")
        ig = ig_["R0029"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0038",ix_)
        arrset(pb.xnames,iv,"C0038")
        ig = ig_["R0002"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0038",ix_)
        arrset(pb.xnames,iv,"C0038")
        ig = ig_["R0029"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0038",ix_)
        arrset(pb.xnames,iv,"C0038")
        ig = ig_["R0030"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0039",ix_)
        arrset(pb.xnames,iv,"C0039")
        ig = ig_["R0003"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0039",ix_)
        arrset(pb.xnames,iv,"C0039")
        ig = ig_["R0030"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0040",ix_)
        arrset(pb.xnames,iv,"C0040")
        ig = ig_["R0004"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0040",ix_)
        arrset(pb.xnames,iv,"C0040")
        ig = ig_["R0030"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0041",ix_)
        arrset(pb.xnames,iv,"C0041")
        ig = ig_["R0005"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0041",ix_)
        arrset(pb.xnames,iv,"C0041")
        ig = ig_["R0030"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0042",ix_)
        arrset(pb.xnames,iv,"C0042")
        ig = ig_["R0006"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0042",ix_)
        arrset(pb.xnames,iv,"C0042")
        ig = ig_["R0029"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0042",ix_)
        arrset(pb.xnames,iv,"C0042")
        ig = ig_["R0031"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0043",ix_)
        arrset(pb.xnames,iv,"C0043")
        ig = ig_["R0007"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0043",ix_)
        arrset(pb.xnames,iv,"C0043")
        ig = ig_["R0031"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0044",ix_)
        arrset(pb.xnames,iv,"C0044")
        ig = ig_["R0008"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0044",ix_)
        arrset(pb.xnames,iv,"C0044")
        ig = ig_["R0031"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0045",ix_)
        arrset(pb.xnames,iv,"C0045")
        ig = ig_["R0009"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0045",ix_)
        arrset(pb.xnames,iv,"C0045")
        ig = ig_["R0031"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0046",ix_)
        arrset(pb.xnames,iv,"C0046")
        ig = ig_["R0010"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0046",ix_)
        arrset(pb.xnames,iv,"C0046")
        ig = ig_["R0031"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0047",ix_)
        arrset(pb.xnames,iv,"C0047")
        ig = ig_["R0011"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0047",ix_)
        arrset(pb.xnames,iv,"C0047")
        ig = ig_["R0029"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0047",ix_)
        arrset(pb.xnames,iv,"C0047")
        ig = ig_["R0032"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0048",ix_)
        arrset(pb.xnames,iv,"C0048")
        ig = ig_["R0012"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0048",ix_)
        arrset(pb.xnames,iv,"C0048")
        ig = ig_["R0032"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0049",ix_)
        arrset(pb.xnames,iv,"C0049")
        ig = ig_["R0013"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0049",ix_)
        arrset(pb.xnames,iv,"C0049")
        ig = ig_["R0032"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0050",ix_)
        arrset(pb.xnames,iv,"C0050")
        ig = ig_["R0014"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0050",ix_)
        arrset(pb.xnames,iv,"C0050")
        ig = ig_["R0032"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0051",ix_)
        arrset(pb.xnames,iv,"C0051")
        ig = ig_["R0015"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0051",ix_)
        arrset(pb.xnames,iv,"C0051")
        ig = ig_["R0029"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0051",ix_)
        arrset(pb.xnames,iv,"C0051")
        ig = ig_["R0033"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0052",ix_)
        arrset(pb.xnames,iv,"C0052")
        ig = ig_["R0016"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0052",ix_)
        arrset(pb.xnames,iv,"C0052")
        ig = ig_["R0033"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0053",ix_)
        arrset(pb.xnames,iv,"C0053")
        ig = ig_["R0017"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0053",ix_)
        arrset(pb.xnames,iv,"C0053")
        ig = ig_["R0033"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0054",ix_)
        arrset(pb.xnames,iv,"C0054")
        ig = ig_["R0018"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0054",ix_)
        arrset(pb.xnames,iv,"C0054")
        ig = ig_["R0029"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0055",ix_)
        arrset(pb.xnames,iv,"C0055")
        ig = ig_["R0001"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0055",ix_)
        arrset(pb.xnames,iv,"C0055")
        ig = ig_["R0034"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0056",ix_)
        arrset(pb.xnames,iv,"C0056")
        ig = ig_["R0002"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0056",ix_)
        arrset(pb.xnames,iv,"C0056")
        ig = ig_["R0034"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0056",ix_)
        arrset(pb.xnames,iv,"C0056")
        ig = ig_["R0035"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0057",ix_)
        arrset(pb.xnames,iv,"C0057")
        ig = ig_["R0003"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0057",ix_)
        arrset(pb.xnames,iv,"C0057")
        ig = ig_["R0035"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0058",ix_)
        arrset(pb.xnames,iv,"C0058")
        ig = ig_["R0004"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0058",ix_)
        arrset(pb.xnames,iv,"C0058")
        ig = ig_["R0035"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0059",ix_)
        arrset(pb.xnames,iv,"C0059")
        ig = ig_["R0005"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0059",ix_)
        arrset(pb.xnames,iv,"C0059")
        ig = ig_["R0035"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0060",ix_)
        arrset(pb.xnames,iv,"C0060")
        ig = ig_["R0006"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0060",ix_)
        arrset(pb.xnames,iv,"C0060")
        ig = ig_["R0034"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0060",ix_)
        arrset(pb.xnames,iv,"C0060")
        ig = ig_["R0036"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0061",ix_)
        arrset(pb.xnames,iv,"C0061")
        ig = ig_["R0007"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0061",ix_)
        arrset(pb.xnames,iv,"C0061")
        ig = ig_["R0036"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0062",ix_)
        arrset(pb.xnames,iv,"C0062")
        ig = ig_["R0008"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0062",ix_)
        arrset(pb.xnames,iv,"C0062")
        ig = ig_["R0036"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0063",ix_)
        arrset(pb.xnames,iv,"C0063")
        ig = ig_["R0009"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0063",ix_)
        arrset(pb.xnames,iv,"C0063")
        ig = ig_["R0036"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0064",ix_)
        arrset(pb.xnames,iv,"C0064")
        ig = ig_["R0010"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0064",ix_)
        arrset(pb.xnames,iv,"C0064")
        ig = ig_["R0036"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0065",ix_)
        arrset(pb.xnames,iv,"C0065")
        ig = ig_["R0011"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0065",ix_)
        arrset(pb.xnames,iv,"C0065")
        ig = ig_["R0034"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0065",ix_)
        arrset(pb.xnames,iv,"C0065")
        ig = ig_["R0037"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0066",ix_)
        arrset(pb.xnames,iv,"C0066")
        ig = ig_["R0012"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0066",ix_)
        arrset(pb.xnames,iv,"C0066")
        ig = ig_["R0037"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0067",ix_)
        arrset(pb.xnames,iv,"C0067")
        ig = ig_["R0013"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0067",ix_)
        arrset(pb.xnames,iv,"C0067")
        ig = ig_["R0037"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0068",ix_)
        arrset(pb.xnames,iv,"C0068")
        ig = ig_["R0014"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0068",ix_)
        arrset(pb.xnames,iv,"C0068")
        ig = ig_["R0037"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0069",ix_)
        arrset(pb.xnames,iv,"C0069")
        ig = ig_["R0015"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0069",ix_)
        arrset(pb.xnames,iv,"C0069")
        ig = ig_["R0034"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0069",ix_)
        arrset(pb.xnames,iv,"C0069")
        ig = ig_["R0038"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0070",ix_)
        arrset(pb.xnames,iv,"C0070")
        ig = ig_["R0016"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0070",ix_)
        arrset(pb.xnames,iv,"C0070")
        ig = ig_["R0038"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0071",ix_)
        arrset(pb.xnames,iv,"C0071")
        ig = ig_["R0017"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0071",ix_)
        arrset(pb.xnames,iv,"C0071")
        ig = ig_["R0038"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0072",ix_)
        arrset(pb.xnames,iv,"C0072")
        ig = ig_["R0018"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0072",ix_)
        arrset(pb.xnames,iv,"C0072")
        ig = ig_["R0034"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0073",ix_)
        arrset(pb.xnames,iv,"C0073")
        ig = ig_["R0001"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0073",ix_)
        arrset(pb.xnames,iv,"C0073")
        ig = ig_["R0039"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0074",ix_)
        arrset(pb.xnames,iv,"C0074")
        ig = ig_["R0002"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0074",ix_)
        arrset(pb.xnames,iv,"C0074")
        ig = ig_["R0039"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0074",ix_)
        arrset(pb.xnames,iv,"C0074")
        ig = ig_["R0040"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0075",ix_)
        arrset(pb.xnames,iv,"C0075")
        ig = ig_["R0003"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0075",ix_)
        arrset(pb.xnames,iv,"C0075")
        ig = ig_["R0040"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0076",ix_)
        arrset(pb.xnames,iv,"C0076")
        ig = ig_["R0004"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0076",ix_)
        arrset(pb.xnames,iv,"C0076")
        ig = ig_["R0040"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0077",ix_)
        arrset(pb.xnames,iv,"C0077")
        ig = ig_["R0005"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0077",ix_)
        arrset(pb.xnames,iv,"C0077")
        ig = ig_["R0040"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0078",ix_)
        arrset(pb.xnames,iv,"C0078")
        ig = ig_["R0006"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0078",ix_)
        arrset(pb.xnames,iv,"C0078")
        ig = ig_["R0039"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0078",ix_)
        arrset(pb.xnames,iv,"C0078")
        ig = ig_["R0041"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0079",ix_)
        arrset(pb.xnames,iv,"C0079")
        ig = ig_["R0007"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0079",ix_)
        arrset(pb.xnames,iv,"C0079")
        ig = ig_["R0041"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0080",ix_)
        arrset(pb.xnames,iv,"C0080")
        ig = ig_["R0008"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0080",ix_)
        arrset(pb.xnames,iv,"C0080")
        ig = ig_["R0041"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0081",ix_)
        arrset(pb.xnames,iv,"C0081")
        ig = ig_["R0009"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0081",ix_)
        arrset(pb.xnames,iv,"C0081")
        ig = ig_["R0041"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0082",ix_)
        arrset(pb.xnames,iv,"C0082")
        ig = ig_["R0010"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0082",ix_)
        arrset(pb.xnames,iv,"C0082")
        ig = ig_["R0041"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0083",ix_)
        arrset(pb.xnames,iv,"C0083")
        ig = ig_["R0011"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0083",ix_)
        arrset(pb.xnames,iv,"C0083")
        ig = ig_["R0039"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0083",ix_)
        arrset(pb.xnames,iv,"C0083")
        ig = ig_["R0042"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0084",ix_)
        arrset(pb.xnames,iv,"C0084")
        ig = ig_["R0012"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0084",ix_)
        arrset(pb.xnames,iv,"C0084")
        ig = ig_["R0042"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0085",ix_)
        arrset(pb.xnames,iv,"C0085")
        ig = ig_["R0013"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0085",ix_)
        arrset(pb.xnames,iv,"C0085")
        ig = ig_["R0042"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0086",ix_)
        arrset(pb.xnames,iv,"C0086")
        ig = ig_["R0014"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0086",ix_)
        arrset(pb.xnames,iv,"C0086")
        ig = ig_["R0042"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0087",ix_)
        arrset(pb.xnames,iv,"C0087")
        ig = ig_["R0015"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0087",ix_)
        arrset(pb.xnames,iv,"C0087")
        ig = ig_["R0039"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0087",ix_)
        arrset(pb.xnames,iv,"C0087")
        ig = ig_["R0043"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0088",ix_)
        arrset(pb.xnames,iv,"C0088")
        ig = ig_["R0016"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0088",ix_)
        arrset(pb.xnames,iv,"C0088")
        ig = ig_["R0043"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0089",ix_)
        arrset(pb.xnames,iv,"C0089")
        ig = ig_["R0017"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0089",ix_)
        arrset(pb.xnames,iv,"C0089")
        ig = ig_["R0043"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0090",ix_)
        arrset(pb.xnames,iv,"C0090")
        ig = ig_["R0018"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0090",ix_)
        arrset(pb.xnames,iv,"C0090")
        ig = ig_["R0039"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0091",ix_)
        arrset(pb.xnames,iv,"C0091")
        ig = ig_["R0001"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0091",ix_)
        arrset(pb.xnames,iv,"C0091")
        ig = ig_["R0044"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0092",ix_)
        arrset(pb.xnames,iv,"C0092")
        ig = ig_["R0002"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0092",ix_)
        arrset(pb.xnames,iv,"C0092")
        ig = ig_["R0044"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0092",ix_)
        arrset(pb.xnames,iv,"C0092")
        ig = ig_["R0045"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0093",ix_)
        arrset(pb.xnames,iv,"C0093")
        ig = ig_["R0003"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0093",ix_)
        arrset(pb.xnames,iv,"C0093")
        ig = ig_["R0045"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0094",ix_)
        arrset(pb.xnames,iv,"C0094")
        ig = ig_["R0004"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0094",ix_)
        arrset(pb.xnames,iv,"C0094")
        ig = ig_["R0045"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0095",ix_)
        arrset(pb.xnames,iv,"C0095")
        ig = ig_["R0005"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0095",ix_)
        arrset(pb.xnames,iv,"C0095")
        ig = ig_["R0045"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0096",ix_)
        arrset(pb.xnames,iv,"C0096")
        ig = ig_["R0006"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0096",ix_)
        arrset(pb.xnames,iv,"C0096")
        ig = ig_["R0044"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0096",ix_)
        arrset(pb.xnames,iv,"C0096")
        ig = ig_["R0046"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0097",ix_)
        arrset(pb.xnames,iv,"C0097")
        ig = ig_["R0007"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0097",ix_)
        arrset(pb.xnames,iv,"C0097")
        ig = ig_["R0046"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0098",ix_)
        arrset(pb.xnames,iv,"C0098")
        ig = ig_["R0008"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0098",ix_)
        arrset(pb.xnames,iv,"C0098")
        ig = ig_["R0046"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0099",ix_)
        arrset(pb.xnames,iv,"C0099")
        ig = ig_["R0009"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0099",ix_)
        arrset(pb.xnames,iv,"C0099")
        ig = ig_["R0046"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0100",ix_)
        arrset(pb.xnames,iv,"C0100")
        ig = ig_["R0010"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0100",ix_)
        arrset(pb.xnames,iv,"C0100")
        ig = ig_["R0046"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0101",ix_)
        arrset(pb.xnames,iv,"C0101")
        ig = ig_["R0011"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0101",ix_)
        arrset(pb.xnames,iv,"C0101")
        ig = ig_["R0044"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0101",ix_)
        arrset(pb.xnames,iv,"C0101")
        ig = ig_["R0047"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0102",ix_)
        arrset(pb.xnames,iv,"C0102")
        ig = ig_["R0012"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0102",ix_)
        arrset(pb.xnames,iv,"C0102")
        ig = ig_["R0047"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0103",ix_)
        arrset(pb.xnames,iv,"C0103")
        ig = ig_["R0013"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0103",ix_)
        arrset(pb.xnames,iv,"C0103")
        ig = ig_["R0047"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0104",ix_)
        arrset(pb.xnames,iv,"C0104")
        ig = ig_["R0014"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0104",ix_)
        arrset(pb.xnames,iv,"C0104")
        ig = ig_["R0047"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0105",ix_)
        arrset(pb.xnames,iv,"C0105")
        ig = ig_["R0015"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0105",ix_)
        arrset(pb.xnames,iv,"C0105")
        ig = ig_["R0044"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0105",ix_)
        arrset(pb.xnames,iv,"C0105")
        ig = ig_["R0048"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0106",ix_)
        arrset(pb.xnames,iv,"C0106")
        ig = ig_["R0016"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0106",ix_)
        arrset(pb.xnames,iv,"C0106")
        ig = ig_["R0048"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0107",ix_)
        arrset(pb.xnames,iv,"C0107")
        ig = ig_["R0017"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0107",ix_)
        arrset(pb.xnames,iv,"C0107")
        ig = ig_["R0048"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0108",ix_)
        arrset(pb.xnames,iv,"C0108")
        ig = ig_["R0018"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0108",ix_)
        arrset(pb.xnames,iv,"C0108")
        ig = ig_["R0044"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0109",ix_)
        arrset(pb.xnames,iv,"C0109")
        ig = ig_["R0001"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0109",ix_)
        arrset(pb.xnames,iv,"C0109")
        ig = ig_["R0049"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0110",ix_)
        arrset(pb.xnames,iv,"C0110")
        ig = ig_["R0002"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0110",ix_)
        arrset(pb.xnames,iv,"C0110")
        ig = ig_["R0049"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0110",ix_)
        arrset(pb.xnames,iv,"C0110")
        ig = ig_["R0050"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0111",ix_)
        arrset(pb.xnames,iv,"C0111")
        ig = ig_["R0003"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0111",ix_)
        arrset(pb.xnames,iv,"C0111")
        ig = ig_["R0050"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0112",ix_)
        arrset(pb.xnames,iv,"C0112")
        ig = ig_["R0004"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0112",ix_)
        arrset(pb.xnames,iv,"C0112")
        ig = ig_["R0050"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0113",ix_)
        arrset(pb.xnames,iv,"C0113")
        ig = ig_["R0005"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0113",ix_)
        arrset(pb.xnames,iv,"C0113")
        ig = ig_["R0050"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0114",ix_)
        arrset(pb.xnames,iv,"C0114")
        ig = ig_["R0006"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0114",ix_)
        arrset(pb.xnames,iv,"C0114")
        ig = ig_["R0049"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0114",ix_)
        arrset(pb.xnames,iv,"C0114")
        ig = ig_["R0051"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0115",ix_)
        arrset(pb.xnames,iv,"C0115")
        ig = ig_["R0007"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0115",ix_)
        arrset(pb.xnames,iv,"C0115")
        ig = ig_["R0051"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0116",ix_)
        arrset(pb.xnames,iv,"C0116")
        ig = ig_["R0008"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0116",ix_)
        arrset(pb.xnames,iv,"C0116")
        ig = ig_["R0051"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0117",ix_)
        arrset(pb.xnames,iv,"C0117")
        ig = ig_["R0009"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0117",ix_)
        arrset(pb.xnames,iv,"C0117")
        ig = ig_["R0051"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0118",ix_)
        arrset(pb.xnames,iv,"C0118")
        ig = ig_["R0010"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0118",ix_)
        arrset(pb.xnames,iv,"C0118")
        ig = ig_["R0051"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0119",ix_)
        arrset(pb.xnames,iv,"C0119")
        ig = ig_["R0011"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0119",ix_)
        arrset(pb.xnames,iv,"C0119")
        ig = ig_["R0049"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0119",ix_)
        arrset(pb.xnames,iv,"C0119")
        ig = ig_["R0052"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0120",ix_)
        arrset(pb.xnames,iv,"C0120")
        ig = ig_["R0012"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0120",ix_)
        arrset(pb.xnames,iv,"C0120")
        ig = ig_["R0052"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0121",ix_)
        arrset(pb.xnames,iv,"C0121")
        ig = ig_["R0013"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0121",ix_)
        arrset(pb.xnames,iv,"C0121")
        ig = ig_["R0052"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0122",ix_)
        arrset(pb.xnames,iv,"C0122")
        ig = ig_["R0014"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0122",ix_)
        arrset(pb.xnames,iv,"C0122")
        ig = ig_["R0052"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0123",ix_)
        arrset(pb.xnames,iv,"C0123")
        ig = ig_["R0015"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0123",ix_)
        arrset(pb.xnames,iv,"C0123")
        ig = ig_["R0049"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0123",ix_)
        arrset(pb.xnames,iv,"C0123")
        ig = ig_["R0053"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0124",ix_)
        arrset(pb.xnames,iv,"C0124")
        ig = ig_["R0016"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0124",ix_)
        arrset(pb.xnames,iv,"C0124")
        ig = ig_["R0053"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0125",ix_)
        arrset(pb.xnames,iv,"C0125")
        ig = ig_["R0017"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0125",ix_)
        arrset(pb.xnames,iv,"C0125")
        ig = ig_["R0053"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0126",ix_)
        arrset(pb.xnames,iv,"C0126")
        ig = ig_["R0018"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0126",ix_)
        arrset(pb.xnames,iv,"C0126")
        ig = ig_["R0049"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0127",ix_)
        arrset(pb.xnames,iv,"C0127")
        ig = ig_["R0001"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0127",ix_)
        arrset(pb.xnames,iv,"C0127")
        ig = ig_["R0054"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0128",ix_)
        arrset(pb.xnames,iv,"C0128")
        ig = ig_["R0002"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0128",ix_)
        arrset(pb.xnames,iv,"C0128")
        ig = ig_["R0054"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0128",ix_)
        arrset(pb.xnames,iv,"C0128")
        ig = ig_["R0055"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0129",ix_)
        arrset(pb.xnames,iv,"C0129")
        ig = ig_["R0003"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0129",ix_)
        arrset(pb.xnames,iv,"C0129")
        ig = ig_["R0055"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0130",ix_)
        arrset(pb.xnames,iv,"C0130")
        ig = ig_["R0004"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0130",ix_)
        arrset(pb.xnames,iv,"C0130")
        ig = ig_["R0055"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0131",ix_)
        arrset(pb.xnames,iv,"C0131")
        ig = ig_["R0005"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0131",ix_)
        arrset(pb.xnames,iv,"C0131")
        ig = ig_["R0055"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0132",ix_)
        arrset(pb.xnames,iv,"C0132")
        ig = ig_["R0006"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0132",ix_)
        arrset(pb.xnames,iv,"C0132")
        ig = ig_["R0054"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0132",ix_)
        arrset(pb.xnames,iv,"C0132")
        ig = ig_["R0056"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0133",ix_)
        arrset(pb.xnames,iv,"C0133")
        ig = ig_["R0007"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0133",ix_)
        arrset(pb.xnames,iv,"C0133")
        ig = ig_["R0056"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0134",ix_)
        arrset(pb.xnames,iv,"C0134")
        ig = ig_["R0008"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0134",ix_)
        arrset(pb.xnames,iv,"C0134")
        ig = ig_["R0056"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0135",ix_)
        arrset(pb.xnames,iv,"C0135")
        ig = ig_["R0009"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0135",ix_)
        arrset(pb.xnames,iv,"C0135")
        ig = ig_["R0056"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0136",ix_)
        arrset(pb.xnames,iv,"C0136")
        ig = ig_["R0010"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0136",ix_)
        arrset(pb.xnames,iv,"C0136")
        ig = ig_["R0056"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0137",ix_)
        arrset(pb.xnames,iv,"C0137")
        ig = ig_["R0011"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0137",ix_)
        arrset(pb.xnames,iv,"C0137")
        ig = ig_["R0054"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0137",ix_)
        arrset(pb.xnames,iv,"C0137")
        ig = ig_["R0057"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0138",ix_)
        arrset(pb.xnames,iv,"C0138")
        ig = ig_["R0012"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0138",ix_)
        arrset(pb.xnames,iv,"C0138")
        ig = ig_["R0057"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0139",ix_)
        arrset(pb.xnames,iv,"C0139")
        ig = ig_["R0013"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0139",ix_)
        arrset(pb.xnames,iv,"C0139")
        ig = ig_["R0057"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0140",ix_)
        arrset(pb.xnames,iv,"C0140")
        ig = ig_["R0014"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0140",ix_)
        arrset(pb.xnames,iv,"C0140")
        ig = ig_["R0057"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0141",ix_)
        arrset(pb.xnames,iv,"C0141")
        ig = ig_["R0015"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0141",ix_)
        arrset(pb.xnames,iv,"C0141")
        ig = ig_["R0054"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0141",ix_)
        arrset(pb.xnames,iv,"C0141")
        ig = ig_["R0058"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0142",ix_)
        arrset(pb.xnames,iv,"C0142")
        ig = ig_["R0016"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0142",ix_)
        arrset(pb.xnames,iv,"C0142")
        ig = ig_["R0058"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0143",ix_)
        arrset(pb.xnames,iv,"C0143")
        ig = ig_["R0017"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0143",ix_)
        arrset(pb.xnames,iv,"C0143")
        ig = ig_["R0058"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0144",ix_)
        arrset(pb.xnames,iv,"C0144")
        ig = ig_["R0018"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0144",ix_)
        arrset(pb.xnames,iv,"C0144")
        ig = ig_["R0054"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0145",ix_)
        arrset(pb.xnames,iv,"C0145")
        ig = ig_["R0001"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0145",ix_)
        arrset(pb.xnames,iv,"C0145")
        ig = ig_["R0059"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0146",ix_)
        arrset(pb.xnames,iv,"C0146")
        ig = ig_["R0002"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0146",ix_)
        arrset(pb.xnames,iv,"C0146")
        ig = ig_["R0059"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0146",ix_)
        arrset(pb.xnames,iv,"C0146")
        ig = ig_["R0060"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0147",ix_)
        arrset(pb.xnames,iv,"C0147")
        ig = ig_["R0003"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0147",ix_)
        arrset(pb.xnames,iv,"C0147")
        ig = ig_["R0060"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0148",ix_)
        arrset(pb.xnames,iv,"C0148")
        ig = ig_["R0004"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0148",ix_)
        arrset(pb.xnames,iv,"C0148")
        ig = ig_["R0060"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0149",ix_)
        arrset(pb.xnames,iv,"C0149")
        ig = ig_["R0005"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0149",ix_)
        arrset(pb.xnames,iv,"C0149")
        ig = ig_["R0060"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0150",ix_)
        arrset(pb.xnames,iv,"C0150")
        ig = ig_["R0006"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0150",ix_)
        arrset(pb.xnames,iv,"C0150")
        ig = ig_["R0059"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0150",ix_)
        arrset(pb.xnames,iv,"C0150")
        ig = ig_["R0061"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0151",ix_)
        arrset(pb.xnames,iv,"C0151")
        ig = ig_["R0007"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0151",ix_)
        arrset(pb.xnames,iv,"C0151")
        ig = ig_["R0061"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0152",ix_)
        arrset(pb.xnames,iv,"C0152")
        ig = ig_["R0008"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0152",ix_)
        arrset(pb.xnames,iv,"C0152")
        ig = ig_["R0061"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0153",ix_)
        arrset(pb.xnames,iv,"C0153")
        ig = ig_["R0009"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0153",ix_)
        arrset(pb.xnames,iv,"C0153")
        ig = ig_["R0061"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0154",ix_)
        arrset(pb.xnames,iv,"C0154")
        ig = ig_["R0010"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0154",ix_)
        arrset(pb.xnames,iv,"C0154")
        ig = ig_["R0061"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0155",ix_)
        arrset(pb.xnames,iv,"C0155")
        ig = ig_["R0011"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0155",ix_)
        arrset(pb.xnames,iv,"C0155")
        ig = ig_["R0059"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0155",ix_)
        arrset(pb.xnames,iv,"C0155")
        ig = ig_["R0062"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0156",ix_)
        arrset(pb.xnames,iv,"C0156")
        ig = ig_["R0012"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0156",ix_)
        arrset(pb.xnames,iv,"C0156")
        ig = ig_["R0062"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0157",ix_)
        arrset(pb.xnames,iv,"C0157")
        ig = ig_["R0013"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0157",ix_)
        arrset(pb.xnames,iv,"C0157")
        ig = ig_["R0062"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0158",ix_)
        arrset(pb.xnames,iv,"C0158")
        ig = ig_["R0014"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0158",ix_)
        arrset(pb.xnames,iv,"C0158")
        ig = ig_["R0062"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0159",ix_)
        arrset(pb.xnames,iv,"C0159")
        ig = ig_["R0015"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0159",ix_)
        arrset(pb.xnames,iv,"C0159")
        ig = ig_["R0059"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0159",ix_)
        arrset(pb.xnames,iv,"C0159")
        ig = ig_["R0063"]
        pbm.A[ig,iv] += Float64(-1)
        iv,ix_,_ = s2mpj_ii("C0160",ix_)
        arrset(pb.xnames,iv,"C0160")
        ig = ig_["R0016"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0160",ix_)
        arrset(pb.xnames,iv,"C0160")
        ig = ig_["R0063"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0161",ix_)
        arrset(pb.xnames,iv,"C0161")
        ig = ig_["R0017"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0161",ix_)
        arrset(pb.xnames,iv,"C0161")
        ig = ig_["R0063"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0162",ix_)
        arrset(pb.xnames,iv,"C0162")
        ig = ig_["R0018"]
        pbm.A[ig,iv] += Float64(1)
        iv,ix_,_ = s2mpj_ii("C0162",ix_)
        arrset(pb.xnames,iv,"C0162")
        ig = ig_["R0059"]
        pbm.A[ig,iv] += Float64(1)
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
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower[ix_["C0001"]] = -8423670
        pb.xupper[ix_["C0001"]] = 8423600
        pb.xlower[ix_["C0002"]] = -2186640
        pb.xupper[ix_["C0002"]] = 2186640
        pb.xlower[ix_["C0003"]] = -993066
        pb.xupper[ix_["C0003"]] = 993060
        pb.xlower[ix_["C0004"]] = -904430
        pb.xupper[ix_["C0004"]] = 904430
        pb.xlower[ix_["C0005"]] = -289145
        pb.xupper[ix_["C0005"]] = 289145
        pb.xlower[ix_["C0006"]] = -1851950
        pb.xupper[ix_["C0006"]] = 1851940
        pb.xlower[ix_["C0007"]] = -62168
        pb.xupper[ix_["C0007"]] = 62168
        pb.xlower[ix_["C0008"]] = -263139
        pb.xupper[ix_["C0008"]] = 263139
        pb.xlower[ix_["C0009"]] = -1117500
        pb.xupper[ix_["C0009"]] = 1117490
        pb.xlower[ix_["C0010"]] = -409143
        pb.xupper[ix_["C0010"]] = 409144
        pb.xlower[ix_["C0011"]] = -2288060
        pb.xupper[ix_["C0011"]] = 2288050
        pb.xlower[ix_["C0012"]] = -242663
        pb.xupper[ix_["C0012"]] = 242663
        pb.xlower[ix_["C0013"]] = -1832280
        pb.xupper[ix_["C0013"]] = 1832280
        pb.xlower[ix_["C0014"]] = -213115
        pb.xupper[ix_["C0014"]] = 213115
        pb.xlower[ix_["C0015"]] = -2096980
        pb.xupper[ix_["C0015"]] = 2096990
        pb.xlower[ix_["C0016"]] = -1376370
        pb.xupper[ix_["C0016"]] = 1376370
        pb.xlower[ix_["C0017"]] = -720616
        pb.xupper[ix_["C0017"]] = 720610
        pb.xupper[ix_["C0018"]] = 0
        pb.xlower[ix_["C0019"]] = 6
        pb.xupper[ix_["C0019"]] = 10
        pb.xlower[ix_["C0020"]] = 1.67
        pb.xupper[ix_["C0020"]] = 12.5
        pb.xlower[ix_["C0021"]] = 1.67
        pb.xupper[ix_["C0021"]] = 2.5
        pb.xupper[ix_["C0022"]] = 0
        pb.xupper[ix_["C0023"]] = 0
        pb.xlower[ix_["C0024"]] = 5
        pb.xupper[ix_["C0024"]] = 7.5
        pb.xlower[ix_["C0025"]] = 1.67
        pb.xupper[ix_["C0025"]] = 2.5
        pb.xupper[ix_["C0026"]] = 0
        pb.xlower[ix_["C0027"]] = 3.33
        pb.xupper[ix_["C0027"]] = 5
        pb.xupper[ix_["C0028"]] = 0
        pb.xupper[ix_["C0029"]] = 0
        pb.xupper[ix_["C0030"]] = 0
        pb.xupper[ix_["C0031"]] = 0
        pb.xupper[ix_["C0032"]] = 0
        pb.xupper[ix_["C0033"]] = 0
        pb.xupper[ix_["C0034"]] = 0
        pb.xupper[ix_["C0035"]] = 0
        pb.xupper[ix_["C0036"]] = 0
        pb.xlower[ix_["C0037"]] = -12.5
        pb.xupper[ix_["C0037"]] = 12.5
        pb.xlower[ix_["C0038"]] = 1.67
        pb.xupper[ix_["C0038"]] = 2.5
        pb.xlower[ix_["C0039"]] = 1.67
        pb.xupper[ix_["C0039"]] = 2.5
        pb.xupper[ix_["C0040"]] = 0
        pb.xupper[ix_["C0041"]] = 0
        pb.xlower[ix_["C0042"]] = 1.67
        pb.xupper[ix_["C0042"]] = 2.5
        pb.xupper[ix_["C0043"]] = 0
        pb.xupper[ix_["C0044"]] = 0
        pb.xlower[ix_["C0045"]] = 1.67
        pb.xupper[ix_["C0045"]] = 2.5
        pb.xupper[ix_["C0046"]] = 0
        pb.xupper[ix_["C0047"]] = 0
        pb.xupper[ix_["C0048"]] = 0
        pb.xupper[ix_["C0049"]] = 0
        pb.xupper[ix_["C0050"]] = 0
        pb.xlower[ix_["C0051"]] = 5
        pb.xupper[ix_["C0051"]] = 7.5
        pb.xlower[ix_["C0052"]] = 5
        pb.xupper[ix_["C0052"]] = 7.5
        pb.xupper[ix_["C0053"]] = 0
        pb.xupper[ix_["C0054"]] = 0
        pb.xlower[ix_["C0055"]] = -1355910
        pb.xupper[ix_["C0055"]] = 1355900
        pb.xlower[ix_["C0056"]] = -359525
        pb.xupper[ix_["C0056"]] = 359521
        pb.xlower[ix_["C0057"]] = -199031
        pb.xupper[ix_["C0057"]] = 199031
        pb.xlower[ix_["C0058"]] = -111995
        pb.xupper[ix_["C0058"]] = 111995
        pb.xlower[ix_["C0059"]] = -48498.5
        pb.xupper[ix_["C0059"]] = 48499
        pb.xlower[ix_["C0060"]] = -321119
        pb.xupper[ix_["C0060"]] = 321119
        pb.xlower[ix_["C0061"]] = -18155.5
        pb.xupper[ix_["C0061"]] = 18155.5
        pb.xlower[ix_["C0062"]] = -46794.5
        pb.xupper[ix_["C0062"]] = 46795
        pb.xlower[ix_["C0063"]] = -172901
        pb.xupper[ix_["C0063"]] = 172901
        pb.xlower[ix_["C0064"]] = -83267.5
        pb.xupper[ix_["C0064"]] = 83267
        pb.xlower[ix_["C0065"]] = -324486
        pb.xupper[ix_["C0065"]] = 324486
        pb.xlower[ix_["C0066"]] = -31883.5
        pb.xupper[ix_["C0066"]] = 31883.5
        pb.xlower[ix_["C0067"]] = -268955
        pb.xupper[ix_["C0067"]] = 268955
        pb.xlower[ix_["C0068"]] = -23647
        pb.xupper[ix_["C0068"]] = 23647
        pb.xlower[ix_["C0069"]] = -350775
        pb.xupper[ix_["C0069"]] = 350771
        pb.xlower[ix_["C0070"]] = -244307
        pb.xupper[ix_["C0070"]] = 244307
        pb.xlower[ix_["C0071"]] = -106468
        pb.xupper[ix_["C0071"]] = 106468
        pb.xupper[ix_["C0072"]] = 0
        pb.xlower[ix_["C0073"]] = -1160260
        pb.xupper[ix_["C0073"]] = 1160270
        pb.xlower[ix_["C0074"]] = -329840
        pb.xupper[ix_["C0074"]] = 329840
        pb.xlower[ix_["C0075"]] = -174019
        pb.xupper[ix_["C0075"]] = 174019
        pb.xlower[ix_["C0076"]] = -110666
        pb.xupper[ix_["C0076"]] = 110666
        pb.xlower[ix_["C0077"]] = -45154.5
        pb.xupper[ix_["C0077"]] = 45155
        pb.xlower[ix_["C0078"]] = -257501
        pb.xupper[ix_["C0078"]] = 257501
        pb.xlower[ix_["C0079"]] = -16066
        pb.xupper[ix_["C0079"]] = 16066
        pb.xlower[ix_["C0080"]] = -47478.5
        pb.xupper[ix_["C0080"]] = 47479
        pb.xlower[ix_["C0081"]] = -125679
        pb.xupper[ix_["C0081"]] = 125679
        pb.xlower[ix_["C0082"]] = -68278
        pb.xupper[ix_["C0082"]] = 68278
        pb.xlower[ix_["C0083"]] = -271785
        pb.xupper[ix_["C0083"]] = 271785
        pb.xlower[ix_["C0084"]] = -37721
        pb.xupper[ix_["C0084"]] = 37721
        pb.xlower[ix_["C0085"]] = -215425
        pb.xupper[ix_["C0085"]] = 215425
        pb.xlower[ix_["C0086"]] = -18638.5
        pb.xupper[ix_["C0086"]] = 18638.5
        pb.xlower[ix_["C0087"]] = -301141
        pb.xupper[ix_["C0087"]] = 301141
        pb.xlower[ix_["C0088"]] = -196197
        pb.xupper[ix_["C0088"]] = 196197
        pb.xlower[ix_["C0089"]] = -104943
        pb.xupper[ix_["C0089"]] = 104943
        pb.xupper[ix_["C0090"]] = 0
        pb.xlower[ix_["C0091"]] = -1252520
        pb.xupper[ix_["C0091"]] = 1252520
        pb.xlower[ix_["C0092"]] = -344481
        pb.xupper[ix_["C0092"]] = 344478
        pb.xlower[ix_["C0093"]] = -177355
        pb.xupper[ix_["C0093"]] = 177355
        pb.xlower[ix_["C0094"]] = -120957
        pb.xupper[ix_["C0094"]] = 120957
        pb.xlower[ix_["C0095"]] = -46169
        pb.xupper[ix_["C0095"]] = 46169
        pb.xlower[ix_["C0096"]] = -267073
        pb.xupper[ix_["C0096"]] = 267073
        pb.xlower[ix_["C0097"]] = -12885
        pb.xupper[ix_["C0097"]] = 12885
        pb.xlower[ix_["C0098"]] = -55465
        pb.xupper[ix_["C0098"]] = 55465
        pb.xlower[ix_["C0099"]] = -125594
        pb.xupper[ix_["C0099"]] = 125594
        pb.xlower[ix_["C0100"]] = -73129.5
        pb.xupper[ix_["C0100"]] = 73129
        pb.xlower[ix_["C0101"]] = -331949
        pb.xupper[ix_["C0101"]] = 331948
        pb.xlower[ix_["C0102"]] = -43652.5
        pb.xupper[ix_["C0102"]] = 43653
        pb.xlower[ix_["C0103"]] = -257510
        pb.xupper[ix_["C0103"]] = 257509
        pb.xlower[ix_["C0104"]] = -30786
        pb.xupper[ix_["C0104"]] = 30786
        pb.xlower[ix_["C0105"]] = -309019
        pb.xupper[ix_["C0105"]] = 309019
        pb.xlower[ix_["C0106"]] = -181745
        pb.xupper[ix_["C0106"]] = 181745
        pb.xlower[ix_["C0107"]] = -127273
        pb.xupper[ix_["C0107"]] = 127273
        pb.xupper[ix_["C0108"]] = 0
        pb.xlower[ix_["C0109"]] = -1399530
        pb.xupper[ix_["C0109"]] = 1399540
        pb.xlower[ix_["C0110"]] = -378265
        pb.xupper[ix_["C0110"]] = 378261
        pb.xlower[ix_["C0111"]] = -209389
        pb.xupper[ix_["C0111"]] = 209389
        pb.xlower[ix_["C0112"]] = -129117
        pb.xupper[ix_["C0112"]] = 129117
        pb.xlower[ix_["C0113"]] = -39759
        pb.xupper[ix_["C0113"]] = 39759
        pb.xlower[ix_["C0114"]] = -310196
        pb.xupper[ix_["C0114"]] = 310196
        pb.xlower[ix_["C0115"]] = -9075
        pb.xupper[ix_["C0115"]] = 9075
        pb.xlower[ix_["C0116"]] = -40899.5
        pb.xupper[ix_["C0116"]] = 40899
        pb.xlower[ix_["C0117"]] = -151689
        pb.xupper[ix_["C0117"]] = 151689
        pb.xlower[ix_["C0118"]] = -108533
        pb.xupper[ix_["C0118"]] = 108533
        pb.xlower[ix_["C0119"]] = -387566
        pb.xupper[ix_["C0119"]] = 387568
        pb.xlower[ix_["C0120"]] = -29976.5
        pb.xupper[ix_["C0120"]] = 29976.5
        pb.xlower[ix_["C0121"]] = -321881
        pb.xupper[ix_["C0121"]] = 321881
        pb.xlower[ix_["C0122"]] = -35708.5
        pb.xupper[ix_["C0122"]] = 35709
        pb.xlower[ix_["C0123"]] = -323511
        pb.xupper[ix_["C0123"]] = 323511
        pb.xlower[ix_["C0124"]] = -201463
        pb.xupper[ix_["C0124"]] = 201463
        pb.xlower[ix_["C0125"]] = -122048
        pb.xupper[ix_["C0125"]] = 122048
        pb.xupper[ix_["C0126"]] = 0
        pb.xlower[ix_["C0127"]] = -3255380
        pb.xupper[ix_["C0127"]] = 3255380
        pb.xlower[ix_["C0128"]] = -774526
        pb.xupper[ix_["C0128"]] = 774520
        pb.xlower[ix_["C0129"]] = -233265
        pb.xupper[ix_["C0129"]] = 233265
        pb.xlower[ix_["C0130"]] = -431697
        pb.xupper[ix_["C0130"]] = 431697
        pb.xlower[ix_["C0131"]] = -109563
        pb.xupper[ix_["C0131"]] = 109563
        pb.xlower[ix_["C0132"]] = -696052
        pb.xupper[ix_["C0132"]] = 696040
        pb.xlower[ix_["C0133"]] = 3480
        pb.xupper[ix_["C0133"]] = 5984
        pb.xlower[ix_["C0134"]] = -72502
        pb.xupper[ix_["C0134"]] = 72502
        pb.xlower[ix_["C0135"]] = -541623
        pb.xupper[ix_["C0135"]] = 541630
        pb.xlower[ix_["C0136"]] = -75935
        pb.xupper[ix_["C0136"]] = 75935
        pb.xlower[ix_["C0137"]] = -972268
        pb.xupper[ix_["C0137"]] = 972280
        pb.xlower[ix_["C0138"]] = -99429.5
        pb.xupper[ix_["C0138"]] = 99429
        pb.xlower[ix_["C0139"]] = -768512
        pb.xupper[ix_["C0139"]] = 768500
        pb.xlower[ix_["C0140"]] = -104335
        pb.xupper[ix_["C0140"]] = 104335
        pb.xlower[ix_["C0141"]] = -812536
        pb.xupper[ix_["C0141"]] = 812530
        pb.xlower[ix_["C0142"]] = -552648
        pb.xupper[ix_["C0142"]] = 552660
        pb.xlower[ix_["C0143"]] = -259881
        pb.xupper[ix_["C0143"]] = 259881
        pb.xupper[ix_["C0144"]] = 0
        pb.xupper[ix_["C0145"]] = 0
        pb.xupper[ix_["C0146"]] = 0
        pb.xupper[ix_["C0147"]] = 0
        pb.xupper[ix_["C0148"]] = 0
        pb.xupper[ix_["C0149"]] = 0
        pb.xupper[ix_["C0150"]] = 0
        pb.xupper[ix_["C0151"]] = 0
        pb.xupper[ix_["C0152"]] = 0
        pb.xupper[ix_["C0153"]] = 0
        pb.xupper[ix_["C0154"]] = 0
        pb.xupper[ix_["C0155"]] = 0
        pb.xupper[ix_["C0156"]] = 0
        pb.xupper[ix_["C0157"]] = 0
        pb.xupper[ix_["C0158"]] = 0
        pb.xupper[ix_["C0159"]] = 0
        pb.xupper[ix_["C0160"]] = 0
        pb.xupper[ix_["C0161"]] = 0
        pb.xupper[ix_["C0162"]] = 0
        #%%%%%%%%%%%%%%%%%%%% QUADRATIC %%%%%%%%%%%%%%%%%%%
        ix1 = ix_["C0001"]
        ix2 = ix_["C0001"]
        pbm.H[ix1,ix2] = Float64(1.187134e-07)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0002"]
        ix2 = ix_["C0002"]
        pbm.H[ix1,ix2] = Float64(4.573226e-07)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0003"]
        ix2 = ix_["C0003"]
        pbm.H[ix1,ix2] = Float64(1.006983e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0004"]
        ix2 = ix_["C0004"]
        pbm.H[ix1,ix2] = Float64(1.105669e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0005"]
        ix2 = ix_["C0005"]
        pbm.H[ix1,ix2] = Float64(3.458478e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0006"]
        ix2 = ix_["C0006"]
        pbm.H[ix1,ix2] = Float64(5.399714e-07)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0007"]
        ix2 = ix_["C0007"]
        pbm.H[ix1,ix2] = Float64(1.608545e-05)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0008"]
        ix2 = ix_["C0008"]
        pbm.H[ix1,ix2] = Float64(3.800266e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0009"]
        ix2 = ix_["C0009"]
        pbm.H[ix1,ix2] = Float64(8.948546e-07)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0010"]
        ix2 = ix_["C0010"]
        pbm.H[ix1,ix2] = Float64(2.444134e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0011"]
        ix2 = ix_["C0011"]
        pbm.H[ix1,ix2] = Float64(4.370514e-07)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0012"]
        ix2 = ix_["C0012"]
        pbm.H[ix1,ix2] = Float64(4.120942e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0013"]
        ix2 = ix_["C0013"]
        pbm.H[ix1,ix2] = Float64(5.457682e-07)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0014"]
        ix2 = ix_["C0014"]
        pbm.H[ix1,ix2] = Float64(4.692302e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0015"]
        ix2 = ix_["C0015"]
        pbm.H[ix1,ix2] = Float64(4.768752e-07)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0016"]
        ix2 = ix_["C0016"]
        pbm.H[ix1,ix2] = Float64(7.265488e-07)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0017"]
        ix2 = ix_["C0017"]
        pbm.H[ix1,ix2] = Float64(1.387704e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0018"]
        ix2 = ix_["C0018"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0019"]
        ix2 = ix_["C0019"]
        pbm.H[ix1,ix2] = Float64(1.000000e-01)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0020"]
        ix2 = ix_["C0020"]
        pbm.H[ix1,ix2] = Float64(4.000000e-01)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0021"]
        ix2 = ix_["C0021"]
        pbm.H[ix1,ix2] = Float64(4.000000e-01)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0022"]
        ix2 = ix_["C0022"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0023"]
        ix2 = ix_["C0023"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0024"]
        ix2 = ix_["C0024"]
        pbm.H[ix1,ix2] = Float64(1.333333e-01)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0025"]
        ix2 = ix_["C0025"]
        pbm.H[ix1,ix2] = Float64(4.000000e-01)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0026"]
        ix2 = ix_["C0026"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0027"]
        ix2 = ix_["C0027"]
        pbm.H[ix1,ix2] = Float64(2.000000e-01)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0028"]
        ix2 = ix_["C0028"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0029"]
        ix2 = ix_["C0029"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0030"]
        ix2 = ix_["C0030"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0031"]
        ix2 = ix_["C0031"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0032"]
        ix2 = ix_["C0032"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0033"]
        ix2 = ix_["C0033"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0034"]
        ix2 = ix_["C0034"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0035"]
        ix2 = ix_["C0035"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0036"]
        ix2 = ix_["C0036"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0037"]
        ix2 = ix_["C0037"]
        pbm.H[ix1,ix2] = Float64(8.000000e-02)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0038"]
        ix2 = ix_["C0038"]
        pbm.H[ix1,ix2] = Float64(4.000000e-01)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0039"]
        ix2 = ix_["C0039"]
        pbm.H[ix1,ix2] = Float64(4.000000e-01)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0040"]
        ix2 = ix_["C0040"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0041"]
        ix2 = ix_["C0041"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0042"]
        ix2 = ix_["C0042"]
        pbm.H[ix1,ix2] = Float64(4.000000e-01)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0043"]
        ix2 = ix_["C0043"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0044"]
        ix2 = ix_["C0044"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0045"]
        ix2 = ix_["C0045"]
        pbm.H[ix1,ix2] = Float64(4.000000e-01)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0046"]
        ix2 = ix_["C0046"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0047"]
        ix2 = ix_["C0047"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0048"]
        ix2 = ix_["C0048"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0049"]
        ix2 = ix_["C0049"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0050"]
        ix2 = ix_["C0050"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0051"]
        ix2 = ix_["C0051"]
        pbm.H[ix1,ix2] = Float64(1.333333e-01)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0052"]
        ix2 = ix_["C0052"]
        pbm.H[ix1,ix2] = Float64(1.333333e-01)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0053"]
        ix2 = ix_["C0053"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0054"]
        ix2 = ix_["C0054"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0055"]
        ix2 = ix_["C0055"]
        pbm.H[ix1,ix2] = Float64(7.375148e-07)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0056"]
        ix2 = ix_["C0056"]
        pbm.H[ix1,ix2] = Float64(2.781452e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0057"]
        ix2 = ix_["C0057"]
        pbm.H[ix1,ix2] = Float64(5.024342e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0058"]
        ix2 = ix_["C0058"]
        pbm.H[ix1,ix2] = Float64(8.928970e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0059"]
        ix2 = ix_["C0059"]
        pbm.H[ix1,ix2] = Float64(2.061920e-05)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0060"]
        ix2 = ix_["C0060"]
        pbm.H[ix1,ix2] = Float64(3.114110e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0061"]
        ix2 = ix_["C0061"]
        pbm.H[ix1,ix2] = Float64(5.507972e-05)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0062"]
        ix2 = ix_["C0062"]
        pbm.H[ix1,ix2] = Float64(2.137004e-05)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0063"]
        ix2 = ix_["C0063"]
        pbm.H[ix1,ix2] = Float64(5.783640e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0064"]
        ix2 = ix_["C0064"]
        pbm.H[ix1,ix2] = Float64(1.200949e-05)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0065"]
        ix2 = ix_["C0065"]
        pbm.H[ix1,ix2] = Float64(3.081798e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0066"]
        ix2 = ix_["C0066"]
        pbm.H[ix1,ix2] = Float64(3.136418e-05)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0067"]
        ix2 = ix_["C0067"]
        pbm.H[ix1,ix2] = Float64(3.718088e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0068"]
        ix2 = ix_["C0068"]
        pbm.H[ix1,ix2] = Float64(4.228866e-05)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0069"]
        ix2 = ix_["C0069"]
        pbm.H[ix1,ix2] = Float64(2.850834e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0070"]
        ix2 = ix_["C0070"]
        pbm.H[ix1,ix2] = Float64(4.093218e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0071"]
        ix2 = ix_["C0071"]
        pbm.H[ix1,ix2] = Float64(9.392494e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0072"]
        ix2 = ix_["C0072"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0073"]
        ix2 = ix_["C0073"]
        pbm.H[ix1,ix2] = Float64(8.618720e-07)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0074"]
        ix2 = ix_["C0074"]
        pbm.H[ix1,ix2] = Float64(3.031772e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0075"]
        ix2 = ix_["C0075"]
        pbm.H[ix1,ix2] = Float64(5.746482e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0076"]
        ix2 = ix_["C0076"]
        pbm.H[ix1,ix2] = Float64(9.036200e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0077"]
        ix2 = ix_["C0077"]
        pbm.H[ix1,ix2] = Float64(2.214618e-05)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0078"]
        ix2 = ix_["C0078"]
        pbm.H[ix1,ix2] = Float64(3.883472e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0079"]
        ix2 = ix_["C0079"]
        pbm.H[ix1,ix2] = Float64(6.224324e-05)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0080"]
        ix2 = ix_["C0080"]
        pbm.H[ix1,ix2] = Float64(2.106216e-05)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0081"]
        ix2 = ix_["C0081"]
        pbm.H[ix1,ix2] = Float64(7.956778e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0082"]
        ix2 = ix_["C0082"]
        pbm.H[ix1,ix2] = Float64(1.464601e-05)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0083"]
        ix2 = ix_["C0083"]
        pbm.H[ix1,ix2] = Float64(3.679378e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0084"]
        ix2 = ix_["C0084"]
        pbm.H[ix1,ix2] = Float64(2.651044e-05)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0085"]
        ix2 = ix_["C0085"]
        pbm.H[ix1,ix2] = Float64(4.641976e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0086"]
        ix2 = ix_["C0086"]
        pbm.H[ix1,ix2] = Float64(5.365238e-05)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0087"]
        ix2 = ix_["C0087"]
        pbm.H[ix1,ix2] = Float64(3.320710e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0088"]
        ix2 = ix_["C0088"]
        pbm.H[ix1,ix2] = Float64(5.096904e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0089"]
        ix2 = ix_["C0089"]
        pbm.H[ix1,ix2] = Float64(9.528982e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0090"]
        ix2 = ix_["C0090"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0091"]
        ix2 = ix_["C0091"]
        pbm.H[ix1,ix2] = Float64(7.983904e-07)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0092"]
        ix2 = ix_["C0092"]
        pbm.H[ix1,ix2] = Float64(2.902918e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0093"]
        ix2 = ix_["C0093"]
        pbm.H[ix1,ix2] = Float64(5.638392e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0094"]
        ix2 = ix_["C0094"]
        pbm.H[ix1,ix2] = Float64(8.267434e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0095"]
        ix2 = ix_["C0095"]
        pbm.H[ix1,ix2] = Float64(2.165956e-05)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0096"]
        ix2 = ix_["C0096"]
        pbm.H[ix1,ix2] = Float64(3.744288e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0097"]
        ix2 = ix_["C0097"]
        pbm.H[ix1,ix2] = Float64(7.760962e-05)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0098"]
        ix2 = ix_["C0098"]
        pbm.H[ix1,ix2] = Float64(1.802939e-05)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0099"]
        ix2 = ix_["C0099"]
        pbm.H[ix1,ix2] = Float64(7.962164e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0100"]
        ix2 = ix_["C0100"]
        pbm.H[ix1,ix2] = Float64(1.367437e-05)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0101"]
        ix2 = ix_["C0101"]
        pbm.H[ix1,ix2] = Float64(3.012516e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0102"]
        ix2 = ix_["C0102"]
        pbm.H[ix1,ix2] = Float64(2.290820e-05)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0103"]
        ix2 = ix_["C0103"]
        pbm.H[ix1,ix2] = Float64(3.883344e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0104"]
        ix2 = ix_["C0104"]
        pbm.H[ix1,ix2] = Float64(3.248230e-05)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0105"]
        ix2 = ix_["C0105"]
        pbm.H[ix1,ix2] = Float64(3.236052e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0106"]
        ix2 = ix_["C0106"]
        pbm.H[ix1,ix2] = Float64(5.502214e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0107"]
        ix2 = ix_["C0107"]
        pbm.H[ix1,ix2] = Float64(7.857096e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0108"]
        ix2 = ix_["C0108"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0109"]
        ix2 = ix_["C0109"]
        pbm.H[ix1,ix2] = Float64(7.145230e-07)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0110"]
        ix2 = ix_["C0110"]
        pbm.H[ix1,ix2] = Float64(2.643652e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0111"]
        ix2 = ix_["C0111"]
        pbm.H[ix1,ix2] = Float64(4.775800e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0112"]
        ix2 = ix_["C0112"]
        pbm.H[ix1,ix2] = Float64(7.744944e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0113"]
        ix2 = ix_["C0113"]
        pbm.H[ix1,ix2] = Float64(2.515154e-05)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0114"]
        ix2 = ix_["C0114"]
        pbm.H[ix1,ix2] = Float64(3.223768e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0115"]
        ix2 = ix_["C0115"]
        pbm.H[ix1,ix2] = Float64(1.101928e-04)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0116"]
        ix2 = ix_["C0116"]
        pbm.H[ix1,ix2] = Float64(2.445018e-05)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0117"]
        ix2 = ix_["C0117"]
        pbm.H[ix1,ix2] = Float64(6.592458e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0118"]
        ix2 = ix_["C0118"]
        pbm.H[ix1,ix2] = Float64(9.213788e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0119"]
        ix2 = ix_["C0119"]
        pbm.H[ix1,ix2] = Float64(2.580206e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0120"]
        ix2 = ix_["C0120"]
        pbm.H[ix1,ix2] = Float64(3.335946e-05)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0121"]
        ix2 = ix_["C0121"]
        pbm.H[ix1,ix2] = Float64(3.106738e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0122"]
        ix2 = ix_["C0122"]
        pbm.H[ix1,ix2] = Float64(2.800454e-05)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0123"]
        ix2 = ix_["C0123"]
        pbm.H[ix1,ix2] = Float64(3.091090e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0124"]
        ix2 = ix_["C0124"]
        pbm.H[ix1,ix2] = Float64(4.963702e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0125"]
        ix2 = ix_["C0125"]
        pbm.H[ix1,ix2] = Float64(8.193498e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0126"]
        ix2 = ix_["C0126"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0127"]
        ix2 = ix_["C0127"]
        pbm.H[ix1,ix2] = Float64(3.071838e-07)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0128"]
        ix2 = ix_["C0128"]
        pbm.H[ix1,ix2] = Float64(1.291114e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0129"]
        ix2 = ix_["C0129"]
        pbm.H[ix1,ix2] = Float64(4.286978e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0130"]
        ix2 = ix_["C0130"]
        pbm.H[ix1,ix2] = Float64(2.316442e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0131"]
        ix2 = ix_["C0131"]
        pbm.H[ix1,ix2] = Float64(9.127128e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0132"]
        ix2 = ix_["C0132"]
        pbm.H[ix1,ix2] = Float64(1.436678e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0133"]
        ix2 = ix_["C0133"]
        pbm.H[ix1,ix2] = Float64(1.671123e-04)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0134"]
        ix2 = ix_["C0134"]
        pbm.H[ix1,ix2] = Float64(1.379272e-05)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0135"]
        ix2 = ix_["C0135"]
        pbm.H[ix1,ix2] = Float64(1.846296e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0136"]
        ix2 = ix_["C0136"]
        pbm.H[ix1,ix2] = Float64(1.316916e-05)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0137"]
        ix2 = ix_["C0137"]
        pbm.H[ix1,ix2] = Float64(1.028521e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0138"]
        ix2 = ix_["C0138"]
        pbm.H[ix1,ix2] = Float64(1.005738e-05)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0139"]
        ix2 = ix_["C0139"]
        pbm.H[ix1,ix2] = Float64(1.301219e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0140"]
        ix2 = ix_["C0140"]
        pbm.H[ix1,ix2] = Float64(9.584512e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0141"]
        ix2 = ix_["C0141"]
        pbm.H[ix1,ix2] = Float64(1.230716e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0142"]
        ix2 = ix_["C0142"]
        pbm.H[ix1,ix2] = Float64(1.809463e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0143"]
        ix2 = ix_["C0143"]
        pbm.H[ix1,ix2] = Float64(3.847908e-06)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0144"]
        ix2 = ix_["C0144"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0145"]
        ix2 = ix_["C0145"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0146"]
        ix2 = ix_["C0146"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0147"]
        ix2 = ix_["C0147"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0148"]
        ix2 = ix_["C0148"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0149"]
        ix2 = ix_["C0149"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0150"]
        ix2 = ix_["C0150"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0151"]
        ix2 = ix_["C0151"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0152"]
        ix2 = ix_["C0152"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0153"]
        ix2 = ix_["C0153"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0154"]
        ix2 = ix_["C0154"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0155"]
        ix2 = ix_["C0155"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0156"]
        ix2 = ix_["C0156"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0157"]
        ix2 = ix_["C0157"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0158"]
        ix2 = ix_["C0158"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0159"]
        ix2 = ix_["C0159"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0160"]
        ix2 = ix_["C0160"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0161"]
        ix2 = ix_["C0161"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["C0162"]
        ix2 = ix_["C0162"]
        pbm.H[ix1,ix2] = Float64(2.000000e+00)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        Hsave = pbm.H[ 1:pb.n, 1:pb.n ]
        pbm.H = Hsave
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons   = collect(1:length(pbm.congrps))
        pb.pbclass = "C-QLR2-RN-162-63"
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

