function AGG(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    An LP with degeneracy
# 
#    Source:
#    The NETLIB collection of test problems.
# 
#    SIF input: (already in MPS format)
# 
#    classification = "C-LLR2-AN-163-488"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 6 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "AGG"

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
        ig,ig_,_ = s2mpj_ii("CAP00101",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00101")
        ig,ig_,_ = s2mpj_ii("CAP00201",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00201")
        ig,ig_,_ = s2mpj_ii("CAP00301",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00301")
        ig,ig_,_ = s2mpj_ii("CAP00401",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00401")
        ig,ig_,_ = s2mpj_ii("CAP00501",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00501")
        ig,ig_,_ = s2mpj_ii("CAP00601",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00601")
        ig,ig_,_ = s2mpj_ii("CAP00701",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00701")
        ig,ig_,_ = s2mpj_ii("CAP00801",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00801")
        ig,ig_,_ = s2mpj_ii("CAP00901",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00901")
        ig,ig_,_ = s2mpj_ii("CAP01001",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01001")
        ig,ig_,_ = s2mpj_ii("CAP01101",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01101")
        ig,ig_,_ = s2mpj_ii("CAP01201",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01201")
        ig,ig_,_ = s2mpj_ii("CAP01301",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01301")
        ig,ig_,_ = s2mpj_ii("CAP01401",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01401")
        ig,ig_,_ = s2mpj_ii("CAP01501",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01501")
        ig,ig_,_ = s2mpj_ii("CAP01601",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01601")
        ig,ig_,_ = s2mpj_ii("CAP01701",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01701")
        ig,ig_,_ = s2mpj_ii("CAP01801",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01801")
        ig,ig_,_ = s2mpj_ii("CAP02001",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02001")
        ig,ig_,_ = s2mpj_ii("CAP02201",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02201")
        ig,ig_,_ = s2mpj_ii("CAP02501",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02501")
        ig,ig_,_ = s2mpj_ii("CAP02701",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02701")
        ig,ig_,_ = s2mpj_ii("CAP03101",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03101")
        ig,ig_,_ = s2mpj_ii("CAP03301",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03301")
        ig,ig_,_ = s2mpj_ii("CAP03701",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03701")
        ig,ig_,_ = s2mpj_ii("CAP04101",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04101")
        ig,ig_,_ = s2mpj_ii("CAP04301",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04301")
        ig,ig_,_ = s2mpj_ii("CAP04601",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04601")
        ig,ig_,_ = s2mpj_ii("CAP05501",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05501")
        ig,ig_,_ = s2mpj_ii("CAP06101",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP06101")
        ig,ig_,_ = s2mpj_ii("CAP06201",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP06201")
        ig,ig_,_ = s2mpj_ii("CAP06301",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP06301")
        ig,ig_,_ = s2mpj_ii("CAP06501",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP06501")
        ig,ig_,_ = s2mpj_ii("CAP00102",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00102")
        ig,ig_,_ = s2mpj_ii("CAP00202",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00202")
        ig,ig_,_ = s2mpj_ii("CAP00302",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00302")
        ig,ig_,_ = s2mpj_ii("CAP00402",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00402")
        ig,ig_,_ = s2mpj_ii("CAP00502",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00502")
        ig,ig_,_ = s2mpj_ii("CAP00602",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00602")
        ig,ig_,_ = s2mpj_ii("CAP00702",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00702")
        ig,ig_,_ = s2mpj_ii("CAP00802",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00802")
        ig,ig_,_ = s2mpj_ii("CAP00902",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00902")
        ig,ig_,_ = s2mpj_ii("CAP01002",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01002")
        ig,ig_,_ = s2mpj_ii("CAP01102",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01102")
        ig,ig_,_ = s2mpj_ii("CAP01202",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01202")
        ig,ig_,_ = s2mpj_ii("CAP01302",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01302")
        ig,ig_,_ = s2mpj_ii("CAP01402",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01402")
        ig,ig_,_ = s2mpj_ii("CAP01502",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01502")
        ig,ig_,_ = s2mpj_ii("CAP01602",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01602")
        ig,ig_,_ = s2mpj_ii("CAP01702",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01702")
        ig,ig_,_ = s2mpj_ii("CAP01802",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01802")
        ig,ig_,_ = s2mpj_ii("CAP01902",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01902")
        ig,ig_,_ = s2mpj_ii("CAP02002",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02002")
        ig,ig_,_ = s2mpj_ii("CAP02102",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02102")
        ig,ig_,_ = s2mpj_ii("CAP02202",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02202")
        ig,ig_,_ = s2mpj_ii("CAP02302",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02302")
        ig,ig_,_ = s2mpj_ii("CAP02402",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02402")
        ig,ig_,_ = s2mpj_ii("CAP02502",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02502")
        ig,ig_,_ = s2mpj_ii("CAP02602",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02602")
        ig,ig_,_ = s2mpj_ii("CAP02702",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02702")
        ig,ig_,_ = s2mpj_ii("CAP02802",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02802")
        ig,ig_,_ = s2mpj_ii("CAP02902",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02902")
        ig,ig_,_ = s2mpj_ii("CAP03002",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03002")
        ig,ig_,_ = s2mpj_ii("CAP03102",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03102")
        ig,ig_,_ = s2mpj_ii("CAP03202",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03202")
        ig,ig_,_ = s2mpj_ii("CAP03302",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03302")
        ig,ig_,_ = s2mpj_ii("CAP03402",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03402")
        ig,ig_,_ = s2mpj_ii("CAP03502",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03502")
        ig,ig_,_ = s2mpj_ii("CAP03602",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03602")
        ig,ig_,_ = s2mpj_ii("CAP03702",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03702")
        ig,ig_,_ = s2mpj_ii("CAP03802",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03802")
        ig,ig_,_ = s2mpj_ii("CAP03902",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03902")
        ig,ig_,_ = s2mpj_ii("CAP04002",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04002")
        ig,ig_,_ = s2mpj_ii("CAP04102",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04102")
        ig,ig_,_ = s2mpj_ii("CAP04202",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04202")
        ig,ig_,_ = s2mpj_ii("CAP04302",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04302")
        ig,ig_,_ = s2mpj_ii("CAP04402",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04402")
        ig,ig_,_ = s2mpj_ii("CAP04502",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04502")
        ig,ig_,_ = s2mpj_ii("CAP04602",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04602")
        ig,ig_,_ = s2mpj_ii("CAP04702",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04702")
        ig,ig_,_ = s2mpj_ii("CAP04802",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04802")
        ig,ig_,_ = s2mpj_ii("CAP04902",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04902")
        ig,ig_,_ = s2mpj_ii("CAP05002",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05002")
        ig,ig_,_ = s2mpj_ii("CAP05102",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05102")
        ig,ig_,_ = s2mpj_ii("CAP05202",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05202")
        ig,ig_,_ = s2mpj_ii("CAP05302",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05302")
        ig,ig_,_ = s2mpj_ii("CAP05402",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05402")
        ig,ig_,_ = s2mpj_ii("CAP05502",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05502")
        ig,ig_,_ = s2mpj_ii("CAP05602",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05602")
        ig,ig_,_ = s2mpj_ii("CAP05702",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05702")
        ig,ig_,_ = s2mpj_ii("CAP05802",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05802")
        ig,ig_,_ = s2mpj_ii("CAP05902",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05902")
        ig,ig_,_ = s2mpj_ii("CAP06002",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP06002")
        ig,ig_,_ = s2mpj_ii("CAP06102",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP06102")
        ig,ig_,_ = s2mpj_ii("CAP06202",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP06202")
        ig,ig_,_ = s2mpj_ii("CAP06302",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP06302")
        ig,ig_,_ = s2mpj_ii("CAP06402",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP06402")
        ig,ig_,_ = s2mpj_ii("CAP06502",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP06502")
        ig,ig_,_ = s2mpj_ii("CAP00103",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00103")
        ig,ig_,_ = s2mpj_ii("CAP00203",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00203")
        ig,ig_,_ = s2mpj_ii("CAP00303",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00303")
        ig,ig_,_ = s2mpj_ii("CAP00403",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00403")
        ig,ig_,_ = s2mpj_ii("CAP00503",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00503")
        ig,ig_,_ = s2mpj_ii("CAP00603",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00603")
        ig,ig_,_ = s2mpj_ii("CAP00703",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00703")
        ig,ig_,_ = s2mpj_ii("CAP00803",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00803")
        ig,ig_,_ = s2mpj_ii("CAP00903",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00903")
        ig,ig_,_ = s2mpj_ii("CAP01003",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01003")
        ig,ig_,_ = s2mpj_ii("CAP01103",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01103")
        ig,ig_,_ = s2mpj_ii("CAP01203",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01203")
        ig,ig_,_ = s2mpj_ii("CAP01303",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01303")
        ig,ig_,_ = s2mpj_ii("CAP01403",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01403")
        ig,ig_,_ = s2mpj_ii("CAP01503",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01503")
        ig,ig_,_ = s2mpj_ii("CAP01603",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01603")
        ig,ig_,_ = s2mpj_ii("CAP01703",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01703")
        ig,ig_,_ = s2mpj_ii("CAP01803",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01803")
        ig,ig_,_ = s2mpj_ii("CAP01903",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01903")
        ig,ig_,_ = s2mpj_ii("CAP02003",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02003")
        ig,ig_,_ = s2mpj_ii("CAP02103",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02103")
        ig,ig_,_ = s2mpj_ii("CAP02203",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02203")
        ig,ig_,_ = s2mpj_ii("CAP02303",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02303")
        ig,ig_,_ = s2mpj_ii("CAP02403",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02403")
        ig,ig_,_ = s2mpj_ii("CAP02503",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02503")
        ig,ig_,_ = s2mpj_ii("CAP02603",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02603")
        ig,ig_,_ = s2mpj_ii("CAP02703",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02703")
        ig,ig_,_ = s2mpj_ii("CAP02803",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02803")
        ig,ig_,_ = s2mpj_ii("CAP02903",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02903")
        ig,ig_,_ = s2mpj_ii("CAP03003",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03003")
        ig,ig_,_ = s2mpj_ii("CAP03103",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03103")
        ig,ig_,_ = s2mpj_ii("CAP03203",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03203")
        ig,ig_,_ = s2mpj_ii("CAP03303",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03303")
        ig,ig_,_ = s2mpj_ii("CAP03403",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03403")
        ig,ig_,_ = s2mpj_ii("CAP03503",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03503")
        ig,ig_,_ = s2mpj_ii("CAP03603",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03603")
        ig,ig_,_ = s2mpj_ii("CAP03703",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03703")
        ig,ig_,_ = s2mpj_ii("CAP03803",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03803")
        ig,ig_,_ = s2mpj_ii("CAP03903",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03903")
        ig,ig_,_ = s2mpj_ii("CAP04003",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04003")
        ig,ig_,_ = s2mpj_ii("CAP04103",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04103")
        ig,ig_,_ = s2mpj_ii("CAP04203",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04203")
        ig,ig_,_ = s2mpj_ii("CAP04303",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04303")
        ig,ig_,_ = s2mpj_ii("CAP04403",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04403")
        ig,ig_,_ = s2mpj_ii("CAP04503",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04503")
        ig,ig_,_ = s2mpj_ii("CAP04603",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04603")
        ig,ig_,_ = s2mpj_ii("CAP04703",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04703")
        ig,ig_,_ = s2mpj_ii("CAP04803",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04803")
        ig,ig_,_ = s2mpj_ii("CAP04903",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04903")
        ig,ig_,_ = s2mpj_ii("CAP05003",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05003")
        ig,ig_,_ = s2mpj_ii("CAP05103",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05103")
        ig,ig_,_ = s2mpj_ii("CAP05203",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05203")
        ig,ig_,_ = s2mpj_ii("CAP05303",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05303")
        ig,ig_,_ = s2mpj_ii("CAP05403",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05403")
        ig,ig_,_ = s2mpj_ii("CAP05503",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05503")
        ig,ig_,_ = s2mpj_ii("CAP05603",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05603")
        ig,ig_,_ = s2mpj_ii("CAP05703",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05703")
        ig,ig_,_ = s2mpj_ii("CAP05803",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05803")
        ig,ig_,_ = s2mpj_ii("CAP05903",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05903")
        ig,ig_,_ = s2mpj_ii("CAP06003",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP06003")
        ig,ig_,_ = s2mpj_ii("CAP06103",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP06103")
        ig,ig_,_ = s2mpj_ii("CAP06203",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP06203")
        ig,ig_,_ = s2mpj_ii("CAP06303",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP06303")
        ig,ig_,_ = s2mpj_ii("CAP06403",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP06403")
        ig,ig_,_ = s2mpj_ii("CAP06503",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP06503")
        ig,ig_,_ = s2mpj_ii("CAP00104",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00104")
        ig,ig_,_ = s2mpj_ii("CAP00204",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00204")
        ig,ig_,_ = s2mpj_ii("CAP00304",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00304")
        ig,ig_,_ = s2mpj_ii("CAP00404",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00404")
        ig,ig_,_ = s2mpj_ii("CAP00504",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00504")
        ig,ig_,_ = s2mpj_ii("CAP00604",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00604")
        ig,ig_,_ = s2mpj_ii("CAP00704",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00704")
        ig,ig_,_ = s2mpj_ii("CAP00804",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00804")
        ig,ig_,_ = s2mpj_ii("CAP00904",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00904")
        ig,ig_,_ = s2mpj_ii("CAP01004",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01004")
        ig,ig_,_ = s2mpj_ii("CAP01104",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01104")
        ig,ig_,_ = s2mpj_ii("CAP01204",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01204")
        ig,ig_,_ = s2mpj_ii("CAP01304",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01304")
        ig,ig_,_ = s2mpj_ii("CAP01404",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01404")
        ig,ig_,_ = s2mpj_ii("CAP01504",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01504")
        ig,ig_,_ = s2mpj_ii("CAP01604",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01604")
        ig,ig_,_ = s2mpj_ii("CAP01704",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01704")
        ig,ig_,_ = s2mpj_ii("CAP01804",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01804")
        ig,ig_,_ = s2mpj_ii("CAP01904",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01904")
        ig,ig_,_ = s2mpj_ii("CAP02004",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02004")
        ig,ig_,_ = s2mpj_ii("CAP02104",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02104")
        ig,ig_,_ = s2mpj_ii("CAP02204",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02204")
        ig,ig_,_ = s2mpj_ii("CAP02304",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02304")
        ig,ig_,_ = s2mpj_ii("CAP02404",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02404")
        ig,ig_,_ = s2mpj_ii("CAP02504",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02504")
        ig,ig_,_ = s2mpj_ii("CAP02604",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02604")
        ig,ig_,_ = s2mpj_ii("CAP02704",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02704")
        ig,ig_,_ = s2mpj_ii("CAP02804",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02804")
        ig,ig_,_ = s2mpj_ii("CAP02904",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02904")
        ig,ig_,_ = s2mpj_ii("CAP03004",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03004")
        ig,ig_,_ = s2mpj_ii("CAP03104",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03104")
        ig,ig_,_ = s2mpj_ii("CAP03204",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03204")
        ig,ig_,_ = s2mpj_ii("CAP03304",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03304")
        ig,ig_,_ = s2mpj_ii("CAP03404",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03404")
        ig,ig_,_ = s2mpj_ii("CAP03504",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03504")
        ig,ig_,_ = s2mpj_ii("CAP03604",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03604")
        ig,ig_,_ = s2mpj_ii("CAP03704",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03704")
        ig,ig_,_ = s2mpj_ii("CAP03804",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03804")
        ig,ig_,_ = s2mpj_ii("CAP03904",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03904")
        ig,ig_,_ = s2mpj_ii("CAP04004",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04004")
        ig,ig_,_ = s2mpj_ii("CAP04104",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04104")
        ig,ig_,_ = s2mpj_ii("CAP04204",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04204")
        ig,ig_,_ = s2mpj_ii("CAP04304",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04304")
        ig,ig_,_ = s2mpj_ii("CAP04404",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04404")
        ig,ig_,_ = s2mpj_ii("CAP04504",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04504")
        ig,ig_,_ = s2mpj_ii("CAP04604",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04604")
        ig,ig_,_ = s2mpj_ii("CAP04704",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04704")
        ig,ig_,_ = s2mpj_ii("CAP04804",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04804")
        ig,ig_,_ = s2mpj_ii("CAP04904",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04904")
        ig,ig_,_ = s2mpj_ii("CAP05004",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05004")
        ig,ig_,_ = s2mpj_ii("CAP05104",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05104")
        ig,ig_,_ = s2mpj_ii("CAP05204",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05204")
        ig,ig_,_ = s2mpj_ii("CAP05304",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05304")
        ig,ig_,_ = s2mpj_ii("CAP05404",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05404")
        ig,ig_,_ = s2mpj_ii("CAP05504",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05504")
        ig,ig_,_ = s2mpj_ii("CAP05604",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05604")
        ig,ig_,_ = s2mpj_ii("CAP05704",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05704")
        ig,ig_,_ = s2mpj_ii("CAP05804",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05804")
        ig,ig_,_ = s2mpj_ii("CAP05904",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05904")
        ig,ig_,_ = s2mpj_ii("CAP06004",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP06004")
        ig,ig_,_ = s2mpj_ii("CAP06104",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP06104")
        ig,ig_,_ = s2mpj_ii("CAP06204",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP06204")
        ig,ig_,_ = s2mpj_ii("CAP06304",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP06304")
        ig,ig_,_ = s2mpj_ii("CAP06404",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP06404")
        ig,ig_,_ = s2mpj_ii("CAP06504",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP06504")
        ig,ig_,_ = s2mpj_ii("CAP00105",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00105")
        ig,ig_,_ = s2mpj_ii("CAP00205",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00205")
        ig,ig_,_ = s2mpj_ii("CAP00305",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00305")
        ig,ig_,_ = s2mpj_ii("CAP00405",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00405")
        ig,ig_,_ = s2mpj_ii("CAP00505",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00505")
        ig,ig_,_ = s2mpj_ii("CAP00605",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00605")
        ig,ig_,_ = s2mpj_ii("CAP00705",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00705")
        ig,ig_,_ = s2mpj_ii("CAP00805",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00805")
        ig,ig_,_ = s2mpj_ii("CAP00905",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00905")
        ig,ig_,_ = s2mpj_ii("CAP01005",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01005")
        ig,ig_,_ = s2mpj_ii("CAP01105",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01105")
        ig,ig_,_ = s2mpj_ii("CAP01205",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01205")
        ig,ig_,_ = s2mpj_ii("CAP01305",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01305")
        ig,ig_,_ = s2mpj_ii("CAP01405",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01405")
        ig,ig_,_ = s2mpj_ii("CAP01505",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01505")
        ig,ig_,_ = s2mpj_ii("CAP01605",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01605")
        ig,ig_,_ = s2mpj_ii("CAP01705",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01705")
        ig,ig_,_ = s2mpj_ii("CAP01805",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01805")
        ig,ig_,_ = s2mpj_ii("CAP01905",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01905")
        ig,ig_,_ = s2mpj_ii("CAP02005",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02005")
        ig,ig_,_ = s2mpj_ii("CAP02105",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02105")
        ig,ig_,_ = s2mpj_ii("CAP02205",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02205")
        ig,ig_,_ = s2mpj_ii("CAP02305",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02305")
        ig,ig_,_ = s2mpj_ii("CAP02405",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02405")
        ig,ig_,_ = s2mpj_ii("CAP02505",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02505")
        ig,ig_,_ = s2mpj_ii("CAP02605",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02605")
        ig,ig_,_ = s2mpj_ii("CAP02705",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02705")
        ig,ig_,_ = s2mpj_ii("CAP02805",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02805")
        ig,ig_,_ = s2mpj_ii("CAP02905",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02905")
        ig,ig_,_ = s2mpj_ii("CAP03005",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03005")
        ig,ig_,_ = s2mpj_ii("CAP03105",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03105")
        ig,ig_,_ = s2mpj_ii("CAP03205",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03205")
        ig,ig_,_ = s2mpj_ii("CAP03305",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03305")
        ig,ig_,_ = s2mpj_ii("CAP03405",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03405")
        ig,ig_,_ = s2mpj_ii("CAP03505",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03505")
        ig,ig_,_ = s2mpj_ii("CAP03605",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03605")
        ig,ig_,_ = s2mpj_ii("CAP03705",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03705")
        ig,ig_,_ = s2mpj_ii("CAP03805",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03805")
        ig,ig_,_ = s2mpj_ii("CAP03905",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03905")
        ig,ig_,_ = s2mpj_ii("CAP04005",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04005")
        ig,ig_,_ = s2mpj_ii("CAP04105",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04105")
        ig,ig_,_ = s2mpj_ii("CAP04205",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04205")
        ig,ig_,_ = s2mpj_ii("CAP04305",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04305")
        ig,ig_,_ = s2mpj_ii("CAP04405",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04405")
        ig,ig_,_ = s2mpj_ii("CAP04505",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04505")
        ig,ig_,_ = s2mpj_ii("CAP04605",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04605")
        ig,ig_,_ = s2mpj_ii("CAP04705",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04705")
        ig,ig_,_ = s2mpj_ii("CAP04805",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04805")
        ig,ig_,_ = s2mpj_ii("CAP04905",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04905")
        ig,ig_,_ = s2mpj_ii("CAP05005",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05005")
        ig,ig_,_ = s2mpj_ii("CAP05105",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05105")
        ig,ig_,_ = s2mpj_ii("CAP05205",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05205")
        ig,ig_,_ = s2mpj_ii("CAP05305",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05305")
        ig,ig_,_ = s2mpj_ii("CAP05405",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05405")
        ig,ig_,_ = s2mpj_ii("CAP05505",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05505")
        ig,ig_,_ = s2mpj_ii("CAP05605",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05605")
        ig,ig_,_ = s2mpj_ii("CAP05705",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05705")
        ig,ig_,_ = s2mpj_ii("CAP05805",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05805")
        ig,ig_,_ = s2mpj_ii("CAP05905",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05905")
        ig,ig_,_ = s2mpj_ii("CAP06005",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP06005")
        ig,ig_,_ = s2mpj_ii("CAP06105",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP06105")
        ig,ig_,_ = s2mpj_ii("CAP06205",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP06205")
        ig,ig_,_ = s2mpj_ii("CAP06305",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP06305")
        ig,ig_,_ = s2mpj_ii("CAP06405",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP06405")
        ig,ig_,_ = s2mpj_ii("CAP06505",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP06505")
        ig,ig_,_ = s2mpj_ii("CAP00106",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00106")
        ig,ig_,_ = s2mpj_ii("CAP00206",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00206")
        ig,ig_,_ = s2mpj_ii("CAP00306",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00306")
        ig,ig_,_ = s2mpj_ii("CAP00406",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00406")
        ig,ig_,_ = s2mpj_ii("CAP00506",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00506")
        ig,ig_,_ = s2mpj_ii("CAP00606",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00606")
        ig,ig_,_ = s2mpj_ii("CAP00706",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00706")
        ig,ig_,_ = s2mpj_ii("CAP00806",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00806")
        ig,ig_,_ = s2mpj_ii("CAP00906",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP00906")
        ig,ig_,_ = s2mpj_ii("CAP01006",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01006")
        ig,ig_,_ = s2mpj_ii("CAP01106",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01106")
        ig,ig_,_ = s2mpj_ii("CAP01206",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01206")
        ig,ig_,_ = s2mpj_ii("CAP01306",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01306")
        ig,ig_,_ = s2mpj_ii("CAP01406",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01406")
        ig,ig_,_ = s2mpj_ii("CAP01506",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01506")
        ig,ig_,_ = s2mpj_ii("CAP01606",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01606")
        ig,ig_,_ = s2mpj_ii("CAP01706",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01706")
        ig,ig_,_ = s2mpj_ii("CAP01806",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01806")
        ig,ig_,_ = s2mpj_ii("CAP01906",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP01906")
        ig,ig_,_ = s2mpj_ii("CAP02006",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02006")
        ig,ig_,_ = s2mpj_ii("CAP02106",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02106")
        ig,ig_,_ = s2mpj_ii("CAP02206",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02206")
        ig,ig_,_ = s2mpj_ii("CAP02306",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02306")
        ig,ig_,_ = s2mpj_ii("CAP02406",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02406")
        ig,ig_,_ = s2mpj_ii("CAP02506",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02506")
        ig,ig_,_ = s2mpj_ii("CAP02606",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02606")
        ig,ig_,_ = s2mpj_ii("CAP02706",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02706")
        ig,ig_,_ = s2mpj_ii("CAP02806",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02806")
        ig,ig_,_ = s2mpj_ii("CAP02906",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP02906")
        ig,ig_,_ = s2mpj_ii("CAP03006",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03006")
        ig,ig_,_ = s2mpj_ii("CAP03106",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03106")
        ig,ig_,_ = s2mpj_ii("CAP03206",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03206")
        ig,ig_,_ = s2mpj_ii("CAP03306",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03306")
        ig,ig_,_ = s2mpj_ii("CAP03406",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03406")
        ig,ig_,_ = s2mpj_ii("CAP03506",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03506")
        ig,ig_,_ = s2mpj_ii("CAP03606",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03606")
        ig,ig_,_ = s2mpj_ii("CAP03706",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03706")
        ig,ig_,_ = s2mpj_ii("CAP03806",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03806")
        ig,ig_,_ = s2mpj_ii("CAP03906",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP03906")
        ig,ig_,_ = s2mpj_ii("CAP04006",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04006")
        ig,ig_,_ = s2mpj_ii("CAP04106",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04106")
        ig,ig_,_ = s2mpj_ii("CAP04206",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04206")
        ig,ig_,_ = s2mpj_ii("CAP04306",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04306")
        ig,ig_,_ = s2mpj_ii("CAP04406",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04406")
        ig,ig_,_ = s2mpj_ii("CAP04506",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04506")
        ig,ig_,_ = s2mpj_ii("CAP04606",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04606")
        ig,ig_,_ = s2mpj_ii("CAP04706",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04706")
        ig,ig_,_ = s2mpj_ii("CAP04806",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04806")
        ig,ig_,_ = s2mpj_ii("CAP04906",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP04906")
        ig,ig_,_ = s2mpj_ii("CAP05006",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05006")
        ig,ig_,_ = s2mpj_ii("CAP05106",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05106")
        ig,ig_,_ = s2mpj_ii("CAP05206",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05206")
        ig,ig_,_ = s2mpj_ii("CAP05306",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05306")
        ig,ig_,_ = s2mpj_ii("CAP05406",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05406")
        ig,ig_,_ = s2mpj_ii("CAP05506",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05506")
        ig,ig_,_ = s2mpj_ii("CAP05606",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05606")
        ig,ig_,_ = s2mpj_ii("CAP05706",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05706")
        ig,ig_,_ = s2mpj_ii("CAP05806",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05806")
        ig,ig_,_ = s2mpj_ii("CAP05906",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP05906")
        ig,ig_,_ = s2mpj_ii("CAP06006",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP06006")
        ig,ig_,_ = s2mpj_ii("CAP06106",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP06106")
        ig,ig_,_ = s2mpj_ii("CAP06206",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP06206")
        ig,ig_,_ = s2mpj_ii("CAP06306",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP06306")
        ig,ig_,_ = s2mpj_ii("CAP06406",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP06406")
        ig,ig_,_ = s2mpj_ii("CAP06506",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CAP06506")
        ig,ig_,_ = s2mpj_ii("MND00102",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND00102")
        ig,ig_,_ = s2mpj_ii("MXD00102",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD00102")
        ig,ig_,_ = s2mpj_ii("MND00202",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND00202")
        ig,ig_,_ = s2mpj_ii("MXD00202",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD00202")
        ig,ig_,_ = s2mpj_ii("MND00502",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND00502")
        ig,ig_,_ = s2mpj_ii("MXD00502",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD00502")
        ig,ig_,_ = s2mpj_ii("MND00702",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND00702")
        ig,ig_,_ = s2mpj_ii("MXD00702",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD00702")
        ig,ig_,_ = s2mpj_ii("MND00802",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND00802")
        ig,ig_,_ = s2mpj_ii("MXD00802",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD00802")
        ig,ig_,_ = s2mpj_ii("MND00902",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND00902")
        ig,ig_,_ = s2mpj_ii("MXD00902",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD00902")
        ig,ig_,_ = s2mpj_ii("MND01002",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND01002")
        ig,ig_,_ = s2mpj_ii("MXD01002",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD01002")
        ig,ig_,_ = s2mpj_ii("MND00103",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND00103")
        ig,ig_,_ = s2mpj_ii("MXD00103",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD00103")
        ig,ig_,_ = s2mpj_ii("MND00203",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND00203")
        ig,ig_,_ = s2mpj_ii("MXD00203",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD00203")
        ig,ig_,_ = s2mpj_ii("MND00303",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND00303")
        ig,ig_,_ = s2mpj_ii("MXD00303",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD00303")
        ig,ig_,_ = s2mpj_ii("MND00403",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND00403")
        ig,ig_,_ = s2mpj_ii("MXD00403",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD00403")
        ig,ig_,_ = s2mpj_ii("MND00503",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND00503")
        ig,ig_,_ = s2mpj_ii("MXD00503",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD00503")
        ig,ig_,_ = s2mpj_ii("MND00603",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND00603")
        ig,ig_,_ = s2mpj_ii("MXD00603",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD00603")
        ig,ig_,_ = s2mpj_ii("MND00703",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND00703")
        ig,ig_,_ = s2mpj_ii("MXD00703",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD00703")
        ig,ig_,_ = s2mpj_ii("MND00803",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND00803")
        ig,ig_,_ = s2mpj_ii("MXD00803",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD00803")
        ig,ig_,_ = s2mpj_ii("MND00903",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND00903")
        ig,ig_,_ = s2mpj_ii("MXD00903",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD00903")
        ig,ig_,_ = s2mpj_ii("MND01003",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND01003")
        ig,ig_,_ = s2mpj_ii("MXD01003",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD01003")
        ig,ig_,_ = s2mpj_ii("MND00104",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND00104")
        ig,ig_,_ = s2mpj_ii("MXD00104",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD00104")
        ig,ig_,_ = s2mpj_ii("MND00204",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND00204")
        ig,ig_,_ = s2mpj_ii("MXD00204",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD00204")
        ig,ig_,_ = s2mpj_ii("MND00304",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND00304")
        ig,ig_,_ = s2mpj_ii("MXD00304",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD00304")
        ig,ig_,_ = s2mpj_ii("MND00404",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND00404")
        ig,ig_,_ = s2mpj_ii("MXD00404",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD00404")
        ig,ig_,_ = s2mpj_ii("MND00504",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND00504")
        ig,ig_,_ = s2mpj_ii("MXD00504",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD00504")
        ig,ig_,_ = s2mpj_ii("MND00604",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND00604")
        ig,ig_,_ = s2mpj_ii("MXD00604",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD00604")
        ig,ig_,_ = s2mpj_ii("MND00704",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND00704")
        ig,ig_,_ = s2mpj_ii("MXD00704",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD00704")
        ig,ig_,_ = s2mpj_ii("MND00804",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND00804")
        ig,ig_,_ = s2mpj_ii("MXD00804",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD00804")
        ig,ig_,_ = s2mpj_ii("MND00904",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND00904")
        ig,ig_,_ = s2mpj_ii("MXD00904",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD00904")
        ig,ig_,_ = s2mpj_ii("MND01004",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND01004")
        ig,ig_,_ = s2mpj_ii("MXD01004",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD01004")
        ig,ig_,_ = s2mpj_ii("MND00105",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND00105")
        ig,ig_,_ = s2mpj_ii("MXD00105",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD00105")
        ig,ig_,_ = s2mpj_ii("MND00205",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND00205")
        ig,ig_,_ = s2mpj_ii("MXD00205",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD00205")
        ig,ig_,_ = s2mpj_ii("MND00305",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND00305")
        ig,ig_,_ = s2mpj_ii("MXD00305",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD00305")
        ig,ig_,_ = s2mpj_ii("MND00405",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND00405")
        ig,ig_,_ = s2mpj_ii("MXD00405",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD00405")
        ig,ig_,_ = s2mpj_ii("MND00505",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND00505")
        ig,ig_,_ = s2mpj_ii("MXD00505",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD00505")
        ig,ig_,_ = s2mpj_ii("MND00605",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND00605")
        ig,ig_,_ = s2mpj_ii("MXD00605",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD00605")
        ig,ig_,_ = s2mpj_ii("MND00705",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND00705")
        ig,ig_,_ = s2mpj_ii("MXD00705",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD00705")
        ig,ig_,_ = s2mpj_ii("MND00805",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND00805")
        ig,ig_,_ = s2mpj_ii("MXD00805",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD00805")
        ig,ig_,_ = s2mpj_ii("MND00905",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND00905")
        ig,ig_,_ = s2mpj_ii("MXD00905",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD00905")
        ig,ig_,_ = s2mpj_ii("MND01005",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND01005")
        ig,ig_,_ = s2mpj_ii("MXD01005",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD01005")
        ig,ig_,_ = s2mpj_ii("MND00106",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND00106")
        ig,ig_,_ = s2mpj_ii("MXD00106",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD00106")
        ig,ig_,_ = s2mpj_ii("MND00206",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND00206")
        ig,ig_,_ = s2mpj_ii("MXD00206",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD00206")
        ig,ig_,_ = s2mpj_ii("MND00306",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND00306")
        ig,ig_,_ = s2mpj_ii("MXD00306",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD00306")
        ig,ig_,_ = s2mpj_ii("MND00406",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND00406")
        ig,ig_,_ = s2mpj_ii("MXD00406",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD00406")
        ig,ig_,_ = s2mpj_ii("MND00506",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND00506")
        ig,ig_,_ = s2mpj_ii("MXD00506",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD00506")
        ig,ig_,_ = s2mpj_ii("MND00606",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND00606")
        ig,ig_,_ = s2mpj_ii("MXD00606",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD00606")
        ig,ig_,_ = s2mpj_ii("MND00706",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND00706")
        ig,ig_,_ = s2mpj_ii("MXD00706",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD00706")
        ig,ig_,_ = s2mpj_ii("MND00806",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND00806")
        ig,ig_,_ = s2mpj_ii("MXD00806",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD00806")
        ig,ig_,_ = s2mpj_ii("MND00906",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND00906")
        ig,ig_,_ = s2mpj_ii("MXD00906",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD00906")
        ig,ig_,_ = s2mpj_ii("MND01006",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"MND01006")
        ig,ig_,_ = s2mpj_ii("MXD01006",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MXD01006")
        ig,ig_,_ = s2mpj_ii("INV00101",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"INV00101")
        ig,ig_,_ = s2mpj_ii("INV00201",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"INV00201")
        ig,ig_,_ = s2mpj_ii("INV00301",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"INV00301")
        ig,ig_,_ = s2mpj_ii("INV00401",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"INV00401")
        ig,ig_,_ = s2mpj_ii("INV00501",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"INV00501")
        ig,ig_,_ = s2mpj_ii("INV00601",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"INV00601")
        ig,ig_,_ = s2mpj_ii("INV00102",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"INV00102")
        ig,ig_,_ = s2mpj_ii("INV00202",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"INV00202")
        ig,ig_,_ = s2mpj_ii("INV00302",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"INV00302")
        ig,ig_,_ = s2mpj_ii("INV00402",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"INV00402")
        ig,ig_,_ = s2mpj_ii("INV00502",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"INV00502")
        ig,ig_,_ = s2mpj_ii("INV00602",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"INV00602")
        ig,ig_,_ = s2mpj_ii("INV00103",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"INV00103")
        ig,ig_,_ = s2mpj_ii("INV00203",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"INV00203")
        ig,ig_,_ = s2mpj_ii("INV00303",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"INV00303")
        ig,ig_,_ = s2mpj_ii("INV00403",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"INV00403")
        ig,ig_,_ = s2mpj_ii("INV00503",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"INV00503")
        ig,ig_,_ = s2mpj_ii("INV00603",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"INV00603")
        ig,ig_,_ = s2mpj_ii("INV00104",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"INV00104")
        ig,ig_,_ = s2mpj_ii("INV00204",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"INV00204")
        ig,ig_,_ = s2mpj_ii("INV00304",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"INV00304")
        ig,ig_,_ = s2mpj_ii("INV00404",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"INV00404")
        ig,ig_,_ = s2mpj_ii("INV00504",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"INV00504")
        ig,ig_,_ = s2mpj_ii("INV00604",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"INV00604")
        ig,ig_,_ = s2mpj_ii("INV00105",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"INV00105")
        ig,ig_,_ = s2mpj_ii("INV00205",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"INV00205")
        ig,ig_,_ = s2mpj_ii("INV00305",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"INV00305")
        ig,ig_,_ = s2mpj_ii("INV00405",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"INV00405")
        ig,ig_,_ = s2mpj_ii("INV00505",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"INV00505")
        ig,ig_,_ = s2mpj_ii("INV00605",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"INV00605")
        ig,ig_,_ = s2mpj_ii("INV00106",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"INV00106")
        ig,ig_,_ = s2mpj_ii("INV00206",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"INV00206")
        ig,ig_,_ = s2mpj_ii("INV00306",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"INV00306")
        ig,ig_,_ = s2mpj_ii("INV00406",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"INV00406")
        ig,ig_,_ = s2mpj_ii("INV00506",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"INV00506")
        ig,ig_,_ = s2mpj_ii("INV00606",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"INV00606")
        ig,ig_,_ = s2mpj_ii("OBJECTIV",ig_)
        arrset(gtype,ig,"<>")
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        ngrp   = length(ig_)
        iv,ix_,_ = s2mpj_ii("Y00102",ix_)
        arrset(pb.xnames,iv,"Y00102")
        ig = ig_["MND00104"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP01401"]
        pbm.A[ig,iv] += Float64(.01773)
        iv,ix_,_ = s2mpj_ii("Y00102",ix_)
        arrset(pb.xnames,iv,"Y00102")
        ig = ig_["CAP01501"]
        pbm.A[ig,iv] += Float64(.01775)
        ig = ig_["CAP06302"]
        pbm.A[ig,iv] += Float64(.00132)
        iv,ix_,_ = s2mpj_ii("Y00102",ix_)
        arrset(pb.xnames,iv,"Y00102")
        ig = ig_["MND00103"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["INV00302"]
        pbm.A[ig,iv] += Float64(-.76829)
        iv,ix_,_ = s2mpj_ii("Y00102",ix_)
        arrset(pb.xnames,iv,"Y00102")
        ig = ig_["MXD00104"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["INV00301"]
        pbm.A[ig,iv] += Float64(-1.1405)
        iv,ix_,_ = s2mpj_ii("Y00102",ix_)
        arrset(pb.xnames,iv,"Y00102")
        ig = ig_["MXD00103"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP01502"]
        pbm.A[ig,iv] += Float64(.01248)
        iv,ix_,_ = s2mpj_ii("Y00102",ix_)
        arrset(pb.xnames,iv,"Y00102")
        ig = ig_["MXD00106"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP01402"]
        pbm.A[ig,iv] += Float64(.01301)
        iv,ix_,_ = s2mpj_ii("Y00102",ix_)
        arrset(pb.xnames,iv,"Y00102")
        ig = ig_["CAP06301"]
        pbm.A[ig,iv] += Float64(.00141)
        ig = ig_["MND00106"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00102",ix_)
        arrset(pb.xnames,iv,"Y00102")
        ig = ig_["MXD00102"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00102"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00102",ix_)
        arrset(pb.xnames,iv,"Y00102")
        ig = ig_["MND00105"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00105"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00102",ix_)
        arrset(pb.xnames,iv,"Y00102")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-31.09)
        iv,ix_,_ = s2mpj_ii("Y00103",ix_)
        arrset(pb.xnames,iv,"Y00103")
        ig = ig_["INV00302"]
        pbm.A[ig,iv] += Float64(-1.05277)
        ig = ig_["CAP06302"]
        pbm.A[ig,iv] += Float64(.0013)
        iv,ix_,_ = s2mpj_ii("Y00103",ix_)
        arrset(pb.xnames,iv,"Y00103")
        ig = ig_["MXD00105"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-29.52)
        iv,ix_,_ = s2mpj_ii("Y00103",ix_)
        arrset(pb.xnames,iv,"Y00103")
        ig = ig_["MND00104"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP01503"]
        pbm.A[ig,iv] += Float64(.01385)
        iv,ix_,_ = s2mpj_ii("Y00103",ix_)
        arrset(pb.xnames,iv,"Y00103")
        ig = ig_["CAP06303"]
        pbm.A[ig,iv] += Float64(.00143)
        ig = ig_["INV00303"]
        pbm.A[ig,iv] += Float64(-.85602)
        iv,ix_,_ = s2mpj_ii("Y00103",ix_)
        arrset(pb.xnames,iv,"Y00103")
        ig = ig_["MND00105"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP01402"]
        pbm.A[ig,iv] += Float64(.01636)
        iv,ix_,_ = s2mpj_ii("Y00103",ix_)
        arrset(pb.xnames,iv,"Y00103")
        ig = ig_["MXD00103"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP01403"]
        pbm.A[ig,iv] += Float64(.01438)
        iv,ix_,_ = s2mpj_ii("Y00103",ix_)
        arrset(pb.xnames,iv,"Y00103")
        ig = ig_["MND00103"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00104"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00103",ix_)
        arrset(pb.xnames,iv,"Y00103")
        ig = ig_["CAP01502"]
        pbm.A[ig,iv] += Float64(.01639)
        ig = ig_["MXD00106"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00103",ix_)
        arrset(pb.xnames,iv,"Y00103")
        ig = ig_["MND00106"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00104",ix_)
        arrset(pb.xnames,iv,"Y00104")
        ig = ig_["MXD00104"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00106"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00104",ix_)
        arrset(pb.xnames,iv,"Y00104")
        ig = ig_["MND00104"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP01504"]
        pbm.A[ig,iv] += Float64(.01385)
        iv,ix_,_ = s2mpj_ii("Y00104",ix_)
        arrset(pb.xnames,iv,"Y00104")
        ig = ig_["MXD00106"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00105"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00104",ix_)
        arrset(pb.xnames,iv,"Y00104")
        ig = ig_["CAP01404"]
        pbm.A[ig,iv] += Float64(.01438)
        ig = ig_["MND00105"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00104",ix_)
        arrset(pb.xnames,iv,"Y00104")
        ig = ig_["CAP01403"]
        pbm.A[ig,iv] += Float64(.01636)
        ig = ig_["CAP01503"]
        pbm.A[ig,iv] += Float64(.01639)
        iv,ix_,_ = s2mpj_ii("Y00104",ix_)
        arrset(pb.xnames,iv,"Y00104")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-28.03)
        ig = ig_["CAP06304"]
        pbm.A[ig,iv] += Float64(.00143)
        iv,ix_,_ = s2mpj_ii("Y00104",ix_)
        arrset(pb.xnames,iv,"Y00104")
        ig = ig_["INV00304"]
        pbm.A[ig,iv] += Float64(-.85602)
        ig = ig_["CAP06303"]
        pbm.A[ig,iv] += Float64(.0013)
        iv,ix_,_ = s2mpj_ii("Y00104",ix_)
        arrset(pb.xnames,iv,"Y00104")
        ig = ig_["INV00303"]
        pbm.A[ig,iv] += Float64(-1.05277)
        iv,ix_,_ = s2mpj_ii("Y00105",ix_)
        arrset(pb.xnames,iv,"Y00105")
        ig = ig_["MXD00106"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00106"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00105",ix_)
        arrset(pb.xnames,iv,"Y00105")
        ig = ig_["CAP01504"]
        pbm.A[ig,iv] += Float64(.01578)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-26.61)
        iv,ix_,_ = s2mpj_ii("Y00105",ix_)
        arrset(pb.xnames,iv,"Y00105")
        ig = ig_["CAP01404"]
        pbm.A[ig,iv] += Float64(.01576)
        ig = ig_["MND00105"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00105",ix_)
        arrset(pb.xnames,iv,"Y00105")
        ig = ig_["INV00305"]
        pbm.A[ig,iv] += Float64(-.89501)
        ig = ig_["INV00304"]
        pbm.A[ig,iv] += Float64(-1.01378)
        iv,ix_,_ = s2mpj_ii("Y00105",ix_)
        arrset(pb.xnames,iv,"Y00105")
        ig = ig_["CAP06305"]
        pbm.A[ig,iv] += Float64(.00147)
        ig = ig_["CAP06304"]
        pbm.A[ig,iv] += Float64(.00125)
        iv,ix_,_ = s2mpj_ii("Y00105",ix_)
        arrset(pb.xnames,iv,"Y00105")
        ig = ig_["MXD00105"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP01405"]
        pbm.A[ig,iv] += Float64(.01498)
        iv,ix_,_ = s2mpj_ii("Y00105",ix_)
        arrset(pb.xnames,iv,"Y00105")
        ig = ig_["CAP01505"]
        pbm.A[ig,iv] += Float64(.01446)
        iv,ix_,_ = s2mpj_ii("Y00106",ix_)
        arrset(pb.xnames,iv,"Y00106")
        ig = ig_["CAP01506"]
        pbm.A[ig,iv] += Float64(.03024)
        ig = ig_["CAP06306"]
        pbm.A[ig,iv] += Float64(.00272)
        iv,ix_,_ = s2mpj_ii("Y00106",ix_)
        arrset(pb.xnames,iv,"Y00106")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-25.26)
        ig = ig_["INV00306"]
        pbm.A[ig,iv] += Float64(-1.90879)
        iv,ix_,_ = s2mpj_ii("Y00106",ix_)
        arrset(pb.xnames,iv,"Y00106")
        ig = ig_["INV00305"]
        pbm.A[ig,iv] += Float64(-1.09488)
        ig = ig_["CAP01405"]
        pbm.A[ig,iv] += Float64(.01702)
        iv,ix_,_ = s2mpj_ii("Y00106",ix_)
        arrset(pb.xnames,iv,"Y00106")
        ig = ig_["CAP01505"]
        pbm.A[ig,iv] += Float64(.01704)
        ig = ig_["MND00106"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00106",ix_)
        arrset(pb.xnames,iv,"Y00106")
        ig = ig_["CAP01406"]
        pbm.A[ig,iv] += Float64(.03074)
        ig = ig_["CAP06305"]
        pbm.A[ig,iv] += Float64(.00135)
        iv,ix_,_ = s2mpj_ii("Y00106",ix_)
        arrset(pb.xnames,iv,"Y00106")
        ig = ig_["MXD00106"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00202",ix_)
        arrset(pb.xnames,iv,"Y00202")
        ig = ig_["CAP06302"]
        pbm.A[ig,iv] += Float64(.00132)
        ig = ig_["CAP00301"]
        pbm.A[ig,iv] += Float64(.02059)
        iv,ix_,_ = s2mpj_ii("Y00202",ix_)
        arrset(pb.xnames,iv,"Y00202")
        ig = ig_["MXD00104"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-31.15)
        iv,ix_,_ = s2mpj_ii("Y00202",ix_)
        arrset(pb.xnames,iv,"Y00202")
        ig = ig_["CAP00102"]
        pbm.A[ig,iv] += Float64(.01176)
        ig = ig_["MND00106"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00202",ix_)
        arrset(pb.xnames,iv,"Y00202")
        ig = ig_["MND00104"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00105"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00202",ix_)
        arrset(pb.xnames,iv,"Y00202")
        ig = ig_["MXD00105"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP06301"]
        pbm.A[ig,iv] += Float64(.00141)
        iv,ix_,_ = s2mpj_ii("Y00202",ix_)
        arrset(pb.xnames,iv,"Y00202")
        ig = ig_["INV00302"]
        pbm.A[ig,iv] += Float64(-.52969)
        ig = ig_["CAP00302"]
        pbm.A[ig,iv] += Float64(.00965)
        iv,ix_,_ = s2mpj_ii("Y00202",ix_)
        arrset(pb.xnames,iv,"Y00202")
        ig = ig_["INV00301"]
        pbm.A[ig,iv] += Float64(-1.3791)
        ig = ig_["MXD00106"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00202",ix_)
        arrset(pb.xnames,iv,"Y00202")
        ig = ig_["MXD00102"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00102"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00202",ix_)
        arrset(pb.xnames,iv,"Y00202")
        ig = ig_["CAP00101"]
        pbm.A[ig,iv] += Float64(.02083)
        ig = ig_["MXD00103"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00202",ix_)
        arrset(pb.xnames,iv,"Y00202")
        ig = ig_["MND00103"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00203",ix_)
        arrset(pb.xnames,iv,"Y00203")
        ig = ig_["MND00104"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00104"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00203",ix_)
        arrset(pb.xnames,iv,"Y00203")
        ig = ig_["MXD00103"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00103"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00203",ix_)
        arrset(pb.xnames,iv,"Y00203")
        ig = ig_["CAP06302"]
        pbm.A[ig,iv] += Float64(.0013)
        ig = ig_["CAP06303"]
        pbm.A[ig,iv] += Float64(.00143)
        iv,ix_,_ = s2mpj_ii("Y00203",ix_)
        arrset(pb.xnames,iv,"Y00203")
        ig = ig_["CAP00103"]
        pbm.A[ig,iv] += Float64(.01336)
        ig = ig_["CAP00303"]
        pbm.A[ig,iv] += Float64(.01123)
        iv,ix_,_ = s2mpj_ii("Y00203",ix_)
        arrset(pb.xnames,iv,"Y00203")
        ig = ig_["MXD00105"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP00302"]
        pbm.A[ig,iv] += Float64(.019)
        iv,ix_,_ = s2mpj_ii("Y00203",ix_)
        arrset(pb.xnames,iv,"Y00203")
        ig = ig_["MND00106"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00106"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00203",ix_)
        arrset(pb.xnames,iv,"Y00203")
        ig = ig_["INV00302"]
        pbm.A[ig,iv] += Float64(-1.27302)
        ig = ig_["INV00303"]
        pbm.A[ig,iv] += Float64(-.63578)
        iv,ix_,_ = s2mpj_ii("Y00203",ix_)
        arrset(pb.xnames,iv,"Y00203")
        ig = ig_["CAP00102"]
        pbm.A[ig,iv] += Float64(.01923)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-29.57)
        iv,ix_,_ = s2mpj_ii("Y00203",ix_)
        arrset(pb.xnames,iv,"Y00203")
        ig = ig_["MND00105"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00204",ix_)
        arrset(pb.xnames,iv,"Y00204")
        ig = ig_["MND00104"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP00103"]
        pbm.A[ig,iv] += Float64(.01923)
        iv,ix_,_ = s2mpj_ii("Y00204",ix_)
        arrset(pb.xnames,iv,"Y00204")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-28.08)
        ig = ig_["MND00106"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00204",ix_)
        arrset(pb.xnames,iv,"Y00204")
        ig = ig_["MXD00106"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP06304"]
        pbm.A[ig,iv] += Float64(.00143)
        iv,ix_,_ = s2mpj_ii("Y00204",ix_)
        arrset(pb.xnames,iv,"Y00204")
        ig = ig_["CAP06303"]
        pbm.A[ig,iv] += Float64(.0013)
        ig = ig_["INV00303"]
        pbm.A[ig,iv] += Float64(-1.27302)
        iv,ix_,_ = s2mpj_ii("Y00204",ix_)
        arrset(pb.xnames,iv,"Y00204")
        ig = ig_["INV00304"]
        pbm.A[ig,iv] += Float64(-.63578)
        ig = ig_["CAP00304"]
        pbm.A[ig,iv] += Float64(.01123)
        iv,ix_,_ = s2mpj_ii("Y00204",ix_)
        arrset(pb.xnames,iv,"Y00204")
        ig = ig_["CAP00104"]
        pbm.A[ig,iv] += Float64(.01336)
        ig = ig_["CAP00303"]
        pbm.A[ig,iv] += Float64(.019)
        iv,ix_,_ = s2mpj_ii("Y00204",ix_)
        arrset(pb.xnames,iv,"Y00204")
        ig = ig_["MXD00104"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00105"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00204",ix_)
        arrset(pb.xnames,iv,"Y00204")
        ig = ig_["MXD00105"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00205",ix_)
        arrset(pb.xnames,iv,"Y00205")
        ig = ig_["INV00304"]
        pbm.A[ig,iv] += Float64(-1.22587)
        ig = ig_["MXD00105"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00205",ix_)
        arrset(pb.xnames,iv,"Y00205")
        ig = ig_["CAP00305"]
        pbm.A[ig,iv] += Float64(.01194)
        ig = ig_["CAP00105"]
        pbm.A[ig,iv] += Float64(.01408)
        iv,ix_,_ = s2mpj_ii("Y00205",ix_)
        arrset(pb.xnames,iv,"Y00205")
        ig = ig_["CAP06304"]
        pbm.A[ig,iv] += Float64(.00125)
        ig = ig_["CAP00304"]
        pbm.A[ig,iv] += Float64(.0183)
        iv,ix_,_ = s2mpj_ii("Y00205",ix_)
        arrset(pb.xnames,iv,"Y00205")
        ig = ig_["INV00305"]
        pbm.A[ig,iv] += Float64(-.68292)
        ig = ig_["CAP06305"]
        pbm.A[ig,iv] += Float64(.00147)
        iv,ix_,_ = s2mpj_ii("Y00205",ix_)
        arrset(pb.xnames,iv,"Y00205")
        ig = ig_["CAP00104"]
        pbm.A[ig,iv] += Float64(.01852)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-26.66)
        iv,ix_,_ = s2mpj_ii("Y00205",ix_)
        arrset(pb.xnames,iv,"Y00205")
        ig = ig_["MND00105"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00106"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00205",ix_)
        arrset(pb.xnames,iv,"Y00205")
        ig = ig_["MXD00106"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00206",ix_)
        arrset(pb.xnames,iv,"Y00206")
        ig = ig_["MND00106"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-25.31)
        iv,ix_,_ = s2mpj_ii("Y00206",ix_)
        arrset(pb.xnames,iv,"Y00206")
        ig = ig_["INV00305"]
        pbm.A[ig,iv] += Float64(-1.32394)
        ig = ig_["INV00306"]
        pbm.A[ig,iv] += Float64(-1.90879)
        iv,ix_,_ = s2mpj_ii("Y00206",ix_)
        arrset(pb.xnames,iv,"Y00206")
        ig = ig_["CAP06306"]
        pbm.A[ig,iv] += Float64(.00272)
        ig = ig_["CAP00105"]
        pbm.A[ig,iv] += Float64(.02)
        iv,ix_,_ = s2mpj_ii("Y00206",ix_)
        arrset(pb.xnames,iv,"Y00206")
        ig = ig_["CAP00306"]
        pbm.A[ig,iv] += Float64(.03024)
        ig = ig_["CAP00305"]
        pbm.A[ig,iv] += Float64(.01976)
        iv,ix_,_ = s2mpj_ii("Y00206",ix_)
        arrset(pb.xnames,iv,"Y00206")
        ig = ig_["CAP00106"]
        pbm.A[ig,iv] += Float64(.03259)
        ig = ig_["CAP06305"]
        pbm.A[ig,iv] += Float64(.00135)
        iv,ix_,_ = s2mpj_ii("Y00206",ix_)
        arrset(pb.xnames,iv,"Y00206")
        ig = ig_["MXD00106"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00302",ix_)
        arrset(pb.xnames,iv,"Y00302")
        ig = ig_["MXD00205"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00205"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00302",ix_)
        arrset(pb.xnames,iv,"Y00302")
        ig = ig_["CAP00501"]
        pbm.A[ig,iv] += Float64(.01455)
        ig = ig_["MXD00202"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00302",ix_)
        arrset(pb.xnames,iv,"Y00302")
        ig = ig_["MND00204"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-25.46)
        iv,ix_,_ = s2mpj_ii("Y00302",ix_)
        arrset(pb.xnames,iv,"Y00302")
        ig = ig_["CAP00401"]
        pbm.A[ig,iv] += Float64(.01426)
        ig = ig_["CAP00502"]
        pbm.A[ig,iv] += Float64(.06023)
        iv,ix_,_ = s2mpj_ii("Y00302",ix_)
        arrset(pb.xnames,iv,"Y00302")
        ig = ig_["MND00202"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00204"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00302",ix_)
        arrset(pb.xnames,iv,"Y00302")
        ig = ig_["CAP06302"]
        pbm.A[ig,iv] += Float64(.00165)
        ig = ig_["MND00206"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00302",ix_)
        arrset(pb.xnames,iv,"Y00302")
        ig = ig_["MXD00206"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["INV00301"]
        pbm.A[ig,iv] += Float64(-.67787)
        iv,ix_,_ = s2mpj_ii("Y00302",ix_)
        arrset(pb.xnames,iv,"Y00302")
        ig = ig_["CAP06301"]
        pbm.A[ig,iv] += Float64(.00015)
        ig = ig_["CAP00402"]
        pbm.A[ig,iv] += Float64(.04611)
        iv,ix_,_ = s2mpj_ii("Y00302",ix_)
        arrset(pb.xnames,iv,"Y00302")
        ig = ig_["INV00302"]
        pbm.A[ig,iv] += Float64(-1.76124)
        ig = ig_["MND00203"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00302",ix_)
        arrset(pb.xnames,iv,"Y00302")
        ig = ig_["MXD00203"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00303",ix_)
        arrset(pb.xnames,iv,"Y00303")
        ig = ig_["CAP00503"]
        pbm.A[ig,iv] += Float64(.06135)
        ig = ig_["CAP00403"]
        pbm.A[ig,iv] += Float64(.0472)
        iv,ix_,_ = s2mpj_ii("Y00303",ix_)
        arrset(pb.xnames,iv,"Y00303")
        ig = ig_["MXD00205"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00203"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00303",ix_)
        arrset(pb.xnames,iv,"Y00303")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-30.04)
        ig = ig_["MND00204"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00303",ix_)
        arrset(pb.xnames,iv,"Y00303")
        ig = ig_["MND00205"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP06303"]
        pbm.A[ig,iv] += Float64(.00166)
        iv,ix_,_ = s2mpj_ii("Y00303",ix_)
        arrset(pb.xnames,iv,"Y00303")
        ig = ig_["MXD00203"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["INV00303"]
        pbm.A[ig,iv] += Float64(-1.81339)
        iv,ix_,_ = s2mpj_ii("Y00303",ix_)
        arrset(pb.xnames,iv,"Y00303")
        ig = ig_["INV00302"]
        pbm.A[ig,iv] += Float64(-.62573)
        ig = ig_["MXD00206"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00303",ix_)
        arrset(pb.xnames,iv,"Y00303")
        ig = ig_["MND00206"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00204"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00303",ix_)
        arrset(pb.xnames,iv,"Y00303")
        ig = ig_["CAP06302"]
        pbm.A[ig,iv] += Float64(.00014)
        ig = ig_["CAP00502"]
        pbm.A[ig,iv] += Float64(.01343)
        iv,ix_,_ = s2mpj_ii("Y00303",ix_)
        arrset(pb.xnames,iv,"Y00303")
        ig = ig_["CAP00402"]
        pbm.A[ig,iv] += Float64(.01316)
        iv,ix_,_ = s2mpj_ii("Y00304",ix_)
        arrset(pb.xnames,iv,"Y00304")
        ig = ig_["CAP00503"]
        pbm.A[ig,iv] += Float64(.01343)
        ig = ig_["MXD00204"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00304",ix_)
        arrset(pb.xnames,iv,"Y00304")
        ig = ig_["MND00204"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP06304"]
        pbm.A[ig,iv] += Float64(.00166)
        iv,ix_,_ = s2mpj_ii("Y00304",ix_)
        arrset(pb.xnames,iv,"Y00304")
        ig = ig_["INV00303"]
        pbm.A[ig,iv] += Float64(-.62573)
        ig = ig_["MND00205"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00304",ix_)
        arrset(pb.xnames,iv,"Y00304")
        ig = ig_["MXD00205"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00206"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00304",ix_)
        arrset(pb.xnames,iv,"Y00304")
        ig = ig_["INV00304"]
        pbm.A[ig,iv] += Float64(-1.81339)
        ig = ig_["CAP00404"]
        pbm.A[ig,iv] += Float64(.0472)
        iv,ix_,_ = s2mpj_ii("Y00304",ix_)
        arrset(pb.xnames,iv,"Y00304")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-28.52)
        ig = ig_["CAP00504"]
        pbm.A[ig,iv] += Float64(.06135)
        iv,ix_,_ = s2mpj_ii("Y00304",ix_)
        arrset(pb.xnames,iv,"Y00304")
        ig = ig_["CAP06303"]
        pbm.A[ig,iv] += Float64(.00014)
        ig = ig_["CAP00403"]
        pbm.A[ig,iv] += Float64(.01316)
        iv,ix_,_ = s2mpj_ii("Y00304",ix_)
        arrset(pb.xnames,iv,"Y00304")
        ig = ig_["MND00206"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00305",ix_)
        arrset(pb.xnames,iv,"Y00305")
        ig = ig_["MND00206"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP06305"]
        pbm.A[ig,iv] += Float64(.00167)
        iv,ix_,_ = s2mpj_ii("Y00305",ix_)
        arrset(pb.xnames,iv,"Y00305")
        ig = ig_["MXD00206"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["INV00305"]
        pbm.A[ig,iv] += Float64(-1.83656)
        iv,ix_,_ = s2mpj_ii("Y00305",ix_)
        arrset(pb.xnames,iv,"Y00305")
        ig = ig_["INV00304"]
        pbm.A[ig,iv] += Float64(-.60255)
        ig = ig_["CAP00505"]
        pbm.A[ig,iv] += Float64(.06184)
        iv,ix_,_ = s2mpj_ii("Y00305",ix_)
        arrset(pb.xnames,iv,"Y00305")
        ig = ig_["CAP00404"]
        pbm.A[ig,iv] += Float64(.01268)
        ig = ig_["MXD00205"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00305",ix_)
        arrset(pb.xnames,iv,"Y00305")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-27.08)
        ig = ig_["CAP00504"]
        pbm.A[ig,iv] += Float64(.01293)
        iv,ix_,_ = s2mpj_ii("Y00305",ix_)
        arrset(pb.xnames,iv,"Y00305")
        ig = ig_["MND00205"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP06304"]
        pbm.A[ig,iv] += Float64(.00013)
        iv,ix_,_ = s2mpj_ii("Y00305",ix_)
        arrset(pb.xnames,iv,"Y00305")
        ig = ig_["CAP00405"]
        pbm.A[ig,iv] += Float64(.04769)
        iv,ix_,_ = s2mpj_ii("Y00306",ix_)
        arrset(pb.xnames,iv,"Y00306")
        ig = ig_["CAP06306"]
        pbm.A[ig,iv] += Float64(.0018)
        ig = ig_["CAP06305"]
        pbm.A[ig,iv] += Float64(.00014)
        iv,ix_,_ = s2mpj_ii("Y00306",ix_)
        arrset(pb.xnames,iv,"Y00306")
        ig = ig_["INV00306"]
        pbm.A[ig,iv] += Float64(-2.43911)
        ig = ig_["CAP00406"]
        pbm.A[ig,iv] += Float64(.06037)
        iv,ix_,_ = s2mpj_ii("Y00306",ix_)
        arrset(pb.xnames,iv,"Y00306")
        ig = ig_["CAP00506"]
        pbm.A[ig,iv] += Float64(.07478)
        ig = ig_["CAP00505"]
        pbm.A[ig,iv] += Float64(.01397)
        iv,ix_,_ = s2mpj_ii("Y00306",ix_)
        arrset(pb.xnames,iv,"Y00306")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-25.71)
        ig = ig_["MND00206"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00306",ix_)
        arrset(pb.xnames,iv,"Y00306")
        ig = ig_["MXD00206"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP00405"]
        pbm.A[ig,iv] += Float64(.01369)
        iv,ix_,_ = s2mpj_ii("Y00306",ix_)
        arrset(pb.xnames,iv,"Y00306")
        ig = ig_["INV00305"]
        pbm.A[ig,iv] += Float64(-.65075)
        iv,ix_,_ = s2mpj_ii("Y00402",ix_)
        arrset(pb.xnames,iv,"Y00402")
        ig = ig_["CAP06302"]
        pbm.A[ig,iv] += Float64(.00165)
        ig = ig_["CAP01701"]
        pbm.A[ig,iv] += Float64(.00928)
        iv,ix_,_ = s2mpj_ii("Y00402",ix_)
        arrset(pb.xnames,iv,"Y00402")
        ig = ig_["CAP06301"]
        pbm.A[ig,iv] += Float64(.00015)
        ig = ig_["MND00206"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00402",ix_)
        arrset(pb.xnames,iv,"Y00402")
        ig = ig_["MXD00206"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-25.42)
        iv,ix_,_ = s2mpj_ii("Y00402",ix_)
        arrset(pb.xnames,iv,"Y00402")
        ig = ig_["CAP01602"]
        pbm.A[ig,iv] += Float64(.05177)
        ig = ig_["MXD00205"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00402",ix_)
        arrset(pb.xnames,iv,"Y00402")
        ig = ig_["MND00205"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00204"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00402",ix_)
        arrset(pb.xnames,iv,"Y00402")
        ig = ig_["CAP01601"]
        pbm.A[ig,iv] += Float64(.0086)
        ig = ig_["INV00301"]
        pbm.A[ig,iv] += Float64(-.37298)
        iv,ix_,_ = s2mpj_ii("Y00402",ix_)
        arrset(pb.xnames,iv,"Y00402")
        ig = ig_["MXD00203"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00203"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00402",ix_)
        arrset(pb.xnames,iv,"Y00402")
        ig = ig_["MND00204"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00202"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00402",ix_)
        arrset(pb.xnames,iv,"Y00402")
        ig = ig_["INV00302"]
        pbm.A[ig,iv] += Float64(-2.06613)
        ig = ig_["MND00202"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00402",ix_)
        arrset(pb.xnames,iv,"Y00402")
        ig = ig_["CAP01702"]
        pbm.A[ig,iv] += Float64(.061)
        iv,ix_,_ = s2mpj_ii("Y00403",ix_)
        arrset(pb.xnames,iv,"Y00403")
        ig = ig_["MXD00206"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["INV00303"]
        pbm.A[ig,iv] += Float64(-2.09482)
        iv,ix_,_ = s2mpj_ii("Y00403",ix_)
        arrset(pb.xnames,iv,"Y00403")
        ig = ig_["MXD00205"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00205"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00403",ix_)
        arrset(pb.xnames,iv,"Y00403")
        ig = ig_["MND00206"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP06303"]
        pbm.A[ig,iv] += Float64(.00166)
        iv,ix_,_ = s2mpj_ii("Y00403",ix_)
        arrset(pb.xnames,iv,"Y00403")
        ig = ig_["CAP01702"]
        pbm.A[ig,iv] += Float64(.00857)
        ig = ig_["MND00204"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00403",ix_)
        arrset(pb.xnames,iv,"Y00403")
        ig = ig_["CAP06302"]
        pbm.A[ig,iv] += Float64(.00014)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-30.0)
        iv,ix_,_ = s2mpj_ii("Y00403",ix_)
        arrset(pb.xnames,iv,"Y00403")
        ig = ig_["CAP01603"]
        pbm.A[ig,iv] += Float64(.05243)
        ig = ig_["CAP01602"]
        pbm.A[ig,iv] += Float64(.00794)
        iv,ix_,_ = s2mpj_ii("Y00403",ix_)
        arrset(pb.xnames,iv,"Y00403")
        ig = ig_["CAP01703"]
        pbm.A[ig,iv] += Float64(.06171)
        ig = ig_["MXD00204"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00403",ix_)
        arrset(pb.xnames,iv,"Y00403")
        ig = ig_["INV00302"]
        pbm.A[ig,iv] += Float64(-.34429)
        ig = ig_["MND00203"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00403",ix_)
        arrset(pb.xnames,iv,"Y00403")
        ig = ig_["MXD00203"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00404",ix_)
        arrset(pb.xnames,iv,"Y00404")
        ig = ig_["MXD00204"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP06304"]
        pbm.A[ig,iv] += Float64(.00166)
        iv,ix_,_ = s2mpj_ii("Y00404",ix_)
        arrset(pb.xnames,iv,"Y00404")
        ig = ig_["INV00303"]
        pbm.A[ig,iv] += Float64(-.34429)
        ig = ig_["MND00204"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00404",ix_)
        arrset(pb.xnames,iv,"Y00404")
        ig = ig_["CAP01604"]
        pbm.A[ig,iv] += Float64(.05243)
        ig = ig_["CAP01704"]
        pbm.A[ig,iv] += Float64(.06171)
        iv,ix_,_ = s2mpj_ii("Y00404",ix_)
        arrset(pb.xnames,iv,"Y00404")
        ig = ig_["INV00304"]
        pbm.A[ig,iv] += Float64(-2.09482)
        ig = ig_["CAP01603"]
        pbm.A[ig,iv] += Float64(.00794)
        iv,ix_,_ = s2mpj_ii("Y00404",ix_)
        arrset(pb.xnames,iv,"Y00404")
        ig = ig_["CAP01703"]
        pbm.A[ig,iv] += Float64(.00857)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-28.48)
        iv,ix_,_ = s2mpj_ii("Y00404",ix_)
        arrset(pb.xnames,iv,"Y00404")
        ig = ig_["MND00206"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00205"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00404",ix_)
        arrset(pb.xnames,iv,"Y00404")
        ig = ig_["MND00205"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP06303"]
        pbm.A[ig,iv] += Float64(.00014)
        iv,ix_,_ = s2mpj_ii("Y00404",ix_)
        arrset(pb.xnames,iv,"Y00404")
        ig = ig_["MXD00206"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00405",ix_)
        arrset(pb.xnames,iv,"Y00405")
        ig = ig_["MND00205"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP01604"]
        pbm.A[ig,iv] += Float64(.00765)
        iv,ix_,_ = s2mpj_ii("Y00405",ix_)
        arrset(pb.xnames,iv,"Y00405")
        ig = ig_["MXD00206"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP01704"]
        pbm.A[ig,iv] += Float64(.00825)
        iv,ix_,_ = s2mpj_ii("Y00405",ix_)
        arrset(pb.xnames,iv,"Y00405")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-27.04)
        ig = ig_["CAP06304"]
        pbm.A[ig,iv] += Float64(.00013)
        iv,ix_,_ = s2mpj_ii("Y00405",ix_)
        arrset(pb.xnames,iv,"Y00405")
        ig = ig_["MND00206"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["INV00305"]
        pbm.A[ig,iv] += Float64(-2.10758)
        iv,ix_,_ = s2mpj_ii("Y00405",ix_)
        arrset(pb.xnames,iv,"Y00405")
        ig = ig_["CAP01605"]
        pbm.A[ig,iv] += Float64(.05272)
        ig = ig_["CAP01705"]
        pbm.A[ig,iv] += Float64(.06203)
        iv,ix_,_ = s2mpj_ii("Y00405",ix_)
        arrset(pb.xnames,iv,"Y00405")
        ig = ig_["CAP06305"]
        pbm.A[ig,iv] += Float64(.00167)
        ig = ig_["INV00304"]
        pbm.A[ig,iv] += Float64(-.33154)
        iv,ix_,_ = s2mpj_ii("Y00405",ix_)
        arrset(pb.xnames,iv,"Y00405")
        ig = ig_["MXD00205"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00406",ix_)
        arrset(pb.xnames,iv,"Y00406")
        ig = ig_["CAP01706"]
        pbm.A[ig,iv] += Float64(.07028)
        ig = ig_["CAP01605"]
        pbm.A[ig,iv] += Float64(.00826)
        iv,ix_,_ = s2mpj_ii("Y00406",ix_)
        arrset(pb.xnames,iv,"Y00406")
        ig = ig_["INV00305"]
        pbm.A[ig,iv] += Float64(-.35806)
        ig = ig_["MND00206"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00406",ix_)
        arrset(pb.xnames,iv,"Y00406")
        ig = ig_["MXD00206"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP06305"]
        pbm.A[ig,iv] += Float64(.00014)
        iv,ix_,_ = s2mpj_ii("Y00406",ix_)
        arrset(pb.xnames,iv,"Y00406")
        ig = ig_["INV00306"]
        pbm.A[ig,iv] += Float64(-2.43911)
        ig = ig_["CAP01606"]
        pbm.A[ig,iv] += Float64(.06037)
        iv,ix_,_ = s2mpj_ii("Y00406",ix_)
        arrset(pb.xnames,iv,"Y00406")
        ig = ig_["CAP01705"]
        pbm.A[ig,iv] += Float64(.00891)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-25.67)
        iv,ix_,_ = s2mpj_ii("Y00406",ix_)
        arrset(pb.xnames,iv,"Y00406")
        ig = ig_["CAP06306"]
        pbm.A[ig,iv] += Float64(.0018)
        iv,ix_,_ = s2mpj_ii("Y00503",ix_)
        arrset(pb.xnames,iv,"Y00503")
        ig = ig_["CAP00401"]
        pbm.A[ig,iv] += Float64(.00401)
        ig = ig_["MND00305"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00503",ix_)
        arrset(pb.xnames,iv,"Y00503")
        ig = ig_["INV00301"]
        pbm.A[ig,iv] += Float64(-.216)
        ig = ig_["MND00304"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00503",ix_)
        arrset(pb.xnames,iv,"Y00503")
        ig = ig_["CAP06002"]
        pbm.A[ig,iv] += Float64(.00225)
        ig = ig_["CAP06003"]
        pbm.A[ig,iv] += Float64(.00124)
        iv,ix_,_ = s2mpj_ii("Y00503",ix_)
        arrset(pb.xnames,iv,"Y00503")
        ig = ig_["CAP00502"]
        pbm.A[ig,iv] += Float64(.03973)
        ig = ig_["MXD00306"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00503",ix_)
        arrset(pb.xnames,iv,"Y00503")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-33.21)
        ig = ig_["MXD00304"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00503",ix_)
        arrset(pb.xnames,iv,"Y00503")
        ig = ig_["MND00306"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00305"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00503",ix_)
        arrset(pb.xnames,iv,"Y00503")
        ig = ig_["CAP00501"]
        pbm.A[ig,iv] += Float64(.00331)
        ig = ig_["INV00302"]
        pbm.A[ig,iv] += Float64(-1.18801)
        iv,ix_,_ = s2mpj_ii("Y00503",ix_)
        arrset(pb.xnames,iv,"Y00503")
        ig = ig_["MXD00303"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP00402"]
        pbm.A[ig,iv] += Float64(.03074)
        iv,ix_,_ = s2mpj_ii("Y00503",ix_)
        arrset(pb.xnames,iv,"Y00503")
        ig = ig_["MND00303"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00504",ix_)
        arrset(pb.xnames,iv,"Y00504")
        ig = ig_["MND00306"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP06004"]
        pbm.A[ig,iv] += Float64(.00124)
        iv,ix_,_ = s2mpj_ii("Y00504",ix_)
        arrset(pb.xnames,iv,"Y00504")
        ig = ig_["MXD00306"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00304"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00504",ix_)
        arrset(pb.xnames,iv,"Y00504")
        ig = ig_["MND00304"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP00403"]
        pbm.A[ig,iv] += Float64(.03341)
        iv,ix_,_ = s2mpj_ii("Y00504",ix_)
        arrset(pb.xnames,iv,"Y00504")
        ig = ig_["CAP00402"]
        pbm.A[ig,iv] += Float64(.00134)
        ig = ig_["CAP06003"]
        pbm.A[ig,iv] += Float64(.00225)
        iv,ix_,_ = s2mpj_ii("Y00504",ix_)
        arrset(pb.xnames,iv,"Y00504")
        ig = ig_["INV00303"]
        pbm.A[ig,iv] += Float64(-1.29601)
        ig = ig_["MXD00305"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00504",ix_)
        arrset(pb.xnames,iv,"Y00504")
        ig = ig_["MND00305"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["INV00302"]
        pbm.A[ig,iv] += Float64(-.108)
        iv,ix_,_ = s2mpj_ii("Y00504",ix_)
        arrset(pb.xnames,iv,"Y00504")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-28.74)
        ig = ig_["CAP00503"]
        pbm.A[ig,iv] += Float64(.04304)
        iv,ix_,_ = s2mpj_ii("Y00505",ix_)
        arrset(pb.xnames,iv,"Y00505")
        ig = ig_["CAP06005"]
        pbm.A[ig,iv] += Float64(.00133)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-27.28)
        iv,ix_,_ = s2mpj_ii("Y00505",ix_)
        arrset(pb.xnames,iv,"Y00505")
        ig = ig_["INV00303"]
        pbm.A[ig,iv] += Float64(-.104)
        ig = ig_["CAP06004"]
        pbm.A[ig,iv] += Float64(.00216)
        iv,ix_,_ = s2mpj_ii("Y00505",ix_)
        arrset(pb.xnames,iv,"Y00505")
        ig = ig_["MXD00306"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP00505"]
        pbm.A[ig,iv] += Float64(.00159)
        iv,ix_,_ = s2mpj_ii("Y00505",ix_)
        arrset(pb.xnames,iv,"Y00505")
        ig = ig_["CAP00403"]
        pbm.A[ig,iv] += Float64(.00129)
        ig = ig_["CAP00404"]
        pbm.A[ig,iv] += Float64(.03346)
        iv,ix_,_ = s2mpj_ii("Y00505",ix_)
        arrset(pb.xnames,iv,"Y00505")
        ig = ig_["CAP00504"]
        pbm.A[ig,iv] += Float64(.04145)
        ig = ig_["INV00304"]
        pbm.A[ig,iv] += Float64(-1.30001)
        iv,ix_,_ = s2mpj_ii("Y00505",ix_)
        arrset(pb.xnames,iv,"Y00505")
        ig = ig_["MND00305"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00305"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00505",ix_)
        arrset(pb.xnames,iv,"Y00505")
        ig = ig_["MND00306"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00506",ix_)
        arrset(pb.xnames,iv,"Y00506")
        ig = ig_["MND00306"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["INV00305"]
        pbm.A[ig,iv] += Float64(-1.51633)
        iv,ix_,_ = s2mpj_ii("Y00506",ix_)
        arrset(pb.xnames,iv,"Y00506")
        ig = ig_["MXD00306"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP00405"]
        pbm.A[ig,iv] += Float64(.03753)
        iv,ix_,_ = s2mpj_ii("Y00506",ix_)
        arrset(pb.xnames,iv,"Y00506")
        ig = ig_["INV00306"]
        pbm.A[ig,iv] += Float64(-1.40401)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-25.9)
        iv,ix_,_ = s2mpj_ii("Y00506",ix_)
        arrset(pb.xnames,iv,"Y00506")
        ig = ig_["CAP06005"]
        pbm.A[ig,iv] += Float64(.00234)
        ig = ig_["CAP06006"]
        pbm.A[ig,iv] += Float64(.00349)
        iv,ix_,_ = s2mpj_ii("Y00506",ix_)
        arrset(pb.xnames,iv,"Y00506")
        ig = ig_["INV00304"]
        pbm.A[ig,iv] += Float64(-.05616)
        ig = ig_["CAP00406"]
        pbm.A[ig,iv] += Float64(.03475)
        iv,ix_,_ = s2mpj_ii("Y00506",ix_)
        arrset(pb.xnames,iv,"Y00506")
        ig = ig_["CAP00506"]
        pbm.A[ig,iv] += Float64(.04304)
        iv,ix_,_ = s2mpj_ii("Y00603",ix_)
        arrset(pb.xnames,iv,"Y00603")
        ig = ig_["MXD00303"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["INV00302"]
        pbm.A[ig,iv] += Float64(-1.29601)
        iv,ix_,_ = s2mpj_ii("Y00603",ix_)
        arrset(pb.xnames,iv,"Y00603")
        ig = ig_["MND00305"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP01602"]
        pbm.A[ig,iv] += Float64(.03208)
        iv,ix_,_ = s2mpj_ii("Y00603",ix_)
        arrset(pb.xnames,iv,"Y00603")
        ig = ig_["MND00303"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00304"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00603",ix_)
        arrset(pb.xnames,iv,"Y00603")
        ig = ig_["MND00306"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP06002"]
        pbm.A[ig,iv] += Float64(.00225)
        iv,ix_,_ = s2mpj_ii("Y00603",ix_)
        arrset(pb.xnames,iv,"Y00603")
        ig = ig_["MXD00306"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP01702"]
        pbm.A[ig,iv] += Float64(.03734)
        iv,ix_,_ = s2mpj_ii("Y00603",ix_)
        arrset(pb.xnames,iv,"Y00603")
        ig = ig_["INV00301"]
        pbm.A[ig,iv] += Float64(-.054)
        ig = ig_["CAP06003"]
        pbm.A[ig,iv] += Float64(.00124)
        iv,ix_,_ = s2mpj_ii("Y00603",ix_)
        arrset(pb.xnames,iv,"Y00603")
        ig = ig_["INV00303"]
        pbm.A[ig,iv] += Float64(-.054)
        ig = ig_["CAP01703"]
        pbm.A[ig,iv] += Float64(.00233)
        iv,ix_,_ = s2mpj_ii("Y00603",ix_)
        arrset(pb.xnames,iv,"Y00603")
        ig = ig_["MND00304"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP01603"]
        pbm.A[ig,iv] += Float64(.00167)
        iv,ix_,_ = s2mpj_ii("Y00603",ix_)
        arrset(pb.xnames,iv,"Y00603")
        ig = ig_["CAP01701"]
        pbm.A[ig,iv] += Float64(.00078)
        ig = ig_["MXD00305"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00603",ix_)
        arrset(pb.xnames,iv,"Y00603")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-33.21)
        ig = ig_["CAP01601"]
        pbm.A[ig,iv] += Float64(.001)
        iv,ix_,_ = s2mpj_ii("Y00604",ix_)
        arrset(pb.xnames,iv,"Y00604")
        ig = ig_["MXD00306"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-28.74)
        iv,ix_,_ = s2mpj_ii("Y00604",ix_)
        arrset(pb.xnames,iv,"Y00604")
        ig = ig_["INV00304"]
        pbm.A[ig,iv] += Float64(-.054)
        ig = ig_["MXD00304"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00604",ix_)
        arrset(pb.xnames,iv,"Y00604")
        ig = ig_["CAP06003"]
        pbm.A[ig,iv] += Float64(.00225)
        ig = ig_["MND00304"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00604",ix_)
        arrset(pb.xnames,iv,"Y00604")
        ig = ig_["CAP01603"]
        pbm.A[ig,iv] += Float64(.03308)
        ig = ig_["CAP01703"]
        pbm.A[ig,iv] += Float64(.03812)
        iv,ix_,_ = s2mpj_ii("Y00604",ix_)
        arrset(pb.xnames,iv,"Y00604")
        ig = ig_["MND00306"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP01604"]
        pbm.A[ig,iv] += Float64(.00167)
        iv,ix_,_ = s2mpj_ii("Y00604",ix_)
        arrset(pb.xnames,iv,"Y00604")
        ig = ig_["INV00303"]
        pbm.A[ig,iv] += Float64(-1.35001)
        ig = ig_["MND00305"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00604",ix_)
        arrset(pb.xnames,iv,"Y00604")
        ig = ig_["CAP06004"]
        pbm.A[ig,iv] += Float64(.00124)
        ig = ig_["CAP01704"]
        pbm.A[ig,iv] += Float64(.00233)
        iv,ix_,_ = s2mpj_ii("Y00604",ix_)
        arrset(pb.xnames,iv,"Y00604")
        ig = ig_["MXD00305"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00605",ix_)
        arrset(pb.xnames,iv,"Y00605")
        ig = ig_["CAP06005"]
        pbm.A[ig,iv] += Float64(.00133)
        ig = ig_["MND00306"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00605",ix_)
        arrset(pb.xnames,iv,"Y00605")
        ig = ig_["INV00305"]
        pbm.A[ig,iv] += Float64(-1.40401)
        ig = ig_["CAP01705"]
        pbm.A[ig,iv] += Float64(.00375)
        iv,ix_,_ = s2mpj_ii("Y00605",ix_)
        arrset(pb.xnames,iv,"Y00605")
        ig = ig_["MXD00305"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP01704"]
        pbm.A[ig,iv] += Float64(.03671)
        iv,ix_,_ = s2mpj_ii("Y00605",ix_)
        arrset(pb.xnames,iv,"Y00605")
        ig = ig_["MND00305"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00306"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00605",ix_)
        arrset(pb.xnames,iv,"Y00605")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-27.28)
        ig = ig_["INV00304"]
        pbm.A[ig,iv] += Float64(-1.30001)
        iv,ix_,_ = s2mpj_ii("Y00605",ix_)
        arrset(pb.xnames,iv,"Y00605")
        ig = ig_["CAP01604"]
        pbm.A[ig,iv] += Float64(.03185)
        ig = ig_["CAP01605"]
        pbm.A[ig,iv] += Float64(.0029)
        iv,ix_,_ = s2mpj_ii("Y00606",ix_)
        arrset(pb.xnames,iv,"Y00606")
        ig = ig_["CAP06005"]
        pbm.A[ig,iv] += Float64(.00234)
        ig = ig_["MXD00306"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00606",ix_)
        arrset(pb.xnames,iv,"Y00606")
        ig = ig_["CAP06006"]
        pbm.A[ig,iv] += Float64(.00349)
        ig = ig_["MND00306"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00606",ix_)
        arrset(pb.xnames,iv,"Y00606")
        ig = ig_["CAP01706"]
        pbm.A[ig,iv] += Float64(.04046)
        ig = ig_["CAP01606"]
        pbm.A[ig,iv] += Float64(.03475)
        iv,ix_,_ = s2mpj_ii("Y00606",ix_)
        arrset(pb.xnames,iv,"Y00606")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-25.9)
        ig = ig_["CAP01705"]
        pbm.A[ig,iv] += Float64(.03965)
        iv,ix_,_ = s2mpj_ii("Y00606",ix_)
        arrset(pb.xnames,iv,"Y00606")
        ig = ig_["CAP01605"]
        pbm.A[ig,iv] += Float64(.0344)
        iv,ix_,_ = s2mpj_ii("Y00703",ix_)
        arrset(pb.xnames,iv,"Y00703")
        ig = ig_["CAP00101"]
        pbm.A[ig,iv] += Float64(.00217)
        ig = ig_["CAP06402"]
        pbm.A[ig,iv] += Float64(.001)
        iv,ix_,_ = s2mpj_ii("Y00703",ix_)
        arrset(pb.xnames,iv,"Y00703")
        ig = ig_["MND00404"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00405"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00703",ix_)
        arrset(pb.xnames,iv,"Y00703")
        ig = ig_["MND00405"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00406"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00703",ix_)
        arrset(pb.xnames,iv,"Y00703")
        ig = ig_["CAP06403"]
        pbm.A[ig,iv] += Float64(.0009)
        ig = ig_["MXD00404"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00703",ix_)
        arrset(pb.xnames,iv,"Y00703")
        ig = ig_["MXD00406"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["INV00602"]
        pbm.A[ig,iv] += Float64(-1.24142)
        iv,ix_,_ = s2mpj_ii("Y00703",ix_)
        arrset(pb.xnames,iv,"Y00703")
        ig = ig_["MND00403"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP00202"]
        pbm.A[ig,iv] += Float64(.0111)
        iv,ix_,_ = s2mpj_ii("Y00703",ix_)
        arrset(pb.xnames,iv,"Y00703")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-87.07)
        ig = ig_["MXD00403"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00703",ix_)
        arrset(pb.xnames,iv,"Y00703")
        ig = ig_["INV00601"]
        pbm.A[ig,iv] += Float64(-.37243)
        ig = ig_["CAP00201"]
        pbm.A[ig,iv] += Float64(.00264)
        iv,ix_,_ = s2mpj_ii("Y00703",ix_)
        arrset(pb.xnames,iv,"Y00703")
        ig = ig_["CAP00102"]
        pbm.A[ig,iv] += Float64(.01192)
        iv,ix_,_ = s2mpj_ii("Y00704",ix_)
        arrset(pb.xnames,iv,"Y00704")
        ig = ig_["CAP00202"]
        pbm.A[ig,iv] += Float64(.00159)
        ig = ig_["MND00404"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00704",ix_)
        arrset(pb.xnames,iv,"Y00704")
        ig = ig_["MXD00405"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP00203"]
        pbm.A[ig,iv] += Float64(.01215)
        iv,ix_,_ = s2mpj_ii("Y00704",ix_)
        arrset(pb.xnames,iv,"Y00704")
        ig = ig_["MND00405"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP00103"]
        pbm.A[ig,iv] += Float64(.01301)
        iv,ix_,_ = s2mpj_ii("Y00704",ix_)
        arrset(pb.xnames,iv,"Y00704")
        ig = ig_["CAP06403"]
        pbm.A[ig,iv] += Float64(.001)
        ig = ig_["INV00602"]
        pbm.A[ig,iv] += Float64(-.24828)
        iv,ix_,_ = s2mpj_ii("Y00704",ix_)
        arrset(pb.xnames,iv,"Y00704")
        ig = ig_["CAP00102"]
        pbm.A[ig,iv] += Float64(.00108)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-82.66)
        iv,ix_,_ = s2mpj_ii("Y00704",ix_)
        arrset(pb.xnames,iv,"Y00704")
        ig = ig_["MND00406"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00404"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00704",ix_)
        arrset(pb.xnames,iv,"Y00704")
        ig = ig_["CAP06404"]
        pbm.A[ig,iv] += Float64(.0009)
        ig = ig_["MXD00406"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00704",ix_)
        arrset(pb.xnames,iv,"Y00704")
        ig = ig_["INV00603"]
        pbm.A[ig,iv] += Float64(-1.36556)
        iv,ix_,_ = s2mpj_ii("Y00705",ix_)
        arrset(pb.xnames,iv,"Y00705")
        ig = ig_["INV00604"]
        pbm.A[ig,iv] += Float64(-1.37476)
        ig = ig_["CAP00103"]
        pbm.A[ig,iv] += Float64(.00104)
        iv,ix_,_ = s2mpj_ii("Y00705",ix_)
        arrset(pb.xnames,iv,"Y00705")
        ig = ig_["MND00406"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-78.48)
        iv,ix_,_ = s2mpj_ii("Y00705",ix_)
        arrset(pb.xnames,iv,"Y00705")
        ig = ig_["INV00603"]
        pbm.A[ig,iv] += Float64(-.23909)
        ig = ig_["CAP06405"]
        pbm.A[ig,iv] += Float64(.00094)
        iv,ix_,_ = s2mpj_ii("Y00705",ix_)
        arrset(pb.xnames,iv,"Y00705")
        ig = ig_["CAP06404"]
        pbm.A[ig,iv] += Float64(.00097)
        ig = ig_["MND00405"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00705",ix_)
        arrset(pb.xnames,iv,"Y00705")
        ig = ig_["MXD00405"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP00104"]
        pbm.A[ig,iv] += Float64(.01305)
        iv,ix_,_ = s2mpj_ii("Y00705",ix_)
        arrset(pb.xnames,iv,"Y00705")
        ig = ig_["CAP00204"]
        pbm.A[ig,iv] += Float64(.01221)
        ig = ig_["CAP00203"]
        pbm.A[ig,iv] += Float64(.00153)
        iv,ix_,_ = s2mpj_ii("Y00705",ix_)
        arrset(pb.xnames,iv,"Y00705")
        ig = ig_["MXD00406"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00706",ix_)
        arrset(pb.xnames,iv,"Y00706")
        ig = ig_["CAP00106"]
        pbm.A[ig,iv] += Float64(.01409)
        ig = ig_["CAP00104"]
        pbm.A[ig,iv] += Float64(.00056)
        iv,ix_,_ = s2mpj_ii("Y00706",ix_)
        arrset(pb.xnames,iv,"Y00706")
        ig = ig_["MND00406"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-74.51)
        iv,ix_,_ = s2mpj_ii("Y00706",ix_)
        arrset(pb.xnames,iv,"Y00706")
        ig = ig_["CAP00105"]
        pbm.A[ig,iv] += Float64(.01522)
        ig = ig_["CAP06406"]
        pbm.A[ig,iv] += Float64(.00191)
        iv,ix_,_ = s2mpj_ii("Y00706",ix_)
        arrset(pb.xnames,iv,"Y00706")
        ig = ig_["CAP06405"]
        pbm.A[ig,iv] += Float64(.00104)
        ig = ig_["CAP00206"]
        pbm.A[ig,iv] += Float64(.01374)
        iv,ix_,_ = s2mpj_ii("Y00706",ix_)
        arrset(pb.xnames,iv,"Y00706")
        ig = ig_["INV00604"]
        pbm.A[ig,iv] += Float64(-.19366)
        ig = ig_["CAP00205"]
        pbm.A[ig,iv] += Float64(.01484)
        iv,ix_,_ = s2mpj_ii("Y00706",ix_)
        arrset(pb.xnames,iv,"Y00706")
        ig = ig_["INV00606"]
        pbm.A[ig,iv] += Float64(-1.61385)
        ig = ig_["CAP00204"]
        pbm.A[ig,iv] += Float64(.0011)
        iv,ix_,_ = s2mpj_ii("Y00706",ix_)
        arrset(pb.xnames,iv,"Y00706")
        ig = ig_["MXD00406"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["INV00605"]
        pbm.A[ig,iv] += Float64(-1.74296)
        iv,ix_,_ = s2mpj_ii("Y00802",ix_)
        arrset(pb.xnames,iv,"Y00802")
        ig = ig_["MXD00503"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP00602"]
        pbm.A[ig,iv] += Float64(.01774)
        iv,ix_,_ = s2mpj_ii("Y00802",ix_)
        arrset(pb.xnames,iv,"Y00802")
        ig = ig_["CAP06501"]
        pbm.A[ig,iv] += Float64(.00033)
        ig = ig_["CAP00201"]
        pbm.A[ig,iv] += Float64(.00241)
        iv,ix_,_ = s2mpj_ii("Y00802",ix_)
        arrset(pb.xnames,iv,"Y00802")
        ig = ig_["INV00102"]
        pbm.A[ig,iv] += Float64(-.64004)
        ig = ig_["MND00503"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00802",ix_)
        arrset(pb.xnames,iv,"Y00802")
        ig = ig_["CAP00202"]
        pbm.A[ig,iv] += Float64(.00401)
        ig = ig_["CAP00601"]
        pbm.A[ig,iv] += Float64(.01267)
        iv,ix_,_ = s2mpj_ii("Y00802",ix_)
        arrset(pb.xnames,iv,"Y00802")
        ig = ig_["MXD00505"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00506"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00802",ix_)
        arrset(pb.xnames,iv,"Y00802")
        ig = ig_["MND00504"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["INV00101"]
        pbm.A[ig,iv] += Float64(-.54157)
        iv,ix_,_ = s2mpj_ii("Y00802",ix_)
        arrset(pb.xnames,iv,"Y00802")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-4.06)
        ig = ig_["MXD00502"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00802",ix_)
        arrset(pb.xnames,iv,"Y00802")
        ig = ig_["MXD00504"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP06502"]
        pbm.A[ig,iv] += Float64(.0023)
        iv,ix_,_ = s2mpj_ii("Y00802",ix_)
        arrset(pb.xnames,iv,"Y00802")
        ig = ig_["MND00502"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00506"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00802",ix_)
        arrset(pb.xnames,iv,"Y00802")
        ig = ig_["MND00505"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00803",ix_)
        arrset(pb.xnames,iv,"Y00803")
        ig = ig_["MXD00506"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00506"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00803",ix_)
        arrset(pb.xnames,iv,"Y00803")
        ig = ig_["CAP06502"]
        pbm.A[ig,iv] += Float64(.0003)
        ig = ig_["MXD00504"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00803",ix_)
        arrset(pb.xnames,iv,"Y00803")
        ig = ig_["MND00504"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["INV00103"]
        pbm.A[ig,iv] += Float64(-.68169)
        iv,ix_,_ = s2mpj_ii("Y00803",ix_)
        arrset(pb.xnames,iv,"Y00803")
        ig = ig_["CAP00202"]
        pbm.A[ig,iv] += Float64(.00222)
        ig = ig_["MND00503"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00803",ix_)
        arrset(pb.xnames,iv,"Y00803")
        ig = ig_["INV00102"]
        pbm.A[ig,iv] += Float64(-.49991)
        ig = ig_["CAP00602"]
        pbm.A[ig,iv] += Float64(.0117)
        iv,ix_,_ = s2mpj_ii("Y00803",ix_)
        arrset(pb.xnames,iv,"Y00803")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-3.86)
        ig = ig_["MXD00503"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00803",ix_)
        arrset(pb.xnames,iv,"Y00803")
        ig = ig_["CAP06503"]
        pbm.A[ig,iv] += Float64(.00232)
        ig = ig_["CAP00203"]
        pbm.A[ig,iv] += Float64(.0042)
        iv,ix_,_ = s2mpj_ii("Y00803",ix_)
        arrset(pb.xnames,iv,"Y00803")
        ig = ig_["MXD00505"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP00603"]
        pbm.A[ig,iv] += Float64(.01872)
        iv,ix_,_ = s2mpj_ii("Y00803",ix_)
        arrset(pb.xnames,iv,"Y00803")
        ig = ig_["MND00505"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00804",ix_)
        arrset(pb.xnames,iv,"Y00804")
        ig = ig_["CAP00603"]
        pbm.A[ig,iv] += Float64(.0117)
        ig = ig_["MND00504"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00804",ix_)
        arrset(pb.xnames,iv,"Y00804")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-3.66)
        ig = ig_["MXD00504"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00804",ix_)
        arrset(pb.xnames,iv,"Y00804")
        ig = ig_["MND00506"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["INV00104"]
        pbm.A[ig,iv] += Float64(-.68169)
        iv,ix_,_ = s2mpj_ii("Y00804",ix_)
        arrset(pb.xnames,iv,"Y00804")
        ig = ig_["CAP00604"]
        pbm.A[ig,iv] += Float64(.01872)
        ig = ig_["CAP06503"]
        pbm.A[ig,iv] += Float64(.0003)
        iv,ix_,_ = s2mpj_ii("Y00804",ix_)
        arrset(pb.xnames,iv,"Y00804")
        ig = ig_["MND00505"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00505"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00804",ix_)
        arrset(pb.xnames,iv,"Y00804")
        ig = ig_["CAP06504"]
        pbm.A[ig,iv] += Float64(.00232)
        ig = ig_["CAP00204"]
        pbm.A[ig,iv] += Float64(.0042)
        iv,ix_,_ = s2mpj_ii("Y00804",ix_)
        arrset(pb.xnames,iv,"Y00804")
        ig = ig_["INV00103"]
        pbm.A[ig,iv] += Float64(-.49991)
        ig = ig_["MXD00506"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00804",ix_)
        arrset(pb.xnames,iv,"Y00804")
        ig = ig_["CAP00203"]
        pbm.A[ig,iv] += Float64(.00222)
        iv,ix_,_ = s2mpj_ii("Y00805",ix_)
        arrset(pb.xnames,iv,"Y00805")
        ig = ig_["MXD00505"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00506"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00805",ix_)
        arrset(pb.xnames,iv,"Y00805")
        ig = ig_["CAP00604"]
        pbm.A[ig,iv] += Float64(.01126)
        ig = ig_["CAP06504"]
        pbm.A[ig,iv] += Float64(.00029)
        iv,ix_,_ = s2mpj_ii("Y00805",ix_)
        arrset(pb.xnames,iv,"Y00805")
        ig = ig_["MND00505"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP06505"]
        pbm.A[ig,iv] += Float64(.00233)
        iv,ix_,_ = s2mpj_ii("Y00805",ix_)
        arrset(pb.xnames,iv,"Y00805")
        ig = ig_["INV00104"]
        pbm.A[ig,iv] += Float64(-.48139)
        ig = ig_["MND00506"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00805",ix_)
        arrset(pb.xnames,iv,"Y00805")
        ig = ig_["CAP00605"]
        pbm.A[ig,iv] += Float64(.01915)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-3.48)
        iv,ix_,_ = s2mpj_ii("Y00805",ix_)
        arrset(pb.xnames,iv,"Y00805")
        ig = ig_["CAP00204"]
        pbm.A[ig,iv] += Float64(.00214)
        ig = ig_["CAP00205"]
        pbm.A[ig,iv] += Float64(.00428)
        iv,ix_,_ = s2mpj_ii("Y00805",ix_)
        arrset(pb.xnames,iv,"Y00805")
        ig = ig_["INV00105"]
        pbm.A[ig,iv] += Float64(-.70021)
        iv,ix_,_ = s2mpj_ii("Y00806",ix_)
        arrset(pb.xnames,iv,"Y00806")
        ig = ig_["CAP00206"]
        pbm.A[ig,iv] += Float64(.00642)
        ig = ig_["MXD00506"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00806",ix_)
        arrset(pb.xnames,iv,"Y00806")
        ig = ig_["CAP06505"]
        pbm.A[ig,iv] += Float64(.00031)
        ig = ig_["MND00506"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00806",ix_)
        arrset(pb.xnames,iv,"Y00806")
        ig = ig_["INV00106"]
        pbm.A[ig,iv] += Float64(-1.1816)
        ig = ig_["CAP06506"]
        pbm.A[ig,iv] += Float64(.00262)
        iv,ix_,_ = s2mpj_ii("Y00806",ix_)
        arrset(pb.xnames,iv,"Y00806")
        ig = ig_["CAP00606"]
        pbm.A[ig,iv] += Float64(.03041)
        ig = ig_["INV00105"]
        pbm.A[ig,iv] += Float64(-.51991)
        iv,ix_,_ = s2mpj_ii("Y00806",ix_)
        arrset(pb.xnames,iv,"Y00806")
        ig = ig_["CAP00205"]
        pbm.A[ig,iv] += Float64(.00231)
        ig = ig_["CAP00605"]
        pbm.A[ig,iv] += Float64(.01217)
        iv,ix_,_ = s2mpj_ii("Y00806",ix_)
        arrset(pb.xnames,iv,"Y00806")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-3.3)
        iv,ix_,_ = s2mpj_ii("Y00902",ix_)
        arrset(pb.xnames,iv,"Y00902")
        ig = ig_["MND00502"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-4.06)
        iv,ix_,_ = s2mpj_ii("Y00902",ix_)
        arrset(pb.xnames,iv,"Y00902")
        ig = ig_["MXD00502"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP01101"]
        pbm.A[ig,iv] += Float64(.01267)
        iv,ix_,_ = s2mpj_ii("Y00902",ix_)
        arrset(pb.xnames,iv,"Y00902")
        ig = ig_["MND00506"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00506"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00902",ix_)
        arrset(pb.xnames,iv,"Y00902")
        ig = ig_["MND00505"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["INV00101"]
        pbm.A[ig,iv] += Float64(-.54157)
        iv,ix_,_ = s2mpj_ii("Y00902",ix_)
        arrset(pb.xnames,iv,"Y00902")
        ig = ig_["CAP06501"]
        pbm.A[ig,iv] += Float64(.00033)
        ig = ig_["INV00102"]
        pbm.A[ig,iv] += Float64(-.64004)
        iv,ix_,_ = s2mpj_ii("Y00902",ix_)
        arrset(pb.xnames,iv,"Y00902")
        ig = ig_["CAP01002"]
        pbm.A[ig,iv] += Float64(.00401)
        ig = ig_["CAP01102"]
        pbm.A[ig,iv] += Float64(.01774)
        iv,ix_,_ = s2mpj_ii("Y00902",ix_)
        arrset(pb.xnames,iv,"Y00902")
        ig = ig_["MND00504"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP01001"]
        pbm.A[ig,iv] += Float64(.00241)
        iv,ix_,_ = s2mpj_ii("Y00902",ix_)
        arrset(pb.xnames,iv,"Y00902")
        ig = ig_["MXD00504"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00505"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00902",ix_)
        arrset(pb.xnames,iv,"Y00902")
        ig = ig_["CAP06502"]
        pbm.A[ig,iv] += Float64(.0023)
        ig = ig_["MXD00503"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00902",ix_)
        arrset(pb.xnames,iv,"Y00902")
        ig = ig_["MND00503"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00903",ix_)
        arrset(pb.xnames,iv,"Y00903")
        ig = ig_["CAP06503"]
        pbm.A[ig,iv] += Float64(.00232)
        ig = ig_["MND00505"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00903",ix_)
        arrset(pb.xnames,iv,"Y00903")
        ig = ig_["MXD00505"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["INV00103"]
        pbm.A[ig,iv] += Float64(-.68169)
        iv,ix_,_ = s2mpj_ii("Y00903",ix_)
        arrset(pb.xnames,iv,"Y00903")
        ig = ig_["CAP06502"]
        pbm.A[ig,iv] += Float64(.0003)
        ig = ig_["MXD00506"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00903",ix_)
        arrset(pb.xnames,iv,"Y00903")
        ig = ig_["INV00102"]
        pbm.A[ig,iv] += Float64(-.49991)
        ig = ig_["CAP01103"]
        pbm.A[ig,iv] += Float64(.01872)
        iv,ix_,_ = s2mpj_ii("Y00903",ix_)
        arrset(pb.xnames,iv,"Y00903")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-3.86)
        ig = ig_["MND00506"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00903",ix_)
        arrset(pb.xnames,iv,"Y00903")
        ig = ig_["MND00503"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP01102"]
        pbm.A[ig,iv] += Float64(.0117)
        iv,ix_,_ = s2mpj_ii("Y00903",ix_)
        arrset(pb.xnames,iv,"Y00903")
        ig = ig_["CAP01003"]
        pbm.A[ig,iv] += Float64(.0042)
        ig = ig_["CAP01002"]
        pbm.A[ig,iv] += Float64(.00222)
        iv,ix_,_ = s2mpj_ii("Y00903",ix_)
        arrset(pb.xnames,iv,"Y00903")
        ig = ig_["MXD00503"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00504"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00903",ix_)
        arrset(pb.xnames,iv,"Y00903")
        ig = ig_["MXD00504"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00904",ix_)
        arrset(pb.xnames,iv,"Y00904")
        ig = ig_["MND00504"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00504"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00904",ix_)
        arrset(pb.xnames,iv,"Y00904")
        ig = ig_["INV00103"]
        pbm.A[ig,iv] += Float64(-.49991)
        ig = ig_["MND00505"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00904",ix_)
        arrset(pb.xnames,iv,"Y00904")
        ig = ig_["CAP06504"]
        pbm.A[ig,iv] += Float64(.00232)
        ig = ig_["MXD00505"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00904",ix_)
        arrset(pb.xnames,iv,"Y00904")
        ig = ig_["MXD00506"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-3.66)
        iv,ix_,_ = s2mpj_ii("Y00904",ix_)
        arrset(pb.xnames,iv,"Y00904")
        ig = ig_["CAP01004"]
        pbm.A[ig,iv] += Float64(.0042)
        ig = ig_["CAP01104"]
        pbm.A[ig,iv] += Float64(.01872)
        iv,ix_,_ = s2mpj_ii("Y00904",ix_)
        arrset(pb.xnames,iv,"Y00904")
        ig = ig_["CAP01003"]
        pbm.A[ig,iv] += Float64(.00222)
        ig = ig_["INV00104"]
        pbm.A[ig,iv] += Float64(-.68169)
        iv,ix_,_ = s2mpj_ii("Y00904",ix_)
        arrset(pb.xnames,iv,"Y00904")
        ig = ig_["CAP06503"]
        pbm.A[ig,iv] += Float64(.0003)
        ig = ig_["CAP01103"]
        pbm.A[ig,iv] += Float64(.0117)
        iv,ix_,_ = s2mpj_ii("Y00904",ix_)
        arrset(pb.xnames,iv,"Y00904")
        ig = ig_["MND00506"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00905",ix_)
        arrset(pb.xnames,iv,"Y00905")
        ig = ig_["INV00104"]
        pbm.A[ig,iv] += Float64(-.48139)
        ig = ig_["INV00105"]
        pbm.A[ig,iv] += Float64(-.70021)
        iv,ix_,_ = s2mpj_ii("Y00905",ix_)
        arrset(pb.xnames,iv,"Y00905")
        ig = ig_["MND00505"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00505"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00905",ix_)
        arrset(pb.xnames,iv,"Y00905")
        ig = ig_["CAP06505"]
        pbm.A[ig,iv] += Float64(.00233)
        ig = ig_["CAP01104"]
        pbm.A[ig,iv] += Float64(.01126)
        iv,ix_,_ = s2mpj_ii("Y00905",ix_)
        arrset(pb.xnames,iv,"Y00905")
        ig = ig_["CAP01004"]
        pbm.A[ig,iv] += Float64(.00214)
        ig = ig_["CAP01005"]
        pbm.A[ig,iv] += Float64(.00428)
        iv,ix_,_ = s2mpj_ii("Y00905",ix_)
        arrset(pb.xnames,iv,"Y00905")
        ig = ig_["MXD00506"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00506"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00905",ix_)
        arrset(pb.xnames,iv,"Y00905")
        ig = ig_["CAP06504"]
        pbm.A[ig,iv] += Float64(.00029)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-3.48)
        iv,ix_,_ = s2mpj_ii("Y00905",ix_)
        arrset(pb.xnames,iv,"Y00905")
        ig = ig_["CAP01105"]
        pbm.A[ig,iv] += Float64(.01915)
        iv,ix_,_ = s2mpj_ii("Y00906",ix_)
        arrset(pb.xnames,iv,"Y00906")
        ig = ig_["CAP01106"]
        pbm.A[ig,iv] += Float64(.03041)
        ig = ig_["CAP06505"]
        pbm.A[ig,iv] += Float64(.00031)
        iv,ix_,_ = s2mpj_ii("Y00906",ix_)
        arrset(pb.xnames,iv,"Y00906")
        ig = ig_["INV00106"]
        pbm.A[ig,iv] += Float64(-1.1816)
        ig = ig_["CAP01006"]
        pbm.A[ig,iv] += Float64(.00642)
        iv,ix_,_ = s2mpj_ii("Y00906",ix_)
        arrset(pb.xnames,iv,"Y00906")
        ig = ig_["MND00506"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-3.3)
        iv,ix_,_ = s2mpj_ii("Y00906",ix_)
        arrset(pb.xnames,iv,"Y00906")
        ig = ig_["CAP06506"]
        pbm.A[ig,iv] += Float64(.00262)
        ig = ig_["CAP01105"]
        pbm.A[ig,iv] += Float64(.01217)
        iv,ix_,_ = s2mpj_ii("Y00906",ix_)
        arrset(pb.xnames,iv,"Y00906")
        ig = ig_["MXD00506"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["INV00105"]
        pbm.A[ig,iv] += Float64(-.51991)
        iv,ix_,_ = s2mpj_ii("Y00906",ix_)
        arrset(pb.xnames,iv,"Y00906")
        ig = ig_["CAP01005"]
        pbm.A[ig,iv] += Float64(.00231)
        iv,ix_,_ = s2mpj_ii("Y01003",ix_)
        arrset(pb.xnames,iv,"Y01003")
        ig = ig_["MXD00606"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00606"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01003",ix_)
        arrset(pb.xnames,iv,"Y01003")
        ig = ig_["CAP06002"]
        pbm.A[ig,iv] += Float64(.00699)
        ig = ig_["MND00603"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01003",ix_)
        arrset(pb.xnames,iv,"Y01003")
        ig = ig_["CAP01002"]
        pbm.A[ig,iv] += Float64(.01383)
        ig = ig_["MXD00605"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01003",ix_)
        arrset(pb.xnames,iv,"Y01003")
        ig = ig_["INV00202"]
        pbm.A[ig,iv] += Float64(-1.2326)
        ig = ig_["INV00201"]
        pbm.A[ig,iv] += Float64(-.22411)
        iv,ix_,_ = s2mpj_ii("Y01003",ix_)
        arrset(pb.xnames,iv,"Y01003")
        ig = ig_["CAP01201"]
        pbm.A[ig,iv] += Float64(.0013)
        ig = ig_["MND00604"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01003",ix_)
        arrset(pb.xnames,iv,"Y01003")
        ig = ig_["MXD00603"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP01001"]
        pbm.A[ig,iv] += Float64(.00115)
        iv,ix_,_ = s2mpj_ii("Y01003",ix_)
        arrset(pb.xnames,iv,"Y01003")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-33.1)
        ig = ig_["MXD00604"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01003",ix_)
        arrset(pb.xnames,iv,"Y01003")
        ig = ig_["MND00605"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP01202"]
        pbm.A[ig,iv] += Float64(.00995)
        iv,ix_,_ = s2mpj_ii("Y01003",ix_)
        arrset(pb.xnames,iv,"Y01003")
        ig = ig_["CAP06003"]
        pbm.A[ig,iv] += Float64(.00384)
        iv,ix_,_ = s2mpj_ii("Y01004",ix_)
        arrset(pb.xnames,iv,"Y01004")
        ig = ig_["CAP01202"]
        pbm.A[ig,iv] += Float64(.00043)
        ig = ig_["CAP01003"]
        pbm.A[ig,iv] += Float64(.01498)
        iv,ix_,_ = s2mpj_ii("Y01004",ix_)
        arrset(pb.xnames,iv,"Y01004")
        ig = ig_["MXD00605"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00605"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01004",ix_)
        arrset(pb.xnames,iv,"Y01004")
        ig = ig_["MXD00606"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP06003"]
        pbm.A[ig,iv] += Float64(.00699)
        iv,ix_,_ = s2mpj_ii("Y01004",ix_)
        arrset(pb.xnames,iv,"Y01004")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-31.42)
        ig = ig_["MXD00604"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01004",ix_)
        arrset(pb.xnames,iv,"Y01004")
        ig = ig_["MND00606"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["INV00202"]
        pbm.A[ig,iv] += Float64(-.11205)
        iv,ix_,_ = s2mpj_ii("Y01004",ix_)
        arrset(pb.xnames,iv,"Y01004")
        ig = ig_["CAP01203"]
        pbm.A[ig,iv] += Float64(.01082)
        ig = ig_["MND00604"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01004",ix_)
        arrset(pb.xnames,iv,"Y01004")
        ig = ig_["CAP06004"]
        pbm.A[ig,iv] += Float64(.00384)
        ig = ig_["INV00203"]
        pbm.A[ig,iv] += Float64(-1.34466)
        iv,ix_,_ = s2mpj_ii("Y01005",ix_)
        arrset(pb.xnames,iv,"Y01005")
        ig = ig_["MXD00605"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP06004"]
        pbm.A[ig,iv] += Float64(.00674)
        iv,ix_,_ = s2mpj_ii("Y01005",ix_)
        arrset(pb.xnames,iv,"Y01005")
        ig = ig_["CAP06005"]
        pbm.A[ig,iv] += Float64(.0041)
        ig = ig_["CAP01203"]
        pbm.A[ig,iv] += Float64(.00042)
        iv,ix_,_ = s2mpj_ii("Y01005",ix_)
        arrset(pb.xnames,iv,"Y01005")
        ig = ig_["INV00204"]
        pbm.A[ig,iv] += Float64(-1.34881)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-29.83)
        iv,ix_,_ = s2mpj_ii("Y01005",ix_)
        arrset(pb.xnames,iv,"Y01005")
        ig = ig_["INV00203"]
        pbm.A[ig,iv] += Float64(-.1079)
        ig = ig_["CAP01204"]
        pbm.A[ig,iv] += Float64(.01083)
        iv,ix_,_ = s2mpj_ii("Y01005",ix_)
        arrset(pb.xnames,iv,"Y01005")
        ig = ig_["MND00605"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP01005"]
        pbm.A[ig,iv] += Float64(.00055)
        iv,ix_,_ = s2mpj_ii("Y01005",ix_)
        arrset(pb.xnames,iv,"Y01005")
        ig = ig_["MXD00606"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00606"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01005",ix_)
        arrset(pb.xnames,iv,"Y01005")
        ig = ig_["CAP01004"]
        pbm.A[ig,iv] += Float64(.01443)
        iv,ix_,_ = s2mpj_ii("Y01006",ix_)
        arrset(pb.xnames,iv,"Y01006")
        ig = ig_["MND00606"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00606"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01006",ix_)
        arrset(pb.xnames,iv,"Y01006")
        ig = ig_["INV00206"]
        pbm.A[ig,iv] += Float64(-1.45671)
        ig = ig_["INV00205"]
        pbm.A[ig,iv] += Float64(-1.57325)
        iv,ix_,_ = s2mpj_ii("Y01006",ix_)
        arrset(pb.xnames,iv,"Y01006")
        ig = ig_["CAP06005"]
        pbm.A[ig,iv] += Float64(.00727)
        ig = ig_["CAP01205"]
        pbm.A[ig,iv] += Float64(.01215)
        iv,ix_,_ = s2mpj_ii("Y01006",ix_)
        arrset(pb.xnames,iv,"Y01006")
        ig = ig_["CAP01005"]
        pbm.A[ig,iv] += Float64(.01558)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-28.32)
        iv,ix_,_ = s2mpj_ii("Y01006",ix_)
        arrset(pb.xnames,iv,"Y01006")
        ig = ig_["INV00204"]
        pbm.A[ig,iv] += Float64(-.05827)
        ig = ig_["CAP01006"]
        pbm.A[ig,iv] += Float64(.01498)
        iv,ix_,_ = s2mpj_ii("Y01006",ix_)
        arrset(pb.xnames,iv,"Y01006")
        ig = ig_["CAP01206"]
        pbm.A[ig,iv] += Float64(.01125)
        ig = ig_["CAP06006"]
        pbm.A[ig,iv] += Float64(.01084)
        iv,ix_,_ = s2mpj_ii("Y01103",ix_)
        arrset(pb.xnames,iv,"Y01103")
        ig = ig_["MXD00604"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00604"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01103",ix_)
        arrset(pb.xnames,iv,"Y01103")
        ig = ig_["CAP06002"]
        pbm.A[ig,iv] += Float64(.00699)
        ig = ig_["MND00605"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01103",ix_)
        arrset(pb.xnames,iv,"Y01103")
        ig = ig_["MXD00605"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP06003"]
        pbm.A[ig,iv] += Float64(.00384)
        iv,ix_,_ = s2mpj_ii("Y01103",ix_)
        arrset(pb.xnames,iv,"Y01103")
        ig = ig_["INV00202"]
        pbm.A[ig,iv] += Float64(-1.2326)
        ig = ig_["CAP00502"]
        pbm.A[ig,iv] += Float64(.01383)
        iv,ix_,_ = s2mpj_ii("Y01103",ix_)
        arrset(pb.xnames,iv,"Y01103")
        ig = ig_["CAP00501"]
        pbm.A[ig,iv] += Float64(.00115)
        ig = ig_["CAP00702"]
        pbm.A[ig,iv] += Float64(.00995)
        iv,ix_,_ = s2mpj_ii("Y01103",ix_)
        arrset(pb.xnames,iv,"Y01103")
        ig = ig_["MND00603"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00606"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01103",ix_)
        arrset(pb.xnames,iv,"Y01103")
        ig = ig_["CAP00701"]
        pbm.A[ig,iv] += Float64(.0013)
        ig = ig_["MXD00606"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01103",ix_)
        arrset(pb.xnames,iv,"Y01103")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-33.1)
        ig = ig_["MXD00603"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01103",ix_)
        arrset(pb.xnames,iv,"Y01103")
        ig = ig_["INV00201"]
        pbm.A[ig,iv] += Float64(-.22411)
        iv,ix_,_ = s2mpj_ii("Y01104",ix_)
        arrset(pb.xnames,iv,"Y01104")
        ig = ig_["CAP00702"]
        pbm.A[ig,iv] += Float64(.00043)
        ig = ig_["CAP00503"]
        pbm.A[ig,iv] += Float64(.01498)
        iv,ix_,_ = s2mpj_ii("Y01104",ix_)
        arrset(pb.xnames,iv,"Y01104")
        ig = ig_["CAP06003"]
        pbm.A[ig,iv] += Float64(.00699)
        ig = ig_["MND00605"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01104",ix_)
        arrset(pb.xnames,iv,"Y01104")
        ig = ig_["MXD00606"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["INV00203"]
        pbm.A[ig,iv] += Float64(-1.34466)
        iv,ix_,_ = s2mpj_ii("Y01104",ix_)
        arrset(pb.xnames,iv,"Y01104")
        ig = ig_["CAP06004"]
        pbm.A[ig,iv] += Float64(.00384)
        ig = ig_["MXD00604"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01104",ix_)
        arrset(pb.xnames,iv,"Y01104")
        ig = ig_["MND00604"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-31.42)
        iv,ix_,_ = s2mpj_ii("Y01104",ix_)
        arrset(pb.xnames,iv,"Y01104")
        ig = ig_["CAP00703"]
        pbm.A[ig,iv] += Float64(.01082)
        ig = ig_["MXD00605"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01104",ix_)
        arrset(pb.xnames,iv,"Y01104")
        ig = ig_["INV00202"]
        pbm.A[ig,iv] += Float64(-.11205)
        ig = ig_["MND00606"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01105",ix_)
        arrset(pb.xnames,iv,"Y01105")
        ig = ig_["INV00204"]
        pbm.A[ig,iv] += Float64(-1.34881)
        ig = ig_["MND00606"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01105",ix_)
        arrset(pb.xnames,iv,"Y01105")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-29.83)
        ig = ig_["CAP00504"]
        pbm.A[ig,iv] += Float64(.01443)
        iv,ix_,_ = s2mpj_ii("Y01105",ix_)
        arrset(pb.xnames,iv,"Y01105")
        ig = ig_["CAP00703"]
        pbm.A[ig,iv] += Float64(.00042)
        ig = ig_["MND00605"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01105",ix_)
        arrset(pb.xnames,iv,"Y01105")
        ig = ig_["CAP00704"]
        pbm.A[ig,iv] += Float64(.01083)
        ig = ig_["CAP00505"]
        pbm.A[ig,iv] += Float64(.00055)
        iv,ix_,_ = s2mpj_ii("Y01105",ix_)
        arrset(pb.xnames,iv,"Y01105")
        ig = ig_["CAP06004"]
        pbm.A[ig,iv] += Float64(.00674)
        ig = ig_["MXD00606"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01105",ix_)
        arrset(pb.xnames,iv,"Y01105")
        ig = ig_["CAP06005"]
        pbm.A[ig,iv] += Float64(.0041)
        ig = ig_["INV00203"]
        pbm.A[ig,iv] += Float64(-.1079)
        iv,ix_,_ = s2mpj_ii("Y01105",ix_)
        arrset(pb.xnames,iv,"Y01105")
        ig = ig_["MXD00605"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01106",ix_)
        arrset(pb.xnames,iv,"Y01106")
        ig = ig_["CAP00705"]
        pbm.A[ig,iv] += Float64(.01215)
        ig = ig_["MXD00606"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01106",ix_)
        arrset(pb.xnames,iv,"Y01106")
        ig = ig_["MND00606"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-28.32)
        iv,ix_,_ = s2mpj_ii("Y01106",ix_)
        arrset(pb.xnames,iv,"Y01106")
        ig = ig_["CAP06006"]
        pbm.A[ig,iv] += Float64(.01084)
        ig = ig_["CAP00505"]
        pbm.A[ig,iv] += Float64(.01558)
        iv,ix_,_ = s2mpj_ii("Y01106",ix_)
        arrset(pb.xnames,iv,"Y01106")
        ig = ig_["INV00206"]
        pbm.A[ig,iv] += Float64(-1.45671)
        ig = ig_["INV00205"]
        pbm.A[ig,iv] += Float64(-1.57325)
        iv,ix_,_ = s2mpj_ii("Y01106",ix_)
        arrset(pb.xnames,iv,"Y01106")
        ig = ig_["CAP00506"]
        pbm.A[ig,iv] += Float64(.01498)
        ig = ig_["CAP06005"]
        pbm.A[ig,iv] += Float64(.00727)
        iv,ix_,_ = s2mpj_ii("Y01106",ix_)
        arrset(pb.xnames,iv,"Y01106")
        ig = ig_["CAP00706"]
        pbm.A[ig,iv] += Float64(.01125)
        ig = ig_["INV00204"]
        pbm.A[ig,iv] += Float64(-.05827)
        iv,ix_,_ = s2mpj_ii("Y01202",ix_)
        arrset(pb.xnames,iv,"Y01202")
        ig = ig_["MND00705"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00703"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01202",ix_)
        arrset(pb.xnames,iv,"Y01202")
        ig = ig_["INV00402"]
        pbm.A[ig,iv] += Float64(-.86606)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-5.6)
        iv,ix_,_ = s2mpj_ii("Y01202",ix_)
        arrset(pb.xnames,iv,"Y01202")
        ig = ig_["CAP06301"]
        pbm.A[ig,iv] += Float64(.00002)
        ig = ig_["CAP06302"]
        pbm.A[ig,iv] += Float64(.0004)
        iv,ix_,_ = s2mpj_ii("Y01202",ix_)
        arrset(pb.xnames,iv,"Y01202")
        ig = ig_["MXD00702"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00706"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01202",ix_)
        arrset(pb.xnames,iv,"Y01202")
        ig = ig_["MND00702"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00705"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01202",ix_)
        arrset(pb.xnames,iv,"Y01202")
        ig = ig_["MND00703"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00704"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01202",ix_)
        arrset(pb.xnames,iv,"Y01202")
        ig = ig_["MXD00704"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP01201"]
        pbm.A[ig,iv] += Float64(.00186)
        iv,ix_,_ = s2mpj_ii("Y01202",ix_)
        arrset(pb.xnames,iv,"Y01202")
        ig = ig_["CAP01202"]
        pbm.A[ig,iv] += Float64(.00706)
        ig = ig_["CAP01002"]
        pbm.A[ig,iv] += Float64(.0099)
        iv,ix_,_ = s2mpj_ii("Y01202",ix_)
        arrset(pb.xnames,iv,"Y01202")
        ig = ig_["MXD00706"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["INV00401"]
        pbm.A[ig,iv] += Float64(-.28869)
        iv,ix_,_ = s2mpj_ii("Y01202",ix_)
        arrset(pb.xnames,iv,"Y01202")
        ig = ig_["CAP01001"]
        pbm.A[ig,iv] += Float64(.00198)
        iv,ix_,_ = s2mpj_ii("Y01203",ix_)
        arrset(pb.xnames,iv,"Y01203")
        ig = ig_["MXD00706"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00706"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01203",ix_)
        arrset(pb.xnames,iv,"Y01203")
        ig = ig_["CAP01203"]
        pbm.A[ig,iv] += Float64(.0072)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-5.32)
        iv,ix_,_ = s2mpj_ii("Y01203",ix_)
        arrset(pb.xnames,iv,"Y01203")
        ig = ig_["MXD00704"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00704"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01203",ix_)
        arrset(pb.xnames,iv,"Y01203")
        ig = ig_["MND00703"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP01003"]
        pbm.A[ig,iv] += Float64(.01005)
        iv,ix_,_ = s2mpj_ii("Y01203",ix_)
        arrset(pb.xnames,iv,"Y01203")
        ig = ig_["MXD00705"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00703"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01203",ix_)
        arrset(pb.xnames,iv,"Y01203")
        ig = ig_["INV00403"]
        pbm.A[ig,iv] += Float64(-.88827)
        ig = ig_["MND00705"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01203",ix_)
        arrset(pb.xnames,iv,"Y01203")
        ig = ig_["CAP01002"]
        pbm.A[ig,iv] += Float64(.00183)
        ig = ig_["INV00402"]
        pbm.A[ig,iv] += Float64(-.26648)
        iv,ix_,_ = s2mpj_ii("Y01203",ix_)
        arrset(pb.xnames,iv,"Y01203")
        ig = ig_["CAP06302"]
        pbm.A[ig,iv] += Float64(.00002)
        ig = ig_["CAP06303"]
        pbm.A[ig,iv] += Float64(.0004)
        iv,ix_,_ = s2mpj_ii("Y01203",ix_)
        arrset(pb.xnames,iv,"Y01203")
        ig = ig_["CAP01202"]
        pbm.A[ig,iv] += Float64(.00171)
        iv,ix_,_ = s2mpj_ii("Y01204",ix_)
        arrset(pb.xnames,iv,"Y01204")
        ig = ig_["CAP01204"]
        pbm.A[ig,iv] += Float64(.0072)
        ig = ig_["INV00403"]
        pbm.A[ig,iv] += Float64(-.26648)
        iv,ix_,_ = s2mpj_ii("Y01204",ix_)
        arrset(pb.xnames,iv,"Y01204")
        ig = ig_["MXD00706"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP01004"]
        pbm.A[ig,iv] += Float64(.01005)
        iv,ix_,_ = s2mpj_ii("Y01204",ix_)
        arrset(pb.xnames,iv,"Y01204")
        ig = ig_["MND00706"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-5.05)
        iv,ix_,_ = s2mpj_ii("Y01204",ix_)
        arrset(pb.xnames,iv,"Y01204")
        ig = ig_["CAP06303"]
        pbm.A[ig,iv] += Float64(.00002)
        ig = ig_["CAP01203"]
        pbm.A[ig,iv] += Float64(.00171)
        iv,ix_,_ = s2mpj_ii("Y01204",ix_)
        arrset(pb.xnames,iv,"Y01204")
        ig = ig_["MXD00704"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP01003"]
        pbm.A[ig,iv] += Float64(.00183)
        iv,ix_,_ = s2mpj_ii("Y01204",ix_)
        arrset(pb.xnames,iv,"Y01204")
        ig = ig_["MXD00705"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00705"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01204",ix_)
        arrset(pb.xnames,iv,"Y01204")
        ig = ig_["CAP06304"]
        pbm.A[ig,iv] += Float64(.0004)
        ig = ig_["MND00704"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01204",ix_)
        arrset(pb.xnames,iv,"Y01204")
        ig = ig_["INV00404"]
        pbm.A[ig,iv] += Float64(-.88827)
        iv,ix_,_ = s2mpj_ii("Y01205",ix_)
        arrset(pb.xnames,iv,"Y01205")
        ig = ig_["CAP06305"]
        pbm.A[ig,iv] += Float64(.0004)
        ig = ig_["CAP01005"]
        pbm.A[ig,iv] += Float64(.01012)
        iv,ix_,_ = s2mpj_ii("Y01205",ix_)
        arrset(pb.xnames,iv,"Y01205")
        ig = ig_["CAP01205"]
        pbm.A[ig,iv] += Float64(.00727)
        ig = ig_["INV00405"]
        pbm.A[ig,iv] += Float64(-.89814)
        iv,ix_,_ = s2mpj_ii("Y01205",ix_)
        arrset(pb.xnames,iv,"Y01205")
        ig = ig_["MND00705"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP06304"]
        pbm.A[ig,iv] += Float64(.00002)
        iv,ix_,_ = s2mpj_ii("Y01205",ix_)
        arrset(pb.xnames,iv,"Y01205")
        ig = ig_["CAP01204"]
        pbm.A[ig,iv] += Float64(.00165)
        ig = ig_["MXD00705"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01205",ix_)
        arrset(pb.xnames,iv,"Y01205")
        ig = ig_["CAP01004"]
        pbm.A[ig,iv] += Float64(.00176)
        ig = ig_["MND00706"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01205",ix_)
        arrset(pb.xnames,iv,"Y01205")
        ig = ig_["MXD00706"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-4.79)
        iv,ix_,_ = s2mpj_ii("Y01205",ix_)
        arrset(pb.xnames,iv,"Y01205")
        ig = ig_["INV00404"]
        pbm.A[ig,iv] += Float64(-.25661)
        iv,ix_,_ = s2mpj_ii("Y01206",ix_)
        arrset(pb.xnames,iv,"Y01206")
        ig = ig_["CAP01205"]
        pbm.A[ig,iv] += Float64(.00178)
        ig = ig_["CAP01005"]
        pbm.A[ig,iv] += Float64(.0019)
        iv,ix_,_ = s2mpj_ii("Y01206",ix_)
        arrset(pb.xnames,iv,"Y01206")
        ig = ig_["CAP06306"]
        pbm.A[ig,iv] += Float64(.00042)
        ig = ig_["INV00405"]
        pbm.A[ig,iv] += Float64(-.27714)
        iv,ix_,_ = s2mpj_ii("Y01206",ix_)
        arrset(pb.xnames,iv,"Y01206")
        ig = ig_["CAP06305"]
        pbm.A[ig,iv] += Float64(.00002)
        ig = ig_["CAP01206"]
        pbm.A[ig,iv] += Float64(.00892)
        iv,ix_,_ = s2mpj_ii("Y01206",ix_)
        arrset(pb.xnames,iv,"Y01206")
        ig = ig_["CAP01006"]
        pbm.A[ig,iv] += Float64(.01188)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-4.55)
        iv,ix_,_ = s2mpj_ii("Y01206",ix_)
        arrset(pb.xnames,iv,"Y01206")
        ig = ig_["INV00406"]
        pbm.A[ig,iv] += Float64(-1.15475)
        ig = ig_["MND00706"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01206",ix_)
        arrset(pb.xnames,iv,"Y01206")
        ig = ig_["MXD00706"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01302",ix_)
        arrset(pb.xnames,iv,"Y01302")
        ig = ig_["MND00706"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP00702"]
        pbm.A[ig,iv] += Float64(.00706)
        iv,ix_,_ = s2mpj_ii("Y01302",ix_)
        arrset(pb.xnames,iv,"Y01302")
        ig = ig_["CAP06301"]
        pbm.A[ig,iv] += Float64(.00002)
        ig = ig_["CAP00501"]
        pbm.A[ig,iv] += Float64(.00198)
        iv,ix_,_ = s2mpj_ii("Y01302",ix_)
        arrset(pb.xnames,iv,"Y01302")
        ig = ig_["CAP00502"]
        pbm.A[ig,iv] += Float64(.0099)
        ig = ig_["CAP06302"]
        pbm.A[ig,iv] += Float64(.0004)
        iv,ix_,_ = s2mpj_ii("Y01302",ix_)
        arrset(pb.xnames,iv,"Y01302")
        ig = ig_["INV00401"]
        pbm.A[ig,iv] += Float64(-.28869)
        ig = ig_["INV00402"]
        pbm.A[ig,iv] += Float64(-.86606)
        iv,ix_,_ = s2mpj_ii("Y01302",ix_)
        arrset(pb.xnames,iv,"Y01302")
        ig = ig_["MXD00702"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00702"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01302",ix_)
        arrset(pb.xnames,iv,"Y01302")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-5.62)
        ig = ig_["MND00704"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01302",ix_)
        arrset(pb.xnames,iv,"Y01302")
        ig = ig_["MXD00706"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00703"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01302",ix_)
        arrset(pb.xnames,iv,"Y01302")
        ig = ig_["MXD00703"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00705"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01302",ix_)
        arrset(pb.xnames,iv,"Y01302")
        ig = ig_["CAP00701"]
        pbm.A[ig,iv] += Float64(.00186)
        ig = ig_["MXD00704"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01302",ix_)
        arrset(pb.xnames,iv,"Y01302")
        ig = ig_["MXD00705"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01303",ix_)
        arrset(pb.xnames,iv,"Y01303")
        ig = ig_["MXD00706"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP06302"]
        pbm.A[ig,iv] += Float64(.00002)
        iv,ix_,_ = s2mpj_ii("Y01303",ix_)
        arrset(pb.xnames,iv,"Y01303")
        ig = ig_["CAP00503"]
        pbm.A[ig,iv] += Float64(.01005)
        ig = ig_["MND00704"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01303",ix_)
        arrset(pb.xnames,iv,"Y01303")
        ig = ig_["MXD00704"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["INV00403"]
        pbm.A[ig,iv] += Float64(-.88827)
        iv,ix_,_ = s2mpj_ii("Y01303",ix_)
        arrset(pb.xnames,iv,"Y01303")
        ig = ig_["CAP00703"]
        pbm.A[ig,iv] += Float64(.0072)
        ig = ig_["INV00402"]
        pbm.A[ig,iv] += Float64(-.26648)
        iv,ix_,_ = s2mpj_ii("Y01303",ix_)
        arrset(pb.xnames,iv,"Y01303")
        ig = ig_["MND00706"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00705"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01303",ix_)
        arrset(pb.xnames,iv,"Y01303")
        ig = ig_["CAP00502"]
        pbm.A[ig,iv] += Float64(.00183)
        ig = ig_["CAP06303"]
        pbm.A[ig,iv] += Float64(.0004)
        iv,ix_,_ = s2mpj_ii("Y01303",ix_)
        arrset(pb.xnames,iv,"Y01303")
        ig = ig_["MND00705"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-5.34)
        iv,ix_,_ = s2mpj_ii("Y01303",ix_)
        arrset(pb.xnames,iv,"Y01303")
        ig = ig_["MXD00703"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00703"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01303",ix_)
        arrset(pb.xnames,iv,"Y01303")
        ig = ig_["CAP00702"]
        pbm.A[ig,iv] += Float64(.00171)
        iv,ix_,_ = s2mpj_ii("Y01304",ix_)
        arrset(pb.xnames,iv,"Y01304")
        ig = ig_["MXD00706"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00705"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01304",ix_)
        arrset(pb.xnames,iv,"Y01304")
        ig = ig_["MND00706"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00705"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01304",ix_)
        arrset(pb.xnames,iv,"Y01304")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-5.07)
        ig = ig_["MND00704"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01304",ix_)
        arrset(pb.xnames,iv,"Y01304")
        ig = ig_["CAP06304"]
        pbm.A[ig,iv] += Float64(.0004)
        ig = ig_["CAP00504"]
        pbm.A[ig,iv] += Float64(.01005)
        iv,ix_,_ = s2mpj_ii("Y01304",ix_)
        arrset(pb.xnames,iv,"Y01304")
        ig = ig_["INV00403"]
        pbm.A[ig,iv] += Float64(-.26648)
        ig = ig_["MXD00704"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01304",ix_)
        arrset(pb.xnames,iv,"Y01304")
        ig = ig_["CAP00704"]
        pbm.A[ig,iv] += Float64(.0072)
        ig = ig_["INV00404"]
        pbm.A[ig,iv] += Float64(-.88827)
        iv,ix_,_ = s2mpj_ii("Y01304",ix_)
        arrset(pb.xnames,iv,"Y01304")
        ig = ig_["CAP06303"]
        pbm.A[ig,iv] += Float64(.00002)
        ig = ig_["CAP00703"]
        pbm.A[ig,iv] += Float64(.00171)
        iv,ix_,_ = s2mpj_ii("Y01304",ix_)
        arrset(pb.xnames,iv,"Y01304")
        ig = ig_["CAP00503"]
        pbm.A[ig,iv] += Float64(.00183)
        iv,ix_,_ = s2mpj_ii("Y01305",ix_)
        arrset(pb.xnames,iv,"Y01305")
        ig = ig_["CAP06304"]
        pbm.A[ig,iv] += Float64(.00002)
        ig = ig_["MND00706"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01305",ix_)
        arrset(pb.xnames,iv,"Y01305")
        ig = ig_["CAP00704"]
        pbm.A[ig,iv] += Float64(.00165)
        ig = ig_["MXD00705"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01305",ix_)
        arrset(pb.xnames,iv,"Y01305")
        ig = ig_["CAP00504"]
        pbm.A[ig,iv] += Float64(.00176)
        ig = ig_["CAP00705"]
        pbm.A[ig,iv] += Float64(.00727)
        iv,ix_,_ = s2mpj_ii("Y01305",ix_)
        arrset(pb.xnames,iv,"Y01305")
        ig = ig_["MND00705"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00706"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01305",ix_)
        arrset(pb.xnames,iv,"Y01305")
        ig = ig_["CAP00505"]
        pbm.A[ig,iv] += Float64(.01012)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-4.81)
        iv,ix_,_ = s2mpj_ii("Y01305",ix_)
        arrset(pb.xnames,iv,"Y01305")
        ig = ig_["CAP06305"]
        pbm.A[ig,iv] += Float64(.0004)
        ig = ig_["INV00405"]
        pbm.A[ig,iv] += Float64(-.89814)
        iv,ix_,_ = s2mpj_ii("Y01305",ix_)
        arrset(pb.xnames,iv,"Y01305")
        ig = ig_["INV00404"]
        pbm.A[ig,iv] += Float64(-.25661)
        iv,ix_,_ = s2mpj_ii("Y01306",ix_)
        arrset(pb.xnames,iv,"Y01306")
        ig = ig_["CAP00705"]
        pbm.A[ig,iv] += Float64(.00178)
        ig = ig_["MND00706"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01306",ix_)
        arrset(pb.xnames,iv,"Y01306")
        ig = ig_["MXD00706"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-4.57)
        iv,ix_,_ = s2mpj_ii("Y01306",ix_)
        arrset(pb.xnames,iv,"Y01306")
        ig = ig_["CAP00506"]
        pbm.A[ig,iv] += Float64(.01188)
        ig = ig_["CAP00706"]
        pbm.A[ig,iv] += Float64(.00892)
        iv,ix_,_ = s2mpj_ii("Y01306",ix_)
        arrset(pb.xnames,iv,"Y01306")
        ig = ig_["CAP06305"]
        pbm.A[ig,iv] += Float64(.00002)
        ig = ig_["CAP00505"]
        pbm.A[ig,iv] += Float64(.0019)
        iv,ix_,_ = s2mpj_ii("Y01306",ix_)
        arrset(pb.xnames,iv,"Y01306")
        ig = ig_["CAP06306"]
        pbm.A[ig,iv] += Float64(.00042)
        ig = ig_["INV00406"]
        pbm.A[ig,iv] += Float64(-1.15475)
        iv,ix_,_ = s2mpj_ii("Y01306",ix_)
        arrset(pb.xnames,iv,"Y01306")
        ig = ig_["INV00405"]
        pbm.A[ig,iv] += Float64(-.27714)
        iv,ix_,_ = s2mpj_ii("Y01402",ix_)
        arrset(pb.xnames,iv,"Y01402")
        ig = ig_["MND00803"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00802"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01402",ix_)
        arrset(pb.xnames,iv,"Y01402")
        ig = ig_["MND00802"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00803"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01402",ix_)
        arrset(pb.xnames,iv,"Y01402")
        ig = ig_["CAP01001"]
        pbm.A[ig,iv] += Float64(.00412)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-5.46)
        iv,ix_,_ = s2mpj_ii("Y01402",ix_)
        arrset(pb.xnames,iv,"Y01402")
        ig = ig_["MXD00805"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["INV00502"]
        pbm.A[ig,iv] += Float64(-.7014)
        iv,ix_,_ = s2mpj_ii("Y01402",ix_)
        arrset(pb.xnames,iv,"Y01402")
        ig = ig_["CAP06201"]
        pbm.A[ig,iv] += Float64(.00083)
        ig = ig_["CAP06202"]
        pbm.A[ig,iv] += Float64(.00442)
        iv,ix_,_ = s2mpj_ii("Y01402",ix_)
        arrset(pb.xnames,iv,"Y01402")
        ig = ig_["MND00805"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP06102"]
        pbm.A[ig,iv] += Float64(.00066)
        iv,ix_,_ = s2mpj_ii("Y01402",ix_)
        arrset(pb.xnames,iv,"Y01402")
        ig = ig_["INV00501"]
        pbm.A[ig,iv] += Float64(-.501)
        ig = ig_["CAP01002"]
        pbm.A[ig,iv] += Float64(.00824)
        iv,ix_,_ = s2mpj_ii("Y01402",ix_)
        arrset(pb.xnames,iv,"Y01402")
        ig = ig_["CAP01202"]
        pbm.A[ig,iv] += Float64(.0058)
        ig = ig_["MXD00806"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01402",ix_)
        arrset(pb.xnames,iv,"Y01402")
        ig = ig_["CAP01201"]
        pbm.A[ig,iv] += Float64(.00348)
        ig = ig_["CAP06101"]
        pbm.A[ig,iv] += Float64(.00006)
        iv,ix_,_ = s2mpj_ii("Y01402",ix_)
        arrset(pb.xnames,iv,"Y01402")
        ig = ig_["MND00806"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00804"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01402",ix_)
        arrset(pb.xnames,iv,"Y01402")
        ig = ig_["MXD00804"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01403",ix_)
        arrset(pb.xnames,iv,"Y01403")
        ig = ig_["CAP01003"]
        pbm.A[ig,iv] += Float64(.00856)
        ig = ig_["MXD00805"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01403",ix_)
        arrset(pb.xnames,iv,"Y01403")
        ig = ig_["MND00805"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00804"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01403",ix_)
        arrset(pb.xnames,iv,"Y01403")
        ig = ig_["MND00806"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-5.18)
        iv,ix_,_ = s2mpj_ii("Y01403",ix_)
        arrset(pb.xnames,iv,"Y01403")
        ig = ig_["CAP06202"]
        pbm.A[ig,iv] += Float64(.00077)
        ig = ig_["CAP01203"]
        pbm.A[ig,iv] += Float64(.00607)
        iv,ix_,_ = s2mpj_ii("Y01403",ix_)
        arrset(pb.xnames,iv,"Y01403")
        ig = ig_["CAP06102"]
        pbm.A[ig,iv] += Float64(.00006)
        ig = ig_["MXD00806"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01403",ix_)
        arrset(pb.xnames,iv,"Y01403")
        ig = ig_["MND00803"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00804"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01403",ix_)
        arrset(pb.xnames,iv,"Y01403")
        ig = ig_["CAP06203"]
        pbm.A[ig,iv] += Float64(.00449)
        ig = ig_["INV00502"]
        pbm.A[ig,iv] += Float64(-.46246)
        iv,ix_,_ = s2mpj_ii("Y01403",ix_)
        arrset(pb.xnames,iv,"Y01403")
        ig = ig_["MXD00803"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP06103"]
        pbm.A[ig,iv] += Float64(.00067)
        iv,ix_,_ = s2mpj_ii("Y01403",ix_)
        arrset(pb.xnames,iv,"Y01403")
        ig = ig_["CAP01202"]
        pbm.A[ig,iv] += Float64(.00321)
        ig = ig_["CAP01002"]
        pbm.A[ig,iv] += Float64(.0038)
        iv,ix_,_ = s2mpj_ii("Y01403",ix_)
        arrset(pb.xnames,iv,"Y01403")
        ig = ig_["INV00503"]
        pbm.A[ig,iv] += Float64(-.73994)
        iv,ix_,_ = s2mpj_ii("Y01404",ix_)
        arrset(pb.xnames,iv,"Y01404")
        ig = ig_["MND00806"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP06204"]
        pbm.A[ig,iv] += Float64(.00449)
        iv,ix_,_ = s2mpj_ii("Y01404",ix_)
        arrset(pb.xnames,iv,"Y01404")
        ig = ig_["MXD00806"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP01203"]
        pbm.A[ig,iv] += Float64(.00321)
        iv,ix_,_ = s2mpj_ii("Y01404",ix_)
        arrset(pb.xnames,iv,"Y01404")
        ig = ig_["CAP06104"]
        pbm.A[ig,iv] += Float64(.00067)
        ig = ig_["CAP06103"]
        pbm.A[ig,iv] += Float64(.00006)
        iv,ix_,_ = s2mpj_ii("Y01404",ix_)
        arrset(pb.xnames,iv,"Y01404")
        ig = ig_["CAP01003"]
        pbm.A[ig,iv] += Float64(.0038)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-4.92)
        iv,ix_,_ = s2mpj_ii("Y01404",ix_)
        arrset(pb.xnames,iv,"Y01404")
        ig = ig_["CAP06203"]
        pbm.A[ig,iv] += Float64(.00077)
        ig = ig_["INV00503"]
        pbm.A[ig,iv] += Float64(-.46246)
        iv,ix_,_ = s2mpj_ii("Y01404",ix_)
        arrset(pb.xnames,iv,"Y01404")
        ig = ig_["MND00805"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["INV00504"]
        pbm.A[ig,iv] += Float64(-.73994)
        iv,ix_,_ = s2mpj_ii("Y01404",ix_)
        arrset(pb.xnames,iv,"Y01404")
        ig = ig_["MXD00805"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00804"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01404",ix_)
        arrset(pb.xnames,iv,"Y01404")
        ig = ig_["MXD00804"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP01004"]
        pbm.A[ig,iv] += Float64(.00856)
        iv,ix_,_ = s2mpj_ii("Y01404",ix_)
        arrset(pb.xnames,iv,"Y01404")
        ig = ig_["CAP01204"]
        pbm.A[ig,iv] += Float64(.00607)
        iv,ix_,_ = s2mpj_ii("Y01405",ix_)
        arrset(pb.xnames,iv,"Y01405")
        ig = ig_["INV00504"]
        pbm.A[ig,iv] += Float64(-.44534)
        ig = ig_["CAP01205"]
        pbm.A[ig,iv] += Float64(.00619)
        iv,ix_,_ = s2mpj_ii("Y01405",ix_)
        arrset(pb.xnames,iv,"Y01405")
        ig = ig_["CAP01005"]
        pbm.A[ig,iv] += Float64(.0087)
        ig = ig_["INV00505"]
        pbm.A[ig,iv] += Float64(-.75707)
        iv,ix_,_ = s2mpj_ii("Y01405",ix_)
        arrset(pb.xnames,iv,"Y01405")
        ig = ig_["CAP06105"]
        pbm.A[ig,iv] += Float64(.00067)
        ig = ig_["CAP06204"]
        pbm.A[ig,iv] += Float64(.00074)
        iv,ix_,_ = s2mpj_ii("Y01405",ix_)
        arrset(pb.xnames,iv,"Y01405")
        ig = ig_["MND00805"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00805"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01405",ix_)
        arrset(pb.xnames,iv,"Y01405")
        ig = ig_["CAP01004"]
        pbm.A[ig,iv] += Float64(.00366)
        ig = ig_["CAP01204"]
        pbm.A[ig,iv] += Float64(.0031)
        iv,ix_,_ = s2mpj_ii("Y01405",ix_)
        arrset(pb.xnames,iv,"Y01405")
        ig = ig_["CAP06104"]
        pbm.A[ig,iv] += Float64(.00005)
        ig = ig_["MND00806"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01405",ix_)
        arrset(pb.xnames,iv,"Y01405")
        ig = ig_["CAP06205"]
        pbm.A[ig,iv] += Float64(.00452)
        ig = ig_["MXD00806"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01405",ix_)
        arrset(pb.xnames,iv,"Y01405")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-4.67)
        iv,ix_,_ = s2mpj_ii("Y01406",ix_)
        arrset(pb.xnames,iv,"Y01406")
        ig = ig_["CAP01206"]
        pbm.A[ig,iv] += Float64(.00929)
        ig = ig_["CAP06105"]
        pbm.A[ig,iv] += Float64(.00006)
        iv,ix_,_ = s2mpj_ii("Y01406",ix_)
        arrset(pb.xnames,iv,"Y01406")
        ig = ig_["CAP01006"]
        pbm.A[ig,iv] += Float64(.01237)
        ig = ig_["CAP06205"]
        pbm.A[ig,iv] += Float64(.0008)
        iv,ix_,_ = s2mpj_ii("Y01406",ix_)
        arrset(pb.xnames,iv,"Y01406")
        ig = ig_["INV00505"]
        pbm.A[ig,iv] += Float64(-.48096)
        ig = ig_["MXD00806"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01406",ix_)
        arrset(pb.xnames,iv,"Y01406")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-4.44)
        ig = ig_["INV00506"]
        pbm.A[ig,iv] += Float64(-1.20241)
        iv,ix_,_ = s2mpj_ii("Y01406",ix_)
        arrset(pb.xnames,iv,"Y01406")
        ig = ig_["CAP06106"]
        pbm.A[ig,iv] += Float64(.00072)
        ig = ig_["MND00806"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01406",ix_)
        arrset(pb.xnames,iv,"Y01406")
        ig = ig_["CAP06206"]
        pbm.A[ig,iv] += Float64(.00525)
        ig = ig_["CAP01005"]
        pbm.A[ig,iv] += Float64(.00396)
        iv,ix_,_ = s2mpj_ii("Y01406",ix_)
        arrset(pb.xnames,iv,"Y01406")
        ig = ig_["CAP01205"]
        pbm.A[ig,iv] += Float64(.00334)
        iv,ix_,_ = s2mpj_ii("Y01502",ix_)
        arrset(pb.xnames,iv,"Y01502")
        ig = ig_["CAP06102"]
        pbm.A[ig,iv] += Float64(.00066)
        ig = ig_["CAP06101"]
        pbm.A[ig,iv] += Float64(.00006)
        iv,ix_,_ = s2mpj_ii("Y01502",ix_)
        arrset(pb.xnames,iv,"Y01502")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-5.46)
        ig = ig_["INV00501"]
        pbm.A[ig,iv] += Float64(-.501)
        iv,ix_,_ = s2mpj_ii("Y01502",ix_)
        arrset(pb.xnames,iv,"Y01502")
        ig = ig_["MXD00806"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP00701"]
        pbm.A[ig,iv] += Float64(.00348)
        iv,ix_,_ = s2mpj_ii("Y01502",ix_)
        arrset(pb.xnames,iv,"Y01502")
        ig = ig_["CAP06201"]
        pbm.A[ig,iv] += Float64(.00083)
        ig = ig_["CAP00501"]
        pbm.A[ig,iv] += Float64(.00412)
        iv,ix_,_ = s2mpj_ii("Y01502",ix_)
        arrset(pb.xnames,iv,"Y01502")
        ig = ig_["CAP06202"]
        pbm.A[ig,iv] += Float64(.00442)
        ig = ig_["MND00806"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01502",ix_)
        arrset(pb.xnames,iv,"Y01502")
        ig = ig_["CAP00502"]
        pbm.A[ig,iv] += Float64(.00824)
        ig = ig_["INV00502"]
        pbm.A[ig,iv] += Float64(-.7014)
        iv,ix_,_ = s2mpj_ii("Y01502",ix_)
        arrset(pb.xnames,iv,"Y01502")
        ig = ig_["MND00805"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00803"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01502",ix_)
        arrset(pb.xnames,iv,"Y01502")
        ig = ig_["MXD00803"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00804"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01502",ix_)
        arrset(pb.xnames,iv,"Y01502")
        ig = ig_["MXD00804"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00805"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01502",ix_)
        arrset(pb.xnames,iv,"Y01502")
        ig = ig_["MND00802"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00802"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01502",ix_)
        arrset(pb.xnames,iv,"Y01502")
        ig = ig_["CAP00702"]
        pbm.A[ig,iv] += Float64(.0058)
        iv,ix_,_ = s2mpj_ii("Y01503",ix_)
        arrset(pb.xnames,iv,"Y01503")
        ig = ig_["CAP00703"]
        pbm.A[ig,iv] += Float64(.00607)
        ig = ig_["MXD00803"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01503",ix_)
        arrset(pb.xnames,iv,"Y01503")
        ig = ig_["INV00503"]
        pbm.A[ig,iv] += Float64(-.73994)
        ig = ig_["MND00804"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01503",ix_)
        arrset(pb.xnames,iv,"Y01503")
        ig = ig_["INV00502"]
        pbm.A[ig,iv] += Float64(-.46246)
        ig = ig_["MXD00804"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01503",ix_)
        arrset(pb.xnames,iv,"Y01503")
        ig = ig_["CAP06203"]
        pbm.A[ig,iv] += Float64(.00449)
        ig = ig_["CAP06202"]
        pbm.A[ig,iv] += Float64(.00077)
        iv,ix_,_ = s2mpj_ii("Y01503",ix_)
        arrset(pb.xnames,iv,"Y01503")
        ig = ig_["MXD00806"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00806"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01503",ix_)
        arrset(pb.xnames,iv,"Y01503")
        ig = ig_["CAP00702"]
        pbm.A[ig,iv] += Float64(.00321)
        ig = ig_["MND00803"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01503",ix_)
        arrset(pb.xnames,iv,"Y01503")
        ig = ig_["CAP06103"]
        pbm.A[ig,iv] += Float64(.00067)
        ig = ig_["MND00805"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01503",ix_)
        arrset(pb.xnames,iv,"Y01503")
        ig = ig_["MXD00805"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP00502"]
        pbm.A[ig,iv] += Float64(.0038)
        iv,ix_,_ = s2mpj_ii("Y01503",ix_)
        arrset(pb.xnames,iv,"Y01503")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-5.18)
        ig = ig_["CAP00503"]
        pbm.A[ig,iv] += Float64(.00856)
        iv,ix_,_ = s2mpj_ii("Y01504",ix_)
        arrset(pb.xnames,iv,"Y01504")
        ig = ig_["CAP00704"]
        pbm.A[ig,iv] += Float64(.00607)
        ig = ig_["CAP00504"]
        pbm.A[ig,iv] += Float64(.00856)
        iv,ix_,_ = s2mpj_ii("Y01504",ix_)
        arrset(pb.xnames,iv,"Y01504")
        ig = ig_["CAP00703"]
        pbm.A[ig,iv] += Float64(.00321)
        ig = ig_["CAP06203"]
        pbm.A[ig,iv] += Float64(.00077)
        iv,ix_,_ = s2mpj_ii("Y01504",ix_)
        arrset(pb.xnames,iv,"Y01504")
        ig = ig_["CAP06104"]
        pbm.A[ig,iv] += Float64(.00067)
        ig = ig_["INV00503"]
        pbm.A[ig,iv] += Float64(-.46246)
        iv,ix_,_ = s2mpj_ii("Y01504",ix_)
        arrset(pb.xnames,iv,"Y01504")
        ig = ig_["MND00804"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00805"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01504",ix_)
        arrset(pb.xnames,iv,"Y01504")
        ig = ig_["CAP06204"]
        pbm.A[ig,iv] += Float64(.00449)
        ig = ig_["CAP00503"]
        pbm.A[ig,iv] += Float64(.0038)
        iv,ix_,_ = s2mpj_ii("Y01504",ix_)
        arrset(pb.xnames,iv,"Y01504")
        ig = ig_["MXD00805"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00804"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01504",ix_)
        arrset(pb.xnames,iv,"Y01504")
        ig = ig_["CAP06103"]
        pbm.A[ig,iv] += Float64(.00006)
        ig = ig_["INV00504"]
        pbm.A[ig,iv] += Float64(-.73994)
        iv,ix_,_ = s2mpj_ii("Y01504",ix_)
        arrset(pb.xnames,iv,"Y01504")
        ig = ig_["MND00806"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00806"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01504",ix_)
        arrset(pb.xnames,iv,"Y01504")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-4.92)
        iv,ix_,_ = s2mpj_ii("Y01505",ix_)
        arrset(pb.xnames,iv,"Y01505")
        ig = ig_["MND00805"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP06104"]
        pbm.A[ig,iv] += Float64(.00005)
        iv,ix_,_ = s2mpj_ii("Y01505",ix_)
        arrset(pb.xnames,iv,"Y01505")
        ig = ig_["CAP06105"]
        pbm.A[ig,iv] += Float64(.00067)
        ig = ig_["MXD00805"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01505",ix_)
        arrset(pb.xnames,iv,"Y01505")
        ig = ig_["CAP06204"]
        pbm.A[ig,iv] += Float64(.00074)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-4.67)
        iv,ix_,_ = s2mpj_ii("Y01505",ix_)
        arrset(pb.xnames,iv,"Y01505")
        ig = ig_["CAP00704"]
        pbm.A[ig,iv] += Float64(.0031)
        ig = ig_["CAP06205"]
        pbm.A[ig,iv] += Float64(.00452)
        iv,ix_,_ = s2mpj_ii("Y01505",ix_)
        arrset(pb.xnames,iv,"Y01505")
        ig = ig_["CAP00504"]
        pbm.A[ig,iv] += Float64(.00366)
        ig = ig_["CAP00705"]
        pbm.A[ig,iv] += Float64(.00619)
        iv,ix_,_ = s2mpj_ii("Y01505",ix_)
        arrset(pb.xnames,iv,"Y01505")
        ig = ig_["INV00504"]
        pbm.A[ig,iv] += Float64(-.44534)
        ig = ig_["CAP00505"]
        pbm.A[ig,iv] += Float64(.0087)
        iv,ix_,_ = s2mpj_ii("Y01505",ix_)
        arrset(pb.xnames,iv,"Y01505")
        ig = ig_["MXD00806"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00806"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01505",ix_)
        arrset(pb.xnames,iv,"Y01505")
        ig = ig_["INV00505"]
        pbm.A[ig,iv] += Float64(-.75707)
        iv,ix_,_ = s2mpj_ii("Y01506",ix_)
        arrset(pb.xnames,iv,"Y01506")
        ig = ig_["CAP06105"]
        pbm.A[ig,iv] += Float64(.00006)
        ig = ig_["INV00506"]
        pbm.A[ig,iv] += Float64(-1.20241)
        iv,ix_,_ = s2mpj_ii("Y01506",ix_)
        arrset(pb.xnames,iv,"Y01506")
        ig = ig_["INV00505"]
        pbm.A[ig,iv] += Float64(-.48096)
        ig = ig_["MXD00806"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01506",ix_)
        arrset(pb.xnames,iv,"Y01506")
        ig = ig_["CAP00706"]
        pbm.A[ig,iv] += Float64(.00929)
        ig = ig_["CAP00705"]
        pbm.A[ig,iv] += Float64(.00334)
        iv,ix_,_ = s2mpj_ii("Y01506",ix_)
        arrset(pb.xnames,iv,"Y01506")
        ig = ig_["CAP06205"]
        pbm.A[ig,iv] += Float64(.0008)
        ig = ig_["CAP00505"]
        pbm.A[ig,iv] += Float64(.00396)
        iv,ix_,_ = s2mpj_ii("Y01506",ix_)
        arrset(pb.xnames,iv,"Y01506")
        ig = ig_["CAP00506"]
        pbm.A[ig,iv] += Float64(.01237)
        ig = ig_["CAP06106"]
        pbm.A[ig,iv] += Float64(.00072)
        iv,ix_,_ = s2mpj_ii("Y01506",ix_)
        arrset(pb.xnames,iv,"Y01506")
        ig = ig_["MND00806"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP06206"]
        pbm.A[ig,iv] += Float64(.00525)
        iv,ix_,_ = s2mpj_ii("Y01506",ix_)
        arrset(pb.xnames,iv,"Y01506")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-4.44)
        iv,ix_,_ = s2mpj_ii("Y01602",ix_)
        arrset(pb.xnames,iv,"Y01602")
        ig = ig_["MND00904"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00906"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01602",ix_)
        arrset(pb.xnames,iv,"Y01602")
        ig = ig_["MXD00903"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00904"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01602",ix_)
        arrset(pb.xnames,iv,"Y01602")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-3.99)
        ig = ig_["MND00903"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01602",ix_)
        arrset(pb.xnames,iv,"Y01602")
        ig = ig_["CAP00802"]
        pbm.A[ig,iv] += Float64(.00063)
        ig = ig_["MXD00906"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01602",ix_)
        arrset(pb.xnames,iv,"Y01602")
        ig = ig_["INV00401"]
        pbm.A[ig,iv] += Float64(-.28869)
        ig = ig_["CAP06302"]
        pbm.A[ig,iv] += Float64(.0004)
        iv,ix_,_ = s2mpj_ii("Y01602",ix_)
        arrset(pb.xnames,iv,"Y01602")
        ig = ig_["INV00402"]
        pbm.A[ig,iv] += Float64(-.86606)
        ig = ig_["MND00902"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01602",ix_)
        arrset(pb.xnames,iv,"Y01602")
        ig = ig_["CAP06301"]
        pbm.A[ig,iv] += Float64(.00002)
        ig = ig_["MXD00905"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01602",ix_)
        arrset(pb.xnames,iv,"Y01602")
        ig = ig_["MXD00902"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP00101"]
        pbm.A[ig,iv] += Float64(.00131)
        iv,ix_,_ = s2mpj_ii("Y01602",ix_)
        arrset(pb.xnames,iv,"Y01602")
        ig = ig_["MND00905"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP00102"]
        pbm.A[ig,iv] += Float64(.00654)
        iv,ix_,_ = s2mpj_ii("Y01602",ix_)
        arrset(pb.xnames,iv,"Y01602")
        ig = ig_["CAP00801"]
        pbm.A[ig,iv] += Float64(.00017)
        iv,ix_,_ = s2mpj_ii("Y01603",ix_)
        arrset(pb.xnames,iv,"Y01603")
        ig = ig_["CAP00102"]
        pbm.A[ig,iv] += Float64(.00121)
        ig = ig_["MXD00903"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01603",ix_)
        arrset(pb.xnames,iv,"Y01603")
        ig = ig_["MXD00904"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00904"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01603",ix_)
        arrset(pb.xnames,iv,"Y01603")
        ig = ig_["MND00903"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP00803"]
        pbm.A[ig,iv] += Float64(.00065)
        iv,ix_,_ = s2mpj_ii("Y01603",ix_)
        arrset(pb.xnames,iv,"Y01603")
        ig = ig_["INV00402"]
        pbm.A[ig,iv] += Float64(-.26648)
        ig = ig_["MND00906"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01603",ix_)
        arrset(pb.xnames,iv,"Y01603")
        ig = ig_["MXD00906"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP00802"]
        pbm.A[ig,iv] += Float64(.00015)
        iv,ix_,_ = s2mpj_ii("Y01603",ix_)
        arrset(pb.xnames,iv,"Y01603")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-3.79)
        ig = ig_["INV00403"]
        pbm.A[ig,iv] += Float64(-.88827)
        iv,ix_,_ = s2mpj_ii("Y01603",ix_)
        arrset(pb.xnames,iv,"Y01603")
        ig = ig_["MND00905"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP00103"]
        pbm.A[ig,iv] += Float64(.00664)
        iv,ix_,_ = s2mpj_ii("Y01603",ix_)
        arrset(pb.xnames,iv,"Y01603")
        ig = ig_["MXD00905"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP06302"]
        pbm.A[ig,iv] += Float64(.00002)
        iv,ix_,_ = s2mpj_ii("Y01603",ix_)
        arrset(pb.xnames,iv,"Y01603")
        ig = ig_["CAP06303"]
        pbm.A[ig,iv] += Float64(.0004)
        iv,ix_,_ = s2mpj_ii("Y01604",ix_)
        arrset(pb.xnames,iv,"Y01604")
        ig = ig_["CAP00103"]
        pbm.A[ig,iv] += Float64(.00121)
        ig = ig_["CAP00803"]
        pbm.A[ig,iv] += Float64(.00015)
        iv,ix_,_ = s2mpj_ii("Y01604",ix_)
        arrset(pb.xnames,iv,"Y01604")
        ig = ig_["CAP00104"]
        pbm.A[ig,iv] += Float64(.00664)
        ig = ig_["MND00906"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01604",ix_)
        arrset(pb.xnames,iv,"Y01604")
        ig = ig_["MXD00906"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP06303"]
        pbm.A[ig,iv] += Float64(.00002)
        iv,ix_,_ = s2mpj_ii("Y01604",ix_)
        arrset(pb.xnames,iv,"Y01604")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-3.59)
        ig = ig_["CAP00804"]
        pbm.A[ig,iv] += Float64(.00065)
        iv,ix_,_ = s2mpj_ii("Y01604",ix_)
        arrset(pb.xnames,iv,"Y01604")
        ig = ig_["MXD00904"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["INV00403"]
        pbm.A[ig,iv] += Float64(-.26648)
        iv,ix_,_ = s2mpj_ii("Y01604",ix_)
        arrset(pb.xnames,iv,"Y01604")
        ig = ig_["MND00905"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP06304"]
        pbm.A[ig,iv] += Float64(.0004)
        iv,ix_,_ = s2mpj_ii("Y01604",ix_)
        arrset(pb.xnames,iv,"Y01604")
        ig = ig_["MXD00905"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00904"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01604",ix_)
        arrset(pb.xnames,iv,"Y01604")
        ig = ig_["INV00404"]
        pbm.A[ig,iv] += Float64(-.88827)
        iv,ix_,_ = s2mpj_ii("Y01605",ix_)
        arrset(pb.xnames,iv,"Y01605")
        ig = ig_["MND00906"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-3.41)
        iv,ix_,_ = s2mpj_ii("Y01605",ix_)
        arrset(pb.xnames,iv,"Y01605")
        ig = ig_["MXD00906"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["INV00404"]
        pbm.A[ig,iv] += Float64(-.25661)
        iv,ix_,_ = s2mpj_ii("Y01605",ix_)
        arrset(pb.xnames,iv,"Y01605")
        ig = ig_["CAP00104"]
        pbm.A[ig,iv] += Float64(.00116)
        ig = ig_["MND00905"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01605",ix_)
        arrset(pb.xnames,iv,"Y01605")
        ig = ig_["CAP00804"]
        pbm.A[ig,iv] += Float64(.00015)
        ig = ig_["INV00405"]
        pbm.A[ig,iv] += Float64(-.89814)
        iv,ix_,_ = s2mpj_ii("Y01605",ix_)
        arrset(pb.xnames,iv,"Y01605")
        ig = ig_["MXD00905"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP06304"]
        pbm.A[ig,iv] += Float64(.00002)
        iv,ix_,_ = s2mpj_ii("Y01605",ix_)
        arrset(pb.xnames,iv,"Y01605")
        ig = ig_["CAP00105"]
        pbm.A[ig,iv] += Float64(.00668)
        ig = ig_["CAP06305"]
        pbm.A[ig,iv] += Float64(.0004)
        iv,ix_,_ = s2mpj_ii("Y01605",ix_)
        arrset(pb.xnames,iv,"Y01605")
        ig = ig_["CAP00805"]
        pbm.A[ig,iv] += Float64(.00065)
        iv,ix_,_ = s2mpj_ii("Y01606",ix_)
        arrset(pb.xnames,iv,"Y01606")
        ig = ig_["INV00405"]
        pbm.A[ig,iv] += Float64(-.27714)
        ig = ig_["CAP00106"]
        pbm.A[ig,iv] += Float64(.00784)
        iv,ix_,_ = s2mpj_ii("Y01606",ix_)
        arrset(pb.xnames,iv,"Y01606")
        ig = ig_["CAP06306"]
        pbm.A[ig,iv] += Float64(.00042)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-3.24)
        iv,ix_,_ = s2mpj_ii("Y01606",ix_)
        arrset(pb.xnames,iv,"Y01606")
        ig = ig_["CAP00805"]
        pbm.A[ig,iv] += Float64(.00016)
        ig = ig_["MND00906"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01606",ix_)
        arrset(pb.xnames,iv,"Y01606")
        ig = ig_["CAP00105"]
        pbm.A[ig,iv] += Float64(.00125)
        ig = ig_["CAP06305"]
        pbm.A[ig,iv] += Float64(.00002)
        iv,ix_,_ = s2mpj_ii("Y01606",ix_)
        arrset(pb.xnames,iv,"Y01606")
        ig = ig_["CAP00806"]
        pbm.A[ig,iv] += Float64(.0008)
        ig = ig_["INV00406"]
        pbm.A[ig,iv] += Float64(-1.15475)
        iv,ix_,_ = s2mpj_ii("Y01606",ix_)
        arrset(pb.xnames,iv,"Y01606")
        ig = ig_["MXD00906"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01702",ix_)
        arrset(pb.xnames,iv,"Y01702")
        ig = ig_["CAP00902"]
        pbm.A[ig,iv] += Float64(.00654)
        ig = ig_["CAP01301"]
        pbm.A[ig,iv] += Float64(.00017)
        iv,ix_,_ = s2mpj_ii("Y01702",ix_)
        arrset(pb.xnames,iv,"Y01702")
        ig = ig_["CAP06301"]
        pbm.A[ig,iv] += Float64(.00002)
        ig = ig_["MND00902"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01702",ix_)
        arrset(pb.xnames,iv,"Y01702")
        ig = ig_["INV00401"]
        pbm.A[ig,iv] += Float64(-.28869)
        ig = ig_["MXD00902"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01702",ix_)
        arrset(pb.xnames,iv,"Y01702")
        ig = ig_["CAP01302"]
        pbm.A[ig,iv] += Float64(.00063)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-3.97)
        iv,ix_,_ = s2mpj_ii("Y01702",ix_)
        arrset(pb.xnames,iv,"Y01702")
        ig = ig_["MND00903"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00903"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01702",ix_)
        arrset(pb.xnames,iv,"Y01702")
        ig = ig_["MND00905"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00905"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01702",ix_)
        arrset(pb.xnames,iv,"Y01702")
        ig = ig_["CAP06302"]
        pbm.A[ig,iv] += Float64(.0004)
        ig = ig_["INV00402"]
        pbm.A[ig,iv] += Float64(-.86606)
        iv,ix_,_ = s2mpj_ii("Y01702",ix_)
        arrset(pb.xnames,iv,"Y01702")
        ig = ig_["MND00906"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00906"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01702",ix_)
        arrset(pb.xnames,iv,"Y01702")
        ig = ig_["CAP00901"]
        pbm.A[ig,iv] += Float64(.00131)
        ig = ig_["MND00904"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01702",ix_)
        arrset(pb.xnames,iv,"Y01702")
        ig = ig_["MXD00904"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01703",ix_)
        arrset(pb.xnames,iv,"Y01703")
        ig = ig_["INV00402"]
        pbm.A[ig,iv] += Float64(-.26648)
        ig = ig_["INV00403"]
        pbm.A[ig,iv] += Float64(-.88827)
        iv,ix_,_ = s2mpj_ii("Y01703",ix_)
        arrset(pb.xnames,iv,"Y01703")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-3.77)
        ig = ig_["CAP01302"]
        pbm.A[ig,iv] += Float64(.00015)
        iv,ix_,_ = s2mpj_ii("Y01703",ix_)
        arrset(pb.xnames,iv,"Y01703")
        ig = ig_["CAP06303"]
        pbm.A[ig,iv] += Float64(.0004)
        ig = ig_["MXD00906"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01703",ix_)
        arrset(pb.xnames,iv,"Y01703")
        ig = ig_["MND00906"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00904"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01703",ix_)
        arrset(pb.xnames,iv,"Y01703")
        ig = ig_["MND00903"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP06302"]
        pbm.A[ig,iv] += Float64(.00002)
        iv,ix_,_ = s2mpj_ii("Y01703",ix_)
        arrset(pb.xnames,iv,"Y01703")
        ig = ig_["MXD00903"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00904"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01703",ix_)
        arrset(pb.xnames,iv,"Y01703")
        ig = ig_["MXD00905"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP00902"]
        pbm.A[ig,iv] += Float64(.00121)
        iv,ix_,_ = s2mpj_ii("Y01703",ix_)
        arrset(pb.xnames,iv,"Y01703")
        ig = ig_["MND00905"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP00903"]
        pbm.A[ig,iv] += Float64(.00664)
        iv,ix_,_ = s2mpj_ii("Y01703",ix_)
        arrset(pb.xnames,iv,"Y01703")
        ig = ig_["CAP01303"]
        pbm.A[ig,iv] += Float64(.00065)
        iv,ix_,_ = s2mpj_ii("Y01704",ix_)
        arrset(pb.xnames,iv,"Y01704")
        ig = ig_["MND00906"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00906"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01704",ix_)
        arrset(pb.xnames,iv,"Y01704")
        ig = ig_["INV00404"]
        pbm.A[ig,iv] += Float64(-.88827)
        ig = ig_["MND00904"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01704",ix_)
        arrset(pb.xnames,iv,"Y01704")
        ig = ig_["CAP01303"]
        pbm.A[ig,iv] += Float64(.00015)
        ig = ig_["INV00403"]
        pbm.A[ig,iv] += Float64(-.26648)
        iv,ix_,_ = s2mpj_ii("Y01704",ix_)
        arrset(pb.xnames,iv,"Y01704")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-3.58)
        ig = ig_["MXD00905"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01704",ix_)
        arrset(pb.xnames,iv,"Y01704")
        ig = ig_["CAP06304"]
        pbm.A[ig,iv] += Float64(.0004)
        ig = ig_["MND00905"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01704",ix_)
        arrset(pb.xnames,iv,"Y01704")
        ig = ig_["CAP00903"]
        pbm.A[ig,iv] += Float64(.00121)
        ig = ig_["CAP01304"]
        pbm.A[ig,iv] += Float64(.00065)
        iv,ix_,_ = s2mpj_ii("Y01704",ix_)
        arrset(pb.xnames,iv,"Y01704")
        ig = ig_["CAP00904"]
        pbm.A[ig,iv] += Float64(.00664)
        ig = ig_["MXD00904"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01704",ix_)
        arrset(pb.xnames,iv,"Y01704")
        ig = ig_["CAP06303"]
        pbm.A[ig,iv] += Float64(.00002)
        iv,ix_,_ = s2mpj_ii("Y01705",ix_)
        arrset(pb.xnames,iv,"Y01705")
        ig = ig_["CAP01305"]
        pbm.A[ig,iv] += Float64(.00065)
        ig = ig_["CAP06305"]
        pbm.A[ig,iv] += Float64(.0004)
        iv,ix_,_ = s2mpj_ii("Y01705",ix_)
        arrset(pb.xnames,iv,"Y01705")
        ig = ig_["MXD00905"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["INV00404"]
        pbm.A[ig,iv] += Float64(-.25661)
        iv,ix_,_ = s2mpj_ii("Y01705",ix_)
        arrset(pb.xnames,iv,"Y01705")
        ig = ig_["INV00405"]
        pbm.A[ig,iv] += Float64(-.89814)
        ig = ig_["MND00905"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01705",ix_)
        arrset(pb.xnames,iv,"Y01705")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-3.4)
        ig = ig_["MND00906"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01705",ix_)
        arrset(pb.xnames,iv,"Y01705")
        ig = ig_["CAP00904"]
        pbm.A[ig,iv] += Float64(.00116)
        ig = ig_["CAP01304"]
        pbm.A[ig,iv] += Float64(.00015)
        iv,ix_,_ = s2mpj_ii("Y01705",ix_)
        arrset(pb.xnames,iv,"Y01705")
        ig = ig_["CAP06304"]
        pbm.A[ig,iv] += Float64(.00002)
        ig = ig_["CAP00905"]
        pbm.A[ig,iv] += Float64(.00668)
        iv,ix_,_ = s2mpj_ii("Y01705",ix_)
        arrset(pb.xnames,iv,"Y01705")
        ig = ig_["MXD00906"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01706",ix_)
        arrset(pb.xnames,iv,"Y01706")
        ig = ig_["INV00406"]
        pbm.A[ig,iv] += Float64(-1.15475)
        ig = ig_["MXD00906"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01706",ix_)
        arrset(pb.xnames,iv,"Y01706")
        ig = ig_["MND00906"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP06306"]
        pbm.A[ig,iv] += Float64(.00042)
        iv,ix_,_ = s2mpj_ii("Y01706",ix_)
        arrset(pb.xnames,iv,"Y01706")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-3.22)
        ig = ig_["CAP01305"]
        pbm.A[ig,iv] += Float64(.00016)
        iv,ix_,_ = s2mpj_ii("Y01706",ix_)
        arrset(pb.xnames,iv,"Y01706")
        ig = ig_["CAP01306"]
        pbm.A[ig,iv] += Float64(.0008)
        ig = ig_["INV00405"]
        pbm.A[ig,iv] += Float64(-.27714)
        iv,ix_,_ = s2mpj_ii("Y01706",ix_)
        arrset(pb.xnames,iv,"Y01706")
        ig = ig_["CAP00905"]
        pbm.A[ig,iv] += Float64(.00125)
        ig = ig_["CAP00906"]
        pbm.A[ig,iv] += Float64(.00784)
        iv,ix_,_ = s2mpj_ii("Y01706",ix_)
        arrset(pb.xnames,iv,"Y01706")
        ig = ig_["CAP06305"]
        pbm.A[ig,iv] += Float64(.00002)
        iv,ix_,_ = s2mpj_ii("Y01802",ix_)
        arrset(pb.xnames,iv,"Y01802")
        ig = ig_["MND00905"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP01401"]
        pbm.A[ig,iv] += Float64(.00077)
        iv,ix_,_ = s2mpj_ii("Y01802",ix_)
        arrset(pb.xnames,iv,"Y01802")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-3.93)
        ig = ig_["MXD00905"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01802",ix_)
        arrset(pb.xnames,iv,"Y01802")
        ig = ig_["INV00401"]
        pbm.A[ig,iv] += Float64(-.14434)
        ig = ig_["MND00906"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01802",ix_)
        arrset(pb.xnames,iv,"Y01802")
        ig = ig_["MND00904"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00906"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01802",ix_)
        arrset(pb.xnames,iv,"Y01802")
        ig = ig_["CAP06301"]
        pbm.A[ig,iv] += Float64(.00002)
        ig = ig_["CAP01802"]
        pbm.A[ig,iv] += Float64(.00071)
        iv,ix_,_ = s2mpj_ii("Y01802",ix_)
        arrset(pb.xnames,iv,"Y01802")
        ig = ig_["MXD00903"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00903"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01802",ix_)
        arrset(pb.xnames,iv,"Y01802")
        ig = ig_["CAP06302"]
        pbm.A[ig,iv] += Float64(.0004)
        ig = ig_["MXD00904"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01802",ix_)
        arrset(pb.xnames,iv,"Y01802")
        ig = ig_["CAP01402"]
        pbm.A[ig,iv] += Float64(.00662)
        ig = ig_["INV00402"]
        pbm.A[ig,iv] += Float64(-1.0104)
        iv,ix_,_ = s2mpj_ii("Y01802",ix_)
        arrset(pb.xnames,iv,"Y01802")
        ig = ig_["CAP01801"]
        pbm.A[ig,iv] += Float64(.00009)
        ig = ig_["MXD00902"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01802",ix_)
        arrset(pb.xnames,iv,"Y01802")
        ig = ig_["MND00902"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01803",ix_)
        arrset(pb.xnames,iv,"Y01803")
        ig = ig_["MND00904"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD00904"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01803",ix_)
        arrset(pb.xnames,iv,"Y01803")
        ig = ig_["INV00403"]
        pbm.A[ig,iv] += Float64(-1.02151)
        ig = ig_["MXD00903"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01803",ix_)
        arrset(pb.xnames,iv,"Y01803")
        ig = ig_["MND00903"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP01402"]
        pbm.A[ig,iv] += Float64(.00071)
        iv,ix_,_ = s2mpj_ii("Y01803",ix_)
        arrset(pb.xnames,iv,"Y01803")
        ig = ig_["CAP01802"]
        pbm.A[ig,iv] += Float64(.00008)
        ig = ig_["CAP06303"]
        pbm.A[ig,iv] += Float64(.0004)
        iv,ix_,_ = s2mpj_ii("Y01803",ix_)
        arrset(pb.xnames,iv,"Y01803")
        ig = ig_["CAP06302"]
        pbm.A[ig,iv] += Float64(.00002)
        ig = ig_["CAP01803"]
        pbm.A[ig,iv] += Float64(.00072)
        iv,ix_,_ = s2mpj_ii("Y01803",ix_)
        arrset(pb.xnames,iv,"Y01803")
        ig = ig_["CAP01403"]
        pbm.A[ig,iv] += Float64(.00668)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-3.73)
        iv,ix_,_ = s2mpj_ii("Y01803",ix_)
        arrset(pb.xnames,iv,"Y01803")
        ig = ig_["MXD00906"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00906"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01803",ix_)
        arrset(pb.xnames,iv,"Y01803")
        ig = ig_["INV00402"]
        pbm.A[ig,iv] += Float64(-.13324)
        ig = ig_["MND00905"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01803",ix_)
        arrset(pb.xnames,iv,"Y01803")
        ig = ig_["MXD00905"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01804",ix_)
        arrset(pb.xnames,iv,"Y01804")
        ig = ig_["CAP01404"]
        pbm.A[ig,iv] += Float64(.00668)
        ig = ig_["MND00906"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01804",ix_)
        arrset(pb.xnames,iv,"Y01804")
        ig = ig_["MXD00906"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP01804"]
        pbm.A[ig,iv] += Float64(.00072)
        iv,ix_,_ = s2mpj_ii("Y01804",ix_)
        arrset(pb.xnames,iv,"Y01804")
        ig = ig_["MND00905"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP06303"]
        pbm.A[ig,iv] += Float64(.00002)
        iv,ix_,_ = s2mpj_ii("Y01804",ix_)
        arrset(pb.xnames,iv,"Y01804")
        ig = ig_["CAP06304"]
        pbm.A[ig,iv] += Float64(.0004)
        ig = ig_["MXD00905"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01804",ix_)
        arrset(pb.xnames,iv,"Y01804")
        ig = ig_["INV00403"]
        pbm.A[ig,iv] += Float64(-.13324)
        ig = ig_["INV00404"]
        pbm.A[ig,iv] += Float64(-1.02151)
        iv,ix_,_ = s2mpj_ii("Y01804",ix_)
        arrset(pb.xnames,iv,"Y01804")
        ig = ig_["MND00904"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-3.54)
        iv,ix_,_ = s2mpj_ii("Y01804",ix_)
        arrset(pb.xnames,iv,"Y01804")
        ig = ig_["MXD00904"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP01403"]
        pbm.A[ig,iv] += Float64(.00071)
        iv,ix_,_ = s2mpj_ii("Y01805",ix_)
        arrset(pb.xnames,iv,"Y01805")
        ig = ig_["MXD00905"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND00905"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01805",ix_)
        arrset(pb.xnames,iv,"Y01805")
        ig = ig_["CAP01805"]
        pbm.A[ig,iv] += Float64(.00072)
        ig = ig_["CAP06305"]
        pbm.A[ig,iv] += Float64(.0004)
        iv,ix_,_ = s2mpj_ii("Y01805",ix_)
        arrset(pb.xnames,iv,"Y01805")
        ig = ig_["CAP01404"]
        pbm.A[ig,iv] += Float64(.00068)
        ig = ig_["CAP06304"]
        pbm.A[ig,iv] += Float64(.00002)
        iv,ix_,_ = s2mpj_ii("Y01805",ix_)
        arrset(pb.xnames,iv,"Y01805")
        ig = ig_["CAP01405"]
        pbm.A[ig,iv] += Float64(.00671)
        ig = ig_["INV00405"]
        pbm.A[ig,iv] += Float64(-1.02644)
        iv,ix_,_ = s2mpj_ii("Y01805",ix_)
        arrset(pb.xnames,iv,"Y01805")
        ig = ig_["CAP01804"]
        pbm.A[ig,iv] += Float64(.00008)
        ig = ig_["INV00404"]
        pbm.A[ig,iv] += Float64(-.12831)
        iv,ix_,_ = s2mpj_ii("Y01805",ix_)
        arrset(pb.xnames,iv,"Y01805")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-3.36)
        ig = ig_["MND00906"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01805",ix_)
        arrset(pb.xnames,iv,"Y01805")
        ig = ig_["MXD00906"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01806",ix_)
        arrset(pb.xnames,iv,"Y01806")
        ig = ig_["INV00405"]
        pbm.A[ig,iv] += Float64(-.13857)
        ig = ig_["MXD00906"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01806",ix_)
        arrset(pb.xnames,iv,"Y01806")
        ig = ig_["CAP01406"]
        pbm.A[ig,iv] += Float64(.00739)
        ig = ig_["CAP01806"]
        pbm.A[ig,iv] += Float64(.0008)
        iv,ix_,_ = s2mpj_ii("Y01806",ix_)
        arrset(pb.xnames,iv,"Y01806")
        ig = ig_["MND00906"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP06306"]
        pbm.A[ig,iv] += Float64(.00042)
        iv,ix_,_ = s2mpj_ii("Y01806",ix_)
        arrset(pb.xnames,iv,"Y01806")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-3.19)
        ig = ig_["CAP06305"]
        pbm.A[ig,iv] += Float64(.00002)
        iv,ix_,_ = s2mpj_ii("Y01806",ix_)
        arrset(pb.xnames,iv,"Y01806")
        ig = ig_["CAP01405"]
        pbm.A[ig,iv] += Float64(.00074)
        ig = ig_["INV00406"]
        pbm.A[ig,iv] += Float64(-1.15475)
        iv,ix_,_ = s2mpj_ii("Y01806",ix_)
        arrset(pb.xnames,iv,"Y01806")
        ig = ig_["CAP01805"]
        pbm.A[ig,iv] += Float64(.00009)
        iv,ix_,_ = s2mpj_ii("Y01902",ix_)
        arrset(pb.xnames,iv,"Y01902")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-3.67)
        ig = ig_["MXD01004"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01902",ix_)
        arrset(pb.xnames,iv,"Y01902")
        ig = ig_["MXD01006"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP06201"]
        pbm.A[ig,iv] += Float64(.00083)
        iv,ix_,_ = s2mpj_ii("Y01902",ix_)
        arrset(pb.xnames,iv,"Y01902")
        ig = ig_["CAP00901"]
        pbm.A[ig,iv] += Float64(.00272)
        ig = ig_["CAP06102"]
        pbm.A[ig,iv] += Float64(.00066)
        iv,ix_,_ = s2mpj_ii("Y01902",ix_)
        arrset(pb.xnames,iv,"Y01902")
        ig = ig_["CAP06101"]
        pbm.A[ig,iv] += Float64(.00006)
        ig = ig_["MND01006"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01902",ix_)
        arrset(pb.xnames,iv,"Y01902")
        ig = ig_["CAP01302"]
        pbm.A[ig,iv] += Float64(.00052)
        ig = ig_["MND01004"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01902",ix_)
        arrset(pb.xnames,iv,"Y01902")
        ig = ig_["MND01005"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP01301"]
        pbm.A[ig,iv] += Float64(.00031)
        iv,ix_,_ = s2mpj_ii("Y01902",ix_)
        arrset(pb.xnames,iv,"Y01902")
        ig = ig_["MND01003"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD01003"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01902",ix_)
        arrset(pb.xnames,iv,"Y01902")
        ig = ig_["MXD01002"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND01002"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01902",ix_)
        arrset(pb.xnames,iv,"Y01902")
        ig = ig_["INV00502"]
        pbm.A[ig,iv] += Float64(-.7014)
        ig = ig_["MXD01005"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01902",ix_)
        arrset(pb.xnames,iv,"Y01902")
        ig = ig_["INV00501"]
        pbm.A[ig,iv] += Float64(-.501)
        ig = ig_["CAP00902"]
        pbm.A[ig,iv] += Float64(.00544)
        iv,ix_,_ = s2mpj_ii("Y01902",ix_)
        arrset(pb.xnames,iv,"Y01902")
        ig = ig_["CAP06202"]
        pbm.A[ig,iv] += Float64(.00442)
        iv,ix_,_ = s2mpj_ii("Y01903",ix_)
        arrset(pb.xnames,iv,"Y01903")
        ig = ig_["MND01006"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP06203"]
        pbm.A[ig,iv] += Float64(.00449)
        iv,ix_,_ = s2mpj_ii("Y01903",ix_)
        arrset(pb.xnames,iv,"Y01903")
        ig = ig_["CAP06103"]
        pbm.A[ig,iv] += Float64(.00067)
        ig = ig_["MXD01006"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01903",ix_)
        arrset(pb.xnames,iv,"Y01903")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-3.49)
        ig = ig_["CAP06202"]
        pbm.A[ig,iv] += Float64(.00077)
        iv,ix_,_ = s2mpj_ii("Y01903",ix_)
        arrset(pb.xnames,iv,"Y01903")
        ig = ig_["CAP06102"]
        pbm.A[ig,iv] += Float64(.00006)
        ig = ig_["CAP00903"]
        pbm.A[ig,iv] += Float64(.00565)
        iv,ix_,_ = s2mpj_ii("Y01903",ix_)
        arrset(pb.xnames,iv,"Y01903")
        ig = ig_["MND01005"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD01005"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01903",ix_)
        arrset(pb.xnames,iv,"Y01903")
        ig = ig_["CAP01302"]
        pbm.A[ig,iv] += Float64(.00029)
        ig = ig_["INV00502"]
        pbm.A[ig,iv] += Float64(-.46246)
        iv,ix_,_ = s2mpj_ii("Y01903",ix_)
        arrset(pb.xnames,iv,"Y01903")
        ig = ig_["MND01003"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD01003"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01903",ix_)
        arrset(pb.xnames,iv,"Y01903")
        ig = ig_["INV00503"]
        pbm.A[ig,iv] += Float64(-.73994)
        ig = ig_["CAP01303"]
        pbm.A[ig,iv] += Float64(.00054)
        iv,ix_,_ = s2mpj_ii("Y01903",ix_)
        arrset(pb.xnames,iv,"Y01903")
        ig = ig_["CAP00902"]
        pbm.A[ig,iv] += Float64(.00251)
        ig = ig_["MND01004"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01903",ix_)
        arrset(pb.xnames,iv,"Y01903")
        ig = ig_["MXD01004"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01904",ix_)
        arrset(pb.xnames,iv,"Y01904")
        ig = ig_["CAP06104"]
        pbm.A[ig,iv] += Float64(.00067)
        ig = ig_["CAP06204"]
        pbm.A[ig,iv] += Float64(.00449)
        iv,ix_,_ = s2mpj_ii("Y01904",ix_)
        arrset(pb.xnames,iv,"Y01904")
        ig = ig_["CAP06103"]
        pbm.A[ig,iv] += Float64(.00006)
        ig = ig_["CAP01303"]
        pbm.A[ig,iv] += Float64(.00029)
        iv,ix_,_ = s2mpj_ii("Y01904",ix_)
        arrset(pb.xnames,iv,"Y01904")
        ig = ig_["MND01005"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP06203"]
        pbm.A[ig,iv] += Float64(.00077)
        iv,ix_,_ = s2mpj_ii("Y01904",ix_)
        arrset(pb.xnames,iv,"Y01904")
        ig = ig_["CAP01304"]
        pbm.A[ig,iv] += Float64(.00054)
        ig = ig_["MXD01005"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01904",ix_)
        arrset(pb.xnames,iv,"Y01904")
        ig = ig_["MND01004"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["INV00504"]
        pbm.A[ig,iv] += Float64(-.73994)
        iv,ix_,_ = s2mpj_ii("Y01904",ix_)
        arrset(pb.xnames,iv,"Y01904")
        ig = ig_["CAP00904"]
        pbm.A[ig,iv] += Float64(.00565)
        ig = ig_["INV00503"]
        pbm.A[ig,iv] += Float64(-.46246)
        iv,ix_,_ = s2mpj_ii("Y01904",ix_)
        arrset(pb.xnames,iv,"Y01904")
        ig = ig_["MXD01006"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD01004"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01904",ix_)
        arrset(pb.xnames,iv,"Y01904")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-3.31)
        ig = ig_["CAP00903"]
        pbm.A[ig,iv] += Float64(.00251)
        iv,ix_,_ = s2mpj_ii("Y01904",ix_)
        arrset(pb.xnames,iv,"Y01904")
        ig = ig_["MND01006"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01905",ix_)
        arrset(pb.xnames,iv,"Y01905")
        ig = ig_["CAP00904"]
        pbm.A[ig,iv] += Float64(.00242)
        ig = ig_["MXD01006"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01905",ix_)
        arrset(pb.xnames,iv,"Y01905")
        ig = ig_["CAP01304"]
        pbm.A[ig,iv] += Float64(.00028)
        ig = ig_["CAP06105"]
        pbm.A[ig,iv] += Float64(.00067)
        iv,ix_,_ = s2mpj_ii("Y01905",ix_)
        arrset(pb.xnames,iv,"Y01905")
        ig = ig_["MXD01005"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP06205"]
        pbm.A[ig,iv] += Float64(.00452)
        iv,ix_,_ = s2mpj_ii("Y01905",ix_)
        arrset(pb.xnames,iv,"Y01905")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-3.14)
        ig = ig_["INV00504"]
        pbm.A[ig,iv] += Float64(-.44534)
        iv,ix_,_ = s2mpj_ii("Y01905",ix_)
        arrset(pb.xnames,iv,"Y01905")
        ig = ig_["CAP01305"]
        pbm.A[ig,iv] += Float64(.00056)
        ig = ig_["CAP06104"]
        pbm.A[ig,iv] += Float64(.00005)
        iv,ix_,_ = s2mpj_ii("Y01905",ix_)
        arrset(pb.xnames,iv,"Y01905")
        ig = ig_["MND01006"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP00905"]
        pbm.A[ig,iv] += Float64(.00575)
        iv,ix_,_ = s2mpj_ii("Y01905",ix_)
        arrset(pb.xnames,iv,"Y01905")
        ig = ig_["CAP06204"]
        pbm.A[ig,iv] += Float64(.00074)
        ig = ig_["MND01005"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01905",ix_)
        arrset(pb.xnames,iv,"Y01905")
        ig = ig_["INV00505"]
        pbm.A[ig,iv] += Float64(-.75707)
        iv,ix_,_ = s2mpj_ii("Y01906",ix_)
        arrset(pb.xnames,iv,"Y01906")
        ig = ig_["CAP01305"]
        pbm.A[ig,iv] += Float64(.0003)
        ig = ig_["CAP06105"]
        pbm.A[ig,iv] += Float64(.00006)
        iv,ix_,_ = s2mpj_ii("Y01906",ix_)
        arrset(pb.xnames,iv,"Y01906")
        ig = ig_["CAP06205"]
        pbm.A[ig,iv] += Float64(.0008)
        ig = ig_["MND01006"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01906",ix_)
        arrset(pb.xnames,iv,"Y01906")
        ig = ig_["INV00505"]
        pbm.A[ig,iv] += Float64(-.48096)
        ig = ig_["MXD01006"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y01906",ix_)
        arrset(pb.xnames,iv,"Y01906")
        ig = ig_["CAP01306"]
        pbm.A[ig,iv] += Float64(.00083)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-2.99)
        iv,ix_,_ = s2mpj_ii("Y01906",ix_)
        arrset(pb.xnames,iv,"Y01906")
        ig = ig_["CAP00905"]
        pbm.A[ig,iv] += Float64(.00261)
        ig = ig_["CAP06106"]
        pbm.A[ig,iv] += Float64(.00072)
        iv,ix_,_ = s2mpj_ii("Y01906",ix_)
        arrset(pb.xnames,iv,"Y01906")
        ig = ig_["INV00506"]
        pbm.A[ig,iv] += Float64(-1.20241)
        ig = ig_["CAP00906"]
        pbm.A[ig,iv] += Float64(.00817)
        iv,ix_,_ = s2mpj_ii("Y02002",ix_)
        arrset(pb.xnames,iv,"Y02002")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-3.65)
        ig = ig_["MND01006"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y02002",ix_)
        arrset(pb.xnames,iv,"Y02002")
        ig = ig_["MXD01006"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP06102"]
        pbm.A[ig,iv] += Float64(.00066)
        iv,ix_,_ = s2mpj_ii("Y02002",ix_)
        arrset(pb.xnames,iv,"Y02002")
        ig = ig_["CAP01802"]
        pbm.A[ig,iv] += Float64(.0006)
        ig = ig_["MND01004"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y02002",ix_)
        arrset(pb.xnames,iv,"Y02002")
        ig = ig_["MXD01004"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD01003"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y02002",ix_)
        arrset(pb.xnames,iv,"Y02002")
        ig = ig_["MND01003"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["INV00502"]
        pbm.A[ig,iv] += Float64(-.85171)
        iv,ix_,_ = s2mpj_ii("Y02002",ix_)
        arrset(pb.xnames,iv,"Y02002")
        ig = ig_["CAP01402"]
        pbm.A[ig,iv] += Float64(.00561)
        ig = ig_["MXD01002"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y02002",ix_)
        arrset(pb.xnames,iv,"Y02002")
        ig = ig_["MND01002"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND01005"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y02002",ix_)
        arrset(pb.xnames,iv,"Y02002")
        ig = ig_["CAP01801"]
        pbm.A[ig,iv] += Float64(.00023)
        ig = ig_["MXD01005"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y02002",ix_)
        arrset(pb.xnames,iv,"Y02002")
        ig = ig_["CAP01401"]
        pbm.A[ig,iv] += Float64(.00209)
        ig = ig_["CAP06201"]
        pbm.A[ig,iv] += Float64(.00083)
        iv,ix_,_ = s2mpj_ii("Y02002",ix_)
        arrset(pb.xnames,iv,"Y02002")
        ig = ig_["CAP06101"]
        pbm.A[ig,iv] += Float64(.00006)
        ig = ig_["INV00501"]
        pbm.A[ig,iv] += Float64(-.3507)
        iv,ix_,_ = s2mpj_ii("Y02002",ix_)
        arrset(pb.xnames,iv,"Y02002")
        ig = ig_["CAP06202"]
        pbm.A[ig,iv] += Float64(.00442)
        iv,ix_,_ = s2mpj_ii("Y02003",ix_)
        arrset(pb.xnames,iv,"Y02003")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-3.46)
        ig = ig_["MXD01004"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y02003",ix_)
        arrset(pb.xnames,iv,"Y02003")
        ig = ig_["MXD01005"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP01803"]
        pbm.A[ig,iv] += Float64(.00062)
        iv,ix_,_ = s2mpj_ii("Y02003",ix_)
        arrset(pb.xnames,iv,"Y02003")
        ig = ig_["MND01006"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD01006"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y02003",ix_)
        arrset(pb.xnames,iv,"Y02003")
        ig = ig_["CAP01402"]
        pbm.A[ig,iv] += Float64(.00192)
        ig = ig_["CAP06103"]
        pbm.A[ig,iv] += Float64(.00067)
        iv,ix_,_ = s2mpj_ii("Y02003",ix_)
        arrset(pb.xnames,iv,"Y02003")
        ig = ig_["MND01003"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP06203"]
        pbm.A[ig,iv] += Float64(.00449)
        iv,ix_,_ = s2mpj_ii("Y02003",ix_)
        arrset(pb.xnames,iv,"Y02003")
        ig = ig_["CAP01403"]
        pbm.A[ig,iv] += Float64(.00577)
        ig = ig_["CAP01802"]
        pbm.A[ig,iv] += Float64(.00022)
        iv,ix_,_ = s2mpj_ii("Y02003",ix_)
        arrset(pb.xnames,iv,"Y02003")
        ig = ig_["INV00502"]
        pbm.A[ig,iv] += Float64(-.32373)
        ig = ig_["INV00503"]
        pbm.A[ig,iv] += Float64(-.87868)
        iv,ix_,_ = s2mpj_ii("Y02003",ix_)
        arrset(pb.xnames,iv,"Y02003")
        ig = ig_["MXD01003"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP06102"]
        pbm.A[ig,iv] += Float64(.00006)
        iv,ix_,_ = s2mpj_ii("Y02003",ix_)
        arrset(pb.xnames,iv,"Y02003")
        ig = ig_["CAP06202"]
        pbm.A[ig,iv] += Float64(.00077)
        ig = ig_["MND01005"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y02003",ix_)
        arrset(pb.xnames,iv,"Y02003")
        ig = ig_["MND01004"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y02004",ix_)
        arrset(pb.xnames,iv,"Y02004")
        ig = ig_["CAP01804"]
        pbm.A[ig,iv] += Float64(.00062)
        ig = ig_["CAP01803"]
        pbm.A[ig,iv] += Float64(.00022)
        iv,ix_,_ = s2mpj_ii("Y02004",ix_)
        arrset(pb.xnames,iv,"Y02004")
        ig = ig_["CAP06103"]
        pbm.A[ig,iv] += Float64(.00006)
        ig = ig_["CAP06203"]
        pbm.A[ig,iv] += Float64(.00077)
        iv,ix_,_ = s2mpj_ii("Y02004",ix_)
        arrset(pb.xnames,iv,"Y02004")
        ig = ig_["CAP01404"]
        pbm.A[ig,iv] += Float64(.00577)
        ig = ig_["CAP06204"]
        pbm.A[ig,iv] += Float64(.00449)
        iv,ix_,_ = s2mpj_ii("Y02004",ix_)
        arrset(pb.xnames,iv,"Y02004")
        ig = ig_["CAP06104"]
        pbm.A[ig,iv] += Float64(.00067)
        ig = ig_["MND01004"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y02004",ix_)
        arrset(pb.xnames,iv,"Y02004")
        ig = ig_["MND01006"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["INV00503"]
        pbm.A[ig,iv] += Float64(-.32373)
        iv,ix_,_ = s2mpj_ii("Y02004",ix_)
        arrset(pb.xnames,iv,"Y02004")
        ig = ig_["MXD01004"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["CAP01403"]
        pbm.A[ig,iv] += Float64(.00192)
        iv,ix_,_ = s2mpj_ii("Y02004",ix_)
        arrset(pb.xnames,iv,"Y02004")
        ig = ig_["MXD01006"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["INV00504"]
        pbm.A[ig,iv] += Float64(-.87868)
        iv,ix_,_ = s2mpj_ii("Y02004",ix_)
        arrset(pb.xnames,iv,"Y02004")
        ig = ig_["MXD01005"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND01005"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y02004",ix_)
        arrset(pb.xnames,iv,"Y02004")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-3.29)
        iv,ix_,_ = s2mpj_ii("Y02005",ix_)
        arrset(pb.xnames,iv,"Y02005")
        ig = ig_["CAP06104"]
        pbm.A[ig,iv] += Float64(.00005)
        ig = ig_["CAP01804"]
        pbm.A[ig,iv] += Float64(.00021)
        iv,ix_,_ = s2mpj_ii("Y02005",ix_)
        arrset(pb.xnames,iv,"Y02005")
        ig = ig_["INV00505"]
        pbm.A[ig,iv] += Float64(-.89067)
        ig = ig_["CAP06204"]
        pbm.A[ig,iv] += Float64(.00074)
        iv,ix_,_ = s2mpj_ii("Y02005",ix_)
        arrset(pb.xnames,iv,"Y02005")
        ig = ig_["INV00504"]
        pbm.A[ig,iv] += Float64(-.31174)
        ig = ig_["CAP01405"]
        pbm.A[ig,iv] += Float64(.00585)
        iv,ix_,_ = s2mpj_ii("Y02005",ix_)
        arrset(pb.xnames,iv,"Y02005")
        ig = ig_["CAP01404"]
        pbm.A[ig,iv] += Float64(.00185)
        ig = ig_["CAP06105"]
        pbm.A[ig,iv] += Float64(.00067)
        iv,ix_,_ = s2mpj_ii("Y02005",ix_)
        arrset(pb.xnames,iv,"Y02005")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-3.12)
        ig = ig_["CAP01805"]
        pbm.A[ig,iv] += Float64(.00062)
        iv,ix_,_ = s2mpj_ii("Y02005",ix_)
        arrset(pb.xnames,iv,"Y02005")
        ig = ig_["CAP06205"]
        pbm.A[ig,iv] += Float64(.00452)
        ig = ig_["MXD01006"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y02005",ix_)
        arrset(pb.xnames,iv,"Y02005")
        ig = ig_["MND01006"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MXD01005"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y02005",ix_)
        arrset(pb.xnames,iv,"Y02005")
        ig = ig_["MND01005"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y02006",ix_)
        arrset(pb.xnames,iv,"Y02006")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(-2.96)
        ig = ig_["CAP06206"]
        pbm.A[ig,iv] += Float64(.00525)
        iv,ix_,_ = s2mpj_ii("Y02006",ix_)
        arrset(pb.xnames,iv,"Y02006")
        ig = ig_["CAP06106"]
        pbm.A[ig,iv] += Float64(.00072)
        ig = ig_["INV00505"]
        pbm.A[ig,iv] += Float64(-.33667)
        iv,ix_,_ = s2mpj_ii("Y02006",ix_)
        arrset(pb.xnames,iv,"Y02006")
        ig = ig_["MXD01006"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["MND01006"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y02006",ix_)
        arrset(pb.xnames,iv,"Y02006")
        ig = ig_["INV00506"]
        pbm.A[ig,iv] += Float64(-1.20241)
        ig = ig_["CAP01806"]
        pbm.A[ig,iv] += Float64(.00083)
        iv,ix_,_ = s2mpj_ii("Y02006",ix_)
        arrset(pb.xnames,iv,"Y02006")
        ig = ig_["CAP01405"]
        pbm.A[ig,iv] += Float64(.002)
        ig = ig_["CAP01406"]
        pbm.A[ig,iv] += Float64(.0077)
        iv,ix_,_ = s2mpj_ii("Y02006",ix_)
        arrset(pb.xnames,iv,"Y02006")
        ig = ig_["CAP06205"]
        pbm.A[ig,iv] += Float64(.0008)
        ig = ig_["CAP06105"]
        pbm.A[ig,iv] += Float64(.00006)
        iv,ix_,_ = s2mpj_ii("Y02006",ix_)
        arrset(pb.xnames,iv,"Y02006")
        ig = ig_["CAP01805"]
        pbm.A[ig,iv] += Float64(.00022)
        iv,ix_,_ = s2mpj_ii("X00103",ix_)
        arrset(pb.xnames,iv,"X00103")
        ig = ig_["CAP02103"]
        pbm.A[ig,iv] += Float64(.00677)
        ig = ig_["CAP01903"]
        pbm.A[ig,iv] += Float64(.02109)
        iv,ix_,_ = s2mpj_ii("X00103",ix_)
        arrset(pb.xnames,iv,"X00103")
        ig = ig_["CAP03101"]
        pbm.A[ig,iv] += Float64(.0033)
        ig = ig_["CAP02003"]
        pbm.A[ig,iv] += Float64(.0013)
        iv,ix_,_ = s2mpj_ii("X00103",ix_)
        arrset(pb.xnames,iv,"X00103")
        ig = ig_["CAP02001"]
        pbm.A[ig,iv] += Float64(.00072)
        ig = ig_["CAP02002"]
        pbm.A[ig,iv] += Float64(.00936)
        iv,ix_,_ = s2mpj_ii("X00103",ix_)
        arrset(pb.xnames,iv,"X00103")
        ig = ig_["CAP02102"]
        pbm.A[ig,iv] += Float64(.00546)
        ig = ig_["CAP03103"]
        pbm.A[ig,iv] += Float64(.0395)
        iv,ix_,_ = s2mpj_ii("X00103",ix_)
        arrset(pb.xnames,iv,"X00103")
        ig = ig_["CAP01902"]
        pbm.A[ig,iv] += Float64(.01285)
        ig = ig_["CAP02803"]
        pbm.A[ig,iv] += Float64(.01142)
        iv,ix_,_ = s2mpj_ii("X00103",ix_)
        arrset(pb.xnames,iv,"X00103")
        ig = ig_["CAP02703"]
        pbm.A[ig,iv] += Float64(.04825)
        ig = ig_["CAP02501"]
        pbm.A[ig,iv] += Float64(.00294)
        iv,ix_,_ = s2mpj_ii("X00103",ix_)
        arrset(pb.xnames,iv,"X00103")
        ig = ig_["CAP02701"]
        pbm.A[ig,iv] += Float64(.00356)
        ig = ig_["CAP02702"]
        pbm.A[ig,iv] += Float64(.10013)
        iv,ix_,_ = s2mpj_ii("X00103",ix_)
        arrset(pb.xnames,iv,"X00103")
        ig = ig_["CAP02802"]
        pbm.A[ig,iv] += Float64(.06479)
        ig = ig_["CAP02502"]
        pbm.A[ig,iv] += Float64(.01241)
        iv,ix_,_ = s2mpj_ii("X00103",ix_)
        arrset(pb.xnames,iv,"X00103")
        ig = ig_["INV00103"]
        pbm.A[ig,iv] += Float64(280.0)
        ig = ig_["CAP03102"]
        pbm.A[ig,iv] += Float64(.08414)
        iv,ix_,_ = s2mpj_ii("X00104",ix_)
        arrset(pb.xnames,iv,"X00104")
        ig = ig_["CAP02003"]
        pbm.A[ig,iv] += Float64(.00947)
        ig = ig_["CAP02002"]
        pbm.A[ig,iv] += Float64(.00042)
        iv,ix_,_ = s2mpj_ii("X00104",ix_)
        arrset(pb.xnames,iv,"X00104")
        ig = ig_["CAP02103"]
        pbm.A[ig,iv] += Float64(.00521)
        ig = ig_["CAP02702"]
        pbm.A[ig,iv] += Float64(.00226)
        iv,ix_,_ = s2mpj_ii("X00104",ix_)
        arrset(pb.xnames,iv,"X00104")
        ig = ig_["CAP02503"]
        pbm.A[ig,iv] += Float64(.01324)
        ig = ig_["CAP03102"]
        pbm.A[ig,iv] += Float64(.0021)
        iv,ix_,_ = s2mpj_ii("X00104",ix_)
        arrset(pb.xnames,iv,"X00104")
        ig = ig_["INV00104"]
        pbm.A[ig,iv] += Float64(280.0)
        ig = ig_["CAP01903"]
        pbm.A[ig,iv] += Float64(.01226)
        iv,ix_,_ = s2mpj_ii("X00104",ix_)
        arrset(pb.xnames,iv,"X00104")
        ig = ig_["CAP02502"]
        pbm.A[ig,iv] += Float64(.00211)
        ig = ig_["CAP03103"]
        pbm.A[ig,iv] += Float64(.08241)
        iv,ix_,_ = s2mpj_ii("X00104",ix_)
        arrset(pb.xnames,iv,"X00104")
        ig = ig_["CAP02703"]
        pbm.A[ig,iv] += Float64(.09784)
        ig = ig_["CAP01904"]
        pbm.A[ig,iv] += Float64(.02167)
        iv,ix_,_ = s2mpj_ii("X00104",ix_)
        arrset(pb.xnames,iv,"X00104")
        ig = ig_["CAP02104"]
        pbm.A[ig,iv] += Float64(.00702)
        ig = ig_["CAP02704"]
        pbm.A[ig,iv] += Float64(.05183)
        iv,ix_,_ = s2mpj_ii("X00104",ix_)
        arrset(pb.xnames,iv,"X00104")
        ig = ig_["CAP02804"]
        pbm.A[ig,iv] += Float64(.01437)
        ig = ig_["CAP03104"]
        pbm.A[ig,iv] += Float64(.04242)
        iv,ix_,_ = s2mpj_ii("X00104",ix_)
        arrset(pb.xnames,iv,"X00104")
        ig = ig_["CAP02004"]
        pbm.A[ig,iv] += Float64(.00149)
        ig = ig_["CAP02803"]
        pbm.A[ig,iv] += Float64(.06185)
        iv,ix_,_ = s2mpj_ii("X00105",ix_)
        arrset(pb.xnames,iv,"X00105")
        ig = ig_["CAP02704"]
        pbm.A[ig,iv] += Float64(.09897)
        ig = ig_["CAP02703"]
        pbm.A[ig,iv] += Float64(.00113)
        iv,ix_,_ = s2mpj_ii("X00105",ix_)
        arrset(pb.xnames,iv,"X00105")
        ig = ig_["CAP01904"]
        pbm.A[ig,iv] += Float64(.01226)
        ig = ig_["CAP02705"]
        pbm.A[ig,iv] += Float64(.05183)
        iv,ix_,_ = s2mpj_ii("X00105",ix_)
        arrset(pb.xnames,iv,"X00105")
        ig = ig_["CAP02503"]
        pbm.A[ig,iv] += Float64(.00141)
        ig = ig_["CAP02104"]
        pbm.A[ig,iv] += Float64(.00521)
        iv,ix_,_ = s2mpj_ii("X00105",ix_)
        arrset(pb.xnames,iv,"X00105")
        ig = ig_["CAP02504"]
        pbm.A[ig,iv] += Float64(.01394)
        ig = ig_["CAP03103"]
        pbm.A[ig,iv] += Float64(.00105)
        iv,ix_,_ = s2mpj_ii("X00105",ix_)
        arrset(pb.xnames,iv,"X00105")
        ig = ig_["CAP03105"]
        pbm.A[ig,iv] += Float64(.04242)
        ig = ig_["CAP02804"]
        pbm.A[ig,iv] += Float64(.06185)
        iv,ix_,_ = s2mpj_ii("X00105",ix_)
        arrset(pb.xnames,iv,"X00105")
        ig = ig_["CAP02805"]
        pbm.A[ig,iv] += Float64(.01437)
        ig = ig_["CAP03104"]
        pbm.A[ig,iv] += Float64(.08346)
        iv,ix_,_ = s2mpj_ii("X00105",ix_)
        arrset(pb.xnames,iv,"X00105")
        ig = ig_["CAP01905"]
        pbm.A[ig,iv] += Float64(.02167)
        ig = ig_["CAP02003"]
        pbm.A[ig,iv] += Float64(.00016)
        iv,ix_,_ = s2mpj_ii("X00105",ix_)
        arrset(pb.xnames,iv,"X00105")
        ig = ig_["CAP02005"]
        pbm.A[ig,iv] += Float64(.00149)
        ig = ig_["CAP02105"]
        pbm.A[ig,iv] += Float64(.00702)
        iv,ix_,_ = s2mpj_ii("X00105",ix_)
        arrset(pb.xnames,iv,"X00105")
        ig = ig_["CAP02004"]
        pbm.A[ig,iv] += Float64(.00973)
        ig = ig_["INV00105"]
        pbm.A[ig,iv] += Float64(280.0)
        iv,ix_,_ = s2mpj_ii("X00106",ix_)
        arrset(pb.xnames,iv,"X00106")
        ig = ig_["CAP03105"]
        pbm.A[ig,iv] += Float64(.08854)
        ig = ig_["CAP02705"]
        pbm.A[ig,iv] += Float64(.10487)
        iv,ix_,_ = s2mpj_ii("X00106",ix_)
        arrset(pb.xnames,iv,"X00106")
        ig = ig_["INV00106"]
        pbm.A[ig,iv] += Float64(280.0)
        ig = ig_["CAP02805"]
        pbm.A[ig,iv] += Float64(.06479)
        iv,ix_,_ = s2mpj_ii("X00106",ix_)
        arrset(pb.xnames,iv,"X00106")
        ig = ig_["CAP02505"]
        pbm.A[ig,iv] += Float64(.01608)
        ig = ig_["CAP01906"]
        pbm.A[ig,iv] += Float64(.03394)
        iv,ix_,_ = s2mpj_ii("X00106",ix_)
        arrset(pb.xnames,iv,"X00106")
        ig = ig_["CAP02006"]
        pbm.A[ig,iv] += Float64(.01139)
        ig = ig_["CAP03106"]
        pbm.A[ig,iv] += Float64(.12694)
        iv,ix_,_ = s2mpj_ii("X00106",ix_)
        arrset(pb.xnames,iv,"X00106")
        ig = ig_["CAP02106"]
        pbm.A[ig,iv] += Float64(.01223)
        ig = ig_["CAP02506"]
        pbm.A[ig,iv] += Float64(.01535)
        iv,ix_,_ = s2mpj_ii("X00106",ix_)
        arrset(pb.xnames,iv,"X00106")
        ig = ig_["CAP02806"]
        pbm.A[ig,iv] += Float64(.07621)
        ig = ig_["CAP02706"]
        pbm.A[ig,iv] += Float64(.15193)
        iv,ix_,_ = s2mpj_ii("X00106",ix_)
        arrset(pb.xnames,iv,"X00106")
        ig = ig_["CAP02105"]
        pbm.A[ig,iv] += Float64(.00546)
        ig = ig_["CAP02005"]
        pbm.A[ig,iv] += Float64(.01036)
        iv,ix_,_ = s2mpj_ii("X00106",ix_)
        arrset(pb.xnames,iv,"X00106")
        ig = ig_["CAP02004"]
        pbm.A[ig,iv] += Float64(.00016)
        ig = ig_["CAP03104"]
        pbm.A[ig,iv] += Float64(.0011)
        iv,ix_,_ = s2mpj_ii("X00106",ix_)
        arrset(pb.xnames,iv,"X00106")
        ig = ig_["CAP02504"]
        pbm.A[ig,iv] += Float64(.00148)
        ig = ig_["CAP02704"]
        pbm.A[ig,iv] += Float64(.00119)
        iv,ix_,_ = s2mpj_ii("X00106",ix_)
        arrset(pb.xnames,iv,"X00106")
        ig = ig_["CAP01905"]
        pbm.A[ig,iv] += Float64(.01285)
        iv,ix_,_ = s2mpj_ii("X00203",ix_)
        arrset(pb.xnames,iv,"X00203")
        ig = ig_["CAP02102"]
        pbm.A[ig,iv] += Float64(.00766)
        ig = ig_["CAP01903"]
        pbm.A[ig,iv] += Float64(.01492)
        iv,ix_,_ = s2mpj_ii("X00203",ix_)
        arrset(pb.xnames,iv,"X00203")
        ig = ig_["CAP02003"]
        pbm.A[ig,iv] += Float64(.00051)
        ig = ig_["CAP02701"]
        pbm.A[ig,iv] += Float64(.00333)
        iv,ix_,_ = s2mpj_ii("X00203",ix_)
        arrset(pb.xnames,iv,"X00203")
        ig = ig_["CAP02501"]
        pbm.A[ig,iv] += Float64(.00296)
        ig = ig_["CAP02502"]
        pbm.A[ig,iv] += Float64(.01161)
        iv,ix_,_ = s2mpj_ii("X00203",ix_)
        arrset(pb.xnames,iv,"X00203")
        ig = ig_["CAP03203"]
        pbm.A[ig,iv] += Float64(.01487)
        ig = ig_["CAP03101"]
        pbm.A[ig,iv] += Float64(.00008)
        iv,ix_,_ = s2mpj_ii("X00203",ix_)
        arrset(pb.xnames,iv,"X00203")
        ig = ig_["CAP02002"]
        pbm.A[ig,iv] += Float64(.0098)
        ig = ig_["CAP02103"]
        pbm.A[ig,iv] += Float64(.00459)
        iv,ix_,_ = s2mpj_ii("X00203",ix_)
        arrset(pb.xnames,iv,"X00203")
        ig = ig_["CAP02303"]
        pbm.A[ig,iv] += Float64(.0206)
        ig = ig_["CAP01902"]
        pbm.A[ig,iv] += Float64(.01768)
        iv,ix_,_ = s2mpj_ii("X00203",ix_)
        arrset(pb.xnames,iv,"X00203")
        ig = ig_["CAP02603"]
        pbm.A[ig,iv] += Float64(.02157)
        ig = ig_["CAP02703"]
        pbm.A[ig,iv] += Float64(.02094)
        iv,ix_,_ = s2mpj_ii("X00203",ix_)
        arrset(pb.xnames,iv,"X00203")
        ig = ig_["CAP02803"]
        pbm.A[ig,iv] += Float64(.00342)
        ig = ig_["CAP03103"]
        pbm.A[ig,iv] += Float64(.03082)
        iv,ix_,_ = s2mpj_ii("X00203",ix_)
        arrset(pb.xnames,iv,"X00203")
        ig = ig_["CAP02001"]
        pbm.A[ig,iv] += Float64(.00079)
        ig = ig_["CAP03102"]
        pbm.A[ig,iv] += Float64(.06879)
        iv,ix_,_ = s2mpj_ii("X00203",ix_)
        arrset(pb.xnames,iv,"X00203")
        ig = ig_["INV00203"]
        pbm.A[ig,iv] += Float64(75.0)
        ig = ig_["CAP03202"]
        pbm.A[ig,iv] += Float64(.02505)
        iv,ix_,_ = s2mpj_ii("X00203",ix_)
        arrset(pb.xnames,iv,"X00203")
        ig = ig_["CAP02802"]
        pbm.A[ig,iv] += Float64(.01722)
        ig = ig_["CAP02602"]
        pbm.A[ig,iv] += Float64(.02373)
        iv,ix_,_ = s2mpj_ii("X00203",ix_)
        arrset(pb.xnames,iv,"X00203")
        ig = ig_["CAP02702"]
        pbm.A[ig,iv] += Float64(.08608)
        ig = ig_["CAP02302"]
        pbm.A[ig,iv] += Float64(.01551)
        iv,ix_,_ = s2mpj_ii("X00204",ix_)
        arrset(pb.xnames,iv,"X00204")
        ig = ig_["CAP03203"]
        pbm.A[ig,iv] += Float64(.02391)
        ig = ig_["CAP03103"]
        pbm.A[ig,iv] += Float64(.06574)
        iv,ix_,_ = s2mpj_ii("X00204",ix_)
        arrset(pb.xnames,iv,"X00204")
        ig = ig_["CAP03204"]
        pbm.A[ig,iv] += Float64(.01601)
        ig = ig_["CAP02803"]
        pbm.A[ig,iv] += Float64(.01644)
        iv,ix_,_ = s2mpj_ii("X00204",ix_)
        arrset(pb.xnames,iv,"X00204")
        ig = ig_["CAP02703"]
        pbm.A[ig,iv] += Float64(.08418)
        ig = ig_["CAP03104"]
        pbm.A[ig,iv] += Float64(.03395)
        iv,ix_,_ = s2mpj_ii("X00204",ix_)
        arrset(pb.xnames,iv,"X00204")
        ig = ig_["CAP02603"]
        pbm.A[ig,iv] += Float64(.02265)
        ig = ig_["CAP02503"]
        pbm.A[ig,iv] += Float64(.01241)
        iv,ix_,_ = s2mpj_ii("X00204",ix_)
        arrset(pb.xnames,iv,"X00204")
        ig = ig_["CAP02303"]
        pbm.A[ig,iv] += Float64(.01481)
        ig = ig_["CAP02103"]
        pbm.A[ig,iv] += Float64(.00731)
        iv,ix_,_ = s2mpj_ii("X00204",ix_)
        arrset(pb.xnames,iv,"X00204")
        ig = ig_["CAP02003"]
        pbm.A[ig,iv] += Float64(.00987)
        ig = ig_["CAP02002"]
        pbm.A[ig,iv] += Float64(.0005)
        iv,ix_,_ = s2mpj_ii("X00204",ix_)
        arrset(pb.xnames,iv,"X00204")
        ig = ig_["CAP01903"]
        pbm.A[ig,iv] += Float64(.01688)
        ig = ig_["CAP02702"]
        pbm.A[ig,iv] += Float64(.00217)
        iv,ix_,_ = s2mpj_ii("X00204",ix_)
        arrset(pb.xnames,iv,"X00204")
        ig = ig_["CAP02604"]
        pbm.A[ig,iv] += Float64(.02265)
        ig = ig_["CAP02804"]
        pbm.A[ig,iv] += Float64(.0042)
        iv,ix_,_ = s2mpj_ii("X00204",ix_)
        arrset(pb.xnames,iv,"X00204")
        ig = ig_["CAP02502"]
        pbm.A[ig,iv] += Float64(.00216)
        ig = ig_["INV00204"]
        pbm.A[ig,iv] += Float64(75.0)
        iv,ix_,_ = s2mpj_ii("X00204",ix_)
        arrset(pb.xnames,iv,"X00204")
        ig = ig_["CAP01904"]
        pbm.A[ig,iv] += Float64(.01572)
        ig = ig_["CAP02304"]
        pbm.A[ig,iv] += Float64(.02131)
        iv,ix_,_ = s2mpj_ii("X00204",ix_)
        arrset(pb.xnames,iv,"X00204")
        ig = ig_["CAP02104"]
        pbm.A[ig,iv] += Float64(.00494)
        ig = ig_["CAP02004"]
        pbm.A[ig,iv] += Float64(.00073)
        iv,ix_,_ = s2mpj_ii("X00204",ix_)
        arrset(pb.xnames,iv,"X00204")
        ig = ig_["CAP02704"]
        pbm.A[ig,iv] += Float64(.024)
        iv,ix_,_ = s2mpj_ii("X00205",ix_)
        arrset(pb.xnames,iv,"X00205")
        ig = ig_["CAP02503"]
        pbm.A[ig,iv] += Float64(.0015)
        ig = ig_["CAP02705"]
        pbm.A[ig,iv] += Float64(.024)
        iv,ix_,_ = s2mpj_ii("X00205",ix_)
        arrset(pb.xnames,iv,"X00205")
        ig = ig_["CAP02104"]
        pbm.A[ig,iv] += Float64(.00731)
        ig = ig_["CAP02703"]
        pbm.A[ig,iv] += Float64(.00117)
        iv,ix_,_ = s2mpj_ii("X00205",ix_)
        arrset(pb.xnames,iv,"X00205")
        ig = ig_["INV00205"]
        pbm.A[ig,iv] += Float64(75.0)
        ig = ig_["CAP03205"]
        pbm.A[ig,iv] += Float64(.01601)
        iv,ix_,_ = s2mpj_ii("X00205",ix_)
        arrset(pb.xnames,iv,"X00205")
        ig = ig_["CAP02004"]
        pbm.A[ig,iv] += Float64(.01012)
        ig = ig_["CAP02605"]
        pbm.A[ig,iv] += Float64(.02265)
        iv,ix_,_ = s2mpj_ii("X00205",ix_)
        arrset(pb.xnames,iv,"X00205")
        ig = ig_["CAP03105"]
        pbm.A[ig,iv] += Float64(.03395)
        ig = ig_["CAP02805"]
        pbm.A[ig,iv] += Float64(.0042)
        iv,ix_,_ = s2mpj_ii("X00205",ix_)
        arrset(pb.xnames,iv,"X00205")
        ig = ig_["CAP02003"]
        pbm.A[ig,iv] += Float64(.00025)
        ig = ig_["CAP01904"]
        pbm.A[ig,iv] += Float64(.01688)
        iv,ix_,_ = s2mpj_ii("X00205",ix_)
        arrset(pb.xnames,iv,"X00205")
        ig = ig_["CAP02804"]
        pbm.A[ig,iv] += Float64(.01644)
        ig = ig_["CAP02005"]
        pbm.A[ig,iv] += Float64(.00073)
        iv,ix_,_ = s2mpj_ii("X00205",ix_)
        arrset(pb.xnames,iv,"X00205")
        ig = ig_["CAP02105"]
        pbm.A[ig,iv] += Float64(.00494)
        ig = ig_["CAP01905"]
        pbm.A[ig,iv] += Float64(.01572)
        iv,ix_,_ = s2mpj_ii("X00205",ix_)
        arrset(pb.xnames,iv,"X00205")
        ig = ig_["CAP02704"]
        pbm.A[ig,iv] += Float64(.08519)
        ig = ig_["CAP03204"]
        pbm.A[ig,iv] += Float64(.02391)
        iv,ix_,_ = s2mpj_ii("X00205",ix_)
        arrset(pb.xnames,iv,"X00205")
        ig = ig_["CAP02604"]
        pbm.A[ig,iv] += Float64(.02265)
        ig = ig_["CAP02305"]
        pbm.A[ig,iv] += Float64(.02131)
        iv,ix_,_ = s2mpj_ii("X00205",ix_)
        arrset(pb.xnames,iv,"X00205")
        ig = ig_["CAP03104"]
        pbm.A[ig,iv] += Float64(.06574)
        ig = ig_["CAP02504"]
        pbm.A[ig,iv] += Float64(.01307)
        iv,ix_,_ = s2mpj_ii("X00205",ix_)
        arrset(pb.xnames,iv,"X00205")
        ig = ig_["CAP02304"]
        pbm.A[ig,iv] += Float64(.01481)
        iv,ix_,_ = s2mpj_ii("X00206",ix_)
        arrset(pb.xnames,iv,"X00206")
        ig = ig_["CAP03206"]
        pbm.A[ig,iv] += Float64(.03992)
        ig = ig_["CAP02704"]
        pbm.A[ig,iv] += Float64(.00122)
        iv,ix_,_ = s2mpj_ii("X00206",ix_)
        arrset(pb.xnames,iv,"X00206")
        ig = ig_["CAP02805"]
        pbm.A[ig,iv] += Float64(.01722)
        ig = ig_["CAP01905"]
        pbm.A[ig,iv] += Float64(.01768)
        iv,ix_,_ = s2mpj_ii("X00206",ix_)
        arrset(pb.xnames,iv,"X00206")
        ig = ig_["CAP02505"]
        pbm.A[ig,iv] += Float64(.01526)
        ig = ig_["CAP02005"]
        pbm.A[ig,iv] += Float64(.01086)
        iv,ix_,_ = s2mpj_ii("X00206",ix_)
        arrset(pb.xnames,iv,"X00206")
        ig = ig_["CAP02705"]
        pbm.A[ig,iv] += Float64(.09047)
        ig = ig_["CAP02004"]
        pbm.A[ig,iv] += Float64(.00026)
        iv,ix_,_ = s2mpj_ii("X00206",ix_)
        arrset(pb.xnames,iv,"X00206")
        ig = ig_["CAP02105"]
        pbm.A[ig,iv] += Float64(.00766)
        ig = ig_["INV00206"]
        pbm.A[ig,iv] += Float64(75.0)
        iv,ix_,_ = s2mpj_ii("X00206",ix_)
        arrset(pb.xnames,iv,"X00206")
        ig = ig_["CAP02605"]
        pbm.A[ig,iv] += Float64(.02373)
        ig = ig_["CAP02305"]
        pbm.A[ig,iv] += Float64(.01551)
        iv,ix_,_ = s2mpj_ii("X00206",ix_)
        arrset(pb.xnames,iv,"X00206")
        ig = ig_["CAP02504"]
        pbm.A[ig,iv] += Float64(.00157)
        ig = ig_["CAP02806"]
        pbm.A[ig,iv] += Float64(.02064)
        iv,ix_,_ = s2mpj_ii("X00206",ix_)
        arrset(pb.xnames,iv,"X00206")
        ig = ig_["CAP03105"]
        pbm.A[ig,iv] += Float64(.06887)
        ig = ig_["CAP02506"]
        pbm.A[ig,iv] += Float64(.01457)
        iv,ix_,_ = s2mpj_ii("X00206",ix_)
        arrset(pb.xnames,iv,"X00206")
        ig = ig_["CAP01906"]
        pbm.A[ig,iv] += Float64(.0326)
        ig = ig_["CAP02306"]
        pbm.A[ig,iv] += Float64(.03611)
        iv,ix_,_ = s2mpj_ii("X00206",ix_)
        arrset(pb.xnames,iv,"X00206")
        ig = ig_["CAP02706"]
        pbm.A[ig,iv] += Float64(.11036)
        ig = ig_["CAP02606"]
        pbm.A[ig,iv] += Float64(.04529)
        iv,ix_,_ = s2mpj_ii("X00206",ix_)
        arrset(pb.xnames,iv,"X00206")
        ig = ig_["CAP02006"]
        pbm.A[ig,iv] += Float64(.0111)
        ig = ig_["CAP02106"]
        pbm.A[ig,iv] += Float64(.01225)
        iv,ix_,_ = s2mpj_ii("X00206",ix_)
        arrset(pb.xnames,iv,"X00206")
        ig = ig_["CAP03106"]
        pbm.A[ig,iv] += Float64(.09969)
        ig = ig_["CAP03205"]
        pbm.A[ig,iv] += Float64(.02505)
        iv,ix_,_ = s2mpj_ii("X00303",ix_)
        arrset(pb.xnames,iv,"X00303")
        ig = ig_["CAP03002"]
        pbm.A[ig,iv] += Float64(.00998)
        ig = ig_["CAP02501"]
        pbm.A[ig,iv] += Float64(.00285)
        iv,ix_,_ = s2mpj_ii("X00303",ix_)
        arrset(pb.xnames,iv,"X00303")
        ig = ig_["CAP02202"]
        pbm.A[ig,iv] += Float64(.12349)
        ig = ig_["CAP02403"]
        pbm.A[ig,iv] += Float64(.00837)
        iv,ix_,_ = s2mpj_ii("X00303",ix_)
        arrset(pb.xnames,iv,"X00303")
        ig = ig_["INV00303"]
        pbm.A[ig,iv] += Float64(43.0)
        ig = ig_["CAP02903"]
        pbm.A[ig,iv] += Float64(.00719)
        iv,ix_,_ = s2mpj_ii("X00303",ix_)
        arrset(pb.xnames,iv,"X00303")
        ig = ig_["CAP02201"]
        pbm.A[ig,iv] += Float64(.00077)
        ig = ig_["CAP02803"]
        pbm.A[ig,iv] += Float64(.00392)
        iv,ix_,_ = s2mpj_ii("X00303",ix_)
        arrset(pb.xnames,iv,"X00303")
        ig = ig_["CAP02502"]
        pbm.A[ig,iv] += Float64(.01187)
        ig = ig_["CAP02902"]
        pbm.A[ig,iv] += Float64(.01593)
        iv,ix_,_ = s2mpj_ii("X00303",ix_)
        arrset(pb.xnames,iv,"X00303")
        ig = ig_["CAP02603"]
        pbm.A[ig,iv] += Float64(.028)
        ig = ig_["CAP02702"]
        pbm.A[ig,iv] += Float64(.11034)
        iv,ix_,_ = s2mpj_ii("X00303",ix_)
        arrset(pb.xnames,iv,"X00303")
        ig = ig_["CAP02703"]
        pbm.A[ig,iv] += Float64(.02005)
        ig = ig_["CAP02001"]
        pbm.A[ig,iv] += Float64(.00079)
        iv,ix_,_ = s2mpj_ii("X00303",ix_)
        arrset(pb.xnames,iv,"X00303")
        ig = ig_["CAP02802"]
        pbm.A[ig,iv] += Float64(.01664)
        ig = ig_["CAP02203"]
        pbm.A[ig,iv] += Float64(.02513)
        iv,ix_,_ = s2mpj_ii("X00303",ix_)
        arrset(pb.xnames,iv,"X00303")
        ig = ig_["CAP02102"]
        pbm.A[ig,iv] += Float64(.02041)
        ig = ig_["CAP02002"]
        pbm.A[ig,iv] += Float64(.00977)
        iv,ix_,_ = s2mpj_ii("X00303",ix_)
        arrset(pb.xnames,iv,"X00303")
        ig = ig_["CAP02003"]
        pbm.A[ig,iv] += Float64(.00054)
        ig = ig_["CAP03003"]
        pbm.A[ig,iv] += Float64(.01077)
        iv,ix_,_ = s2mpj_ii("X00303",ix_)
        arrset(pb.xnames,iv,"X00303")
        ig = ig_["CAP03202"]
        pbm.A[ig,iv] += Float64(.01831)
        ig = ig_["CAP02701"]
        pbm.A[ig,iv] += Float64(.00576)
        iv,ix_,_ = s2mpj_ii("X00303",ix_)
        arrset(pb.xnames,iv,"X00303")
        ig = ig_["CAP03102"]
        pbm.A[ig,iv] += Float64(.00594)
        ig = ig_["CAP02402"]
        pbm.A[ig,iv] += Float64(.01363)
        iv,ix_,_ = s2mpj_ii("X00303",ix_)
        arrset(pb.xnames,iv,"X00303")
        ig = ig_["CAP02602"]
        pbm.A[ig,iv] += Float64(.01608)
        ig = ig_["CAP02103"]
        pbm.A[ig,iv] += Float64(.02973)
        iv,ix_,_ = s2mpj_ii("X00303",ix_)
        arrset(pb.xnames,iv,"X00303")
        ig = ig_["CAP03203"]
        pbm.A[ig,iv] += Float64(.02077)
        ig = ig_["CAP01902"]
        pbm.A[ig,iv] += Float64(.01255)
        iv,ix_,_ = s2mpj_ii("X00303",ix_)
        arrset(pb.xnames,iv,"X00303")
        ig = ig_["CAP01903"]
        pbm.A[ig,iv] += Float64(.01979)
        ig = ig_["CAP03103"]
        pbm.A[ig,iv] += Float64(.01198)
        iv,ix_,_ = s2mpj_ii("X00304",ix_)
        arrset(pb.xnames,iv,"X00304")
        ig = ig_["CAP02502"]
        pbm.A[ig,iv] += Float64(.00205)
        ig = ig_["CAP02204"]
        pbm.A[ig,iv] += Float64(.03078)
        iv,ix_,_ = s2mpj_ii("X00304",ix_)
        arrset(pb.xnames,iv,"X00304")
        ig = ig_["CAP02703"]
        pbm.A[ig,iv] += Float64(.1083)
        ig = ig_["CAP02104"]
        pbm.A[ig,iv] += Float64(.03066)
        iv,ix_,_ = s2mpj_ii("X00304",ix_)
        arrset(pb.xnames,iv,"X00304")
        ig = ig_["CAP02903"]
        pbm.A[ig,iv] += Float64(.0152)
        ig = ig_["CAP02404"]
        pbm.A[ig,iv] += Float64(.00899)
        iv,ix_,_ = s2mpj_ii("X00304",ix_)
        arrset(pb.xnames,iv,"X00304")
        ig = ig_["CAP02702"]
        pbm.A[ig,iv] += Float64(.00376)
        ig = ig_["CAP03204"]
        pbm.A[ig,iv] += Float64(.0216)
        iv,ix_,_ = s2mpj_ii("X00304",ix_)
        arrset(pb.xnames,iv,"X00304")
        ig = ig_["CAP02004"]
        pbm.A[ig,iv] += Float64(.00077)
        ig = ig_["CAP03003"]
        pbm.A[ig,iv] += Float64(.00953)
        iv,ix_,_ = s2mpj_ii("X00304",ix_)
        arrset(pb.xnames,iv,"X00304")
        ig = ig_["CAP02803"]
        pbm.A[ig,iv] += Float64(.01589)
        ig = ig_["CAP03203"]
        pbm.A[ig,iv] += Float64(.01748)
        iv,ix_,_ = s2mpj_ii("X00304",ix_)
        arrset(pb.xnames,iv,"X00304")
        ig = ig_["CAP03103"]
        pbm.A[ig,iv] += Float64(.00567)
        ig = ig_["CAP02403"]
        pbm.A[ig,iv] += Float64(.01301)
        iv,ix_,_ = s2mpj_ii("X00304",ix_)
        arrset(pb.xnames,iv,"X00304")
        ig = ig_["CAP01903"]
        pbm.A[ig,iv] += Float64(.01198)
        ig = ig_["CAP02002"]
        pbm.A[ig,iv] += Float64(.0005)
        iv,ix_,_ = s2mpj_ii("X00304",ix_)
        arrset(pb.xnames,iv,"X00304")
        ig = ig_["CAP02904"]
        pbm.A[ig,iv] += Float64(.00791)
        ig = ig_["CAP02804"]
        pbm.A[ig,iv] += Float64(.00467)
        iv,ix_,_ = s2mpj_ii("X00304",ix_)
        arrset(pb.xnames,iv,"X00304")
        ig = ig_["CAP02003"]
        pbm.A[ig,iv] += Float64(.00983)
        ig = ig_["CAP03104"]
        pbm.A[ig,iv] += Float64(.01225)
        iv,ix_,_ = s2mpj_ii("X00304",ix_)
        arrset(pb.xnames,iv,"X00304")
        ig = ig_["CAP02103"]
        pbm.A[ig,iv] += Float64(.01948)
        ig = ig_["CAP02203"]
        pbm.A[ig,iv] += Float64(.11861)
        iv,ix_,_ = s2mpj_ii("X00304",ix_)
        arrset(pb.xnames,iv,"X00304")
        ig = ig_["CAP02704"]
        pbm.A[ig,iv] += Float64(.02409)
        ig = ig_["CAP03004"]
        pbm.A[ig,iv] += Float64(.01122)
        iv,ix_,_ = s2mpj_ii("X00304",ix_)
        arrset(pb.xnames,iv,"X00304")
        ig = ig_["CAP01904"]
        pbm.A[ig,iv] += Float64(.02036)
        ig = ig_["CAP02604"]
        pbm.A[ig,iv] += Float64(.02873)
        iv,ix_,_ = s2mpj_ii("X00304",ix_)
        arrset(pb.xnames,iv,"X00304")
        ig = ig_["INV00304"]
        pbm.A[ig,iv] += Float64(43.0)
        ig = ig_["CAP02503"]
        pbm.A[ig,iv] += Float64(.01266)
        iv,ix_,_ = s2mpj_ii("X00304",ix_)
        arrset(pb.xnames,iv,"X00304")
        ig = ig_["CAP02603"]
        pbm.A[ig,iv] += Float64(.01535)
        iv,ix_,_ = s2mpj_ii("X00305",ix_)
        arrset(pb.xnames,iv,"X00305")
        ig = ig_["CAP01905"]
        pbm.A[ig,iv] += Float64(.02036)
        ig = ig_["CAP02504"]
        pbm.A[ig,iv] += Float64(.01333)
        iv,ix_,_ = s2mpj_ii("X00305",ix_)
        arrset(pb.xnames,iv,"X00305")
        ig = ig_["CAP02205"]
        pbm.A[ig,iv] += Float64(.03078)
        ig = ig_["CAP02904"]
        pbm.A[ig,iv] += Float64(.0152)
        iv,ix_,_ = s2mpj_ii("X00305",ix_)
        arrset(pb.xnames,iv,"X00305")
        ig = ig_["CAP02204"]
        pbm.A[ig,iv] += Float64(.11861)
        ig = ig_["CAP02005"]
        pbm.A[ig,iv] += Float64(.00077)
        iv,ix_,_ = s2mpj_ii("X00305",ix_)
        arrset(pb.xnames,iv,"X00305")
        ig = ig_["CAP02404"]
        pbm.A[ig,iv] += Float64(.01301)
        ig = ig_["CAP02804"]
        pbm.A[ig,iv] += Float64(.01589)
        iv,ix_,_ = s2mpj_ii("X00305",ix_)
        arrset(pb.xnames,iv,"X00305")
        ig = ig_["CAP02604"]
        pbm.A[ig,iv] += Float64(.01535)
        ig = ig_["CAP02704"]
        pbm.A[ig,iv] += Float64(.10954)
        iv,ix_,_ = s2mpj_ii("X00305",ix_)
        arrset(pb.xnames,iv,"X00305")
        ig = ig_["CAP02105"]
        pbm.A[ig,iv] += Float64(.03066)
        ig = ig_["CAP03104"]
        pbm.A[ig,iv] += Float64(.00567)
        iv,ix_,_ = s2mpj_ii("X00305",ix_)
        arrset(pb.xnames,iv,"X00305")
        ig = ig_["CAP03004"]
        pbm.A[ig,iv] += Float64(.00953)
        ig = ig_["CAP02503"]
        pbm.A[ig,iv] += Float64(.00138)
        iv,ix_,_ = s2mpj_ii("X00305",ix_)
        arrset(pb.xnames,iv,"X00305")
        ig = ig_["CAP02605"]
        pbm.A[ig,iv] += Float64(.02873)
        ig = ig_["CAP03205"]
        pbm.A[ig,iv] += Float64(.0216)
        iv,ix_,_ = s2mpj_ii("X00305",ix_)
        arrset(pb.xnames,iv,"X00305")
        ig = ig_["CAP03105"]
        pbm.A[ig,iv] += Float64(.01225)
        ig = ig_["CAP02104"]
        pbm.A[ig,iv] += Float64(.01948)
        iv,ix_,_ = s2mpj_ii("X00305",ix_)
        arrset(pb.xnames,iv,"X00305")
        ig = ig_["CAP02703"]
        pbm.A[ig,iv] += Float64(.00253)
        ig = ig_["CAP03204"]
        pbm.A[ig,iv] += Float64(.01748)
        iv,ix_,_ = s2mpj_ii("X00305",ix_)
        arrset(pb.xnames,iv,"X00305")
        ig = ig_["CAP02003"]
        pbm.A[ig,iv] += Float64(.00025)
        ig = ig_["CAP02705"]
        pbm.A[ig,iv] += Float64(.02409)
        iv,ix_,_ = s2mpj_ii("X00305",ix_)
        arrset(pb.xnames,iv,"X00305")
        ig = ig_["CAP02004"]
        pbm.A[ig,iv] += Float64(.01009)
        ig = ig_["CAP01904"]
        pbm.A[ig,iv] += Float64(.01198)
        iv,ix_,_ = s2mpj_ii("X00305",ix_)
        arrset(pb.xnames,iv,"X00305")
        ig = ig_["CAP02905"]
        pbm.A[ig,iv] += Float64(.00791)
        ig = ig_["CAP02405"]
        pbm.A[ig,iv] += Float64(.00899)
        iv,ix_,_ = s2mpj_ii("X00305",ix_)
        arrset(pb.xnames,iv,"X00305")
        ig = ig_["CAP03005"]
        pbm.A[ig,iv] += Float64(.01122)
        ig = ig_["CAP02805"]
        pbm.A[ig,iv] += Float64(.00467)
        iv,ix_,_ = s2mpj_ii("X00305",ix_)
        arrset(pb.xnames,iv,"X00305")
        ig = ig_["INV00305"]
        pbm.A[ig,iv] += Float64(43.0)
        iv,ix_,_ = s2mpj_ii("X00306",ix_)
        arrset(pb.xnames,iv,"X00306")
        ig = ig_["CAP03206"]
        pbm.A[ig,iv] += Float64(.03908)
        ig = ig_["CAP03006"]
        pbm.A[ig,iv] += Float64(.02075)
        iv,ix_,_ = s2mpj_ii("X00306",ix_)
        arrset(pb.xnames,iv,"X00306")
        ig = ig_["CAP02106"]
        pbm.A[ig,iv] += Float64(.05014)
        ig = ig_["CAP02206"]
        pbm.A[ig,iv] += Float64(.14939)
        iv,ix_,_ = s2mpj_ii("X00306",ix_)
        arrset(pb.xnames,iv,"X00306")
        ig = ig_["CAP03106"]
        pbm.A[ig,iv] += Float64(.01793)
        ig = ig_["CAP02705"]
        pbm.A[ig,iv] += Float64(.1174)
        iv,ix_,_ = s2mpj_ii("X00306",ix_)
        arrset(pb.xnames,iv,"X00306")
        ig = ig_["CAP03105"]
        pbm.A[ig,iv] += Float64(.00594)
        ig = ig_["CAP02704"]
        pbm.A[ig,iv] += Float64(.00265)
        iv,ix_,_ = s2mpj_ii("X00306",ix_)
        arrset(pb.xnames,iv,"X00306")
        ig = ig_["CAP02806"]
        pbm.A[ig,iv] += Float64(.02056)
        ig = ig_["CAP02006"]
        pbm.A[ig,iv] += Float64(.0111)
        iv,ix_,_ = s2mpj_ii("X00306",ix_)
        arrset(pb.xnames,iv,"X00306")
        ig = ig_["CAP03005"]
        pbm.A[ig,iv] += Float64(.00998)
        ig = ig_["CAP02905"]
        pbm.A[ig,iv] += Float64(.01593)
        iv,ix_,_ = s2mpj_ii("X00306",ix_)
        arrset(pb.xnames,iv,"X00306")
        ig = ig_["CAP02805"]
        pbm.A[ig,iv] += Float64(.01664)
        ig = ig_["INV00306"]
        pbm.A[ig,iv] += Float64(43.0)
        iv,ix_,_ = s2mpj_ii("X00306",ix_)
        arrset(pb.xnames,iv,"X00306")
        ig = ig_["CAP02105"]
        pbm.A[ig,iv] += Float64(.02041)
        ig = ig_["CAP02405"]
        pbm.A[ig,iv] += Float64(.01363)
        iv,ix_,_ = s2mpj_ii("X00306",ix_)
        arrset(pb.xnames,iv,"X00306")
        ig = ig_["CAP02004"]
        pbm.A[ig,iv] += Float64(.00026)
        ig = ig_["CAP02906"]
        pbm.A[ig,iv] += Float64(.02312)
        iv,ix_,_ = s2mpj_ii("X00306",ix_)
        arrset(pb.xnames,iv,"X00306")
        ig = ig_["CAP02606"]
        pbm.A[ig,iv] += Float64(.04408)
        ig = ig_["CAP02505"]
        pbm.A[ig,iv] += Float64(.01542)
        iv,ix_,_ = s2mpj_ii("X00306",ix_)
        arrset(pb.xnames,iv,"X00306")
        ig = ig_["CAP02005"]
        pbm.A[ig,iv] += Float64(.01083)
        ig = ig_["CAP02506"]
        pbm.A[ig,iv] += Float64(.01472)
        iv,ix_,_ = s2mpj_ii("X00306",ix_)
        arrset(pb.xnames,iv,"X00306")
        ig = ig_["CAP02406"]
        pbm.A[ig,iv] += Float64(.02201)
        ig = ig_["CAP02205"]
        pbm.A[ig,iv] += Float64(.12425)
        iv,ix_,_ = s2mpj_ii("X00306",ix_)
        arrset(pb.xnames,iv,"X00306")
        ig = ig_["CAP02706"]
        pbm.A[ig,iv] += Float64(.13616)
        ig = ig_["CAP02504"]
        pbm.A[ig,iv] += Float64(.00145)
        iv,ix_,_ = s2mpj_ii("X00306",ix_)
        arrset(pb.xnames,iv,"X00306")
        ig = ig_["CAP01905"]
        pbm.A[ig,iv] += Float64(.01255)
        ig = ig_["CAP01906"]
        pbm.A[ig,iv] += Float64(.03234)
        iv,ix_,_ = s2mpj_ii("X00306",ix_)
        arrset(pb.xnames,iv,"X00306")
        ig = ig_["CAP02605"]
        pbm.A[ig,iv] += Float64(.01608)
        ig = ig_["CAP03205"]
        pbm.A[ig,iv] += Float64(.01831)
        iv,ix_,_ = s2mpj_ii("X00403",ix_)
        arrset(pb.xnames,iv,"X00403")
        ig = ig_["CAP04703"]
        pbm.A[ig,iv] += Float64(.01583)
        ig = ig_["CAP04403"]
        pbm.A[ig,iv] += Float64(.01286)
        iv,ix_,_ = s2mpj_ii("X00403",ix_)
        arrset(pb.xnames,iv,"X00403")
        ig = ig_["CAP04603"]
        pbm.A[ig,iv] += Float64(.00104)
        ig = ig_["CAP04503"]
        pbm.A[ig,iv] += Float64(.01688)
        iv,ix_,_ = s2mpj_ii("X00403",ix_)
        arrset(pb.xnames,iv,"X00403")
        ig = ig_["CAP04402"]
        pbm.A[ig,iv] += Float64(.00521)
        ig = ig_["CAP03301"]
        pbm.A[ig,iv] += Float64(.00006)
        iv,ix_,_ = s2mpj_ii("X00403",ix_)
        arrset(pb.xnames,iv,"X00403")
        ig = ig_["CAP04102"]
        pbm.A[ig,iv] += Float64(.0117)
        ig = ig_["CAP03702"]
        pbm.A[ig,iv] += Float64(.03304)
        iv,ix_,_ = s2mpj_ii("X00403",ix_)
        arrset(pb.xnames,iv,"X00403")
        ig = ig_["CAP03402"]
        pbm.A[ig,iv] += Float64(.02558)
        ig = ig_["CAP03802"]
        pbm.A[ig,iv] += Float64(.02733)
        iv,ix_,_ = s2mpj_ii("X00403",ix_)
        arrset(pb.xnames,iv,"X00403")
        ig = ig_["CAP04002"]
        pbm.A[ig,iv] += Float64(.01621)
        ig = ig_["CAP05403"]
        pbm.A[ig,iv] += Float64(.01089)
        iv,ix_,_ = s2mpj_ii("X00403",ix_)
        arrset(pb.xnames,iv,"X00403")
        ig = ig_["CAP05303"]
        pbm.A[ig,iv] += Float64(.00591)
        ig = ig_["CAP04903"]
        pbm.A[ig,iv] += Float64(.00431)
        iv,ix_,_ = s2mpj_ii("X00403",ix_)
        arrset(pb.xnames,iv,"X00403")
        ig = ig_["CAP04202"]
        pbm.A[ig,iv] += Float64(.01126)
        ig = ig_["CAP05103"]
        pbm.A[ig,iv] += Float64(.00439)
        iv,ix_,_ = s2mpj_ii("X00403",ix_)
        arrset(pb.xnames,iv,"X00403")
        ig = ig_["CAP05003"]
        pbm.A[ig,iv] += Float64(.00685)
        ig = ig_["CAP04302"]
        pbm.A[ig,iv] += Float64(.02743)
        iv,ix_,_ = s2mpj_ii("X00403",ix_)
        arrset(pb.xnames,iv,"X00403")
        ig = ig_["CAP03803"]
        pbm.A[ig,iv] += Float64(.01706)
        ig = ig_["CAP04103"]
        pbm.A[ig,iv] += Float64(.0005)
        iv,ix_,_ = s2mpj_ii("X00403",ix_)
        arrset(pb.xnames,iv,"X00403")
        ig = ig_["CAP03703"]
        pbm.A[ig,iv] += Float64(.00362)
        ig = ig_["CAP03302"]
        pbm.A[ig,iv] += Float64(.01548)
        iv,ix_,_ = s2mpj_ii("X00403",ix_)
        arrset(pb.xnames,iv,"X00403")
        ig = ig_["CAP04502"]
        pbm.A[ig,iv] += Float64(.02252)
        ig = ig_["CAP04203"]
        pbm.A[ig,iv] += Float64(.02393)
        iv,ix_,_ = s2mpj_ii("X00403",ix_)
        arrset(pb.xnames,iv,"X00403")
        ig = ig_["CAP04003"]
        pbm.A[ig,iv] += Float64(.02556)
        ig = ig_["CAP03403"]
        pbm.A[ig,iv] += Float64(.03192)
        iv,ix_,_ = s2mpj_ii("X00403",ix_)
        arrset(pb.xnames,iv,"X00403")
        ig = ig_["CAP04602"]
        pbm.A[ig,iv] += Float64(.01566)
        ig = ig_["CAP03303"]
        pbm.A[ig,iv] += Float64(.00255)
        iv,ix_,_ = s2mpj_ii("X00403",ix_)
        arrset(pb.xnames,iv,"X00403")
        ig = ig_["CAP04702"]
        pbm.A[ig,iv] += Float64(.01004)
        ig = ig_["CAP04303"]
        pbm.A[ig,iv] += Float64(.00576)
        iv,ix_,_ = s2mpj_ii("X00403",ix_)
        arrset(pb.xnames,iv,"X00403")
        ig = ig_["CAP05402"]
        pbm.A[ig,iv] += Float64(.00786)
        ig = ig_["CAP05302"]
        pbm.A[ig,iv] += Float64(.01148)
        iv,ix_,_ = s2mpj_ii("X00403",ix_)
        arrset(pb.xnames,iv,"X00403")
        ig = ig_["CAP04601"]
        pbm.A[ig,iv] += Float64(.00027)
        ig = ig_["CAP05102"]
        pbm.A[ig,iv] += Float64(.0132)
        iv,ix_,_ = s2mpj_ii("X00403",ix_)
        arrset(pb.xnames,iv,"X00403")
        ig = ig_["CAP05002"]
        pbm.A[ig,iv] += Float64(.01633)
        ig = ig_["INV00303"]
        pbm.A[ig,iv] += Float64(153.0)
        iv,ix_,_ = s2mpj_ii("X00403",ix_)
        arrset(pb.xnames,iv,"X00403")
        ig = ig_["CAP04902"]
        pbm.A[ig,iv] += Float64(.01393)
        ig = ig_["CAP04101"]
        pbm.A[ig,iv] += Float64(.00047)
        iv,ix_,_ = s2mpj_ii("X00404",ix_)
        arrset(pb.xnames,iv,"X00404")
        ig = ig_["CAP05004"]
        pbm.A[ig,iv] += Float64(.00685)
        ig = ig_["CAP05103"]
        pbm.A[ig,iv] += Float64(.0132)
        iv,ix_,_ = s2mpj_ii("X00404",ix_)
        arrset(pb.xnames,iv,"X00404")
        ig = ig_["CAP03803"]
        pbm.A[ig,iv] += Float64(.02733)
        ig = ig_["CAP05104"]
        pbm.A[ig,iv] += Float64(.00439)
        iv,ix_,_ = s2mpj_ii("X00404",ix_)
        arrset(pb.xnames,iv,"X00404")
        ig = ig_["CAP04904"]
        pbm.A[ig,iv] += Float64(.00431)
        ig = ig_["CAP04704"]
        pbm.A[ig,iv] += Float64(.01583)
        iv,ix_,_ = s2mpj_ii("X00404",ix_)
        arrset(pb.xnames,iv,"X00404")
        ig = ig_["CAP05003"]
        pbm.A[ig,iv] += Float64(.01633)
        ig = ig_["CAP03403"]
        pbm.A[ig,iv] += Float64(.02558)
        iv,ix_,_ = s2mpj_ii("X00404",ix_)
        arrset(pb.xnames,iv,"X00404")
        ig = ig_["CAP05304"]
        pbm.A[ig,iv] += Float64(.00591)
        ig = ig_["CAP05404"]
        pbm.A[ig,iv] += Float64(.01089)
        iv,ix_,_ = s2mpj_ii("X00404",ix_)
        arrset(pb.xnames,iv,"X00404")
        ig = ig_["CAP05303"]
        pbm.A[ig,iv] += Float64(.01148)
        ig = ig_["CAP03703"]
        pbm.A[ig,iv] += Float64(.03304)
        iv,ix_,_ = s2mpj_ii("X00404",ix_)
        arrset(pb.xnames,iv,"X00404")
        ig = ig_["CAP05403"]
        pbm.A[ig,iv] += Float64(.00786)
        ig = ig_["CAP03303"]
        pbm.A[ig,iv] += Float64(.01554)
        iv,ix_,_ = s2mpj_ii("X00404",ix_)
        arrset(pb.xnames,iv,"X00404")
        ig = ig_["CAP03804"]
        pbm.A[ig,iv] += Float64(.01706)
        ig = ig_["CAP04203"]
        pbm.A[ig,iv] += Float64(.01126)
        iv,ix_,_ = s2mpj_ii("X00404",ix_)
        arrset(pb.xnames,iv,"X00404")
        ig = ig_["CAP04303"]
        pbm.A[ig,iv] += Float64(.02743)
        ig = ig_["CAP04404"]
        pbm.A[ig,iv] += Float64(.01286)
        iv,ix_,_ = s2mpj_ii("X00404",ix_)
        arrset(pb.xnames,iv,"X00404")
        ig = ig_["INV00304"]
        pbm.A[ig,iv] += Float64(153.0)
        ig = ig_["CAP03404"]
        pbm.A[ig,iv] += Float64(.03192)
        iv,ix_,_ = s2mpj_ii("X00404",ix_)
        arrset(pb.xnames,iv,"X00404")
        ig = ig_["CAP04304"]
        pbm.A[ig,iv] += Float64(.00576)
        ig = ig_["CAP04603"]
        pbm.A[ig,iv] += Float64(.01593)
        iv,ix_,_ = s2mpj_ii("X00404",ix_)
        arrset(pb.xnames,iv,"X00404")
        ig = ig_["CAP04403"]
        pbm.A[ig,iv] += Float64(.00521)
        ig = ig_["CAP04204"]
        pbm.A[ig,iv] += Float64(.02393)
        iv,ix_,_ = s2mpj_ii("X00404",ix_)
        arrset(pb.xnames,iv,"X00404")
        ig = ig_["CAP04003"]
        pbm.A[ig,iv] += Float64(.01621)
        ig = ig_["CAP03304"]
        pbm.A[ig,iv] += Float64(.00255)
        iv,ix_,_ = s2mpj_ii("X00404",ix_)
        arrset(pb.xnames,iv,"X00404")
        ig = ig_["CAP04104"]
        pbm.A[ig,iv] += Float64(.0005)
        ig = ig_["CAP04503"]
        pbm.A[ig,iv] += Float64(.02252)
        iv,ix_,_ = s2mpj_ii("X00404",ix_)
        arrset(pb.xnames,iv,"X00404")
        ig = ig_["CAP04703"]
        pbm.A[ig,iv] += Float64(.01004)
        ig = ig_["CAP04504"]
        pbm.A[ig,iv] += Float64(.01688)
        iv,ix_,_ = s2mpj_ii("X00404",ix_)
        arrset(pb.xnames,iv,"X00404")
        ig = ig_["CAP04103"]
        pbm.A[ig,iv] += Float64(.01217)
        ig = ig_["CAP03704"]
        pbm.A[ig,iv] += Float64(.00362)
        iv,ix_,_ = s2mpj_ii("X00404",ix_)
        arrset(pb.xnames,iv,"X00404")
        ig = ig_["CAP04604"]
        pbm.A[ig,iv] += Float64(.00104)
        ig = ig_["CAP04903"]
        pbm.A[ig,iv] += Float64(.01393)
        iv,ix_,_ = s2mpj_ii("X00404",ix_)
        arrset(pb.xnames,iv,"X00404")
        ig = ig_["CAP04004"]
        pbm.A[ig,iv] += Float64(.02556)
        iv,ix_,_ = s2mpj_ii("X00405",ix_)
        arrset(pb.xnames,iv,"X00405")
        ig = ig_["CAP03305"]
        pbm.A[ig,iv] += Float64(.00313)
        ig = ig_["CAP04505"]
        pbm.A[ig,iv] += Float64(.01772)
        iv,ix_,_ = s2mpj_ii("X00405",ix_)
        arrset(pb.xnames,iv,"X00405")
        ig = ig_["CAP03304"]
        pbm.A[ig,iv] += Float64(.01496)
        ig = ig_["CAP04004"]
        pbm.A[ig,iv] += Float64(.01561)
        iv,ix_,_ = s2mpj_ii("X00405",ix_)
        arrset(pb.xnames,iv,"X00405")
        ig = ig_["CAP03804"]
        pbm.A[ig,iv] += Float64(.02632)
        ig = ig_["CAP04305"]
        pbm.A[ig,iv] += Float64(.00678)
        iv,ix_,_ = s2mpj_ii("X00405",ix_)
        arrset(pb.xnames,iv,"X00405")
        ig = ig_["CAP04405"]
        pbm.A[ig,iv] += Float64(.01306)
        ig = ig_["CAP04105"]
        pbm.A[ig,iv] += Float64(.00095)
        iv,ix_,_ = s2mpj_ii("X00405",ix_)
        arrset(pb.xnames,iv,"X00405")
        ig = ig_["CAP03405"]
        pbm.A[ig,iv] += Float64(.03287)
        ig = ig_["CAP03705"]
        pbm.A[ig,iv] += Float64(.00485)
        iv,ix_,_ = s2mpj_ii("X00405",ix_)
        arrset(pb.xnames,iv,"X00405")
        ig = ig_["CAP03805"]
        pbm.A[ig,iv] += Float64(.01807)
        ig = ig_["CAP03404"]
        pbm.A[ig,iv] += Float64(.02464)
        iv,ix_,_ = s2mpj_ii("X00405",ix_)
        arrset(pb.xnames,iv,"X00405")
        ig = ig_["CAP04005"]
        pbm.A[ig,iv] += Float64(.02616)
        ig = ig_["CAP04205"]
        pbm.A[ig,iv] += Float64(.02435)
        iv,ix_,_ = s2mpj_ii("X00405",ix_)
        arrset(pb.xnames,iv,"X00405")
        ig = ig_["CAP05105"]
        pbm.A[ig,iv] += Float64(.00488)
        ig = ig_["CAP04604"]
        pbm.A[ig,iv] += Float64(.01534)
        iv,ix_,_ = s2mpj_ii("X00405",ix_)
        arrset(pb.xnames,iv,"X00405")
        ig = ig_["CAP04704"]
        pbm.A[ig,iv] += Float64(.00967)
        ig = ig_["CAP04904"]
        pbm.A[ig,iv] += Float64(.01341)
        iv,ix_,_ = s2mpj_ii("X00405",ix_)
        arrset(pb.xnames,iv,"X00405")
        ig = ig_["CAP05004"]
        pbm.A[ig,iv] += Float64(.01573)
        ig = ig_["CAP04504"]
        pbm.A[ig,iv] += Float64(.02169)
        iv,ix_,_ = s2mpj_ii("X00405",ix_)
        arrset(pb.xnames,iv,"X00405")
        ig = ig_["INV00305"]
        pbm.A[ig,iv] += Float64(153.0)
        ig = ig_["CAP05104"]
        pbm.A[ig,iv] += Float64(.01271)
        iv,ix_,_ = s2mpj_ii("X00405",ix_)
        arrset(pb.xnames,iv,"X00405")
        ig = ig_["CAP04404"]
        pbm.A[ig,iv] += Float64(.00502)
        ig = ig_["CAP05304"]
        pbm.A[ig,iv] += Float64(.01106)
        iv,ix_,_ = s2mpj_ii("X00405",ix_)
        arrset(pb.xnames,iv,"X00405")
        ig = ig_["CAP05404"]
        pbm.A[ig,iv] += Float64(.00757)
        ig = ig_["CAP04605"]
        pbm.A[ig,iv] += Float64(.00163)
        iv,ix_,_ = s2mpj_ii("X00405",ix_)
        arrset(pb.xnames,iv,"X00405")
        ig = ig_["CAP04304"]
        pbm.A[ig,iv] += Float64(.02642)
        ig = ig_["CAP05405"]
        pbm.A[ig,iv] += Float64(.01118)
        iv,ix_,_ = s2mpj_ii("X00405",ix_)
        arrset(pb.xnames,iv,"X00405")
        ig = ig_["CAP05305"]
        pbm.A[ig,iv] += Float64(.00634)
        ig = ig_["CAP04204"]
        pbm.A[ig,iv] += Float64(.01084)
        iv,ix_,_ = s2mpj_ii("X00405",ix_)
        arrset(pb.xnames,iv,"X00405")
        ig = ig_["CAP05005"]
        pbm.A[ig,iv] += Float64(.00745)
        ig = ig_["CAP04905"]
        pbm.A[ig,iv] += Float64(.00483)
        iv,ix_,_ = s2mpj_ii("X00405",ix_)
        arrset(pb.xnames,iv,"X00405")
        ig = ig_["CAP04104"]
        pbm.A[ig,iv] += Float64(.01172)
        ig = ig_["CAP04705"]
        pbm.A[ig,iv] += Float64(.0162)
        iv,ix_,_ = s2mpj_ii("X00405",ix_)
        arrset(pb.xnames,iv,"X00405")
        ig = ig_["CAP03704"]
        pbm.A[ig,iv] += Float64(.03181)
        iv,ix_,_ = s2mpj_ii("X00406",ix_)
        arrset(pb.xnames,iv,"X00406")
        ig = ig_["CAP04405"]
        pbm.A[ig,iv] += Float64(.00542)
        ig = ig_["CAP05305"]
        pbm.A[ig,iv] += Float64(.01194)
        iv,ix_,_ = s2mpj_ii("X00406",ix_)
        arrset(pb.xnames,iv,"X00406")
        ig = ig_["CAP05105"]
        pbm.A[ig,iv] += Float64(.01373)
        ig = ig_["CAP05005"]
        pbm.A[ig,iv] += Float64(.01699)
        iv,ix_,_ = s2mpj_ii("X00406",ix_)
        arrset(pb.xnames,iv,"X00406")
        ig = ig_["CAP04905"]
        pbm.A[ig,iv] += Float64(.01448)
        ig = ig_["CAP04705"]
        pbm.A[ig,iv] += Float64(.01044)
        iv,ix_,_ = s2mpj_ii("X00406",ix_)
        arrset(pb.xnames,iv,"X00406")
        ig = ig_["CAP04605"]
        pbm.A[ig,iv] += Float64(.01656)
        ig = ig_["CAP04505"]
        pbm.A[ig,iv] += Float64(.02342)
        iv,ix_,_ = s2mpj_ii("X00406",ix_)
        arrset(pb.xnames,iv,"X00406")
        ig = ig_["CAP04305"]
        pbm.A[ig,iv] += Float64(.02853)
        ig = ig_["CAP05405"]
        pbm.A[ig,iv] += Float64(.00817)
        iv,ix_,_ = s2mpj_ii("X00406",ix_)
        arrset(pb.xnames,iv,"X00406")
        ig = ig_["CAP04205"]
        pbm.A[ig,iv] += Float64(.01171)
        ig = ig_["CAP04105"]
        pbm.A[ig,iv] += Float64(.01266)
        iv,ix_,_ = s2mpj_ii("X00406",ix_)
        arrset(pb.xnames,iv,"X00406")
        ig = ig_["CAP04005"]
        pbm.A[ig,iv] += Float64(.01686)
        ig = ig_["CAP04906"]
        pbm.A[ig,iv] += Float64(.01824)
        iv,ix_,_ = s2mpj_ii("X00406",ix_)
        arrset(pb.xnames,iv,"X00406")
        ig = ig_["CAP03705"]
        pbm.A[ig,iv] += Float64(.03436)
        ig = ig_["CAP03405"]
        pbm.A[ig,iv] += Float64(.02661)
        iv,ix_,_ = s2mpj_ii("X00406",ix_)
        arrset(pb.xnames,iv,"X00406")
        ig = ig_["CAP03306"]
        pbm.A[ig,iv] += Float64(.01809)
        ig = ig_["CAP04406"]
        pbm.A[ig,iv] += Float64(.01808)
        iv,ix_,_ = s2mpj_ii("X00406",ix_)
        arrset(pb.xnames,iv,"X00406")
        ig = ig_["CAP04706"]
        pbm.A[ig,iv] += Float64(.02586)
        ig = ig_["CAP04606"]
        pbm.A[ig,iv] += Float64(.01696)
        iv,ix_,_ = s2mpj_ii("X00406",ix_)
        arrset(pb.xnames,iv,"X00406")
        ig = ig_["CAP04506"]
        pbm.A[ig,iv] += Float64(.0394)
        ig = ig_["CAP05006"]
        pbm.A[ig,iv] += Float64(.02318)
        iv,ix_,_ = s2mpj_ii("X00406",ix_)
        arrset(pb.xnames,iv,"X00406")
        ig = ig_["CAP05106"]
        pbm.A[ig,iv] += Float64(.01759)
        ig = ig_["CAP05306"]
        pbm.A[ig,iv] += Float64(.0174)
        iv,ix_,_ = s2mpj_ii("X00406",ix_)
        arrset(pb.xnames,iv,"X00406")
        ig = ig_["CAP05406"]
        pbm.A[ig,iv] += Float64(.01875)
        ig = ig_["CAP04306"]
        pbm.A[ig,iv] += Float64(.03319)
        iv,ix_,_ = s2mpj_ii("X00406",ix_)
        arrset(pb.xnames,iv,"X00406")
        ig = ig_["CAP03305"]
        pbm.A[ig,iv] += Float64(.01616)
        ig = ig_["CAP04206"]
        pbm.A[ig,iv] += Float64(.03519)
        iv,ix_,_ = s2mpj_ii("X00406",ix_)
        arrset(pb.xnames,iv,"X00406")
        ig = ig_["CAP04106"]
        pbm.A[ig,iv] += Float64(.01267)
        ig = ig_["CAP04006"]
        pbm.A[ig,iv] += Float64(.04177)
        iv,ix_,_ = s2mpj_ii("X00406",ix_)
        arrset(pb.xnames,iv,"X00406")
        ig = ig_["CAP03806"]
        pbm.A[ig,iv] += Float64(.04439)
        ig = ig_["CAP03706"]
        pbm.A[ig,iv] += Float64(.03666)
        iv,ix_,_ = s2mpj_ii("X00406",ix_)
        arrset(pb.xnames,iv,"X00406")
        ig = ig_["CAP03406"]
        pbm.A[ig,iv] += Float64(.05751)
        ig = ig_["CAP03805"]
        pbm.A[ig,iv] += Float64(.02843)
        iv,ix_,_ = s2mpj_ii("X00503",ix_)
        arrset(pb.xnames,iv,"X00503")
        ig = ig_["CAP03002"]
        pbm.A[ig,iv] += Float64(.00998)
        ig = ig_["CAP02802"]
        pbm.A[ig,iv] += Float64(.01664)
        iv,ix_,_ = s2mpj_ii("X00503",ix_)
        arrset(pb.xnames,iv,"X00503")
        ig = ig_["CAP02803"]
        pbm.A[ig,iv] += Float64(.00392)
        ig = ig_["CAP02703"]
        pbm.A[ig,iv] += Float64(.02005)
        iv,ix_,_ = s2mpj_ii("X00503",ix_)
        arrset(pb.xnames,iv,"X00503")
        ig = ig_["CAP02701"]
        pbm.A[ig,iv] += Float64(.00576)
        ig = ig_["CAP03102"]
        pbm.A[ig,iv] += Float64(.00594)
        iv,ix_,_ = s2mpj_ii("X00503",ix_)
        arrset(pb.xnames,iv,"X00503")
        ig = ig_["CAP02702"]
        pbm.A[ig,iv] += Float64(.11034)
        ig = ig_["CAP03202"]
        pbm.A[ig,iv] += Float64(.01831)
        iv,ix_,_ = s2mpj_ii("X00503",ix_)
        arrset(pb.xnames,iv,"X00503")
        ig = ig_["INV00403"]
        pbm.A[ig,iv] += Float64(43.0)
        ig = ig_["CAP02902"]
        pbm.A[ig,iv] += Float64(.01593)
        iv,ix_,_ = s2mpj_ii("X00503",ix_)
        arrset(pb.xnames,iv,"X00503")
        ig = ig_["CAP02603"]
        pbm.A[ig,iv] += Float64(.028)
        ig = ig_["CAP02903"]
        pbm.A[ig,iv] += Float64(.00719)
        iv,ix_,_ = s2mpj_ii("X00503",ix_)
        arrset(pb.xnames,iv,"X00503")
        ig = ig_["CAP03003"]
        pbm.A[ig,iv] += Float64(.01077)
        ig = ig_["CAP03203"]
        pbm.A[ig,iv] += Float64(.02077)
        iv,ix_,_ = s2mpj_ii("X00503",ix_)
        arrset(pb.xnames,iv,"X00503")
        ig = ig_["CAP02501"]
        pbm.A[ig,iv] += Float64(.00285)
        ig = ig_["CAP02201"]
        pbm.A[ig,iv] += Float64(.00077)
        iv,ix_,_ = s2mpj_ii("X00503",ix_)
        arrset(pb.xnames,iv,"X00503")
        ig = ig_["CAP02001"]
        pbm.A[ig,iv] += Float64(.00079)
        ig = ig_["CAP01903"]
        pbm.A[ig,iv] += Float64(.01979)
        iv,ix_,_ = s2mpj_ii("X00503",ix_)
        arrset(pb.xnames,iv,"X00503")
        ig = ig_["CAP03103"]
        pbm.A[ig,iv] += Float64(.01198)
        ig = ig_["CAP02203"]
        pbm.A[ig,iv] += Float64(.02513)
        iv,ix_,_ = s2mpj_ii("X00503",ix_)
        arrset(pb.xnames,iv,"X00503")
        ig = ig_["CAP02003"]
        pbm.A[ig,iv] += Float64(.00054)
        ig = ig_["CAP02202"]
        pbm.A[ig,iv] += Float64(.12349)
        iv,ix_,_ = s2mpj_ii("X00503",ix_)
        arrset(pb.xnames,iv,"X00503")
        ig = ig_["CAP01902"]
        pbm.A[ig,iv] += Float64(.01255)
        ig = ig_["CAP02402"]
        pbm.A[ig,iv] += Float64(.01363)
        iv,ix_,_ = s2mpj_ii("X00503",ix_)
        arrset(pb.xnames,iv,"X00503")
        ig = ig_["CAP02102"]
        pbm.A[ig,iv] += Float64(.02041)
        ig = ig_["CAP02502"]
        pbm.A[ig,iv] += Float64(.01187)
        iv,ix_,_ = s2mpj_ii("X00503",ix_)
        arrset(pb.xnames,iv,"X00503")
        ig = ig_["CAP02103"]
        pbm.A[ig,iv] += Float64(.02973)
        ig = ig_["CAP02002"]
        pbm.A[ig,iv] += Float64(.00977)
        iv,ix_,_ = s2mpj_ii("X00503",ix_)
        arrset(pb.xnames,iv,"X00503")
        ig = ig_["CAP02602"]
        pbm.A[ig,iv] += Float64(.01608)
        ig = ig_["CAP02403"]
        pbm.A[ig,iv] += Float64(.00837)
        iv,ix_,_ = s2mpj_ii("X00504",ix_)
        arrset(pb.xnames,iv,"X00504")
        ig = ig_["CAP03103"]
        pbm.A[ig,iv] += Float64(.00567)
        ig = ig_["CAP03204"]
        pbm.A[ig,iv] += Float64(.0216)
        iv,ix_,_ = s2mpj_ii("X00504",ix_)
        arrset(pb.xnames,iv,"X00504")
        ig = ig_["CAP01903"]
        pbm.A[ig,iv] += Float64(.01198)
        ig = ig_["CAP02003"]
        pbm.A[ig,iv] += Float64(.00983)
        iv,ix_,_ = s2mpj_ii("X00504",ix_)
        arrset(pb.xnames,iv,"X00504")
        ig = ig_["CAP02503"]
        pbm.A[ig,iv] += Float64(.01266)
        ig = ig_["CAP03203"]
        pbm.A[ig,iv] += Float64(.01748)
        iv,ix_,_ = s2mpj_ii("X00504",ix_)
        arrset(pb.xnames,iv,"X00504")
        ig = ig_["CAP02903"]
        pbm.A[ig,iv] += Float64(.0152)
        ig = ig_["CAP02103"]
        pbm.A[ig,iv] += Float64(.01948)
        iv,ix_,_ = s2mpj_ii("X00504",ix_)
        arrset(pb.xnames,iv,"X00504")
        ig = ig_["CAP02803"]
        pbm.A[ig,iv] += Float64(.01589)
        ig = ig_["CAP02703"]
        pbm.A[ig,iv] += Float64(.1083)
        iv,ix_,_ = s2mpj_ii("X00504",ix_)
        arrset(pb.xnames,iv,"X00504")
        ig = ig_["CAP02203"]
        pbm.A[ig,iv] += Float64(.11861)
        ig = ig_["CAP02403"]
        pbm.A[ig,iv] += Float64(.01301)
        iv,ix_,_ = s2mpj_ii("X00504",ix_)
        arrset(pb.xnames,iv,"X00504")
        ig = ig_["CAP02603"]
        pbm.A[ig,iv] += Float64(.01535)
        ig = ig_["CAP03003"]
        pbm.A[ig,iv] += Float64(.00953)
        iv,ix_,_ = s2mpj_ii("X00504",ix_)
        arrset(pb.xnames,iv,"X00504")
        ig = ig_["CAP03004"]
        pbm.A[ig,iv] += Float64(.01122)
        ig = ig_["CAP02804"]
        pbm.A[ig,iv] += Float64(.00467)
        iv,ix_,_ = s2mpj_ii("X00504",ix_)
        arrset(pb.xnames,iv,"X00504")
        ig = ig_["INV00404"]
        pbm.A[ig,iv] += Float64(43.0)
        ig = ig_["CAP02702"]
        pbm.A[ig,iv] += Float64(.00376)
        iv,ix_,_ = s2mpj_ii("X00504",ix_)
        arrset(pb.xnames,iv,"X00504")
        ig = ig_["CAP02704"]
        pbm.A[ig,iv] += Float64(.02409)
        ig = ig_["CAP02904"]
        pbm.A[ig,iv] += Float64(.00791)
        iv,ix_,_ = s2mpj_ii("X00504",ix_)
        arrset(pb.xnames,iv,"X00504")
        ig = ig_["CAP02604"]
        pbm.A[ig,iv] += Float64(.02873)
        ig = ig_["CAP02104"]
        pbm.A[ig,iv] += Float64(.03066)
        iv,ix_,_ = s2mpj_ii("X00504",ix_)
        arrset(pb.xnames,iv,"X00504")
        ig = ig_["CAP02502"]
        pbm.A[ig,iv] += Float64(.00205)
        ig = ig_["CAP03104"]
        pbm.A[ig,iv] += Float64(.01225)
        iv,ix_,_ = s2mpj_ii("X00504",ix_)
        arrset(pb.xnames,iv,"X00504")
        ig = ig_["CAP02004"]
        pbm.A[ig,iv] += Float64(.00077)
        ig = ig_["CAP02002"]
        pbm.A[ig,iv] += Float64(.0005)
        iv,ix_,_ = s2mpj_ii("X00504",ix_)
        arrset(pb.xnames,iv,"X00504")
        ig = ig_["CAP01904"]
        pbm.A[ig,iv] += Float64(.02036)
        ig = ig_["CAP02204"]
        pbm.A[ig,iv] += Float64(.03078)
        iv,ix_,_ = s2mpj_ii("X00504",ix_)
        arrset(pb.xnames,iv,"X00504")
        ig = ig_["CAP02404"]
        pbm.A[ig,iv] += Float64(.00899)
        iv,ix_,_ = s2mpj_ii("X00505",ix_)
        arrset(pb.xnames,iv,"X00505")
        ig = ig_["CAP02703"]
        pbm.A[ig,iv] += Float64(.00253)
        ig = ig_["CAP01905"]
        pbm.A[ig,iv] += Float64(.02036)
        iv,ix_,_ = s2mpj_ii("X00505",ix_)
        arrset(pb.xnames,iv,"X00505")
        ig = ig_["CAP02405"]
        pbm.A[ig,iv] += Float64(.00899)
        ig = ig_["INV00405"]
        pbm.A[ig,iv] += Float64(43.0)
        iv,ix_,_ = s2mpj_ii("X00505",ix_)
        arrset(pb.xnames,iv,"X00505")
        ig = ig_["CAP02805"]
        pbm.A[ig,iv] += Float64(.00467)
        ig = ig_["CAP02704"]
        pbm.A[ig,iv] += Float64(.10954)
        iv,ix_,_ = s2mpj_ii("X00505",ix_)
        arrset(pb.xnames,iv,"X00505")
        ig = ig_["CAP02003"]
        pbm.A[ig,iv] += Float64(.00025)
        ig = ig_["CAP02404"]
        pbm.A[ig,iv] += Float64(.01301)
        iv,ix_,_ = s2mpj_ii("X00505",ix_)
        arrset(pb.xnames,iv,"X00505")
        ig = ig_["CAP02905"]
        pbm.A[ig,iv] += Float64(.00791)
        ig = ig_["CAP02104"]
        pbm.A[ig,iv] += Float64(.01948)
        iv,ix_,_ = s2mpj_ii("X00505",ix_)
        arrset(pb.xnames,iv,"X00505")
        ig = ig_["CAP02503"]
        pbm.A[ig,iv] += Float64(.00138)
        ig = ig_["CAP02705"]
        pbm.A[ig,iv] += Float64(.02409)
        iv,ix_,_ = s2mpj_ii("X00505",ix_)
        arrset(pb.xnames,iv,"X00505")
        ig = ig_["CAP02605"]
        pbm.A[ig,iv] += Float64(.02873)
        ig = ig_["CAP02204"]
        pbm.A[ig,iv] += Float64(.11861)
        iv,ix_,_ = s2mpj_ii("X00505",ix_)
        arrset(pb.xnames,iv,"X00505")
        ig = ig_["CAP02004"]
        pbm.A[ig,iv] += Float64(.01009)
        ig = ig_["CAP02604"]
        pbm.A[ig,iv] += Float64(.01535)
        iv,ix_,_ = s2mpj_ii("X00505",ix_)
        arrset(pb.xnames,iv,"X00505")
        ig = ig_["CAP03004"]
        pbm.A[ig,iv] += Float64(.00953)
        ig = ig_["CAP02504"]
        pbm.A[ig,iv] += Float64(.01333)
        iv,ix_,_ = s2mpj_ii("X00505",ix_)
        arrset(pb.xnames,iv,"X00505")
        ig = ig_["CAP02904"]
        pbm.A[ig,iv] += Float64(.0152)
        ig = ig_["CAP02205"]
        pbm.A[ig,iv] += Float64(.03078)
        iv,ix_,_ = s2mpj_ii("X00505",ix_)
        arrset(pb.xnames,iv,"X00505")
        ig = ig_["CAP01904"]
        pbm.A[ig,iv] += Float64(.01198)
        ig = ig_["CAP02005"]
        pbm.A[ig,iv] += Float64(.00077)
        iv,ix_,_ = s2mpj_ii("X00505",ix_)
        arrset(pb.xnames,iv,"X00505")
        ig = ig_["CAP02804"]
        pbm.A[ig,iv] += Float64(.01589)
        ig = ig_["CAP03204"]
        pbm.A[ig,iv] += Float64(.01748)
        iv,ix_,_ = s2mpj_ii("X00505",ix_)
        arrset(pb.xnames,iv,"X00505")
        ig = ig_["CAP02105"]
        pbm.A[ig,iv] += Float64(.03066)
        ig = ig_["CAP03205"]
        pbm.A[ig,iv] += Float64(.0216)
        iv,ix_,_ = s2mpj_ii("X00505",ix_)
        arrset(pb.xnames,iv,"X00505")
        ig = ig_["CAP03105"]
        pbm.A[ig,iv] += Float64(.01225)
        ig = ig_["CAP03104"]
        pbm.A[ig,iv] += Float64(.00567)
        iv,ix_,_ = s2mpj_ii("X00505",ix_)
        arrset(pb.xnames,iv,"X00505")
        ig = ig_["CAP03005"]
        pbm.A[ig,iv] += Float64(.01122)
        iv,ix_,_ = s2mpj_ii("X00506",ix_)
        arrset(pb.xnames,iv,"X00506")
        ig = ig_["CAP02704"]
        pbm.A[ig,iv] += Float64(.00265)
        ig = ig_["INV00406"]
        pbm.A[ig,iv] += Float64(43.0)
        iv,ix_,_ = s2mpj_ii("X00506",ix_)
        arrset(pb.xnames,iv,"X00506")
        ig = ig_["CAP02905"]
        pbm.A[ig,iv] += Float64(.01593)
        ig = ig_["CAP02504"]
        pbm.A[ig,iv] += Float64(.00145)
        iv,ix_,_ = s2mpj_ii("X00506",ix_)
        arrset(pb.xnames,iv,"X00506")
        ig = ig_["CAP02705"]
        pbm.A[ig,iv] += Float64(.1174)
        ig = ig_["CAP03205"]
        pbm.A[ig,iv] += Float64(.01831)
        iv,ix_,_ = s2mpj_ii("X00506",ix_)
        arrset(pb.xnames,iv,"X00506")
        ig = ig_["CAP03105"]
        pbm.A[ig,iv] += Float64(.00594)
        ig = ig_["CAP02805"]
        pbm.A[ig,iv] += Float64(.01664)
        iv,ix_,_ = s2mpj_ii("X00506",ix_)
        arrset(pb.xnames,iv,"X00506")
        ig = ig_["CAP03005"]
        pbm.A[ig,iv] += Float64(.00998)
        ig = ig_["CAP02106"]
        pbm.A[ig,iv] += Float64(.05014)
        iv,ix_,_ = s2mpj_ii("X00506",ix_)
        arrset(pb.xnames,iv,"X00506")
        ig = ig_["CAP03006"]
        pbm.A[ig,iv] += Float64(.02075)
        ig = ig_["CAP02706"]
        pbm.A[ig,iv] += Float64(.13616)
        iv,ix_,_ = s2mpj_ii("X00506",ix_)
        arrset(pb.xnames,iv,"X00506")
        ig = ig_["CAP01905"]
        pbm.A[ig,iv] += Float64(.01255)
        ig = ig_["CAP02806"]
        pbm.A[ig,iv] += Float64(.02056)
        iv,ix_,_ = s2mpj_ii("X00506",ix_)
        arrset(pb.xnames,iv,"X00506")
        ig = ig_["CAP02906"]
        pbm.A[ig,iv] += Float64(.02312)
        ig = ig_["CAP02606"]
        pbm.A[ig,iv] += Float64(.04408)
        iv,ix_,_ = s2mpj_ii("X00506",ix_)
        arrset(pb.xnames,iv,"X00506")
        ig = ig_["CAP02004"]
        pbm.A[ig,iv] += Float64(.00026)
        ig = ig_["CAP03106"]
        pbm.A[ig,iv] += Float64(.01793)
        iv,ix_,_ = s2mpj_ii("X00506",ix_)
        arrset(pb.xnames,iv,"X00506")
        ig = ig_["CAP02005"]
        pbm.A[ig,iv] += Float64(.01083)
        ig = ig_["CAP02205"]
        pbm.A[ig,iv] += Float64(.12425)
        iv,ix_,_ = s2mpj_ii("X00506",ix_)
        arrset(pb.xnames,iv,"X00506")
        ig = ig_["CAP03206"]
        pbm.A[ig,iv] += Float64(.03908)
        ig = ig_["CAP01906"]
        pbm.A[ig,iv] += Float64(.03234)
        iv,ix_,_ = s2mpj_ii("X00506",ix_)
        arrset(pb.xnames,iv,"X00506")
        ig = ig_["CAP02506"]
        pbm.A[ig,iv] += Float64(.01472)
        ig = ig_["CAP02406"]
        pbm.A[ig,iv] += Float64(.02201)
        iv,ix_,_ = s2mpj_ii("X00506",ix_)
        arrset(pb.xnames,iv,"X00506")
        ig = ig_["CAP02405"]
        pbm.A[ig,iv] += Float64(.01363)
        ig = ig_["CAP02505"]
        pbm.A[ig,iv] += Float64(.01542)
        iv,ix_,_ = s2mpj_ii("X00506",ix_)
        arrset(pb.xnames,iv,"X00506")
        ig = ig_["CAP02105"]
        pbm.A[ig,iv] += Float64(.02041)
        ig = ig_["CAP02206"]
        pbm.A[ig,iv] += Float64(.14939)
        iv,ix_,_ = s2mpj_ii("X00506",ix_)
        arrset(pb.xnames,iv,"X00506")
        ig = ig_["CAP02605"]
        pbm.A[ig,iv] += Float64(.01608)
        ig = ig_["CAP02006"]
        pbm.A[ig,iv] += Float64(.0111)
        iv,ix_,_ = s2mpj_ii("X00603",ix_)
        arrset(pb.xnames,iv,"X00603")
        ig = ig_["CAP05303"]
        pbm.A[ig,iv] += Float64(.00591)
        ig = ig_["CAP05102"]
        pbm.A[ig,iv] += Float64(.0132)
        iv,ix_,_ = s2mpj_ii("X00603",ix_)
        arrset(pb.xnames,iv,"X00603")
        ig = ig_["CAP03303"]
        pbm.A[ig,iv] += Float64(.00255)
        ig = ig_["CAP03803"]
        pbm.A[ig,iv] += Float64(.01706)
        iv,ix_,_ = s2mpj_ii("X00603",ix_)
        arrset(pb.xnames,iv,"X00603")
        ig = ig_["CAP04601"]
        pbm.A[ig,iv] += Float64(.00027)
        ig = ig_["CAP05403"]
        pbm.A[ig,iv] += Float64(.01089)
        iv,ix_,_ = s2mpj_ii("X00603",ix_)
        arrset(pb.xnames,iv,"X00603")
        ig = ig_["CAP05302"]
        pbm.A[ig,iv] += Float64(.01148)
        ig = ig_["CAP05402"]
        pbm.A[ig,iv] += Float64(.00786)
        iv,ix_,_ = s2mpj_ii("X00603",ix_)
        arrset(pb.xnames,iv,"X00603")
        ig = ig_["CAP03403"]
        pbm.A[ig,iv] += Float64(.03192)
        ig = ig_["CAP03703"]
        pbm.A[ig,iv] += Float64(.00362)
        iv,ix_,_ = s2mpj_ii("X00603",ix_)
        arrset(pb.xnames,iv,"X00603")
        ig = ig_["CAP03301"]
        pbm.A[ig,iv] += Float64(.00006)
        ig = ig_["CAP04101"]
        pbm.A[ig,iv] += Float64(.00047)
        iv,ix_,_ = s2mpj_ii("X00603",ix_)
        arrset(pb.xnames,iv,"X00603")
        ig = ig_["CAP04503"]
        pbm.A[ig,iv] += Float64(.01688)
        ig = ig_["CAP04702"]
        pbm.A[ig,iv] += Float64(.01004)
        iv,ix_,_ = s2mpj_ii("X00603",ix_)
        arrset(pb.xnames,iv,"X00603")
        ig = ig_["CAP04202"]
        pbm.A[ig,iv] += Float64(.01126)
        ig = ig_["CAP04903"]
        pbm.A[ig,iv] += Float64(.00431)
        iv,ix_,_ = s2mpj_ii("X00603",ix_)
        arrset(pb.xnames,iv,"X00603")
        ig = ig_["CAP04003"]
        pbm.A[ig,iv] += Float64(.02556)
        ig = ig_["CAP04402"]
        pbm.A[ig,iv] += Float64(.00521)
        iv,ix_,_ = s2mpj_ii("X00603",ix_)
        arrset(pb.xnames,iv,"X00603")
        ig = ig_["CAP04203"]
        pbm.A[ig,iv] += Float64(.02393)
        ig = ig_["CAP04103"]
        pbm.A[ig,iv] += Float64(.0005)
        iv,ix_,_ = s2mpj_ii("X00603",ix_)
        arrset(pb.xnames,iv,"X00603")
        ig = ig_["CAP03302"]
        pbm.A[ig,iv] += Float64(.01548)
        ig = ig_["CAP04902"]
        pbm.A[ig,iv] += Float64(.01393)
        iv,ix_,_ = s2mpj_ii("X00603",ix_)
        arrset(pb.xnames,iv,"X00603")
        ig = ig_["CAP04502"]
        pbm.A[ig,iv] += Float64(.02252)
        ig = ig_["CAP05103"]
        pbm.A[ig,iv] += Float64(.00439)
        iv,ix_,_ = s2mpj_ii("X00603",ix_)
        arrset(pb.xnames,iv,"X00603")
        ig = ig_["CAP04403"]
        pbm.A[ig,iv] += Float64(.01286)
        ig = ig_["CAP05003"]
        pbm.A[ig,iv] += Float64(.00685)
        iv,ix_,_ = s2mpj_ii("X00603",ix_)
        arrset(pb.xnames,iv,"X00603")
        ig = ig_["CAP04602"]
        pbm.A[ig,iv] += Float64(.01566)
        ig = ig_["INV00403"]
        pbm.A[ig,iv] += Float64(143.0)
        iv,ix_,_ = s2mpj_ii("X00603",ix_)
        arrset(pb.xnames,iv,"X00603")
        ig = ig_["CAP04302"]
        pbm.A[ig,iv] += Float64(.02743)
        ig = ig_["CAP03802"]
        pbm.A[ig,iv] += Float64(.02733)
        iv,ix_,_ = s2mpj_ii("X00603",ix_)
        arrset(pb.xnames,iv,"X00603")
        ig = ig_["CAP04002"]
        pbm.A[ig,iv] += Float64(.01621)
        ig = ig_["CAP04603"]
        pbm.A[ig,iv] += Float64(.00104)
        iv,ix_,_ = s2mpj_ii("X00603",ix_)
        arrset(pb.xnames,iv,"X00603")
        ig = ig_["CAP05002"]
        pbm.A[ig,iv] += Float64(.01633)
        ig = ig_["CAP03702"]
        pbm.A[ig,iv] += Float64(.03304)
        iv,ix_,_ = s2mpj_ii("X00603",ix_)
        arrset(pb.xnames,iv,"X00603")
        ig = ig_["CAP04303"]
        pbm.A[ig,iv] += Float64(.00576)
        ig = ig_["CAP03402"]
        pbm.A[ig,iv] += Float64(.02558)
        iv,ix_,_ = s2mpj_ii("X00603",ix_)
        arrset(pb.xnames,iv,"X00603")
        ig = ig_["CAP04102"]
        pbm.A[ig,iv] += Float64(.0117)
        ig = ig_["CAP04703"]
        pbm.A[ig,iv] += Float64(.01583)
        iv,ix_,_ = s2mpj_ii("X00604",ix_)
        arrset(pb.xnames,iv,"X00604")
        ig = ig_["CAP04903"]
        pbm.A[ig,iv] += Float64(.01393)
        ig = ig_["CAP03804"]
        pbm.A[ig,iv] += Float64(.01706)
        iv,ix_,_ = s2mpj_ii("X00604",ix_)
        arrset(pb.xnames,iv,"X00604")
        ig = ig_["CAP04304"]
        pbm.A[ig,iv] += Float64(.00576)
        ig = ig_["CAP05003"]
        pbm.A[ig,iv] += Float64(.01633)
        iv,ix_,_ = s2mpj_ii("X00604",ix_)
        arrset(pb.xnames,iv,"X00604")
        ig = ig_["CAP04104"]
        pbm.A[ig,iv] += Float64(.0005)
        ig = ig_["CAP03704"]
        pbm.A[ig,iv] += Float64(.00362)
        iv,ix_,_ = s2mpj_ii("X00604",ix_)
        arrset(pb.xnames,iv,"X00604")
        ig = ig_["CAP04204"]
        pbm.A[ig,iv] += Float64(.02393)
        ig = ig_["CAP04303"]
        pbm.A[ig,iv] += Float64(.02743)
        iv,ix_,_ = s2mpj_ii("X00604",ix_)
        arrset(pb.xnames,iv,"X00604")
        ig = ig_["CAP04203"]
        pbm.A[ig,iv] += Float64(.01126)
        ig = ig_["CAP04004"]
        pbm.A[ig,iv] += Float64(.02556)
        iv,ix_,_ = s2mpj_ii("X00604",ix_)
        arrset(pb.xnames,iv,"X00604")
        ig = ig_["CAP03403"]
        pbm.A[ig,iv] += Float64(.02558)
        ig = ig_["CAP05004"]
        pbm.A[ig,iv] += Float64(.00685)
        iv,ix_,_ = s2mpj_ii("X00604",ix_)
        arrset(pb.xnames,iv,"X00604")
        ig = ig_["CAP04703"]
        pbm.A[ig,iv] += Float64(.01004)
        ig = ig_["CAP05104"]
        pbm.A[ig,iv] += Float64(.00439)
        iv,ix_,_ = s2mpj_ii("X00604",ix_)
        arrset(pb.xnames,iv,"X00604")
        ig = ig_["INV00404"]
        pbm.A[ig,iv] += Float64(143.0)
        ig = ig_["CAP03304"]
        pbm.A[ig,iv] += Float64(.00255)
        iv,ix_,_ = s2mpj_ii("X00604",ix_)
        arrset(pb.xnames,iv,"X00604")
        ig = ig_["CAP05303"]
        pbm.A[ig,iv] += Float64(.01148)
        ig = ig_["CAP04603"]
        pbm.A[ig,iv] += Float64(.01593)
        iv,ix_,_ = s2mpj_ii("X00604",ix_)
        arrset(pb.xnames,iv,"X00604")
        ig = ig_["CAP03404"]
        pbm.A[ig,iv] += Float64(.03192)
        ig = ig_["CAP03803"]
        pbm.A[ig,iv] += Float64(.02733)
        iv,ix_,_ = s2mpj_ii("X00604",ix_)
        arrset(pb.xnames,iv,"X00604")
        ig = ig_["CAP04704"]
        pbm.A[ig,iv] += Float64(.01583)
        ig = ig_["CAP04904"]
        pbm.A[ig,iv] += Float64(.00431)
        iv,ix_,_ = s2mpj_ii("X00604",ix_)
        arrset(pb.xnames,iv,"X00604")
        ig = ig_["CAP04503"]
        pbm.A[ig,iv] += Float64(.02252)
        ig = ig_["CAP04604"]
        pbm.A[ig,iv] += Float64(.00104)
        iv,ix_,_ = s2mpj_ii("X00604",ix_)
        arrset(pb.xnames,iv,"X00604")
        ig = ig_["CAP05103"]
        pbm.A[ig,iv] += Float64(.0132)
        ig = ig_["CAP05403"]
        pbm.A[ig,iv] += Float64(.00786)
        iv,ix_,_ = s2mpj_ii("X00604",ix_)
        arrset(pb.xnames,iv,"X00604")
        ig = ig_["CAP04404"]
        pbm.A[ig,iv] += Float64(.01286)
        ig = ig_["CAP05404"]
        pbm.A[ig,iv] += Float64(.01089)
        iv,ix_,_ = s2mpj_ii("X00604",ix_)
        arrset(pb.xnames,iv,"X00604")
        ig = ig_["CAP03703"]
        pbm.A[ig,iv] += Float64(.03304)
        ig = ig_["CAP04403"]
        pbm.A[ig,iv] += Float64(.00521)
        iv,ix_,_ = s2mpj_ii("X00604",ix_)
        arrset(pb.xnames,iv,"X00604")
        ig = ig_["CAP05304"]
        pbm.A[ig,iv] += Float64(.00591)
        ig = ig_["CAP04003"]
        pbm.A[ig,iv] += Float64(.01621)
        iv,ix_,_ = s2mpj_ii("X00604",ix_)
        arrset(pb.xnames,iv,"X00604")
        ig = ig_["CAP04103"]
        pbm.A[ig,iv] += Float64(.01217)
        ig = ig_["CAP04504"]
        pbm.A[ig,iv] += Float64(.01688)
        iv,ix_,_ = s2mpj_ii("X00604",ix_)
        arrset(pb.xnames,iv,"X00604")
        ig = ig_["CAP03303"]
        pbm.A[ig,iv] += Float64(.01554)
        iv,ix_,_ = s2mpj_ii("X00605",ix_)
        arrset(pb.xnames,iv,"X00605")
        ig = ig_["CAP04105"]
        pbm.A[ig,iv] += Float64(.00095)
        ig = ig_["CAP04505"]
        pbm.A[ig,iv] += Float64(.01772)
        iv,ix_,_ = s2mpj_ii("X00605",ix_)
        arrset(pb.xnames,iv,"X00605")
        ig = ig_["CAP04005"]
        pbm.A[ig,iv] += Float64(.02616)
        ig = ig_["CAP03805"]
        pbm.A[ig,iv] += Float64(.01807)
        iv,ix_,_ = s2mpj_ii("X00605",ix_)
        arrset(pb.xnames,iv,"X00605")
        ig = ig_["CAP04205"]
        pbm.A[ig,iv] += Float64(.02435)
        ig = ig_["CAP03705"]
        pbm.A[ig,iv] += Float64(.00485)
        iv,ix_,_ = s2mpj_ii("X00605",ix_)
        arrset(pb.xnames,iv,"X00605")
        ig = ig_["CAP03405"]
        pbm.A[ig,iv] += Float64(.03287)
        ig = ig_["CAP04305"]
        pbm.A[ig,iv] += Float64(.00678)
        iv,ix_,_ = s2mpj_ii("X00605",ix_)
        arrset(pb.xnames,iv,"X00605")
        ig = ig_["CAP04405"]
        pbm.A[ig,iv] += Float64(.01306)
        ig = ig_["CAP04605"]
        pbm.A[ig,iv] += Float64(.00163)
        iv,ix_,_ = s2mpj_ii("X00605",ix_)
        arrset(pb.xnames,iv,"X00605")
        ig = ig_["CAP05405"]
        pbm.A[ig,iv] += Float64(.01118)
        ig = ig_["CAP04705"]
        pbm.A[ig,iv] += Float64(.0162)
        iv,ix_,_ = s2mpj_ii("X00605",ix_)
        arrset(pb.xnames,iv,"X00605")
        ig = ig_["CAP04905"]
        pbm.A[ig,iv] += Float64(.00483)
        ig = ig_["INV00405"]
        pbm.A[ig,iv] += Float64(143.0)
        iv,ix_,_ = s2mpj_ii("X00605",ix_)
        arrset(pb.xnames,iv,"X00605")
        ig = ig_["CAP05005"]
        pbm.A[ig,iv] += Float64(.00745)
        ig = ig_["CAP05105"]
        pbm.A[ig,iv] += Float64(.00488)
        iv,ix_,_ = s2mpj_ii("X00605",ix_)
        arrset(pb.xnames,iv,"X00605")
        ig = ig_["CAP05305"]
        pbm.A[ig,iv] += Float64(.00634)
        ig = ig_["CAP03305"]
        pbm.A[ig,iv] += Float64(.00313)
        iv,ix_,_ = s2mpj_ii("X00605",ix_)
        arrset(pb.xnames,iv,"X00605")
        ig = ig_["CAP04304"]
        pbm.A[ig,iv] += Float64(.02642)
        ig = ig_["CAP04204"]
        pbm.A[ig,iv] += Float64(.01084)
        iv,ix_,_ = s2mpj_ii("X00605",ix_)
        arrset(pb.xnames,iv,"X00605")
        ig = ig_["CAP05104"]
        pbm.A[ig,iv] += Float64(.01271)
        ig = ig_["CAP05304"]
        pbm.A[ig,iv] += Float64(.01106)
        iv,ix_,_ = s2mpj_ii("X00605",ix_)
        arrset(pb.xnames,iv,"X00605")
        ig = ig_["CAP03704"]
        pbm.A[ig,iv] += Float64(.03181)
        ig = ig_["CAP03304"]
        pbm.A[ig,iv] += Float64(.01496)
        iv,ix_,_ = s2mpj_ii("X00605",ix_)
        arrset(pb.xnames,iv,"X00605")
        ig = ig_["CAP04704"]
        pbm.A[ig,iv] += Float64(.00967)
        ig = ig_["CAP03804"]
        pbm.A[ig,iv] += Float64(.02632)
        iv,ix_,_ = s2mpj_ii("X00605",ix_)
        arrset(pb.xnames,iv,"X00605")
        ig = ig_["CAP04404"]
        pbm.A[ig,iv] += Float64(.00502)
        ig = ig_["CAP05404"]
        pbm.A[ig,iv] += Float64(.00757)
        iv,ix_,_ = s2mpj_ii("X00605",ix_)
        arrset(pb.xnames,iv,"X00605")
        ig = ig_["CAP03404"]
        pbm.A[ig,iv] += Float64(.02464)
        ig = ig_["CAP04504"]
        pbm.A[ig,iv] += Float64(.02169)
        iv,ix_,_ = s2mpj_ii("X00605",ix_)
        arrset(pb.xnames,iv,"X00605")
        ig = ig_["CAP04604"]
        pbm.A[ig,iv] += Float64(.01534)
        ig = ig_["CAP04904"]
        pbm.A[ig,iv] += Float64(.01341)
        iv,ix_,_ = s2mpj_ii("X00605",ix_)
        arrset(pb.xnames,iv,"X00605")
        ig = ig_["CAP04104"]
        pbm.A[ig,iv] += Float64(.01172)
        ig = ig_["CAP04004"]
        pbm.A[ig,iv] += Float64(.01561)
        iv,ix_,_ = s2mpj_ii("X00605",ix_)
        arrset(pb.xnames,iv,"X00605")
        ig = ig_["CAP05004"]
        pbm.A[ig,iv] += Float64(.01573)
        iv,ix_,_ = s2mpj_ii("X00606",ix_)
        arrset(pb.xnames,iv,"X00606")
        ig = ig_["CAP04606"]
        pbm.A[ig,iv] += Float64(.01696)
        ig = ig_["CAP04206"]
        pbm.A[ig,iv] += Float64(.03519)
        iv,ix_,_ = s2mpj_ii("X00606",ix_)
        arrset(pb.xnames,iv,"X00606")
        ig = ig_["CAP04305"]
        pbm.A[ig,iv] += Float64(.02853)
        ig = ig_["CAP04106"]
        pbm.A[ig,iv] += Float64(.01267)
        iv,ix_,_ = s2mpj_ii("X00606",ix_)
        arrset(pb.xnames,iv,"X00606")
        ig = ig_["CAP04006"]
        pbm.A[ig,iv] += Float64(.04177)
        ig = ig_["CAP04205"]
        pbm.A[ig,iv] += Float64(.01171)
        iv,ix_,_ = s2mpj_ii("X00606",ix_)
        arrset(pb.xnames,iv,"X00606")
        ig = ig_["CAP04906"]
        pbm.A[ig,iv] += Float64(.01824)
        ig = ig_["CAP04405"]
        pbm.A[ig,iv] += Float64(.00542)
        iv,ix_,_ = s2mpj_ii("X00606",ix_)
        arrset(pb.xnames,iv,"X00606")
        ig = ig_["CAP03806"]
        pbm.A[ig,iv] += Float64(.04439)
        ig = ig_["CAP03705"]
        pbm.A[ig,iv] += Float64(.03436)
        iv,ix_,_ = s2mpj_ii("X00606",ix_)
        arrset(pb.xnames,iv,"X00606")
        ig = ig_["CAP04306"]
        pbm.A[ig,iv] += Float64(.03319)
        ig = ig_["CAP04506"]
        pbm.A[ig,iv] += Float64(.0394)
        iv,ix_,_ = s2mpj_ii("X00606",ix_)
        arrset(pb.xnames,iv,"X00606")
        ig = ig_["CAP04406"]
        pbm.A[ig,iv] += Float64(.01808)
        ig = ig_["CAP03805"]
        pbm.A[ig,iv] += Float64(.02843)
        iv,ix_,_ = s2mpj_ii("X00606",ix_)
        arrset(pb.xnames,iv,"X00606")
        ig = ig_["CAP03305"]
        pbm.A[ig,iv] += Float64(.01616)
        ig = ig_["CAP05106"]
        pbm.A[ig,iv] += Float64(.01759)
        iv,ix_,_ = s2mpj_ii("X00606",ix_)
        arrset(pb.xnames,iv,"X00606")
        ig = ig_["CAP05406"]
        pbm.A[ig,iv] += Float64(.01875)
        ig = ig_["CAP05006"]
        pbm.A[ig,iv] += Float64(.02318)
        iv,ix_,_ = s2mpj_ii("X00606",ix_)
        arrset(pb.xnames,iv,"X00606")
        ig = ig_["CAP04005"]
        pbm.A[ig,iv] += Float64(.01686)
        ig = ig_["CAP03405"]
        pbm.A[ig,iv] += Float64(.02661)
        iv,ix_,_ = s2mpj_ii("X00606",ix_)
        arrset(pb.xnames,iv,"X00606")
        ig = ig_["CAP04105"]
        pbm.A[ig,iv] += Float64(.01266)
        ig = ig_["CAP05306"]
        pbm.A[ig,iv] += Float64(.0174)
        iv,ix_,_ = s2mpj_ii("X00606",ix_)
        arrset(pb.xnames,iv,"X00606")
        ig = ig_["CAP04706"]
        pbm.A[ig,iv] += Float64(.02586)
        ig = ig_["CAP03306"]
        pbm.A[ig,iv] += Float64(.01809)
        iv,ix_,_ = s2mpj_ii("X00606",ix_)
        arrset(pb.xnames,iv,"X00606")
        ig = ig_["CAP04505"]
        pbm.A[ig,iv] += Float64(.02342)
        ig = ig_["CAP05105"]
        pbm.A[ig,iv] += Float64(.01373)
        iv,ix_,_ = s2mpj_ii("X00606",ix_)
        arrset(pb.xnames,iv,"X00606")
        ig = ig_["CAP05405"]
        pbm.A[ig,iv] += Float64(.00817)
        ig = ig_["CAP05005"]
        pbm.A[ig,iv] += Float64(.01699)
        iv,ix_,_ = s2mpj_ii("X00606",ix_)
        arrset(pb.xnames,iv,"X00606")
        ig = ig_["CAP05305"]
        pbm.A[ig,iv] += Float64(.01194)
        ig = ig_["CAP04905"]
        pbm.A[ig,iv] += Float64(.01448)
        iv,ix_,_ = s2mpj_ii("X00606",ix_)
        arrset(pb.xnames,iv,"X00606")
        ig = ig_["CAP03406"]
        pbm.A[ig,iv] += Float64(.05751)
        ig = ig_["CAP04705"]
        pbm.A[ig,iv] += Float64(.01044)
        iv,ix_,_ = s2mpj_ii("X00606",ix_)
        arrset(pb.xnames,iv,"X00606")
        ig = ig_["CAP04605"]
        pbm.A[ig,iv] += Float64(.01656)
        ig = ig_["CAP03706"]
        pbm.A[ig,iv] += Float64(.03666)
        iv,ix_,_ = s2mpj_ii("X00703",ix_)
        arrset(pb.xnames,iv,"X00703")
        ig = ig_["CAP03303"]
        pbm.A[ig,iv] += Float64(.00376)
        ig = ig_["CAP03902"]
        pbm.A[ig,iv] += Float64(.00969)
        iv,ix_,_ = s2mpj_ii("X00703",ix_)
        arrset(pb.xnames,iv,"X00703")
        ig = ig_["CAP03301"]
        pbm.A[ig,iv] += Float64(.00093)
        ig = ig_["CAP03802"]
        pbm.A[ig,iv] += Float64(.01945)
        iv,ix_,_ = s2mpj_ii("X00703",ix_)
        arrset(pb.xnames,iv,"X00703")
        ig = ig_["CAP05102"]
        pbm.A[ig,iv] += Float64(.00862)
        ig = ig_["CAP04102"]
        pbm.A[ig,iv] += Float64(.01116)
        iv,ix_,_ = s2mpj_ii("X00703",ix_)
        arrset(pb.xnames,iv,"X00703")
        ig = ig_["CAP04203"]
        pbm.A[ig,iv] += Float64(.0244)
        ig = ig_["CAP04402"]
        pbm.A[ig,iv] += Float64(.00496)
        iv,ix_,_ = s2mpj_ii("X00703",ix_)
        arrset(pb.xnames,iv,"X00703")
        ig = ig_["CAP04202"]
        pbm.A[ig,iv] += Float64(.00978)
        ig = ig_["CAP05703"]
        pbm.A[ig,iv] += Float64(.02434)
        iv,ix_,_ = s2mpj_ii("X00703",ix_)
        arrset(pb.xnames,iv,"X00703")
        ig = ig_["CAP03502"]
        pbm.A[ig,iv] += Float64(.01463)
        ig = ig_["CAP03302"]
        pbm.A[ig,iv] += Float64(.01962)
        iv,ix_,_ = s2mpj_ii("X00703",ix_)
        arrset(pb.xnames,iv,"X00703")
        ig = ig_["CAP05903"]
        pbm.A[ig,iv] += Float64(.02055)
        ig = ig_["CAP03903"]
        pbm.A[ig,iv] += Float64(.0079)
        iv,ix_,_ = s2mpj_ii("X00703",ix_)
        arrset(pb.xnames,iv,"X00703")
        ig = ig_["CAP04602"]
        pbm.A[ig,iv] += Float64(.01478)
        ig = ig_["CAP03701"]
        pbm.A[ig,iv] += Float64(.0026)
        iv,ix_,_ = s2mpj_ii("X00703",ix_)
        arrset(pb.xnames,iv,"X00703")
        ig = ig_["CAP03402"]
        pbm.A[ig,iv] += Float64(.0313)
        ig = ig_["CAP05803"]
        pbm.A[ig,iv] += Float64(.01677)
        iv,ix_,_ = s2mpj_ii("X00703",ix_)
        arrset(pb.xnames,iv,"X00703")
        ig = ig_["CAP03803"]
        pbm.A[ig,iv] += Float64(.0227)
        ig = ig_["CAP04303"]
        pbm.A[ig,iv] += Float64(.00396)
        iv,ix_,_ = s2mpj_ii("X00703",ix_)
        arrset(pb.xnames,iv,"X00703")
        ig = ig_["CAP04403"]
        pbm.A[ig,iv] += Float64(.01299)
        ig = ig_["CAP03702"]
        pbm.A[ig,iv] += Float64(.0251)
        iv,ix_,_ = s2mpj_ii("X00703",ix_)
        arrset(pb.xnames,iv,"X00703")
        ig = ig_["CAP04302"]
        pbm.A[ig,iv] += Float64(.02908)
        ig = ig_["CAP05302"]
        pbm.A[ig,iv] += Float64(.01822)
        iv,ix_,_ = s2mpj_ii("X00703",ix_)
        arrset(pb.xnames,iv,"X00703")
        ig = ig_["CAP04301"]
        pbm.A[ig,iv] += Float64(.00012)
        ig = ig_["CAP05303"]
        pbm.A[ig,iv] += Float64(.01977)
        iv,ix_,_ = s2mpj_ii("X00703",ix_)
        arrset(pb.xnames,iv,"X00703")
        ig = ig_["CAP03503"]
        pbm.A[ig,iv] += Float64(.03477)
        ig = ig_["CAP05902"]
        pbm.A[ig,iv] += Float64(.01455)
        iv,ix_,_ = s2mpj_ii("X00703",ix_)
        arrset(pb.xnames,iv,"X00703")
        ig = ig_["INV00503"]
        pbm.A[ig,iv] += Float64(424.0)
        ig = ig_["CAP05103"]
        pbm.A[ig,iv] += Float64(.00539)
        iv,ix_,_ = s2mpj_ii("X00703",ix_)
        arrset(pb.xnames,iv,"X00703")
        ig = ig_["CAP05402"]
        pbm.A[ig,iv] += Float64(.00542)
        ig = ig_["CAP04601"]
        pbm.A[ig,iv] += Float64(.00258)
        iv,ix_,_ = s2mpj_ii("X00703",ix_)
        arrset(pb.xnames,iv,"X00703")
        ig = ig_["CAP05502"]
        pbm.A[ig,iv] += Float64(.01259)
        ig = ig_["CAP05501"]
        pbm.A[ig,iv] += Float64(.00229)
        iv,ix_,_ = s2mpj_ii("X00703",ix_)
        arrset(pb.xnames,iv,"X00703")
        ig = ig_["CAP05602"]
        pbm.A[ig,iv] += Float64(.02266)
        ig = ig_["CAP04101"]
        pbm.A[ig,iv] += Float64(.00203)
        iv,ix_,_ = s2mpj_ii("X00703",ix_)
        arrset(pb.xnames,iv,"X00703")
        ig = ig_["CAP05702"]
        pbm.A[ig,iv] += Float64(.01024)
        ig = ig_["CAP03403"]
        pbm.A[ig,iv] += Float64(.04854)
        iv,ix_,_ = s2mpj_ii("X00703",ix_)
        arrset(pb.xnames,iv,"X00703")
        ig = ig_["CAP05403"]
        pbm.A[ig,iv] += Float64(.01247)
        ig = ig_["CAP05802"]
        pbm.A[ig,iv] += Float64(.01936)
        iv,ix_,_ = s2mpj_ii("X00703",ix_)
        arrset(pb.xnames,iv,"X00703")
        ig = ig_["CAP05603"]
        pbm.A[ig,iv] += Float64(.01208)
        iv,ix_,_ = s2mpj_ii("X00704",ix_)
        arrset(pb.xnames,iv,"X00704")
        ig = ig_["CAP05304"]
        pbm.A[ig,iv] += Float64(.01977)
        ig = ig_["CAP04204"]
        pbm.A[ig,iv] += Float64(.0244)
        iv,ix_,_ = s2mpj_ii("X00704",ix_)
        arrset(pb.xnames,iv,"X00704")
        ig = ig_["CAP05903"]
        pbm.A[ig,iv] += Float64(.01455)
        ig = ig_["CAP04303"]
        pbm.A[ig,iv] += Float64(.0292)
        iv,ix_,_ = s2mpj_ii("X00704",ix_)
        arrset(pb.xnames,iv,"X00704")
        ig = ig_["CAP03703"]
        pbm.A[ig,iv] += Float64(.02723)
        ig = ig_["CAP05403"]
        pbm.A[ig,iv] += Float64(.00542)
        iv,ix_,_ = s2mpj_ii("X00704",ix_)
        arrset(pb.xnames,iv,"X00704")
        ig = ig_["CAP05303"]
        pbm.A[ig,iv] += Float64(.01822)
        ig = ig_["CAP05803"]
        pbm.A[ig,iv] += Float64(.01936)
        iv,ix_,_ = s2mpj_ii("X00704",ix_)
        arrset(pb.xnames,iv,"X00704")
        ig = ig_["CAP05103"]
        pbm.A[ig,iv] += Float64(.00862)
        ig = ig_["INV00504"]
        pbm.A[ig,iv] += Float64(424.0)
        iv,ix_,_ = s2mpj_ii("X00704",ix_)
        arrset(pb.xnames,iv,"X00704")
        ig = ig_["CAP05704"]
        pbm.A[ig,iv] += Float64(.02434)
        ig = ig_["CAP03504"]
        pbm.A[ig,iv] += Float64(.03477)
        iv,ix_,_ = s2mpj_ii("X00704",ix_)
        arrset(pb.xnames,iv,"X00704")
        ig = ig_["CAP03302"]
        pbm.A[ig,iv] += Float64(.00045)
        ig = ig_["CAP05604"]
        pbm.A[ig,iv] += Float64(.01208)
        iv,ix_,_ = s2mpj_ii("X00704",ix_)
        arrset(pb.xnames,iv,"X00704")
        ig = ig_["CAP04103"]
        pbm.A[ig,iv] += Float64(.01218)
        ig = ig_["CAP05404"]
        pbm.A[ig,iv] += Float64(.01247)
        iv,ix_,_ = s2mpj_ii("X00704",ix_)
        arrset(pb.xnames,iv,"X00704")
        ig = ig_["CAP04602"]
        pbm.A[ig,iv] += Float64(.00125)
        ig = ig_["CAP03404"]
        pbm.A[ig,iv] += Float64(.04854)
        iv,ix_,_ = s2mpj_ii("X00704",ix_)
        arrset(pb.xnames,iv,"X00704")
        ig = ig_["CAP05603"]
        pbm.A[ig,iv] += Float64(.02266)
        ig = ig_["CAP04404"]
        pbm.A[ig,iv] += Float64(.01299)
        iv,ix_,_ = s2mpj_ii("X00704",ix_)
        arrset(pb.xnames,iv,"X00704")
        ig = ig_["CAP05804"]
        pbm.A[ig,iv] += Float64(.01677)
        ig = ig_["CAP04304"]
        pbm.A[ig,iv] += Float64(.00396)
        iv,ix_,_ = s2mpj_ii("X00704",ix_)
        arrset(pb.xnames,iv,"X00704")
        ig = ig_["CAP04203"]
        pbm.A[ig,iv] += Float64(.00978)
        ig = ig_["CAP03304"]
        pbm.A[ig,iv] += Float64(.00376)
        iv,ix_,_ = s2mpj_ii("X00704",ix_)
        arrset(pb.xnames,iv,"X00704")
        ig = ig_["CAP03904"]
        pbm.A[ig,iv] += Float64(.0079)
        ig = ig_["CAP03803"]
        pbm.A[ig,iv] += Float64(.01945)
        iv,ix_,_ = s2mpj_ii("X00704",ix_)
        arrset(pb.xnames,iv,"X00704")
        ig = ig_["CAP03503"]
        pbm.A[ig,iv] += Float64(.01463)
        ig = ig_["CAP05104"]
        pbm.A[ig,iv] += Float64(.00539)
        iv,ix_,_ = s2mpj_ii("X00704",ix_)
        arrset(pb.xnames,iv,"X00704")
        ig = ig_["CAP03403"]
        pbm.A[ig,iv] += Float64(.0313)
        ig = ig_["CAP05503"]
        pbm.A[ig,iv] += Float64(.01374)
        iv,ix_,_ = s2mpj_ii("X00704",ix_)
        arrset(pb.xnames,iv,"X00704")
        ig = ig_["CAP04102"]
        pbm.A[ig,iv] += Float64(.00101)
        ig = ig_["CAP04403"]
        pbm.A[ig,iv] += Float64(.00496)
        iv,ix_,_ = s2mpj_ii("X00704",ix_)
        arrset(pb.xnames,iv,"X00704")
        ig = ig_["CAP04603"]
        pbm.A[ig,iv] += Float64(.01612)
        ig = ig_["CAP05904"]
        pbm.A[ig,iv] += Float64(.02055)
        iv,ix_,_ = s2mpj_ii("X00704",ix_)
        arrset(pb.xnames,iv,"X00704")
        ig = ig_["CAP03303"]
        pbm.A[ig,iv] += Float64(.0201)
        ig = ig_["CAP03804"]
        pbm.A[ig,iv] += Float64(.0227)
        iv,ix_,_ = s2mpj_ii("X00704",ix_)
        arrset(pb.xnames,iv,"X00704")
        ig = ig_["CAP05502"]
        pbm.A[ig,iv] += Float64(.00114)
        ig = ig_["CAP03702"]
        pbm.A[ig,iv] += Float64(.00047)
        iv,ix_,_ = s2mpj_ii("X00704",ix_)
        arrset(pb.xnames,iv,"X00704")
        ig = ig_["CAP05703"]
        pbm.A[ig,iv] += Float64(.01024)
        ig = ig_["CAP03903"]
        pbm.A[ig,iv] += Float64(.00969)
        iv,ix_,_ = s2mpj_ii("X00705",ix_)
        arrset(pb.xnames,iv,"X00705")
        ig = ig_["CAP04104"]
        pbm.A[ig,iv] += Float64(.01221)
        ig = ig_["CAP03303"]
        pbm.A[ig,iv] += Float64(.00043)
        iv,ix_,_ = s2mpj_ii("X00705",ix_)
        arrset(pb.xnames,iv,"X00705")
        ig = ig_["CAP05605"]
        pbm.A[ig,iv] += Float64(.01292)
        ig = ig_["CAP05704"]
        pbm.A[ig,iv] += Float64(.00986)
        iv,ix_,_ = s2mpj_ii("X00705",ix_)
        arrset(pb.xnames,iv,"X00705")
        ig = ig_["INV00505"]
        pbm.A[ig,iv] += Float64(424.0)
        ig = ig_["CAP05105"]
        pbm.A[ig,iv] += Float64(.00571)
        iv,ix_,_ = s2mpj_ii("X00705",ix_)
        arrset(pb.xnames,iv,"X00705")
        ig = ig_["CAP05904"]
        pbm.A[ig,iv] += Float64(.01401)
        ig = ig_["CAP05405"]
        pbm.A[ig,iv] += Float64(.01267)
        iv,ix_,_ = s2mpj_ii("X00705",ix_)
        arrset(pb.xnames,iv,"X00705")
        ig = ig_["CAP05905"]
        pbm.A[ig,iv] += Float64(.02109)
        ig = ig_["CAP03704"]
        pbm.A[ig,iv] += Float64(.02668)
        iv,ix_,_ = s2mpj_ii("X00705",ix_)
        arrset(pb.xnames,iv,"X00705")
        ig = ig_["CAP05705"]
        pbm.A[ig,iv] += Float64(.02472)
        ig = ig_["CAP03904"]
        pbm.A[ig,iv] += Float64(.00933)
        iv,ix_,_ = s2mpj_ii("X00705",ix_)
        arrset(pb.xnames,iv,"X00705")
        ig = ig_["CAP05804"]
        pbm.A[ig,iv] += Float64(.01864)
        ig = ig_["CAP05805"]
        pbm.A[ig,iv] += Float64(.01749)
        iv,ix_,_ = s2mpj_ii("X00705",ix_)
        arrset(pb.xnames,iv,"X00705")
        ig = ig_["CAP04204"]
        pbm.A[ig,iv] += Float64(.00942)
        ig = ig_["CAP05503"]
        pbm.A[ig,iv] += Float64(.0011)
        iv,ix_,_ = s2mpj_ii("X00705",ix_)
        arrset(pb.xnames,iv,"X00705")
        ig = ig_["CAP03804"]
        pbm.A[ig,iv] += Float64(.01873)
        ig = ig_["CAP05305"]
        pbm.A[ig,iv] += Float64(.02044)
        iv,ix_,_ = s2mpj_ii("X00705",ix_)
        arrset(pb.xnames,iv,"X00705")
        ig = ig_["CAP05404"]
        pbm.A[ig,iv] += Float64(.00522)
        ig = ig_["CAP03305"]
        pbm.A[ig,iv] += Float64(.00429)
        iv,ix_,_ = s2mpj_ii("X00705",ix_)
        arrset(pb.xnames,iv,"X00705")
        ig = ig_["CAP04604"]
        pbm.A[ig,iv] += Float64(.01616)
        ig = ig_["CAP04305"]
        pbm.A[ig,iv] += Float64(.00504)
        iv,ix_,_ = s2mpj_ii("X00705",ix_)
        arrset(pb.xnames,iv,"X00705")
        ig = ig_["CAP03304"]
        pbm.A[ig,iv] += Float64(.01958)
        ig = ig_["CAP04205"]
        pbm.A[ig,iv] += Float64(.02476)
        iv,ix_,_ = s2mpj_ii("X00705",ix_)
        arrset(pb.xnames,iv,"X00705")
        ig = ig_["CAP03703"]
        pbm.A[ig,iv] += Float64(.00045)
        ig = ig_["CAP05104"]
        pbm.A[ig,iv] += Float64(.0083)
        iv,ix_,_ = s2mpj_ii("X00705",ix_)
        arrset(pb.xnames,iv,"X00705")
        ig = ig_["CAP03905"]
        pbm.A[ig,iv] += Float64(.00826)
        ig = ig_["CAP04405"]
        pbm.A[ig,iv] += Float64(.01317)
        iv,ix_,_ = s2mpj_ii("X00705",ix_)
        arrset(pb.xnames,iv,"X00705")
        ig = ig_["CAP03805"]
        pbm.A[ig,iv] += Float64(.02342)
        ig = ig_["CAP04603"]
        pbm.A[ig,iv] += Float64(.0012)
        iv,ix_,_ = s2mpj_ii("X00705",ix_)
        arrset(pb.xnames,iv,"X00705")
        ig = ig_["CAP03405"]
        pbm.A[ig,iv] += Float64(.0497)
        ig = ig_["CAP03505"]
        pbm.A[ig,iv] += Float64(.03531)
        iv,ix_,_ = s2mpj_ii("X00705",ix_)
        arrset(pb.xnames,iv,"X00705")
        ig = ig_["CAP04304"]
        pbm.A[ig,iv] += Float64(.02812)
        ig = ig_["CAP03705"]
        pbm.A[ig,iv] += Float64(.00057)
        iv,ix_,_ = s2mpj_ii("X00705",ix_)
        arrset(pb.xnames,iv,"X00705")
        ig = ig_["CAP05604"]
        pbm.A[ig,iv] += Float64(.02182)
        ig = ig_["CAP04103"]
        pbm.A[ig,iv] += Float64(.00098)
        iv,ix_,_ = s2mpj_ii("X00705",ix_)
        arrset(pb.xnames,iv,"X00705")
        ig = ig_["CAP04404"]
        pbm.A[ig,iv] += Float64(.00477)
        ig = ig_["CAP03404"]
        pbm.A[ig,iv] += Float64(.03014)
        iv,ix_,_ = s2mpj_ii("X00705",ix_)
        arrset(pb.xnames,iv,"X00705")
        ig = ig_["CAP05504"]
        pbm.A[ig,iv] += Float64(.01378)
        ig = ig_["CAP05304"]
        pbm.A[ig,iv] += Float64(.01755)
        iv,ix_,_ = s2mpj_ii("X00705",ix_)
        arrset(pb.xnames,iv,"X00705")
        ig = ig_["CAP03504"]
        pbm.A[ig,iv] += Float64(.01409)
        iv,ix_,_ = s2mpj_ii("X00706",ix_)
        arrset(pb.xnames,iv,"X00706")
        ig = ig_["CAP05406"]
        pbm.A[ig,iv] += Float64(.01789)
        ig = ig_["CAP05906"]
        pbm.A[ig,iv] += Float64(.0351)
        iv,ix_,_ = s2mpj_ii("X00706",ix_)
        arrset(pb.xnames,iv,"X00706")
        ig = ig_["CAP04406"]
        pbm.A[ig,iv] += Float64(.01795)
        ig = ig_["CAP03306"]
        pbm.A[ig,iv] += Float64(.0243)
        iv,ix_,_ = s2mpj_ii("X00706",ix_)
        arrset(pb.xnames,iv,"X00706")
        ig = ig_["CAP05306"]
        pbm.A[ig,iv] += Float64(.03799)
        ig = ig_["CAP05506"]
        pbm.A[ig,iv] += Float64(.01488)
        iv,ix_,_ = s2mpj_ii("X00706",ix_)
        arrset(pb.xnames,iv,"X00706")
        ig = ig_["CAP05106"]
        pbm.A[ig,iv] += Float64(.01401)
        ig = ig_["CAP03706"]
        pbm.A[ig,iv] += Float64(.0277)
        iv,ix_,_ = s2mpj_ii("X00706",ix_)
        arrset(pb.xnames,iv,"X00706")
        ig = ig_["CAP04106"]
        pbm.A[ig,iv] += Float64(.01319)
        ig = ig_["CAP03806"]
        pbm.A[ig,iv] += Float64(.04215)
        iv,ix_,_ = s2mpj_ii("X00706",ix_)
        arrset(pb.xnames,iv,"X00706")
        ig = ig_["CAP05806"]
        pbm.A[ig,iv] += Float64(.03613)
        ig = ig_["CAP03906"]
        pbm.A[ig,iv] += Float64(.01759)
        iv,ix_,_ = s2mpj_ii("X00706",ix_)
        arrset(pb.xnames,iv,"X00706")
        ig = ig_["CAP03304"]
        pbm.A[ig,iv] += Float64(.00022)
        ig = ig_["CAP05706"]
        pbm.A[ig,iv] += Float64(.03458)
        iv,ix_,_ = s2mpj_ii("X00706",ix_)
        arrset(pb.xnames,iv,"X00706")
        ig = ig_["CAP03506"]
        pbm.A[ig,iv] += Float64(.0494)
        ig = ig_["CAP04306"]
        pbm.A[ig,iv] += Float64(.03316)
        iv,ix_,_ = s2mpj_ii("X00706",ix_)
        arrset(pb.xnames,iv,"X00706")
        ig = ig_["CAP03406"]
        pbm.A[ig,iv] += Float64(.07984)
        ig = ig_["CAP04606"]
        pbm.A[ig,iv] += Float64(.01736)
        iv,ix_,_ = s2mpj_ii("X00706",ix_)
        arrset(pb.xnames,iv,"X00706")
        ig = ig_["CAP04206"]
        pbm.A[ig,iv] += Float64(.03418)
        ig = ig_["CAP05606"]
        pbm.A[ig,iv] += Float64(.03474)
        iv,ix_,_ = s2mpj_ii("X00706",ix_)
        arrset(pb.xnames,iv,"X00706")
        ig = ig_["INV00506"]
        pbm.A[ig,iv] += Float64(424.0)
        ig = ig_["CAP03705"]
        pbm.A[ig,iv] += Float64(.0293)
        iv,ix_,_ = s2mpj_ii("X00706",ix_)
        arrset(pb.xnames,iv,"X00706")
        ig = ig_["CAP05405"]
        pbm.A[ig,iv] += Float64(.00564)
        ig = ig_["CAP05305"]
        pbm.A[ig,iv] += Float64(.01895)
        iv,ix_,_ = s2mpj_ii("X00706",ix_)
        arrset(pb.xnames,iv,"X00706")
        ig = ig_["CAP03305"]
        pbm.A[ig,iv] += Float64(.02162)
        ig = ig_["CAP05105"]
        pbm.A[ig,iv] += Float64(.00897)
        iv,ix_,_ = s2mpj_ii("X00706",ix_)
        arrset(pb.xnames,iv,"X00706")
        ig = ig_["CAP03405"]
        pbm.A[ig,iv] += Float64(.03255)
        ig = ig_["CAP03505"]
        pbm.A[ig,iv] += Float64(.01522)
        iv,ix_,_ = s2mpj_ii("X00706",ix_)
        arrset(pb.xnames,iv,"X00706")
        ig = ig_["CAP03805"]
        pbm.A[ig,iv] += Float64(.02023)
        ig = ig_["CAP05505"]
        pbm.A[ig,iv] += Float64(.01607)
        iv,ix_,_ = s2mpj_ii("X00706",ix_)
        arrset(pb.xnames,iv,"X00706")
        ig = ig_["CAP03905"]
        pbm.A[ig,iv] += Float64(.01007)
        ig = ig_["CAP04105"]
        pbm.A[ig,iv] += Float64(.01425)
        iv,ix_,_ = s2mpj_ii("X00706",ix_)
        arrset(pb.xnames,iv,"X00706")
        ig = ig_["CAP04205"]
        pbm.A[ig,iv] += Float64(.01017)
        ig = ig_["CAP04605"]
        pbm.A[ig,iv] += Float64(.01875)
        iv,ix_,_ = s2mpj_ii("X00706",ix_)
        arrset(pb.xnames,iv,"X00706")
        ig = ig_["CAP04305"]
        pbm.A[ig,iv] += Float64(.03037)
        ig = ig_["CAP04405"]
        pbm.A[ig,iv] += Float64(.00515)
        iv,ix_,_ = s2mpj_ii("X00706",ix_)
        arrset(pb.xnames,iv,"X00706")
        ig = ig_["CAP04604"]
        pbm.A[ig,iv] += Float64(.0006)
        ig = ig_["CAP04104"]
        pbm.A[ig,iv] += Float64(.00053)
        iv,ix_,_ = s2mpj_ii("X00706",ix_)
        arrset(pb.xnames,iv,"X00706")
        ig = ig_["CAP05605"]
        pbm.A[ig,iv] += Float64(.02356)
        ig = ig_["CAP05705"]
        pbm.A[ig,iv] += Float64(.01065)
        iv,ix_,_ = s2mpj_ii("X00706",ix_)
        arrset(pb.xnames,iv,"X00706")
        ig = ig_["CAP05504"]
        pbm.A[ig,iv] += Float64(.0006)
        ig = ig_["CAP05905"]
        pbm.A[ig,iv] += Float64(.01513)
        iv,ix_,_ = s2mpj_ii("X00706",ix_)
        arrset(pb.xnames,iv,"X00706")
        ig = ig_["CAP05805"]
        pbm.A[ig,iv] += Float64(.02013)
        iv,ix_,_ = s2mpj_ii("X00803",ix_)
        arrset(pb.xnames,iv,"X00803")
        ig = ig_["CAP05103"]
        pbm.A[ig,iv] += Float64(.00692)
        ig = ig_["CAP05102"]
        pbm.A[ig,iv] += Float64(.02245)
        iv,ix_,_ = s2mpj_ii("X00803",ix_)
        arrset(pb.xnames,iv,"X00803")
        ig = ig_["CAP04403"]
        pbm.A[ig,iv] += Float64(.01289)
        ig = ig_["CAP03701"]
        pbm.A[ig,iv] += Float64(.00081)
        iv,ix_,_ = s2mpj_ii("X00803",ix_)
        arrset(pb.xnames,iv,"X00803")
        ig = ig_["CAP05203"]
        pbm.A[ig,iv] += Float64(.0323)
        ig = ig_["CAP05002"]
        pbm.A[ig,iv] += Float64(.01034)
        iv,ix_,_ = s2mpj_ii("X00803",ix_)
        arrset(pb.xnames,iv,"X00803")
        ig = ig_["CAP05403"]
        pbm.A[ig,iv] += Float64(.01151)
        ig = ig_["CAP04902"]
        pbm.A[ig,iv] += Float64(.01364)
        iv,ix_,_ = s2mpj_ii("X00803",ix_)
        arrset(pb.xnames,iv,"X00803")
        ig = ig_["CAP04601"]
        pbm.A[ig,iv] += Float64(.00097)
        ig = ig_["CAP03603"]
        pbm.A[ig,iv] += Float64(.01934)
        iv,ix_,_ = s2mpj_ii("X00803",ix_)
        arrset(pb.xnames,iv,"X00803")
        ig = ig_["CAP05903"]
        pbm.A[ig,iv] += Float64(.02063)
        ig = ig_["CAP05902"]
        pbm.A[ig,iv] += Float64(.01835)
        iv,ix_,_ = s2mpj_ii("X00803",ix_)
        arrset(pb.xnames,iv,"X00803")
        ig = ig_["CAP04903"]
        pbm.A[ig,iv] += Float64(.00453)
        ig = ig_["CAP04803"]
        pbm.A[ig,iv] += Float64(.00494)
        iv,ix_,_ = s2mpj_ii("X00803",ix_)
        arrset(pb.xnames,iv,"X00803")
        ig = ig_["CAP05802"]
        pbm.A[ig,iv] += Float64(.02611)
        ig = ig_["CAP05702"]
        pbm.A[ig,iv] += Float64(.01309)
        iv,ix_,_ = s2mpj_ii("X00803",ix_)
        arrset(pb.xnames,iv,"X00803")
        ig = ig_["CAP05602"]
        pbm.A[ig,iv] += Float64(.03388)
        ig = ig_["CAP05402"]
        pbm.A[ig,iv] += Float64(.00679)
        iv,ix_,_ = s2mpj_ii("X00803",ix_)
        arrset(pb.xnames,iv,"X00803")
        ig = ig_["CAP04101"]
        pbm.A[ig,iv] += Float64(.00098)
        ig = ig_["CAP04603"]
        pbm.A[ig,iv] += Float64(.00037)
        iv,ix_,_ = s2mpj_ii("X00803",ix_)
        arrset(pb.xnames,iv,"X00803")
        ig = ig_["CAP05603"]
        pbm.A[ig,iv] += Float64(.00898)
        ig = ig_["CAP05202"]
        pbm.A[ig,iv] += Float64(.01958)
        iv,ix_,_ = s2mpj_ii("X00803",ix_)
        arrset(pb.xnames,iv,"X00803")
        ig = ig_["CAP03703"]
        pbm.A[ig,iv] += Float64(.00397)
        ig = ig_["CAP05003"]
        pbm.A[ig,iv] += Float64(.00678)
        iv,ix_,_ = s2mpj_ii("X00803",ix_)
        arrset(pb.xnames,iv,"X00803")
        ig = ig_["INV00603"]
        pbm.A[ig,iv] += Float64(26.0)
        ig = ig_["CAP03303"]
        pbm.A[ig,iv] += Float64(.00271)
        iv,ix_,_ = s2mpj_ii("X00803",ix_)
        arrset(pb.xnames,iv,"X00803")
        ig = ig_["CAP03302"]
        pbm.A[ig,iv] += Float64(.01518)
        ig = ig_["CAP04103"]
        pbm.A[ig,iv] += Float64(.00002)
        iv,ix_,_ = s2mpj_ii("X00803",ix_)
        arrset(pb.xnames,iv,"X00803")
        ig = ig_["CAP03602"]
        pbm.A[ig,iv] += Float64(.03149)
        ig = ig_["CAP04602"]
        pbm.A[ig,iv] += Float64(.01598)
        iv,ix_,_ = s2mpj_ii("X00803",ix_)
        arrset(pb.xnames,iv,"X00803")
        ig = ig_["CAP04402"]
        pbm.A[ig,iv] += Float64(.00529)
        ig = ig_["CAP03702"]
        pbm.A[ig,iv] += Float64(.05747)
        iv,ix_,_ = s2mpj_ii("X00803",ix_)
        arrset(pb.xnames,iv,"X00803")
        ig = ig_["CAP03301"]
        pbm.A[ig,iv] += Float64(.00033)
        ig = ig_["CAP05803"]
        pbm.A[ig,iv] += Float64(.01413)
        iv,ix_,_ = s2mpj_ii("X00803",ix_)
        arrset(pb.xnames,iv,"X00803")
        ig = ig_["CAP04102"]
        pbm.A[ig,iv] += Float64(.01196)
        ig = ig_["CAP04302"]
        pbm.A[ig,iv] += Float64(.03242)
        iv,ix_,_ = s2mpj_ii("X00803",ix_)
        arrset(pb.xnames,iv,"X00803")
        ig = ig_["CAP04802"]
        pbm.A[ig,iv] += Float64(.01307)
        ig = ig_["CAP05703"]
        pbm.A[ig,iv] += Float64(.02476)
        iv,ix_,_ = s2mpj_ii("X00803",ix_)
        arrset(pb.xnames,iv,"X00803")
        ig = ig_["CAP04303"]
        pbm.A[ig,iv] += Float64(.00489)
        iv,ix_,_ = s2mpj_ii("X00804",ix_)
        arrset(pb.xnames,iv,"X00804")
        ig = ig_["CAP05404"]
        pbm.A[ig,iv] += Float64(.01151)
        ig = ig_["CAP04103"]
        pbm.A[ig,iv] += Float64(.01294)
        iv,ix_,_ = s2mpj_ii("X00804",ix_)
        arrset(pb.xnames,iv,"X00804")
        ig = ig_["CAP03704"]
        pbm.A[ig,iv] += Float64(.00397)
        ig = ig_["CAP05204"]
        pbm.A[ig,iv] += Float64(.0323)
        iv,ix_,_ = s2mpj_ii("X00804",ix_)
        arrset(pb.xnames,iv,"X00804")
        ig = ig_["CAP05104"]
        pbm.A[ig,iv] += Float64(.00692)
        ig = ig_["CAP04603"]
        pbm.A[ig,iv] += Float64(.01695)
        iv,ix_,_ = s2mpj_ii("X00804",ix_)
        arrset(pb.xnames,iv,"X00804")
        ig = ig_["CAP05603"]
        pbm.A[ig,iv] += Float64(.03388)
        ig = ig_["CAP03703"]
        pbm.A[ig,iv] += Float64(.05828)
        iv,ix_,_ = s2mpj_ii("X00804",ix_)
        arrset(pb.xnames,iv,"X00804")
        ig = ig_["INV00604"]
        pbm.A[ig,iv] += Float64(26.0)
        ig = ig_["CAP03303"]
        pbm.A[ig,iv] += Float64(.0155)
        iv,ix_,_ = s2mpj_ii("X00804",ix_)
        arrset(pb.xnames,iv,"X00804")
        ig = ig_["CAP03603"]
        pbm.A[ig,iv] += Float64(.03149)
        ig = ig_["CAP05904"]
        pbm.A[ig,iv] += Float64(.02063)
        iv,ix_,_ = s2mpj_ii("X00804",ix_)
        arrset(pb.xnames,iv,"X00804")
        ig = ig_["CAP05403"]
        pbm.A[ig,iv] += Float64(.00679)
        ig = ig_["CAP05804"]
        pbm.A[ig,iv] += Float64(.01413)
        iv,ix_,_ = s2mpj_ii("X00804",ix_)
        arrset(pb.xnames,iv,"X00804")
        ig = ig_["CAP05704"]
        pbm.A[ig,iv] += Float64(.02476)
        ig = ig_["CAP05203"]
        pbm.A[ig,iv] += Float64(.01958)
        iv,ix_,_ = s2mpj_ii("X00804",ix_)
        arrset(pb.xnames,iv,"X00804")
        ig = ig_["CAP04803"]
        pbm.A[ig,iv] += Float64(.01307)
        ig = ig_["CAP05103"]
        pbm.A[ig,iv] += Float64(.02245)
        iv,ix_,_ = s2mpj_ii("X00804",ix_)
        arrset(pb.xnames,iv,"X00804")
        ig = ig_["CAP05604"]
        pbm.A[ig,iv] += Float64(.00898)
        ig = ig_["CAP05003"]
        pbm.A[ig,iv] += Float64(.01034)
        iv,ix_,_ = s2mpj_ii("X00804",ix_)
        arrset(pb.xnames,iv,"X00804")
        ig = ig_["CAP05703"]
        pbm.A[ig,iv] += Float64(.01309)
        ig = ig_["CAP04903"]
        pbm.A[ig,iv] += Float64(.01364)
        iv,ix_,_ = s2mpj_ii("X00804",ix_)
        arrset(pb.xnames,iv,"X00804")
        ig = ig_["CAP04104"]
        pbm.A[ig,iv] += Float64(.00002)
        ig = ig_["CAP03604"]
        pbm.A[ig,iv] += Float64(.01934)
        iv,ix_,_ = s2mpj_ii("X00804",ix_)
        arrset(pb.xnames,iv,"X00804")
        ig = ig_["CAP05004"]
        pbm.A[ig,iv] += Float64(.00678)
        ig = ig_["CAP05903"]
        pbm.A[ig,iv] += Float64(.01835)
        iv,ix_,_ = s2mpj_ii("X00804",ix_)
        arrset(pb.xnames,iv,"X00804")
        ig = ig_["CAP04404"]
        pbm.A[ig,iv] += Float64(.01289)
        ig = ig_["CAP04304"]
        pbm.A[ig,iv] += Float64(.00489)
        iv,ix_,_ = s2mpj_ii("X00804",ix_)
        arrset(pb.xnames,iv,"X00804")
        ig = ig_["CAP04904"]
        pbm.A[ig,iv] += Float64(.00453)
        ig = ig_["CAP04804"]
        pbm.A[ig,iv] += Float64(.00494)
        iv,ix_,_ = s2mpj_ii("X00804",ix_)
        arrset(pb.xnames,iv,"X00804")
        ig = ig_["CAP04403"]
        pbm.A[ig,iv] += Float64(.00529)
        ig = ig_["CAP05803"]
        pbm.A[ig,iv] += Float64(.02611)
        iv,ix_,_ = s2mpj_ii("X00804",ix_)
        arrset(pb.xnames,iv,"X00804")
        ig = ig_["CAP04604"]
        pbm.A[ig,iv] += Float64(.00037)
        ig = ig_["CAP03304"]
        pbm.A[ig,iv] += Float64(.00271)
        iv,ix_,_ = s2mpj_ii("X00804",ix_)
        arrset(pb.xnames,iv,"X00804")
        ig = ig_["CAP04303"]
        pbm.A[ig,iv] += Float64(.03242)
        iv,ix_,_ = s2mpj_ii("X00805",ix_)
        arrset(pb.xnames,iv,"X00805")
        ig = ig_["CAP04105"]
        pbm.A[ig,iv] += Float64(.0005)
        ig = ig_["CAP04805"]
        pbm.A[ig,iv] += Float64(.00542)
        iv,ix_,_ = s2mpj_ii("X00805",ix_)
        arrset(pb.xnames,iv,"X00805")
        ig = ig_["CAP05905"]
        pbm.A[ig,iv] += Float64(.02131)
        ig = ig_["CAP05604"]
        pbm.A[ig,iv] += Float64(.03262)
        iv,ix_,_ = s2mpj_ii("X00805",ix_)
        arrset(pb.xnames,iv,"X00805")
        ig = ig_["CAP04605"]
        pbm.A[ig,iv] += Float64(.00099)
        ig = ig_["CAP04104"]
        pbm.A[ig,iv] += Float64(.01246)
        iv,ix_,_ = s2mpj_ii("X00805",ix_)
        arrset(pb.xnames,iv,"X00805")
        ig = ig_["CAP05704"]
        pbm.A[ig,iv] += Float64(.0126)
        ig = ig_["CAP04604"]
        pbm.A[ig,iv] += Float64(.01632)
        iv,ix_,_ = s2mpj_ii("X00805",ix_)
        arrset(pb.xnames,iv,"X00805")
        ig = ig_["CAP05804"]
        pbm.A[ig,iv] += Float64(.02515)
        ig = ig_["CAP04804"]
        pbm.A[ig,iv] += Float64(.01258)
        iv,ix_,_ = s2mpj_ii("X00805",ix_)
        arrset(pb.xnames,iv,"X00805")
        ig = ig_["CAP04305"]
        pbm.A[ig,iv] += Float64(.00609)
        ig = ig_["INV00605"]
        pbm.A[ig,iv] += Float64(26.0)
        iv,ix_,_ = s2mpj_ii("X00805",ix_)
        arrset(pb.xnames,iv,"X00805")
        ig = ig_["CAP05904"]
        pbm.A[ig,iv] += Float64(.01767)
        ig = ig_["CAP04404"]
        pbm.A[ig,iv] += Float64(.0051)
        iv,ix_,_ = s2mpj_ii("X00805",ix_)
        arrset(pb.xnames,iv,"X00805")
        ig = ig_["CAP04405"]
        pbm.A[ig,iv] += Float64(.01309)
        ig = ig_["CAP03304"]
        pbm.A[ig,iv] += Float64(.01493)
        iv,ix_,_ = s2mpj_ii("X00805",ix_)
        arrset(pb.xnames,iv,"X00805")
        ig = ig_["CAP05105"]
        pbm.A[ig,iv] += Float64(.00776)
        ig = ig_["CAP05805"]
        pbm.A[ig,iv] += Float64(.0151)
        iv,ix_,_ = s2mpj_ii("X00805",ix_)
        arrset(pb.xnames,iv,"X00805")
        ig = ig_["CAP03604"]
        pbm.A[ig,iv] += Float64(.03033)
        ig = ig_["CAP05205"]
        pbm.A[ig,iv] += Float64(.03303)
        iv,ix_,_ = s2mpj_ii("X00805",ix_)
        arrset(pb.xnames,iv,"X00805")
        ig = ig_["CAP03605"]
        pbm.A[ig,iv] += Float64(.0205)
        ig = ig_["CAP03305"]
        pbm.A[ig,iv] += Float64(.00328)
        iv,ix_,_ = s2mpj_ii("X00805",ix_)
        arrset(pb.xnames,iv,"X00805")
        ig = ig_["CAP05005"]
        pbm.A[ig,iv] += Float64(.00717)
        ig = ig_["CAP05104"]
        pbm.A[ig,iv] += Float64(.02161)
        iv,ix_,_ = s2mpj_ii("X00805",ix_)
        arrset(pb.xnames,iv,"X00805")
        ig = ig_["CAP03704"]
        pbm.A[ig,iv] += Float64(.05612)
        ig = ig_["CAP05004"]
        pbm.A[ig,iv] += Float64(.00996)
        iv,ix_,_ = s2mpj_ii("X00805",ix_)
        arrset(pb.xnames,iv,"X00805")
        ig = ig_["CAP05204"]
        pbm.A[ig,iv] += Float64(.01886)
        ig = ig_["CAP04304"]
        pbm.A[ig,iv] += Float64(.03122)
        iv,ix_,_ = s2mpj_ii("X00805",ix_)
        arrset(pb.xnames,iv,"X00805")
        ig = ig_["CAP03705"]
        pbm.A[ig,iv] += Float64(.00613)
        ig = ig_["CAP05405"]
        pbm.A[ig,iv] += Float64(.01176)
        iv,ix_,_ = s2mpj_ii("X00805",ix_)
        arrset(pb.xnames,iv,"X00805")
        ig = ig_["CAP04904"]
        pbm.A[ig,iv] += Float64(.01313)
        ig = ig_["CAP04905"]
        pbm.A[ig,iv] += Float64(.00503)
        iv,ix_,_ = s2mpj_ii("X00805",ix_)
        arrset(pb.xnames,iv,"X00805")
        ig = ig_["CAP05605"]
        pbm.A[ig,iv] += Float64(.01024)
        ig = ig_["CAP05404"]
        pbm.A[ig,iv] += Float64(.00654)
        iv,ix_,_ = s2mpj_ii("X00805",ix_)
        arrset(pb.xnames,iv,"X00805")
        ig = ig_["CAP05705"]
        pbm.A[ig,iv] += Float64(.02525)
        iv,ix_,_ = s2mpj_ii("X00806",ix_)
        arrset(pb.xnames,iv,"X00806")
        ig = ig_["CAP05705"]
        pbm.A[ig,iv] += Float64(.01361)
        ig = ig_["CAP04805"]
        pbm.A[ig,iv] += Float64(.01359)
        iv,ix_,_ = s2mpj_ii("X00806",ix_)
        arrset(pb.xnames,iv,"X00806")
        ig = ig_["CAP05805"]
        pbm.A[ig,iv] += Float64(.02716)
        ig = ig_["CAP05605"]
        pbm.A[ig,iv] += Float64(.03523)
        iv,ix_,_ = s2mpj_ii("X00806",ix_)
        arrset(pb.xnames,iv,"X00806")
        ig = ig_["CAP05405"]
        pbm.A[ig,iv] += Float64(.00707)
        ig = ig_["CAP05205"]
        pbm.A[ig,iv] += Float64(.02037)
        iv,ix_,_ = s2mpj_ii("X00806",ix_)
        arrset(pb.xnames,iv,"X00806")
        ig = ig_["CAP05105"]
        pbm.A[ig,iv] += Float64(.02334)
        ig = ig_["CAP05005"]
        pbm.A[ig,iv] += Float64(.01076)
        iv,ix_,_ = s2mpj_ii("X00806",ix_)
        arrset(pb.xnames,iv,"X00806")
        ig = ig_["CAP04905"]
        pbm.A[ig,iv] += Float64(.01418)
        ig = ig_["CAP04605"]
        pbm.A[ig,iv] += Float64(.01763)
        iv,ix_,_ = s2mpj_ii("X00806",ix_)
        arrset(pb.xnames,iv,"X00806")
        ig = ig_["CAP05905"]
        pbm.A[ig,iv] += Float64(.01908)
        ig = ig_["CAP04405"]
        pbm.A[ig,iv] += Float64(.00551)
        iv,ix_,_ = s2mpj_ii("X00806",ix_)
        arrset(pb.xnames,iv,"X00806")
        ig = ig_["CAP04305"]
        pbm.A[ig,iv] += Float64(.03372)
        ig = ig_["CAP04105"]
        pbm.A[ig,iv] += Float64(.01346)
        iv,ix_,_ = s2mpj_ii("X00806",ix_)
        arrset(pb.xnames,iv,"X00806")
        ig = ig_["CAP03705"]
        pbm.A[ig,iv] += Float64(.06061)
        ig = ig_["CAP03605"]
        pbm.A[ig,iv] += Float64(.03275)
        iv,ix_,_ = s2mpj_ii("X00806",ix_)
        arrset(pb.xnames,iv,"X00806")
        ig = ig_["CAP03305"]
        pbm.A[ig,iv] += Float64(.01612)
        ig = ig_["CAP03306"]
        pbm.A[ig,iv] += Float64(.01821)
        iv,ix_,_ = s2mpj_ii("X00806",ix_)
        arrset(pb.xnames,iv,"X00806")
        ig = ig_["CAP04906"]
        pbm.A[ig,iv] += Float64(.01817)
        ig = ig_["CAP05806"]
        pbm.A[ig,iv] += Float64(.04025)
        iv,ix_,_ = s2mpj_ii("X00806",ix_)
        arrset(pb.xnames,iv,"X00806")
        ig = ig_["CAP05706"]
        pbm.A[ig,iv] += Float64(.03785)
        ig = ig_["CAP05606"]
        pbm.A[ig,iv] += Float64(.04286)
        iv,ix_,_ = s2mpj_ii("X00806",ix_)
        arrset(pb.xnames,iv,"X00806")
        ig = ig_["CAP05406"]
        pbm.A[ig,iv] += Float64(.01831)
        ig = ig_["CAP05206"]
        pbm.A[ig,iv] += Float64(.05189)
        iv,ix_,_ = s2mpj_ii("X00806",ix_)
        arrset(pb.xnames,iv,"X00806")
        ig = ig_["CAP05106"]
        pbm.A[ig,iv] += Float64(.02937)
        ig = ig_["CAP05006"]
        pbm.A[ig,iv] += Float64(.01713)
        iv,ix_,_ = s2mpj_ii("X00806",ix_)
        arrset(pb.xnames,iv,"X00806")
        ig = ig_["CAP04806"]
        pbm.A[ig,iv] += Float64(.018)
        ig = ig_["CAP05906"]
        pbm.A[ig,iv] += Float64(.03898)
        iv,ix_,_ = s2mpj_ii("X00806",ix_)
        arrset(pb.xnames,iv,"X00806")
        ig = ig_["CAP04606"]
        pbm.A[ig,iv] += Float64(.01732)
        ig = ig_["CAP04406"]
        pbm.A[ig,iv] += Float64(.01818)
        iv,ix_,_ = s2mpj_ii("X00806",ix_)
        arrset(pb.xnames,iv,"X00806")
        ig = ig_["CAP04306"]
        pbm.A[ig,iv] += Float64(.03731)
        ig = ig_["CAP04106"]
        pbm.A[ig,iv] += Float64(.01296)
        iv,ix_,_ = s2mpj_ii("X00806",ix_)
        arrset(pb.xnames,iv,"X00806")
        ig = ig_["CAP03706"]
        pbm.A[ig,iv] += Float64(.06226)
        ig = ig_["CAP03606"]
        pbm.A[ig,iv] += Float64(.05083)
        iv,ix_,_ = s2mpj_ii("I00101",ix_)
        arrset(pb.xnames,iv,"I00101")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(100.08)
        ig = ig_["INV00102"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("I00101",ix_)
        arrset(pb.xnames,iv,"I00101")
        ig = ig_["INV00101"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("I00102",ix_)
        arrset(pb.xnames,iv,"I00102")
        ig = ig_["INV00102"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["INV00103"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("I00102",ix_)
        arrset(pb.xnames,iv,"I00102")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(100.08)
        iv,ix_,_ = s2mpj_ii("I00103",ix_)
        arrset(pb.xnames,iv,"I00103")
        ig = ig_["INV00103"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["INV00104"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("I00103",ix_)
        arrset(pb.xnames,iv,"I00103")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(100.08)
        iv,ix_,_ = s2mpj_ii("I00104",ix_)
        arrset(pb.xnames,iv,"I00104")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(100.08)
        ig = ig_["INV00104"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("I00104",ix_)
        arrset(pb.xnames,iv,"I00104")
        ig = ig_["INV00105"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("I00105",ix_)
        arrset(pb.xnames,iv,"I00105")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(100.08)
        ig = ig_["INV00106"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("I00105",ix_)
        arrset(pb.xnames,iv,"I00105")
        ig = ig_["INV00105"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("I00106",ix_)
        arrset(pb.xnames,iv,"I00106")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(100.08)
        ig = ig_["INV00106"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("I00201",ix_)
        arrset(pb.xnames,iv,"I00201")
        ig = ig_["INV00201"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["INV00202"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("I00201",ix_)
        arrset(pb.xnames,iv,"I00201")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(100.08)
        iv,ix_,_ = s2mpj_ii("I00202",ix_)
        arrset(pb.xnames,iv,"I00202")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(100.08)
        ig = ig_["INV00203"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("I00202",ix_)
        arrset(pb.xnames,iv,"I00202")
        ig = ig_["INV00202"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("I00203",ix_)
        arrset(pb.xnames,iv,"I00203")
        ig = ig_["INV00203"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["INV00204"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("I00203",ix_)
        arrset(pb.xnames,iv,"I00203")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(100.08)
        iv,ix_,_ = s2mpj_ii("I00204",ix_)
        arrset(pb.xnames,iv,"I00204")
        ig = ig_["INV00205"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["INV00204"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("I00204",ix_)
        arrset(pb.xnames,iv,"I00204")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(100.08)
        iv,ix_,_ = s2mpj_ii("I00205",ix_)
        arrset(pb.xnames,iv,"I00205")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(100.08)
        ig = ig_["INV00205"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("I00205",ix_)
        arrset(pb.xnames,iv,"I00205")
        ig = ig_["INV00206"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("I00206",ix_)
        arrset(pb.xnames,iv,"I00206")
        ig = ig_["INV00206"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(100.08)
        iv,ix_,_ = s2mpj_ii("I00301",ix_)
        arrset(pb.xnames,iv,"I00301")
        ig = ig_["INV00302"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(100.08)
        iv,ix_,_ = s2mpj_ii("I00301",ix_)
        arrset(pb.xnames,iv,"I00301")
        ig = ig_["INV00301"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("I00302",ix_)
        arrset(pb.xnames,iv,"I00302")
        ig = ig_["INV00303"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["INV00302"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("I00302",ix_)
        arrset(pb.xnames,iv,"I00302")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(100.08)
        iv,ix_,_ = s2mpj_ii("I00303",ix_)
        arrset(pb.xnames,iv,"I00303")
        ig = ig_["INV00304"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["INV00303"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("I00303",ix_)
        arrset(pb.xnames,iv,"I00303")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(100.08)
        iv,ix_,_ = s2mpj_ii("I00304",ix_)
        arrset(pb.xnames,iv,"I00304")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(100.08)
        ig = ig_["INV00305"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("I00304",ix_)
        arrset(pb.xnames,iv,"I00304")
        ig = ig_["INV00304"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("I00305",ix_)
        arrset(pb.xnames,iv,"I00305")
        ig = ig_["INV00306"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(100.08)
        iv,ix_,_ = s2mpj_ii("I00305",ix_)
        arrset(pb.xnames,iv,"I00305")
        ig = ig_["INV00305"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("I00306",ix_)
        arrset(pb.xnames,iv,"I00306")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(100.08)
        ig = ig_["INV00306"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("I00401",ix_)
        arrset(pb.xnames,iv,"I00401")
        ig = ig_["INV00401"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["INV00402"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("I00401",ix_)
        arrset(pb.xnames,iv,"I00401")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(100.08)
        iv,ix_,_ = s2mpj_ii("I00402",ix_)
        arrset(pb.xnames,iv,"I00402")
        ig = ig_["INV00402"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["INV00403"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("I00402",ix_)
        arrset(pb.xnames,iv,"I00402")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(100.08)
        iv,ix_,_ = s2mpj_ii("I00403",ix_)
        arrset(pb.xnames,iv,"I00403")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(100.08)
        ig = ig_["INV00404"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("I00403",ix_)
        arrset(pb.xnames,iv,"I00403")
        ig = ig_["INV00403"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("I00404",ix_)
        arrset(pb.xnames,iv,"I00404")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(100.08)
        ig = ig_["INV00404"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("I00404",ix_)
        arrset(pb.xnames,iv,"I00404")
        ig = ig_["INV00405"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("I00405",ix_)
        arrset(pb.xnames,iv,"I00405")
        ig = ig_["INV00406"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["INV00405"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("I00405",ix_)
        arrset(pb.xnames,iv,"I00405")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(100.08)
        iv,ix_,_ = s2mpj_ii("I00406",ix_)
        arrset(pb.xnames,iv,"I00406")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(100.08)
        ig = ig_["INV00406"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("I00501",ix_)
        arrset(pb.xnames,iv,"I00501")
        ig = ig_["INV00502"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(100.08)
        iv,ix_,_ = s2mpj_ii("I00501",ix_)
        arrset(pb.xnames,iv,"I00501")
        ig = ig_["INV00501"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("I00502",ix_)
        arrset(pb.xnames,iv,"I00502")
        ig = ig_["INV00503"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(100.08)
        iv,ix_,_ = s2mpj_ii("I00502",ix_)
        arrset(pb.xnames,iv,"I00502")
        ig = ig_["INV00502"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("I00503",ix_)
        arrset(pb.xnames,iv,"I00503")
        ig = ig_["INV00503"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["INV00504"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("I00503",ix_)
        arrset(pb.xnames,iv,"I00503")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(100.08)
        iv,ix_,_ = s2mpj_ii("I00504",ix_)
        arrset(pb.xnames,iv,"I00504")
        ig = ig_["INV00504"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(100.08)
        iv,ix_,_ = s2mpj_ii("I00504",ix_)
        arrset(pb.xnames,iv,"I00504")
        ig = ig_["INV00505"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("I00505",ix_)
        arrset(pb.xnames,iv,"I00505")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(100.08)
        ig = ig_["INV00505"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("I00505",ix_)
        arrset(pb.xnames,iv,"I00505")
        ig = ig_["INV00506"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("I00506",ix_)
        arrset(pb.xnames,iv,"I00506")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(100.08)
        ig = ig_["INV00506"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("I00601",ix_)
        arrset(pb.xnames,iv,"I00601")
        ig = ig_["INV00602"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(100.08)
        iv,ix_,_ = s2mpj_ii("I00601",ix_)
        arrset(pb.xnames,iv,"I00601")
        ig = ig_["INV00601"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("I00602",ix_)
        arrset(pb.xnames,iv,"I00602")
        ig = ig_["INV00602"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["INV00603"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("I00602",ix_)
        arrset(pb.xnames,iv,"I00602")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(100.08)
        iv,ix_,_ = s2mpj_ii("I00603",ix_)
        arrset(pb.xnames,iv,"I00603")
        ig = ig_["INV00604"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["INV00603"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("I00603",ix_)
        arrset(pb.xnames,iv,"I00603")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(100.08)
        iv,ix_,_ = s2mpj_ii("I00604",ix_)
        arrset(pb.xnames,iv,"I00604")
        ig = ig_["INV00604"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["INV00605"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("I00604",ix_)
        arrset(pb.xnames,iv,"I00604")
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(100.08)
        iv,ix_,_ = s2mpj_ii("I00605",ix_)
        arrset(pb.xnames,iv,"I00605")
        ig = ig_["INV00605"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(100.08)
        iv,ix_,_ = s2mpj_ii("I00605",ix_)
        arrset(pb.xnames,iv,"I00605")
        ig = ig_["INV00606"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("I00606",ix_)
        arrset(pb.xnames,iv,"I00606")
        ig = ig_["INV00606"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["OBJECTIV"]
        pbm.A[ig,iv] += Float64(100.08)
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
        pbm.gconst[ig_["CAP00101"]] = Float64(23995.8)
        pbm.gconst[ig_["CAP00201"]] = Float64(2133.0)
        pbm.gconst[ig_["CAP00301"]] = Float64(6398.9)
        pbm.gconst[ig_["CAP00401"]] = Float64(7998.6)
        pbm.gconst[ig_["CAP00501"]] = Float64(29328.1)
        pbm.gconst[ig_["CAP00601"]] = Float64(7465.3)
        pbm.gconst[ig_["CAP00701"]] = Float64(11731.3)
        pbm.gconst[ig_["CAP00801"]] = Float64(1599.7)
        pbm.gconst[ig_["CAP00901"]] = Float64(40824.0)
        pbm.gconst[ig_["CAP01001"]] = Float64(27216.0)
        pbm.gconst[ig_["CAP01101"]] = Float64(6237.0)
        pbm.gconst[ig_["CAP01201"]] = Float64(2835.0)
        pbm.gconst[ig_["CAP01301"]] = Float64(1701.0)
        pbm.gconst[ig_["CAP01401"]] = Float64(6479.7)
        pbm.gconst[ig_["CAP01501"]] = Float64(540.0)
        pbm.gconst[ig_["CAP01601"]] = Float64(540.0)
        pbm.gconst[ig_["CAP01701"]] = Float64(540.0)
        pbm.gconst[ig_["CAP01801"]] = Float64(1080.0)
        pbm.gconst[ig_["CAP02001"]] = Float64(993.6)
        pbm.gconst[ig_["CAP02201"]] = Float64(441.6)
        pbm.gconst[ig_["CAP02501"]] = Float64(552.0)
        pbm.gconst[ig_["CAP02701"]] = Float64(2550.2)
        pbm.gconst[ig_["CAP03101"]] = Float64(728.6)
        pbm.gconst[ig_["CAP03301"]] = Float64(1296.0)
        pbm.gconst[ig_["CAP03701"]] = Float64(3013.2)
        pbm.gconst[ig_["CAP04101"]] = Float64(751.7)
        pbm.gconst[ig_["CAP04301"]] = Float64(1224.7)
        pbm.gconst[ig_["CAP04601"]] = Float64(1503.4)
        pbm.gconst[ig_["CAP05501"]] = Float64(1769.0)
        pbm.gconst[ig_["CAP06101"]] = Float64(3864.0)
        pbm.gconst[ig_["CAP06201"]] = Float64(1449.0)
        pbm.gconst[ig_["CAP06301"]] = Float64(3732.7)
        pbm.gconst[ig_["CAP06501"]] = Float64(2649.6)
        pbm.gconst[ig_["CAP00102"]] = Float64(21329.6)
        pbm.gconst[ig_["CAP00202"]] = Float64(1896.0)
        pbm.gconst[ig_["CAP00302"]] = Float64(5687.9)
        pbm.gconst[ig_["CAP00402"]] = Float64(7109.9)
        pbm.gconst[ig_["CAP00502"]] = Float64(26069.5)
        pbm.gconst[ig_["CAP00602"]] = Float64(6635.9)
        pbm.gconst[ig_["CAP00702"]] = Float64(10427.8)
        pbm.gconst[ig_["CAP00802"]] = Float64(1422.0)
        pbm.gconst[ig_["CAP00902"]] = Float64(36288.0)
        pbm.gconst[ig_["CAP01002"]] = Float64(24192.0)
        pbm.gconst[ig_["CAP01102"]] = Float64(5544.0)
        pbm.gconst[ig_["CAP01202"]] = Float64(2520.0)
        pbm.gconst[ig_["CAP01302"]] = Float64(1512.0)
        pbm.gconst[ig_["CAP01402"]] = Float64(5759.8)
        pbm.gconst[ig_["CAP01502"]] = Float64(480.0)
        pbm.gconst[ig_["CAP01602"]] = Float64(480.0)
        pbm.gconst[ig_["CAP01702"]] = Float64(480.0)
        pbm.gconst[ig_["CAP01802"]] = Float64(960.0)
        pbm.gconst[ig_["CAP01902"]] = Float64(652.8)
        pbm.gconst[ig_["CAP02002"]] = Float64(864.0)
        pbm.gconst[ig_["CAP02102"]] = Float64(859.7)
        pbm.gconst[ig_["CAP02202"]] = Float64(384.0)
        pbm.gconst[ig_["CAP02302"]] = Float64(573.1)
        pbm.gconst[ig_["CAP02402"]] = Float64(576.0)
        pbm.gconst[ig_["CAP02502"]] = Float64(480.0)
        pbm.gconst[ig_["CAP02602"]] = Float64(1488.0)
        pbm.gconst[ig_["CAP02702"]] = Float64(2217.6)
        pbm.gconst[ig_["CAP02802"]] = Float64(1392.0)
        pbm.gconst[ig_["CAP02902"]] = Float64(432.0)
        pbm.gconst[ig_["CAP03002"]] = Float64(480.0)
        pbm.gconst[ig_["CAP03102"]] = Float64(633.6)
        pbm.gconst[ig_["CAP03202"]] = Float64(1219.2)
        pbm.gconst[ig_["CAP03302"]] = Float64(1152.0)
        pbm.gconst[ig_["CAP03402"]] = Float64(2102.4)
        pbm.gconst[ig_["CAP03502"]] = Float64(115.2)
        pbm.gconst[ig_["CAP03602"]] = Float64(1209.6)
        pbm.gconst[ig_["CAP03702"]] = Float64(2678.4)
        pbm.gconst[ig_["CAP03802"]] = Float64(1157.8)
        pbm.gconst[ig_["CAP03902"]] = Float64(299.5)
        pbm.gconst[ig_["CAP04002"]] = Float64(328.3)
        pbm.gconst[ig_["CAP04102"]] = Float64(668.2)
        pbm.gconst[ig_["CAP04202"]] = Float64(1036.8)
        pbm.gconst[ig_["CAP04302"]] = Float64(1088.6)
        pbm.gconst[ig_["CAP04402"]] = Float64(362.9)
        pbm.gconst[ig_["CAP04502"]] = Float64(1002.2)
        pbm.gconst[ig_["CAP04602"]] = Float64(1336.3)
        pbm.gconst[ig_["CAP04702"]] = Float64(1336.3)
        pbm.gconst[ig_["CAP04802"]] = Float64(668.2)
        pbm.gconst[ig_["CAP04902"]] = Float64(668.2)
        pbm.gconst[ig_["CAP05002"]] = Float64(668.2)
        pbm.gconst[ig_["CAP05102"]] = Float64(679.7)
        pbm.gconst[ig_["CAP05202"]] = Float64(1670.4)
        pbm.gconst[ig_["CAP05302"]] = Float64(576.0)
        pbm.gconst[ig_["CAP05402"]] = Float64(725.8)
        pbm.gconst[ig_["CAP05502"]] = Float64(1572.5)
        pbm.gconst[ig_["CAP05602"]] = Float64(985.0)
        pbm.gconst[ig_["CAP05702"]] = Float64(898.6)
        pbm.gconst[ig_["CAP05802"]] = Float64(985.0)
        pbm.gconst[ig_["CAP05902"]] = Float64(673.9)
        pbm.gconst[ig_["CAP06002"]] = Float64(3960.0)
        pbm.gconst[ig_["CAP06102"]] = Float64(3360.0)
        pbm.gconst[ig_["CAP06202"]] = Float64(1260.0)
        pbm.gconst[ig_["CAP06302"]] = Float64(3317.9)
        pbm.gconst[ig_["CAP06402"]] = Float64(1728.0)
        pbm.gconst[ig_["CAP06502"]] = Float64(2304.0)
        pbm.gconst[ig_["CAP00103"]] = Float64(23107.0)
        pbm.gconst[ig_["CAP00203"]] = Float64(2054.0)
        pbm.gconst[ig_["CAP00303"]] = Float64(6161.9)
        pbm.gconst[ig_["CAP00403"]] = Float64(7702.3)
        pbm.gconst[ig_["CAP00503"]] = Float64(28241.9)
        pbm.gconst[ig_["CAP00603"]] = Float64(7188.9)
        pbm.gconst[ig_["CAP00703"]] = Float64(11296.8)
        pbm.gconst[ig_["CAP00803"]] = Float64(1540.5)
        pbm.gconst[ig_["CAP00903"]] = Float64(39312.0)
        pbm.gconst[ig_["CAP01003"]] = Float64(26208.0)
        pbm.gconst[ig_["CAP01103"]] = Float64(6006.0)
        pbm.gconst[ig_["CAP01203"]] = Float64(2730.0)
        pbm.gconst[ig_["CAP01303"]] = Float64(1638.0)
        pbm.gconst[ig_["CAP01403"]] = Float64(6239.7)
        pbm.gconst[ig_["CAP01503"]] = Float64(520.0)
        pbm.gconst[ig_["CAP01603"]] = Float64(520.0)
        pbm.gconst[ig_["CAP01703"]] = Float64(520.0)
        pbm.gconst[ig_["CAP01803"]] = Float64(1040.0)
        pbm.gconst[ig_["CAP01903"]] = Float64(685.4)
        pbm.gconst[ig_["CAP02003"]] = Float64(907.2)
        pbm.gconst[ig_["CAP02103"]] = Float64(902.7)
        pbm.gconst[ig_["CAP02203"]] = Float64(403.2)
        pbm.gconst[ig_["CAP02303"]] = Float64(601.8)
        pbm.gconst[ig_["CAP02403"]] = Float64(604.8)
        pbm.gconst[ig_["CAP02503"]] = Float64(504.0)
        pbm.gconst[ig_["CAP02603"]] = Float64(1562.4)
        pbm.gconst[ig_["CAP02703"]] = Float64(2328.5)
        pbm.gconst[ig_["CAP02803"]] = Float64(1461.6)
        pbm.gconst[ig_["CAP02903"]] = Float64(453.6)
        pbm.gconst[ig_["CAP03003"]] = Float64(504.0)
        pbm.gconst[ig_["CAP03103"]] = Float64(665.3)
        pbm.gconst[ig_["CAP03203"]] = Float64(1280.2)
        pbm.gconst[ig_["CAP03303"]] = Float64(1248.0)
        pbm.gconst[ig_["CAP03403"]] = Float64(2277.6)
        pbm.gconst[ig_["CAP03503"]] = Float64(124.8)
        pbm.gconst[ig_["CAP03603"]] = Float64(1310.4)
        pbm.gconst[ig_["CAP03703"]] = Float64(2901.6)
        pbm.gconst[ig_["CAP03803"]] = Float64(1254.2)
        pbm.gconst[ig_["CAP03903"]] = Float64(324.5)
        pbm.gconst[ig_["CAP04003"]] = Float64(355.7)
        pbm.gconst[ig_["CAP04103"]] = Float64(723.8)
        pbm.gconst[ig_["CAP04203"]] = Float64(1123.2)
        pbm.gconst[ig_["CAP04303"]] = Float64(1179.4)
        pbm.gconst[ig_["CAP04403"]] = Float64(393.1)
        pbm.gconst[ig_["CAP04503"]] = Float64(1085.8)
        pbm.gconst[ig_["CAP04603"]] = Float64(1447.7)
        pbm.gconst[ig_["CAP04703"]] = Float64(1447.7)
        pbm.gconst[ig_["CAP04803"]] = Float64(723.8)
        pbm.gconst[ig_["CAP04903"]] = Float64(723.8)
        pbm.gconst[ig_["CAP05003"]] = Float64(723.8)
        pbm.gconst[ig_["CAP05103"]] = Float64(736.3)
        pbm.gconst[ig_["CAP05203"]] = Float64(1809.6)
        pbm.gconst[ig_["CAP05303"]] = Float64(624.0)
        pbm.gconst[ig_["CAP05403"]] = Float64(786.2)
        pbm.gconst[ig_["CAP05503"]] = Float64(1703.5)
        pbm.gconst[ig_["CAP05603"]] = Float64(1067.0)
        pbm.gconst[ig_["CAP05703"]] = Float64(973.4)
        pbm.gconst[ig_["CAP05803"]] = Float64(1067.0)
        pbm.gconst[ig_["CAP05903"]] = Float64(730.1)
        pbm.gconst[ig_["CAP06003"]] = Float64(4158.0)
        pbm.gconst[ig_["CAP06103"]] = Float64(3528.0)
        pbm.gconst[ig_["CAP06203"]] = Float64(1323.0)
        pbm.gconst[ig_["CAP06303"]] = Float64(3594.4)
        pbm.gconst[ig_["CAP06403"]] = Float64(1814.4)
        pbm.gconst[ig_["CAP06503"]] = Float64(2419.2)
        pbm.gconst[ig_["CAP00104"]] = Float64(23107.0)
        pbm.gconst[ig_["CAP00204"]] = Float64(2054.0)
        pbm.gconst[ig_["CAP00304"]] = Float64(6161.9)
        pbm.gconst[ig_["CAP00404"]] = Float64(7702.3)
        pbm.gconst[ig_["CAP00504"]] = Float64(28241.9)
        pbm.gconst[ig_["CAP00604"]] = Float64(7188.9)
        pbm.gconst[ig_["CAP00704"]] = Float64(11296.8)
        pbm.gconst[ig_["CAP00804"]] = Float64(1540.5)
        pbm.gconst[ig_["CAP00904"]] = Float64(39312.0)
        pbm.gconst[ig_["CAP01004"]] = Float64(26208.0)
        pbm.gconst[ig_["CAP01104"]] = Float64(6006.0)
        pbm.gconst[ig_["CAP01204"]] = Float64(2730.0)
        pbm.gconst[ig_["CAP01304"]] = Float64(1638.0)
        pbm.gconst[ig_["CAP01404"]] = Float64(6239.7)
        pbm.gconst[ig_["CAP01504"]] = Float64(520.0)
        pbm.gconst[ig_["CAP01604"]] = Float64(520.0)
        pbm.gconst[ig_["CAP01704"]] = Float64(520.0)
        pbm.gconst[ig_["CAP01804"]] = Float64(1040.0)
        pbm.gconst[ig_["CAP01904"]] = Float64(718.1)
        pbm.gconst[ig_["CAP02004"]] = Float64(950.4)
        pbm.gconst[ig_["CAP02104"]] = Float64(945.6)
        pbm.gconst[ig_["CAP02204"]] = Float64(422.4)
        pbm.gconst[ig_["CAP02304"]] = Float64(630.4)
        pbm.gconst[ig_["CAP02404"]] = Float64(633.6)
        pbm.gconst[ig_["CAP02504"]] = Float64(528.0)
        pbm.gconst[ig_["CAP02604"]] = Float64(1636.8)
        pbm.gconst[ig_["CAP02704"]] = Float64(2439.4)
        pbm.gconst[ig_["CAP02804"]] = Float64(1531.2)
        pbm.gconst[ig_["CAP02904"]] = Float64(475.2)
        pbm.gconst[ig_["CAP03004"]] = Float64(528.0)
        pbm.gconst[ig_["CAP03104"]] = Float64(697.0)
        pbm.gconst[ig_["CAP03204"]] = Float64(1341.1)
        pbm.gconst[ig_["CAP03304"]] = Float64(1248.0)
        pbm.gconst[ig_["CAP03404"]] = Float64(2277.6)
        pbm.gconst[ig_["CAP03504"]] = Float64(124.8)
        pbm.gconst[ig_["CAP03604"]] = Float64(1310.4)
        pbm.gconst[ig_["CAP03704"]] = Float64(2901.6)
        pbm.gconst[ig_["CAP03804"]] = Float64(1254.2)
        pbm.gconst[ig_["CAP03904"]] = Float64(324.5)
        pbm.gconst[ig_["CAP04004"]] = Float64(355.7)
        pbm.gconst[ig_["CAP04104"]] = Float64(723.8)
        pbm.gconst[ig_["CAP04204"]] = Float64(1123.2)
        pbm.gconst[ig_["CAP04304"]] = Float64(1179.4)
        pbm.gconst[ig_["CAP04404"]] = Float64(393.1)
        pbm.gconst[ig_["CAP04504"]] = Float64(1085.8)
        pbm.gconst[ig_["CAP04604"]] = Float64(1447.7)
        pbm.gconst[ig_["CAP04704"]] = Float64(1447.7)
        pbm.gconst[ig_["CAP04804"]] = Float64(723.8)
        pbm.gconst[ig_["CAP04904"]] = Float64(723.8)
        pbm.gconst[ig_["CAP05004"]] = Float64(723.8)
        pbm.gconst[ig_["CAP05104"]] = Float64(736.3)
        pbm.gconst[ig_["CAP05204"]] = Float64(1809.6)
        pbm.gconst[ig_["CAP05304"]] = Float64(624.0)
        pbm.gconst[ig_["CAP05404"]] = Float64(786.2)
        pbm.gconst[ig_["CAP05504"]] = Float64(1703.5)
        pbm.gconst[ig_["CAP05604"]] = Float64(1067.0)
        pbm.gconst[ig_["CAP05704"]] = Float64(973.4)
        pbm.gconst[ig_["CAP05804"]] = Float64(1067.0)
        pbm.gconst[ig_["CAP05904"]] = Float64(730.1)
        pbm.gconst[ig_["CAP06004"]] = Float64(4356.0)
        pbm.gconst[ig_["CAP06104"]] = Float64(3696.0)
        pbm.gconst[ig_["CAP06204"]] = Float64(1386.0)
        pbm.gconst[ig_["CAP06304"]] = Float64(3594.4)
        pbm.gconst[ig_["CAP06404"]] = Float64(1900.8)
        pbm.gconst[ig_["CAP06504"]] = Float64(2534.4)
        pbm.gconst[ig_["CAP00105"]] = Float64(23995.8)
        pbm.gconst[ig_["CAP00205"]] = Float64(2133.0)
        pbm.gconst[ig_["CAP00305"]] = Float64(6398.9)
        pbm.gconst[ig_["CAP00405"]] = Float64(7998.6)
        pbm.gconst[ig_["CAP00505"]] = Float64(29328.1)
        pbm.gconst[ig_["CAP00605"]] = Float64(7465.3)
        pbm.gconst[ig_["CAP00705"]] = Float64(11731.3)
        pbm.gconst[ig_["CAP00805"]] = Float64(1599.7)
        pbm.gconst[ig_["CAP00905"]] = Float64(40824.0)
        pbm.gconst[ig_["CAP01005"]] = Float64(27216.0)
        pbm.gconst[ig_["CAP01105"]] = Float64(6237.0)
        pbm.gconst[ig_["CAP01205"]] = Float64(2835.0)
        pbm.gconst[ig_["CAP01305"]] = Float64(1701.0)
        pbm.gconst[ig_["CAP01405"]] = Float64(6479.7)
        pbm.gconst[ig_["CAP01505"]] = Float64(540.0)
        pbm.gconst[ig_["CAP01605"]] = Float64(540.0)
        pbm.gconst[ig_["CAP01705"]] = Float64(540.0)
        pbm.gconst[ig_["CAP01805"]] = Float64(1080.0)
        pbm.gconst[ig_["CAP01905"]] = Float64(718.1)
        pbm.gconst[ig_["CAP02005"]] = Float64(950.4)
        pbm.gconst[ig_["CAP02105"]] = Float64(945.6)
        pbm.gconst[ig_["CAP02205"]] = Float64(422.4)
        pbm.gconst[ig_["CAP02305"]] = Float64(630.4)
        pbm.gconst[ig_["CAP02405"]] = Float64(633.6)
        pbm.gconst[ig_["CAP02505"]] = Float64(528.0)
        pbm.gconst[ig_["CAP02605"]] = Float64(1636.8)
        pbm.gconst[ig_["CAP02705"]] = Float64(2439.4)
        pbm.gconst[ig_["CAP02805"]] = Float64(1531.2)
        pbm.gconst[ig_["CAP02905"]] = Float64(475.2)
        pbm.gconst[ig_["CAP03005"]] = Float64(528.0)
        pbm.gconst[ig_["CAP03105"]] = Float64(697.0)
        pbm.gconst[ig_["CAP03205"]] = Float64(1341.1)
        pbm.gconst[ig_["CAP03305"]] = Float64(1296.0)
        pbm.gconst[ig_["CAP03405"]] = Float64(2365.2)
        pbm.gconst[ig_["CAP03505"]] = Float64(129.6)
        pbm.gconst[ig_["CAP03605"]] = Float64(1360.8)
        pbm.gconst[ig_["CAP03705"]] = Float64(3013.2)
        pbm.gconst[ig_["CAP03805"]] = Float64(1302.5)
        pbm.gconst[ig_["CAP03905"]] = Float64(337.0)
        pbm.gconst[ig_["CAP04005"]] = Float64(369.4)
        pbm.gconst[ig_["CAP04105"]] = Float64(751.7)
        pbm.gconst[ig_["CAP04205"]] = Float64(1166.4)
        pbm.gconst[ig_["CAP04305"]] = Float64(1224.7)
        pbm.gconst[ig_["CAP04405"]] = Float64(408.2)
        pbm.gconst[ig_["CAP04505"]] = Float64(1127.5)
        pbm.gconst[ig_["CAP04605"]] = Float64(1503.4)
        pbm.gconst[ig_["CAP04705"]] = Float64(1503.4)
        pbm.gconst[ig_["CAP04805"]] = Float64(751.7)
        pbm.gconst[ig_["CAP04905"]] = Float64(751.7)
        pbm.gconst[ig_["CAP05005"]] = Float64(751.7)
        pbm.gconst[ig_["CAP05105"]] = Float64(764.6)
        pbm.gconst[ig_["CAP05205"]] = Float64(1879.2)
        pbm.gconst[ig_["CAP05305"]] = Float64(648.0)
        pbm.gconst[ig_["CAP05405"]] = Float64(816.5)
        pbm.gconst[ig_["CAP05505"]] = Float64(1769.0)
        pbm.gconst[ig_["CAP05605"]] = Float64(1108.1)
        pbm.gconst[ig_["CAP05705"]] = Float64(1010.9)
        pbm.gconst[ig_["CAP05805"]] = Float64(1108.1)
        pbm.gconst[ig_["CAP05905"]] = Float64(758.2)
        pbm.gconst[ig_["CAP06005"]] = Float64(4356.0)
        pbm.gconst[ig_["CAP06105"]] = Float64(3696.0)
        pbm.gconst[ig_["CAP06205"]] = Float64(1386.0)
        pbm.gconst[ig_["CAP06305"]] = Float64(3732.7)
        pbm.gconst[ig_["CAP06405"]] = Float64(1900.8)
        pbm.gconst[ig_["CAP06505"]] = Float64(2534.4)
        pbm.gconst[ig_["CAP00106"]] = Float64(22218.3)
        pbm.gconst[ig_["CAP00206"]] = Float64(1975.0)
        pbm.gconst[ig_["CAP00306"]] = Float64(5924.9)
        pbm.gconst[ig_["CAP00406"]] = Float64(7406.1)
        pbm.gconst[ig_["CAP00506"]] = Float64(27155.7)
        pbm.gconst[ig_["CAP00606"]] = Float64(6912.4)
        pbm.gconst[ig_["CAP00706"]] = Float64(10862.3)
        pbm.gconst[ig_["CAP00806"]] = Float64(1481.2)
        pbm.gconst[ig_["CAP00906"]] = Float64(37800.0)
        pbm.gconst[ig_["CAP01006"]] = Float64(25200.0)
        pbm.gconst[ig_["CAP01106"]] = Float64(5775.0)
        pbm.gconst[ig_["CAP01206"]] = Float64(2625.0)
        pbm.gconst[ig_["CAP01306"]] = Float64(1575.0)
        pbm.gconst[ig_["CAP01406"]] = Float64(1999.9)
        pbm.gconst[ig_["CAP01506"]] = Float64(166.7)
        pbm.gconst[ig_["CAP01606"]] = Float64(166.7)
        pbm.gconst[ig_["CAP01706"]] = Float64(166.7)
        pbm.gconst[ig_["CAP01806"]] = Float64(333.3)
        pbm.gconst[ig_["CAP01906"]] = Float64(685.4)
        pbm.gconst[ig_["CAP02006"]] = Float64(907.2)
        pbm.gconst[ig_["CAP02106"]] = Float64(902.7)
        pbm.gconst[ig_["CAP02206"]] = Float64(403.2)
        pbm.gconst[ig_["CAP02306"]] = Float64(601.8)
        pbm.gconst[ig_["CAP02406"]] = Float64(604.8)
        pbm.gconst[ig_["CAP02506"]] = Float64(504.0)
        pbm.gconst[ig_["CAP02606"]] = Float64(1562.4)
        pbm.gconst[ig_["CAP02706"]] = Float64(2328.5)
        pbm.gconst[ig_["CAP02806"]] = Float64(1461.6)
        pbm.gconst[ig_["CAP02906"]] = Float64(453.6)
        pbm.gconst[ig_["CAP03006"]] = Float64(504.0)
        pbm.gconst[ig_["CAP03106"]] = Float64(665.3)
        pbm.gconst[ig_["CAP03206"]] = Float64(1280.2)
        pbm.gconst[ig_["CAP03306"]] = Float64(1200.0)
        pbm.gconst[ig_["CAP03406"]] = Float64(2190.0)
        pbm.gconst[ig_["CAP03506"]] = Float64(120.0)
        pbm.gconst[ig_["CAP03606"]] = Float64(1260.0)
        pbm.gconst[ig_["CAP03706"]] = Float64(2790.0)
        pbm.gconst[ig_["CAP03806"]] = Float64(1206.0)
        pbm.gconst[ig_["CAP03906"]] = Float64(312.0)
        pbm.gconst[ig_["CAP04006"]] = Float64(342.0)
        pbm.gconst[ig_["CAP04106"]] = Float64(696.0)
        pbm.gconst[ig_["CAP04206"]] = Float64(1080.0)
        pbm.gconst[ig_["CAP04306"]] = Float64(1134.0)
        pbm.gconst[ig_["CAP04406"]] = Float64(378.0)
        pbm.gconst[ig_["CAP04506"]] = Float64(1044.0)
        pbm.gconst[ig_["CAP04606"]] = Float64(1392.0)
        pbm.gconst[ig_["CAP04706"]] = Float64(1392.0)
        pbm.gconst[ig_["CAP04806"]] = Float64(696.0)
        pbm.gconst[ig_["CAP04906"]] = Float64(696.0)
        pbm.gconst[ig_["CAP05006"]] = Float64(696.0)
        pbm.gconst[ig_["CAP05106"]] = Float64(708.0)
        pbm.gconst[ig_["CAP05206"]] = Float64(1740.0)
        pbm.gconst[ig_["CAP05306"]] = Float64(600.0)
        pbm.gconst[ig_["CAP05406"]] = Float64(756.0)
        pbm.gconst[ig_["CAP05506"]] = Float64(1638.0)
        pbm.gconst[ig_["CAP05606"]] = Float64(1026.0)
        pbm.gconst[ig_["CAP05706"]] = Float64(936.0)
        pbm.gconst[ig_["CAP05806"]] = Float64(1026.0)
        pbm.gconst[ig_["CAP05906"]] = Float64(702.0)
        pbm.gconst[ig_["CAP06006"]] = Float64(4158.0)
        pbm.gconst[ig_["CAP06106"]] = Float64(3528.0)
        pbm.gconst[ig_["CAP06206"]] = Float64(1323.0)
        pbm.gconst[ig_["CAP06306"]] = Float64(3456.2)
        pbm.gconst[ig_["CAP06406"]] = Float64(1814.4)
        pbm.gconst[ig_["CAP06506"]] = Float64(2419.2)
        pbm.gconst[ig_["MXD00102"]] = Float64(56200.0)
        pbm.gconst[ig_["MXD00202"]] = Float64(245714.0)
        pbm.gconst[ig_["MXD00502"]] = Float64(1201600.0)
        pbm.gconst[ig_["MXD00702"]] = Float64(78000.0)
        pbm.gconst[ig_["MXD00802"]] = Float64(180000.0)
        pbm.gconst[ig_["MXD00902"]] = Float64(1391199.0)
        pbm.gconst[ig_["MXD01002"]] = Float64(450000.0)
        pbm.gconst[ig_["MXD00103"]] = Float64(98200.0)
        pbm.gconst[ig_["MXD00203"]] = Float64(377143.0)
        pbm.gconst[ig_["MXD00303"]] = Float64(142824.0)
        pbm.gconst[ig_["MXD00403"]] = Float64(16000.0)
        pbm.gconst[ig_["MXD00503"]] = Float64(1984799.0)
        pbm.gconst[ig_["MXD00603"]] = Float64(34200.0)
        pbm.gconst[ig_["MXD00703"]] = Float64(155200.0)
        pbm.gconst[ig_["MXD00803"]] = Float64(258000.0)
        pbm.gconst[ig_["MXD00903"]] = Float64(2414999.0)
        pbm.gconst[ig_["MXD01003"]] = Float64(706000.0)
        pbm.gconst[ig_["MND00104"]] = Float64(89091.0)
        pbm.gconst[ig_["MXD00104"]] = Float64(232200.0)
        pbm.gconst[ig_["MXD00204"]] = Float64(491428.0)
        pbm.gconst[ig_["MND00304"]] = Float64(13640.0)
        pbm.gconst[ig_["MXD00304"]] = Float64(288706.0)
        pbm.gconst[ig_["MXD00404"]] = Float64(20400.0)
        pbm.gconst[ig_["MND00504"]] = Float64(433908.0)
        pbm.gconst[ig_["MXD00504"]] = Float64(2669599.0)
        pbm.gconst[ig_["MXD00604"]] = Float64(54200.0)
        pbm.gconst[ig_["MND00704"]] = Float64(99890.0)
        pbm.gconst[ig_["MXD00704"]] = Float64(214400.0)
        pbm.gconst[ig_["MND00804"]] = Float64(149985.0)
        pbm.gconst[ig_["MXD00804"]] = Float64(336000.0)
        pbm.gconst[ig_["MND00904"]] = Float64(68537.0)
        pbm.gconst[ig_["MXD00904"]] = Float64(3635799.0)
        pbm.gconst[ig_["MND01004"]] = Float64(118071.0)
        pbm.gconst[ig_["MXD01004"]] = Float64(966000.0)
        pbm.gconst[ig_["MND00105"]] = Float64(129087.0)
        pbm.gconst[ig_["MXD00105"]] = Float64(366200.0)
        pbm.gconst[ig_["MND00205"]] = Float64(150791.0)
        pbm.gconst[ig_["MXD00205"]] = Float64(492571.0)
        pbm.gconst[ig_["MND00305"]] = Float64(200788.0)
        pbm.gconst[ig_["MXD00305"]] = Float64(486353.0)
        pbm.gconst[ig_["MND00405"]] = Float64(11699.0)
        pbm.gconst[ig_["MXD00405"]] = Float64(25400.0)
        pbm.gconst[ig_["MND00505"]] = Float64(884246.0)
        pbm.gconst[ig_["MXD00505"]] = Float64(3241599.0)
        pbm.gconst[ig_["MND00605"]] = Float64(24096.0)
        pbm.gconst[ig_["MXD00605"]] = Float64(68200.0)
        pbm.gconst[ig_["MND00705"]] = Float64(125187.0)
        pbm.gconst[ig_["MXD00705"]] = Float64(273600.0)
        pbm.gconst[ig_["MND00805"]] = Float64(179982.0)
        pbm.gconst[ig_["MXD00805"]] = Float64(414000.0)
        pbm.gconst[ig_["MND00905"]] = Float64(1849407.0)
        pbm.gconst[ig_["MXD00905"]] = Float64(4798598.0)
        pbm.gconst[ig_["MND01005"]] = Float64(334408.0)
        pbm.gconst[ig_["MXD01005"]] = Float64(1226000.0)
        pbm.gconst[ig_["MND00106"]] = Float64(169083.0)
        pbm.gconst[ig_["MXD00106"]] = Float64(500200.0)
        pbm.gconst[ig_["MND00206"]] = Float64(150791.0)
        pbm.gconst[ig_["MXD00206"]] = Float64(493714.0)
        pbm.gconst[ig_["MND00306"]] = Float64(233184.0)
        pbm.gconst[ig_["MXD00306"]] = Float64(702823.0)
        pbm.gconst[ig_["MND00406"]] = Float64(11699.0)
        pbm.gconst[ig_["MXD00406"]] = Float64(30800.0)
        pbm.gconst[ig_["MND00506"]] = Float64(1301415.0)
        pbm.gconst[ig_["MXD00506"]] = Float64(3817599.0)
        pbm.gconst[ig_["MND00606"]] = Float64(30096.0)
        pbm.gconst[ig_["MXD00606"]] = Float64(80200.0)
        pbm.gconst[ig_["MND00706"]] = Float64(125187.0)
        pbm.gconst[ig_["MXD00706"]] = Float64(345200.0)
        pbm.gconst[ig_["MND00806"]] = Float64(209979.0)
        pbm.gconst[ig_["MXD00806"]] = Float64(492000.0)
        pbm.gconst[ig_["MND00906"]] = Float64(1849407.0)
        pbm.gconst[ig_["MXD00906"]] = Float64(6141396.0)
        pbm.gconst[ig_["MND01006"]] = Float64(644932.0)
        pbm.gconst[ig_["MXD01006"]] = Float64(1486000.0)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN                   ???
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
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
        pb.lincons   = collect(1:length(pbm.congrps))
        pb.pbclass = "C-LLR2-AN-163-488"
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

