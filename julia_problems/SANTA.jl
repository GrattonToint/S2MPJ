function SANTA(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : SANTA
#    --------
# 
#    The Santa problem as suggested in a Christmas competition
#    by Jens Jensen (Scientific Computing, STFC). To quote Jens,
# 
#    Santa and His Elves
# 
#    SCD Christmas programming challenge 2016
# 
#    Christmas has come to the Santa Claus Department – or rather, the SCD
#    is coming to Christmas. Santa is flying around the world, presently
#    presenting presents. Ho, ho, ho! No striking air crew on Santa’s sleigh!
#    No airport strikes on the North Pole.
# 
#    For the purpose of this exercise, the Earth is round as a perfect ball,
#    with radius precisely 6,371,000 metres. However, everything is at the
#    same longitude and latitude as the “real” Earth. So for example, the
#    Greenwich observatory is at 51°28'40"N 0°00'04"W both on the “real”
#    Earth and on Santa’s Earth. (Also ignore rotation of the Earth and
#    anything practical like that.)
# 
#    Santa sets off from the North Pole along 2°6'57.6" E bearing south
#    (obviously), and also bearing presents (obviously). Whenever Santa
#    leaves a location, he leaves an elf behind, in order to help unwrapping
#    presents; the elf does this and then flies out independently to meet up
#    with Santa at the next location - this means that Santa only needs two
#    Elves. Here’s how:
# 
#    1. Santa leaves the North Pole, setting out for location A. Elf 1 is
#    left behind (in this particular case, not to unwrap presents, but to
#    turn the lights off, and ensure the oven is off – it's elf'n'safety,
#    you know.)
# 
#    2. Santa arrives in location A and hands them their present. Now Elf 2
#    is with Santa; Elf 1 is still at the NP.
# 
#    3. Santa leaves location A, leaving behind Elf 2. Santa flies on to
#    location B; Elf 1, who remained at the North Pole, also flies to B and
#    meets Santa there; Elf 2 is left behind at A.
# 
#    4. Santa arrives at location B along with Elf 1, and hands out
#    presents. Santa then leaves location B to fly to C, leaving behind Elf 1
#    at location B. Meanwhile Elf 2, having finished helping at location A,
#    leaves location A to fly on to C, to meet Santa there.
# 
#    5. Santa arrives from B at location C; Elf 2 also arrives into C from
#    location A. Elf 1 remains at B until Santa flies onward to location D.
# 
#    6. At the last hop, Santa needs a rest and flies to 31°46'42.4" S
#    144°46'12.9" W.  The Elves also fly to this location - maps show no land
#    here but it is hidden. Either that or we got the coordinates wrong.
#    In either case Santa and elves do fly to this location.
# 
#    The following table shows the distance of Santa's hops, as well as those
#    of the elves, with the distance given in metres:
# 
#    Who     Hop  Distance travelled
#    Santa   1    5405238
#            2    623852
#            3    1005461
#            4    7470967
#            5    3632559
#            6    10206818
#            7    7967212
#            8    5896361
#            9    8337266
#            10   13019505
#            11   8690818
#            12   8971302
#    Elf1    1    4866724
#            2    6833740
#            3    13489586
#            4    9195575
#            5    9704793
#            6    12498127
#    Elf2    1    1375828
#            2    4917407
#            3    10617953
#            4    10996150
#            5    7901038
#            6    8971302
# 
#    What is Santa’s route?  What sort of presents is he carrying?
# 
#    Bonus question: did you really need to know the starting direction?
# 
#    Added by Nick: the problem has many local minimizers of the sum of squares
#    of infeasibility, but it is only the solution with zero residuals
#    that is of interest.
# 
#    SIF input: Nick Gould, Dec 2016.
# 
#    classification = "C-NOR2-AN-21-23"
# 
#    Number of stops on Santa's path (path goes from index 0 to 12)
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "SANTA"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["S"] = 12
        v_["1"] = 1
        v_["2"] = 2
        v_["3"] = 3
        v_["4"] = 4
        v_["180.0"] = 180.0
        v_["S-1"] = -1+v_["S"]
        v_["RS"] = Float64(v_["S"])
        v_["PI/4"] = atan(1.0)
        v_["PI"] = 4.0*v_["PI/4"]
        v_["PI/180"] = v_["PI"]/v_["180.0"]
        v_["PHI0"] = 90.0
        v_["LAM0"] = 0.0
        v_["PHI12"] = -31.77844444
        v_["LAM12"] = -144.77025
        v_["LAM1"] = 2.116
        v_["PHI0"] = v_["PHI0"]*v_["PI/180"]
        v_["LAM0"] = v_["LAM0"]*v_["PI/180"]
        v_["PHI12"] = v_["PHI12"]*v_["PI/180"]
        v_["LAM12"] = v_["LAM12"]*v_["PI/180"]
        v_["LAM1"] = v_["LAM1"]*v_["PI/180"]
        v_["DPHI"] = v_["PHI12"]-v_["PHI0"]
        v_["DLAM"] = v_["LAM12"]-v_["LAM0"]
        v_["DPHI/S"] = v_["DPHI"]/v_["RS"]
        v_["DLAM/S"] = v_["DLAM"]/v_["RS"]
        v_["RADIUS"] = 6371000.0
        v_["D0,1"] = 5405238.0
        v_["D0,2"] = 4866724.0
        v_["D1,2"] = 623852.0
        v_["D1,3"] = 1375828.0
        v_["D2,3"] = 1005461.0
        v_["D2,4"] = 6833740.0
        v_["D3,4"] = 7470967.0
        v_["D3,5"] = 4917407.0
        v_["D4,5"] = 3632559.0
        v_["D4,6"] = 13489586.0
        v_["D5,6"] = 10206818.0
        v_["D5,7"] = 10617953.0
        v_["D6,7"] = 7967212.0
        v_["D6,8"] = 9195575.0
        v_["D7,8"] = 5896361.0
        v_["D7,9"] = 10996150.0
        v_["D8,9"] = 8337266.0
        v_["D8,10"] = 9704793.0
        v_["D9,10"] = 13019505.0
        v_["D9,11"] = 7901038.0
        v_["D10,11"] = 8690818.0
        v_["D10,12"] = 12498127.0
        v_["D11,12"] = 8971302.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        iv,ix_,_ = s2mpj_ii("PHI1",ix_)
        arrset(pb.xnames,iv,"PHI1")
        for I = Int64(v_["2"]):Int64(v_["S-1"])
            iv,ix_,_ = s2mpj_ii("PHI"*string(I),ix_)
            arrset(pb.xnames,iv,"PHI"*string(I))
            iv,ix_,_ = s2mpj_ii("LAM"*string(I),ix_)
            arrset(pb.xnames,iv,"LAM"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("R0,1",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R0,1")
        for I = Int64(v_["2"]):Int64(v_["S"])
            v_["I1"] = -1+I
            v_["I2"] = -2+I
            ig,ig_,_ = s2mpj_ii("R"*string(Int64(v_["I2"]))*","*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"R"*string(Int64(v_["I2"]))*","*string(I))
            ig,ig_,_ = s2mpj_ii("R"*string(Int64(v_["I1"]))*","*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"R"*string(Int64(v_["I1"]))*","*string(I))
        end
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
        v_["D/RAD"] = v_["D0,1"]/v_["RADIUS"]
        v_["CD/RAD"] = cos(v_["D/RAD"])
        pbm.gconst[ig_["R0,1"]] = Float64(v_["CD/RAD"])
        for I = Int64(v_["2"]):Int64(v_["S"])
            v_["I2"] = -2+I
            v_["D/RAD"] = v_["D"*string(Int64(v_["I2"]))*","*string(I)]/v_["RADIUS"]
            v_["CD/RAD"] = cos(v_["D/RAD"])
            pbm.gconst[ig_["R"*string(Int64(v_["I2"]))*","*string(I)]]  = (
                  Float64(v_["CD/RAD"]))
            v_["I1"] = -1+I
            v_["D/RAD"] = v_["D"*string(Int64(v_["I1"]))*","*string(I)]/v_["RADIUS"]
            v_["CD/RAD"] = cos(v_["D/RAD"])
            pbm.gconst[ig_["R"*string(Int64(v_["I1"]))*","*string(I)]]  = (
                  Float64(v_["CD/RAD"]))
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = fill(-1000.0,pb.n)
        pb.xupper = fill(1000.0,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        pb.x0[ix_["PHI1"]] = Float64(0.7223835215)
        pb.x0[ix_["PHI2"]] = Float64(0.8069093428)
        pb.x0[ix_["LAM2"]] = Float64(-0.031657133)
        pb.x0[ix_["PHI3"]] = Float64(0.9310164154)
        pb.x0[ix_["LAM3"]] = Float64(0.1199353230)
        pb.x0[ix_["PHI4"]] = Float64(6.6067392710)
        pb.x0[ix_["LAM4"]] = Float64(-1.214314477)
        pb.x0[ix_["PHI5"]] = Float64(-3.530946794)
        pb.x0[ix_["LAM5"]] = Float64(2.5329493980)
        pb.x0[ix_["PHI6"]] = Float64(-9.798251905)
        pb.x0[ix_["LAM6"]] = Float64(4.3021328700)
        pb.x0[ix_["PHI7"]] = Float64(14.632267534)
        pb.x0[ix_["LAM7"]] = Float64(-12.96253311)
        pb.x0[ix_["PHI8"]] = Float64(2.0349445303)
        pb.x0[ix_["LAM8"]] = Float64(-4.050000443)
        pb.x0[ix_["PHI9"]] = Float64(-28.45607804)
        pb.x0[ix_["LAM9"]] = Float64(22.430117198)
        pb.x0[ix_["PHI10"]] = Float64(16.034035489)
        pb.x0[ix_["LAM10"]] = Float64(-17.28050167)
        pb.x0[ix_["PHI11"]] = Float64(0.8717052037)
        pb.x0[ix_["LAM11"]] = Float64(-0.833052840)
        v_["PHIS"] = v_["DPHI/S"]
        v_["START"] = v_["PHI0"]+v_["PHIS"]
        for I = Int64(v_["2"]):Int64(v_["S-1"])
            v_["RI"] = Float64(I)
            v_["PHIS"] = v_["DPHI/S"]*v_["RI"]
            v_["START"] = v_["PHI0"]+v_["PHIS"]
            v_["LAMS"] = v_["DLAM/S"]*v_["RI"]
            v_["START"] = v_["LAM0"]+v_["LAMS"]
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eE", iet_)
        loaset(elftv,it,1,"PHI1")
        loaset(elftv,it,2,"PHI2")
        loaset(elftv,it,3,"LAM1")
        loaset(elftv,it,4,"LAM2")
        it,iet_,_ = s2mpj_ii( "eE3", iet_)
        loaset(elftv,it,1,"PHI1")
        loaset(elftv,it,2,"PHI2")
        loaset(elftv,it,3,"LAM1")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"LAMF")
        it,iet_,_ = s2mpj_ii( "eE2", iet_)
        loaset(elftv,it,1,"PHI1")
        loaset(elftv,it,2,"LAM1")
        loaset(elftp,it,1,"PHIF")
        loaset(elftp,it,2,"LAMF")
        it,iet_,_ = s2mpj_ii( "eE1", iet_)
        loaset(elftv,it,1,"PHI1")
        loaset(elftp,it,1,"PHIF")
        loaset(elftp,it,2,"LAMF")
        loaset(elftp,it,3,"LAMS")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "E0,1"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eE1")
        arrset(ielftype,ie,iet_["eE1"])
        vname = "PHI1"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-1000.0),Float64(1000.0),nothing))
        posev = findfirst(x->x=="PHI1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="PHIF",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(v_["PHI0"]))
        posep = findfirst(x->x=="LAMF",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(v_["LAM0"]))
        posep = findfirst(x->x=="LAMS",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(v_["LAM1"]))
        ename = "E0,2"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eE2")
        arrset(ielftype,ie,iet_["eE2"])
        vname = "PHI2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-1000.0),Float64(1000.0),nothing))
        posev = findfirst(x->x=="PHI1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "LAM2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-1000.0),Float64(1000.0),nothing))
        posev = findfirst(x->x=="LAM1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="PHIF",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(v_["PHI0"]))
        posep = findfirst(x->x=="LAMF",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(v_["LAM0"]))
        ename = "E1,2"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eE3")
        arrset(ielftype,ie,iet_["eE3"])
        vname = "PHI1"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-1000.0),Float64(1000.0),nothing))
        posev = findfirst(x->x=="PHI1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "PHI2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-1000.0),Float64(1000.0),nothing))
        posev = findfirst(x->x=="PHI2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "LAM2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-1000.0),Float64(1000.0),nothing))
        posev = findfirst(x->x=="LAM1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="LAMF",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(v_["LAM1"]))
        ename = "E1,3"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eE3")
        arrset(ielftype,ie,iet_["eE3"])
        vname = "PHI1"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-1000.0),Float64(1000.0),nothing))
        posev = findfirst(x->x=="PHI1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "PHI3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-1000.0),Float64(1000.0),nothing))
        posev = findfirst(x->x=="PHI2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "LAM3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-1000.0),Float64(1000.0),nothing))
        posev = findfirst(x->x=="LAM1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="LAMF",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(v_["LAM1"]))
        ename = "E2,3"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eE")
        arrset(ielftype,ie,iet_["eE"])
        vname = "PHI2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-1000.0),Float64(1000.0),nothing))
        posev = findfirst(x->x=="PHI1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "PHI3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-1000.0),Float64(1000.0),nothing))
        posev = findfirst(x->x=="PHI2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "LAM2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-1000.0),Float64(1000.0),nothing))
        posev = findfirst(x->x=="LAM1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "LAM3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-1000.0),Float64(1000.0),nothing))
        posev = findfirst(x->x=="LAM2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        for I = Int64(v_["4"]):Int64(v_["S-1"])
            v_["I2"] = -2+I
            ename = "E"*string(Int64(v_["I2"]))*","*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eE")
            arrset(ielftype,ie,iet_["eE"])
            ename = "E"*string(Int64(v_["I2"]))*","*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "PHI"*string(Int64(v_["I2"]))
            iv,ix_,pb  = (
                  s2mpj_nlx(vname,ix_,pb,1,Float64(-1000.0),Float64(1000.0),nothing))
            posev = findfirst(x->x=="PHI1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "E"*string(Int64(v_["I2"]))*","*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "PHI"*string(I)
            iv,ix_,pb  = (
                  s2mpj_nlx(vname,ix_,pb,1,Float64(-1000.0),Float64(1000.0),nothing))
            posev = findfirst(x->x=="PHI2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "E"*string(Int64(v_["I2"]))*","*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "LAM"*string(Int64(v_["I2"]))
            iv,ix_,pb  = (
                  s2mpj_nlx(vname,ix_,pb,1,Float64(-1000.0),Float64(1000.0),nothing))
            posev = findfirst(x->x=="LAM1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "E"*string(Int64(v_["I2"]))*","*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "LAM"*string(I)
            iv,ix_,pb  = (
                  s2mpj_nlx(vname,ix_,pb,1,Float64(-1000.0),Float64(1000.0),nothing))
            posev = findfirst(x->x=="LAM2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            v_["I1"] = -1+I
            ename = "E"*string(Int64(v_["I1"]))*","*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eE")
            arrset(ielftype,ie,iet_["eE"])
            ename = "E"*string(Int64(v_["I1"]))*","*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "PHI"*string(Int64(v_["I1"]))
            iv,ix_,pb  = (
                  s2mpj_nlx(vname,ix_,pb,1,Float64(-1000.0),Float64(1000.0),nothing))
            posev = findfirst(x->x=="PHI1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "E"*string(Int64(v_["I1"]))*","*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "PHI"*string(I)
            iv,ix_,pb  = (
                  s2mpj_nlx(vname,ix_,pb,1,Float64(-1000.0),Float64(1000.0),nothing))
            posev = findfirst(x->x=="PHI2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "E"*string(Int64(v_["I1"]))*","*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "LAM"*string(Int64(v_["I1"]))
            iv,ix_,pb  = (
                  s2mpj_nlx(vname,ix_,pb,1,Float64(-1000.0),Float64(1000.0),nothing))
            posev = findfirst(x->x=="LAM1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "E"*string(Int64(v_["I1"]))*","*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "LAM"*string(I)
            iv,ix_,pb  = (
                  s2mpj_nlx(vname,ix_,pb,1,Float64(-1000.0),Float64(1000.0),nothing))
            posev = findfirst(x->x=="LAM2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        ename = "E10,12"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eE2")
        arrset(ielftype,ie,iet_["eE2"])
        vname = "PHI10"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-1000.0),Float64(1000.0),nothing))
        posev = findfirst(x->x=="PHI1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "LAM10"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-1000.0),Float64(1000.0),nothing))
        posev = findfirst(x->x=="LAM1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="PHIF",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(v_["PHI12"]))
        posep = findfirst(x->x=="LAMF",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(v_["LAM12"]))
        ename = "E11,12"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eE2")
        arrset(ielftype,ie,iet_["eE2"])
        vname = "PHI11"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-1000.0),Float64(1000.0),nothing))
        posev = findfirst(x->x=="PHI1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "LAM11"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-1000.0),Float64(1000.0),nothing))
        posev = findfirst(x->x=="LAM1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="PHIF",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(v_["PHI12"]))
        posep = findfirst(x->x=="LAMF",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(v_["LAM12"]))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["R0,1"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E0,1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        for I = Int64(v_["2"]):Int64(v_["S"])
            v_["I2"] = -2+I
            ig = ig_["R"*string(Int64(v_["I2"]))*","*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["I2"]))*","*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
            v_["I1"] = -1+I
            ig = ig_["R"*string(Int64(v_["I1"]))*","*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["I1"]))*","*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SANTA               0.0
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
        pb.pbclass = "C-NOR2-AN-21-23"
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

    elseif action == "eE"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        S1 = sin(EV_[1])
        S2 = sin(EV_[2])
        C1 = cos(EV_[1])
        C2 = cos(EV_[2])
        C = cos(EV_[3]-EV_[4])
        S = sin(EV_[3]-EV_[4])
        C1C2S = C1*C2*S
        C1C2C = C1*C2*C
        C1S2S = C1*S2*S
        S1C2S = S1*C2*S
        f_   = S1*S2+C1*C2*C
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = C1*S2-S1*C2*C
            g_[2] = S1*C2-C1*S2*C
            g_[3] = -C1C2S
            g_[4] = C1C2S
            if nargout>2
                H_ = zeros(Float64,4,4)
                H_[1,1] = -S1*S2-C1*C2*C
                H_[2,1] = C1*C2+S1*S2*C
                H_[1,2] = H_[2,1]
                H_[2,2] = -S1*S2-C1*C2*C
                H_[3,1] = S1C2S
                H_[1,3] = H_[3,1]
                H_[3,2] = C1S2S
                H_[2,3] = H_[3,2]
                H_[3,3] = -C1C2C
                H_[1,4] = -S1C2S
                H_[4,1] = H_[1,4]
                H_[2,4] = -C1S2S
                H_[4,2] = H_[2,4]
                H_[3,4] = C1C2C
                H_[4,3] = H_[3,4]
                H_[4,4] = -C1C2C
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eE3"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        S1 = sin(EV_[1])
        S2 = sin(EV_[2])
        C1 = cos(EV_[1])
        C2 = cos(EV_[2])
        C = cos(EV_[3]-pbm.elpar[iel_][1])
        S = sin(EV_[3]-pbm.elpar[iel_][1])
        f_   = S1*S2+C1*C2*C
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = C1*S2-S1*C2*C
            g_[2] = S1*C2-C1*S2*C
            g_[3] = -C1*C2*S
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,1] = -S1*S2-C1*C2*C
                H_[2,1] = C1*C2+S1*S2*C
                H_[1,2] = H_[2,1]
                H_[2,2] = -S1*S2-C1*C2*C
                H_[3,1] = S1*C2*S
                H_[1,3] = H_[3,1]
                H_[3,2] = C1*S2*S
                H_[2,3] = H_[3,2]
                H_[3,3] = -C1*C2*C
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eE2"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        S1 = sin(EV_[1])
        SF = sin(pbm.elpar[iel_][1])
        C1 = cos(EV_[1])
        CF = cos(pbm.elpar[iel_][1])
        C = cos(EV_[2]-pbm.elpar[iel_][2])
        S = sin(EV_[2]-pbm.elpar[iel_][2])
        f_   = S1*SF+C1*CF*C
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = C1*SF-S1*CF*C
            g_[2] = -C1*CF*S
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = -S1*SF-C1*CF*C
                H_[2,1] = S1*CF*S
                H_[1,2] = H_[2,1]
                H_[2,2] = -C1*CF*C
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eE1"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        S1 = sin(EV_[1])
        SF = sin(pbm.elpar[iel_][1])
        C1 = cos(EV_[1])
        CF = cos(pbm.elpar[iel_][1])
        C = cos(pbm.elpar[iel_][3]-pbm.elpar[iel_][2])
        S = sin(pbm.elpar[iel_][3]-pbm.elpar[iel_][2])
        f_   = S1*SF+C1*CF*C
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = C1*SF-S1*CF*C
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = -S1*SF-C1*CF*C
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

