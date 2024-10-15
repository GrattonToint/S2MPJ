function n3PK(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : n3PK
#    *********
# 
#    A problem arising in the estimation of structured O/D matrix
# 
#    Source:  
#    M. Bierlaire, private communication
#    see also
#    M. Bierlaire and Ph. L. Toint,
#    MEUSE: an origin-destination estimator that exploits structure,
#    Transportation Research B, 29, 1, 47--60, 1995.
# 
#    SIF input: Ph. Toint, Dec 1989, Corrected July 1993.
# 
#    classification = "C-SBR2-MN-30-0"
# 
#  Parameters
# 
#  Number of parking columns
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "n3PK"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["NPKC"] = 3
        v_["NPKC-1"] = -1+v_["NPKC"]
        v_["NPKC+1"] = 1+v_["NPKC"]
        v_["NCENT"] = 6
        v_["NCENT-1"] = -1+v_["NCENT"]
        v_["RNCENT-1"] = Float64(v_["NCENT-1"])
        v_["GAMMA"] = 1.0000e+04
        v_["FT0"] = 0.500000
        v_["FT1"] = 0.500000
        v_["FT2"] = 0.500000
        v_["WFT0"] = 1.000000
        v_["WFT1"] = 1.000000
        v_["WFT2"] = 1.000000
        v_["COUNT"] = 9
        v_["COUNT-1"] = -1+v_["COUNT"]
        v_["DEFW"] = 999.999953
        v_["0"] = 0
        v_["1"] = 1
        v_["COU0"] = 910.000000
        v_["COU1"] = 175.000000
        v_["COU2"] = 1915.000000
        v_["COU3"] = 450.000000
        v_["COU4"] = 260.000000
        v_["COU5"] = 80.000000
        v_["COU6"] = 670.000000
        v_["COU7"] = 1450.000000
        v_["COU8"] = 990.000000
        v_["PHI0"] = 1.000000
        v_["PHI1"] = 1.000000
        v_["PHI2"] = 1.000000
        v_["PHI3"] = 1.000000
        v_["PHI4"] = 1.000000
        v_["PHI5"] = 1.000000
        v_["PHI6"] = 1.000000
        v_["PHI7"] = 1.000000
        v_["PHI8"] = 1.000000
        for I = Int64(v_["0"]):Int64(v_["COUNT-1"])
            v_["PHI"*string(I)] = v_["PHI"*string(I)]/v_["GAMMA"]
        end
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        iv,ix_,_ = s2mpj_ii("A1,0",ix_)
        arrset(pb.xnames,iv,"A1,0")
        iv,ix_,_ = s2mpj_ii("A2,0",ix_)
        arrset(pb.xnames,iv,"A2,0")
        iv,ix_,_ = s2mpj_ii("A3,0",ix_)
        arrset(pb.xnames,iv,"A3,0")
        iv,ix_,_ = s2mpj_ii("A4,0",ix_)
        arrset(pb.xnames,iv,"A4,0")
        iv,ix_,_ = s2mpj_ii("A5,0",ix_)
        arrset(pb.xnames,iv,"A5,0")
        iv,ix_,_ = s2mpj_ii("A0,1",ix_)
        arrset(pb.xnames,iv,"A0,1")
        iv,ix_,_ = s2mpj_ii("A2,1",ix_)
        arrset(pb.xnames,iv,"A2,1")
        iv,ix_,_ = s2mpj_ii("A3,1",ix_)
        arrset(pb.xnames,iv,"A3,1")
        iv,ix_,_ = s2mpj_ii("A4,1",ix_)
        arrset(pb.xnames,iv,"A4,1")
        iv,ix_,_ = s2mpj_ii("A5,1",ix_)
        arrset(pb.xnames,iv,"A5,1")
        iv,ix_,_ = s2mpj_ii("A0,2",ix_)
        arrset(pb.xnames,iv,"A0,2")
        iv,ix_,_ = s2mpj_ii("A1,2",ix_)
        arrset(pb.xnames,iv,"A1,2")
        iv,ix_,_ = s2mpj_ii("A3,2",ix_)
        arrset(pb.xnames,iv,"A3,2")
        iv,ix_,_ = s2mpj_ii("A4,2",ix_)
        arrset(pb.xnames,iv,"A4,2")
        iv,ix_,_ = s2mpj_ii("A5,2",ix_)
        arrset(pb.xnames,iv,"A5,2")
        for J = Int64(v_["NPKC"]):Int64(v_["NCENT-1"])
            v_["J+1"] = 1+J
            v_["J-1"] = -1+J
            for I = Int64(v_["0"]):Int64(v_["J-1"])
                iv,ix_,_ = s2mpj_ii("T"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"T"*string(I)*","*string(J))
            end
            for I = Int64(v_["J+1"]):Int64(v_["NCENT-1"])
                iv,ix_,_ = s2mpj_ii("T"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"T"*string(I)*","*string(J))
            end
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("G0,3",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T0,3"]
        pbm.A[ig,iv] += Float64(0.010000)
        ig,ig_,_ = s2mpj_ii("G1,3",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T1,3"]
        pbm.A[ig,iv] += Float64(0.007143)
        ig,ig_,_ = s2mpj_ii("G2,3",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T2,3"]
        pbm.A[ig,iv] += Float64(0.008333)
        ig,ig_,_ = s2mpj_ii("G4,3",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T4,3"]
        pbm.A[ig,iv] += Float64(0.050000)
        ig,ig_,_ = s2mpj_ii("G5,3",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T5,3"]
        pbm.A[ig,iv] += Float64(0.050000)
        ig,ig_,_ = s2mpj_ii("G0,4",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T0,4"]
        pbm.A[ig,iv] += Float64(0.005000)
        ig,ig_,_ = s2mpj_ii("G1,4",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T1,4"]
        pbm.A[ig,iv] += Float64(0.005556)
        ig,ig_,_ = s2mpj_ii("G2,4",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T2,4"]
        pbm.A[ig,iv] += Float64(0.050000)
        ig,ig_,_ = s2mpj_ii("G3,4",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T3,4"]
        pbm.A[ig,iv] += Float64(0.001667)
        ig,ig_,_ = s2mpj_ii("G5,4",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T5,4"]
        pbm.A[ig,iv] += Float64(0.025000)
        ig,ig_,_ = s2mpj_ii("G0,5",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T0,5"]
        pbm.A[ig,iv] += Float64(0.020000)
        ig,ig_,_ = s2mpj_ii("G1,5",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T1,5"]
        pbm.A[ig,iv] += Float64(0.033333)
        ig,ig_,_ = s2mpj_ii("G2,5",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T2,5"]
        pbm.A[ig,iv] += Float64(0.014286)
        ig,ig_,_ = s2mpj_ii("G3,5",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T3,5"]
        pbm.A[ig,iv] += Float64(0.006667)
        ig,ig_,_ = s2mpj_ii("G4,5",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T4,5"]
        pbm.A[ig,iv] += Float64(0.050000)
        v_["TMP"] = 5.000000*v_["FT0"]
        v_["TMP1"] = 1.0/v_["TMP"]
        ig,ig_,_ = s2mpj_ii("H0",ig_)
        arrset(gtype,ig,"<>")
        arrset(pbm.gscale,ig,Float64(v_["WFT0"]))
        iv = ix_["A1,0"]
        pbm.A[ig,iv] += Float64(v_["TMP1"])
        iv = ix_["A2,0"]
        pbm.A[ig,iv] += Float64(v_["TMP1"])
        iv = ix_["A3,0"]
        pbm.A[ig,iv] += Float64(v_["TMP1"])
        iv = ix_["A4,0"]
        pbm.A[ig,iv] += Float64(v_["TMP1"])
        iv = ix_["A5,0"]
        pbm.A[ig,iv] += Float64(v_["TMP1"])
        v_["TMP"] = 5.000000*v_["FT1"]
        v_["TMP1"] = 1.0/v_["TMP"]
        ig,ig_,_ = s2mpj_ii("H1",ig_)
        arrset(gtype,ig,"<>")
        arrset(pbm.gscale,ig,Float64(v_["WFT1"]))
        iv = ix_["A0,1"]
        pbm.A[ig,iv] += Float64(v_["TMP1"])
        iv = ix_["A2,1"]
        pbm.A[ig,iv] += Float64(v_["TMP1"])
        iv = ix_["A3,1"]
        pbm.A[ig,iv] += Float64(v_["TMP1"])
        iv = ix_["A4,1"]
        pbm.A[ig,iv] += Float64(v_["TMP1"])
        iv = ix_["A5,1"]
        pbm.A[ig,iv] += Float64(v_["TMP1"])
        v_["TMP"] = 5.000000*v_["FT2"]
        v_["TMP1"] = 1.0/v_["TMP"]
        ig,ig_,_ = s2mpj_ii("H2",ig_)
        arrset(gtype,ig,"<>")
        arrset(pbm.gscale,ig,Float64(v_["WFT2"]))
        iv = ix_["A0,2"]
        pbm.A[ig,iv] += Float64(v_["TMP1"])
        iv = ix_["A1,2"]
        pbm.A[ig,iv] += Float64(v_["TMP1"])
        iv = ix_["A3,2"]
        pbm.A[ig,iv] += Float64(v_["TMP1"])
        iv = ix_["A4,2"]
        pbm.A[ig,iv] += Float64(v_["TMP1"])
        iv = ix_["A5,2"]
        pbm.A[ig,iv] += Float64(v_["TMP1"])
        for I = Int64(v_["0"]):Int64(v_["COUNT-1"])
            ig,ig_,_ = s2mpj_ii("K"*string(I),ig_)
            arrset(gtype,ig,"<>")
            arrset(pbm.gscale,ig,Float64(v_["PHI"*string(I)]))
        end
        v_["TMP1"] = 200.000000/v_["COU7"]
        v_["TMP"] = 1.000000*v_["TMP1"]
        ig,ig_,_ = s2mpj_ii("K7",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A1,0"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP1"] = 200.000000/v_["COU4"]
        v_["TMP"] = 1.000000*v_["TMP1"]
        ig,ig_,_ = s2mpj_ii("K4",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A1,0"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP1"] = 200.000000/v_["COU2"]
        v_["TMP"] = 1.000000*v_["TMP1"]
        ig,ig_,_ = s2mpj_ii("K2",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A1,0"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP1"] = 480.000000/v_["COU8"]
        v_["TMP"] = 1.000000*v_["TMP1"]
        ig,ig_,_ = s2mpj_ii("K8",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A2,0"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP1"] = 480.000000/v_["COU7"]
        v_["TMP"] = 1.000000*v_["TMP1"]
        ig,ig_,_ = s2mpj_ii("K7",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A2,0"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP1"] = 480.000000/v_["COU6"]
        v_["TMP"] = 1.000000*v_["TMP1"]
        ig,ig_,_ = s2mpj_ii("K6",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A2,0"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP1"] = 480.000000/v_["COU2"]
        v_["TMP"] = 1.000000*v_["TMP1"]
        ig,ig_,_ = s2mpj_ii("K2",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A2,0"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP1"] = 120.000000/v_["COU2"]
        v_["TMP"] = 1.000000*v_["TMP1"]
        iv = ix_["A3,0"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP1"] = 360.000000/v_["COU7"]
        v_["TMP"] = 1.000000*v_["TMP1"]
        ig,ig_,_ = s2mpj_ii("K7",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A4,0"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP1"] = 360.000000/v_["COU2"]
        v_["TMP"] = 1.000000*v_["TMP1"]
        ig,ig_,_ = s2mpj_ii("K2",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A4,0"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP1"] = 560.000000/v_["COU8"]
        v_["TMP"] = 1.000000*v_["TMP1"]
        ig,ig_,_ = s2mpj_ii("K8",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A5,0"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP1"] = 560.000000/v_["COU7"]
        v_["TMP"] = 1.000000*v_["TMP1"]
        ig,ig_,_ = s2mpj_ii("K7",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A5,0"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP1"] = 560.000000/v_["COU2"]
        v_["TMP"] = 1.000000*v_["TMP1"]
        ig,ig_,_ = s2mpj_ii("K2",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A5,0"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP1"] = 240.000000/v_["COU0"]
        v_["TMP"] = 1.000000*v_["TMP1"]
        ig,ig_,_ = s2mpj_ii("K0",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A0,1"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP1"] = 400.000000/v_["COU8"]
        v_["TMP"] = 1.000000*v_["TMP1"]
        ig,ig_,_ = s2mpj_ii("K8",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A2,1"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP1"] = 400.000000/v_["COU7"]
        v_["TMP"] = 1.000000*v_["TMP1"]
        ig,ig_,_ = s2mpj_ii("K7",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A2,1"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP1"] = 400.000000/v_["COU6"]
        v_["TMP"] = 1.000000*v_["TMP1"]
        ig,ig_,_ = s2mpj_ii("K6",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A2,1"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP1"] = 400.000000/v_["COU2"]
        v_["TMP"] = 1.000000*v_["TMP1"]
        ig,ig_,_ = s2mpj_ii("K2",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A2,1"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP1"] = 400.000000/v_["COU0"]
        v_["TMP"] = 1.000000*v_["TMP1"]
        ig,ig_,_ = s2mpj_ii("K0",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A2,1"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP1"] = 420.000000/v_["COU2"]
        v_["TMP"] = 1.000000*v_["TMP1"]
        ig,ig_,_ = s2mpj_ii("K2",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A3,1"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP1"] = 420.000000/v_["COU0"]
        v_["TMP"] = 1.000000*v_["TMP1"]
        ig,ig_,_ = s2mpj_ii("K0",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A3,1"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP1"] = 180.000000/v_["COU7"]
        v_["TMP"] = 1.000000*v_["TMP1"]
        ig,ig_,_ = s2mpj_ii("K7",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A4,1"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP1"] = 180.000000/v_["COU2"]
        v_["TMP"] = 1.000000*v_["TMP1"]
        ig,ig_,_ = s2mpj_ii("K2",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A4,1"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP1"] = 180.000000/v_["COU0"]
        v_["TMP"] = 1.000000*v_["TMP1"]
        ig,ig_,_ = s2mpj_ii("K0",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A4,1"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP1"] = 320.000000/v_["COU8"]
        v_["TMP"] = 1.000000*v_["TMP1"]
        ig,ig_,_ = s2mpj_ii("K8",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A5,1"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP1"] = 320.000000/v_["COU7"]
        v_["TMP"] = 1.000000*v_["TMP1"]
        ig,ig_,_ = s2mpj_ii("K7",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A5,1"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP1"] = 320.000000/v_["COU2"]
        v_["TMP"] = 1.000000*v_["TMP1"]
        ig,ig_,_ = s2mpj_ii("K2",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A5,1"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP1"] = 320.000000/v_["COU0"]
        v_["TMP"] = 1.000000*v_["TMP1"]
        ig,ig_,_ = s2mpj_ii("K0",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A5,1"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP1"] = 20.000000/v_["COU1"]
        v_["TMP"] = 1.000000*v_["TMP1"]
        ig,ig_,_ = s2mpj_ii("K1",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A0,2"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP1"] = 20.000000/v_["COU0"]
        v_["TMP"] = 1.000000*v_["TMP1"]
        ig,ig_,_ = s2mpj_ii("K0",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A0,2"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP1"] = 60.000000/v_["COU1"]
        v_["TMP"] = 1.000000*v_["TMP1"]
        ig,ig_,_ = s2mpj_ii("K1",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A1,2"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP1"] = 40.000000/v_["COU2"]
        v_["TMP"] = 1.000000*v_["TMP1"]
        ig,ig_,_ = s2mpj_ii("K2",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A3,2"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP1"] = 40.000000/v_["COU1"]
        v_["TMP"] = 1.000000*v_["TMP1"]
        ig,ig_,_ = s2mpj_ii("K1",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A3,2"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP1"] = 40.000000/v_["COU0"]
        v_["TMP"] = 1.000000*v_["TMP1"]
        ig,ig_,_ = s2mpj_ii("K0",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A3,2"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP1"] = 120.000000/v_["COU5"]
        v_["TMP"] = 1.000000*v_["TMP1"]
        ig,ig_,_ = s2mpj_ii("K5",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A4,2"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP1"] = 20.000000/v_["COU8"]
        v_["TMP"] = 1.000000*v_["TMP1"]
        ig,ig_,_ = s2mpj_ii("K8",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A5,2"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP1"] = 20.000000/v_["COU5"]
        v_["TMP"] = 1.000000*v_["TMP1"]
        ig,ig_,_ = s2mpj_ii("K5",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A5,2"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP"] = 1.000000/v_["COU7"]
        ig,ig_,_ = s2mpj_ii("K7",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T0,3"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP"] = 1.000000/v_["COU3"]
        ig,ig_,_ = s2mpj_ii("K3",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T0,3"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP"] = 1.000000/v_["COU7"]
        ig,ig_,_ = s2mpj_ii("K7",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T1,3"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP"] = 1.000000/v_["COU4"]
        ig,ig_,_ = s2mpj_ii("K4",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T1,3"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP"] = 1.000000/v_["COU8"]
        ig,ig_,_ = s2mpj_ii("K8",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T2,3"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP"] = 1.000000/v_["COU7"]
        ig,ig_,_ = s2mpj_ii("K7",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T2,3"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP"] = 1.000000/v_["COU6"]
        ig,ig_,_ = s2mpj_ii("K6",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T2,3"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP"] = 1.000000/v_["COU7"]
        ig,ig_,_ = s2mpj_ii("K7",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T4,3"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP"] = 1.000000/v_["COU8"]
        ig,ig_,_ = s2mpj_ii("K8",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T5,3"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP"] = 1.000000/v_["COU7"]
        ig,ig_,_ = s2mpj_ii("K7",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T5,3"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP"] = 1.000000/v_["COU3"]
        ig,ig_,_ = s2mpj_ii("K3",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T0,4"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP"] = 1.000000/v_["COU4"]
        ig,ig_,_ = s2mpj_ii("K4",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T1,4"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP"] = 1.000000/v_["COU8"]
        ig,ig_,_ = s2mpj_ii("K8",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T2,4"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP"] = 1.000000/v_["COU6"]
        ig,ig_,_ = s2mpj_ii("K6",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T2,4"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP"] = 1.000000/v_["COU3"]
        ig,ig_,_ = s2mpj_ii("K3",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T3,4"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP"] = 1.000000/v_["COU2"]
        ig,ig_,_ = s2mpj_ii("K2",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T3,4"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP"] = 1.000000/v_["COU8"]
        ig,ig_,_ = s2mpj_ii("K8",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T5,4"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP"] = 1.000000/v_["COU6"]
        ig,ig_,_ = s2mpj_ii("K6",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T0,5"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP"] = 1.000000/v_["COU1"]
        ig,ig_,_ = s2mpj_ii("K1",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T0,5"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP"] = 1.000000/v_["COU0"]
        ig,ig_,_ = s2mpj_ii("K0",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T0,5"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP"] = 1.000000/v_["COU6"]
        ig,ig_,_ = s2mpj_ii("K6",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T1,5"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP"] = 1.000000/v_["COU1"]
        ig,ig_,_ = s2mpj_ii("K1",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T1,5"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP"] = 1.000000/v_["COU6"]
        ig,ig_,_ = s2mpj_ii("K6",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T2,5"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP"] = 1.000000/v_["COU6"]
        iv = ix_["T3,5"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP"] = 1.000000/v_["COU2"]
        ig,ig_,_ = s2mpj_ii("K2",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T3,5"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP"] = 1.000000/v_["COU1"]
        ig,ig_,_ = s2mpj_ii("K1",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T3,5"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP"] = 1.000000/v_["COU0"]
        ig,ig_,_ = s2mpj_ii("K0",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T3,5"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP"] = 1.000000/v_["COU6"]
        ig,ig_,_ = s2mpj_ii("K6",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T4,5"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        v_["TMP"] = 1.000000/v_["COU5"]
        ig,ig_,_ = s2mpj_ii("K5",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T4,5"]
        pbm.A[ig,iv] += Float64(v_["TMP"])
        ig,ig_,_ = s2mpj_ii("L1,0",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A1,0"]
        pbm.A[ig,iv] += Float64(-0.800000)
        arrset(pbm.gscale,ig,Float64(0.500000))
        iv = ix_["A2,0"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A3,0"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A4,0"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A5,0"]
        pbm.A[ig,iv] += Float64(0.200000)
        ig,ig_,_ = s2mpj_ii("L2,0",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A1,0"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A2,0"]
        pbm.A[ig,iv] += Float64(-0.800000)
        arrset(pbm.gscale,ig,Float64(0.500000))
        iv = ix_["A3,0"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A4,0"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A5,0"]
        pbm.A[ig,iv] += Float64(0.200000)
        ig,ig_,_ = s2mpj_ii("L3,0",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A1,0"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A2,0"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A3,0"]
        pbm.A[ig,iv] += Float64(-0.800000)
        arrset(pbm.gscale,ig,Float64(0.500000))
        iv = ix_["A4,0"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A5,0"]
        pbm.A[ig,iv] += Float64(0.200000)
        ig,ig_,_ = s2mpj_ii("L4,0",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A1,0"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A2,0"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A3,0"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A4,0"]
        pbm.A[ig,iv] += Float64(-0.800000)
        arrset(pbm.gscale,ig,Float64(0.500000))
        iv = ix_["A5,0"]
        pbm.A[ig,iv] += Float64(0.200000)
        ig,ig_,_ = s2mpj_ii("L5,0",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A1,0"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A2,0"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A3,0"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A4,0"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A5,0"]
        pbm.A[ig,iv] += Float64(-0.800000)
        arrset(pbm.gscale,ig,Float64(0.500000))
        ig,ig_,_ = s2mpj_ii("L0,1",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A0,1"]
        pbm.A[ig,iv] += Float64(-0.800000)
        arrset(pbm.gscale,ig,Float64(0.500000))
        iv = ix_["A2,1"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A3,1"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A4,1"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A5,1"]
        pbm.A[ig,iv] += Float64(0.200000)
        ig,ig_,_ = s2mpj_ii("L2,1",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A0,1"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A2,1"]
        pbm.A[ig,iv] += Float64(-0.800000)
        arrset(pbm.gscale,ig,Float64(0.500000))
        iv = ix_["A3,1"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A4,1"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A5,1"]
        pbm.A[ig,iv] += Float64(0.200000)
        ig,ig_,_ = s2mpj_ii("L3,1",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A0,1"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A2,1"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A3,1"]
        pbm.A[ig,iv] += Float64(-0.800000)
        arrset(pbm.gscale,ig,Float64(0.500000))
        iv = ix_["A4,1"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A5,1"]
        pbm.A[ig,iv] += Float64(0.200000)
        ig,ig_,_ = s2mpj_ii("L4,1",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A0,1"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A2,1"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A3,1"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A4,1"]
        pbm.A[ig,iv] += Float64(-0.800000)
        arrset(pbm.gscale,ig,Float64(0.500000))
        iv = ix_["A5,1"]
        pbm.A[ig,iv] += Float64(0.200000)
        ig,ig_,_ = s2mpj_ii("L5,1",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A0,1"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A2,1"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A3,1"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A4,1"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A5,1"]
        pbm.A[ig,iv] += Float64(-0.800000)
        arrset(pbm.gscale,ig,Float64(0.500000))
        ig,ig_,_ = s2mpj_ii("L0,2",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A0,2"]
        pbm.A[ig,iv] += Float64(-0.800000)
        arrset(pbm.gscale,ig,Float64(0.500000))
        iv = ix_["A1,2"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A3,2"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A4,2"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A5,2"]
        pbm.A[ig,iv] += Float64(0.200000)
        ig,ig_,_ = s2mpj_ii("L1,2",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A0,2"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A1,2"]
        pbm.A[ig,iv] += Float64(-0.800000)
        arrset(pbm.gscale,ig,Float64(0.500000))
        iv = ix_["A3,2"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A4,2"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A5,2"]
        pbm.A[ig,iv] += Float64(0.200000)
        ig,ig_,_ = s2mpj_ii("L3,2",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A0,2"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A1,2"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A3,2"]
        pbm.A[ig,iv] += Float64(-0.800000)
        arrset(pbm.gscale,ig,Float64(0.500000))
        iv = ix_["A4,2"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A5,2"]
        pbm.A[ig,iv] += Float64(0.200000)
        ig,ig_,_ = s2mpj_ii("L4,2",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A0,2"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A1,2"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A3,2"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A4,2"]
        pbm.A[ig,iv] += Float64(-0.800000)
        arrset(pbm.gscale,ig,Float64(0.500000))
        iv = ix_["A5,2"]
        pbm.A[ig,iv] += Float64(0.200000)
        ig,ig_,_ = s2mpj_ii("L5,2",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["A0,2"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A1,2"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A3,2"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A4,2"]
        pbm.A[ig,iv] += Float64(0.200000)
        iv = ix_["A5,2"]
        pbm.A[ig,iv] += Float64(-0.800000)
        arrset(pbm.gscale,ig,Float64(0.500000))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        for J = Int64(v_["NPKC"]):Int64(v_["NCENT-1"])
            v_["J+1"] = 1+J
            v_["J-1"] = -1+J
            for I = Int64(v_["0"]):Int64(v_["J-1"])
                pbm.gconst[ig_["G"*string(I)*","*string(J)]] = Float64(1.0)
            end
            for I = Int64(v_["J+1"]):Int64(v_["NCENT-1"])
                pbm.gconst[ig_["G"*string(I)*","*string(J)]] = Float64(1.0)
            end
        end
        for J = Int64(v_["0"]):Int64(v_["NPKC-1"])
            pbm.gconst[ig_["H"*string(J)]] = Float64(1.0)
        end
        for J = Int64(v_["0"]):Int64(v_["COUNT-1"])
            pbm.gconst[ig_["K"*string(J)]] = Float64(1.0)
        end
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.x0[ix_["A1,0"]] = Float64(v_["FT0"])
        pb.x0[ix_["A2,0"]] = Float64(v_["FT0"])
        pb.x0[ix_["A3,0"]] = Float64(v_["FT0"])
        pb.x0[ix_["A4,0"]] = Float64(v_["FT0"])
        pb.x0[ix_["A5,0"]] = Float64(v_["FT0"])
        pb.x0[ix_["A0,1"]] = Float64(v_["FT1"])
        pb.x0[ix_["A2,1"]] = Float64(v_["FT1"])
        pb.x0[ix_["A3,1"]] = Float64(v_["FT1"])
        pb.x0[ix_["A4,1"]] = Float64(v_["FT1"])
        pb.x0[ix_["A5,1"]] = Float64(v_["FT1"])
        pb.x0[ix_["A0,2"]] = Float64(v_["FT2"])
        pb.x0[ix_["A1,2"]] = Float64(v_["FT2"])
        pb.x0[ix_["A3,2"]] = Float64(v_["FT2"])
        pb.x0[ix_["A4,2"]] = Float64(v_["FT2"])
        pb.x0[ix_["A5,2"]] = Float64(v_["FT2"])
        pb.x0[ix_["T0,3"]] = Float64(100.000000)
        pb.x0[ix_["T1,3"]] = Float64(140.000000)
        pb.x0[ix_["T2,3"]] = Float64(120.000000)
        pb.x0[ix_["T4,3"]] = Float64(20.000000)
        pb.x0[ix_["T5,3"]] = Float64(20.000000)
        pb.x0[ix_["T0,4"]] = Float64(200.000000)
        pb.x0[ix_["T1,4"]] = Float64(180.000000)
        pb.x0[ix_["T2,4"]] = Float64(20.000000)
        pb.x0[ix_["T3,4"]] = Float64(600.000000)
        pb.x0[ix_["T5,4"]] = Float64(40.000000)
        pb.x0[ix_["T0,5"]] = Float64(50.000000)
        pb.x0[ix_["T1,5"]] = Float64(30.000000)
        pb.x0[ix_["T2,5"]] = Float64(70.000000)
        pb.x0[ix_["T3,5"]] = Float64(150.000000)
        pb.x0[ix_["T4,5"]] = Float64(20.000000)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gSQUARE",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-SBR2-MN-30-0"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
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

