function YORKNET(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    A problem arising in the modelling of the Yorkshire water system.
# 
#    Source:
#    an problem submitted for the LANCELOT licence.
# 
#    SIF input: B. Ulanicki, Water Software Systems,De Montfort University,
#               The Gateway, Leicester LE1 9BH, UK.
#               e-mail: bul@uk.ac.dmu * Tel no.0533 577070
#               correction by S. Gratton & Ph. Toint, May 2024
# 
#    classification = "C-SOR2-AY-312-256"
# 
# DECLARE CONSTANTS DESCRIBING NETWORK
# 
# STANDARD DECLARATIONS
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "YORKNET"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["0"] = 0
        v_["1"] = 1
        v_["2"] = 2
        v_["3"] = 3
        v_["4"] = 4
        v_["5"] = 5
        v_["6"] = 6
        v_["7"] = 7
        v_["8"] = 8
        v_["9"] = 9
        v_["10"] = 10
        v_["11"] = 11
        v_["12"] = 12
        v_["13"] = 13
        v_["14"] = 14
        v_["15"] = 15
        v_["16"] = 16
        v_["17"] = 17
        v_["18"] = 18
        v_["19"] = 19
        v_["ONE"] = 1
        v_["NSTEP"] = 8
        v_["NSTEP+1"] = 8+v_["NSTEP"]
        v_["ST1"] = 0.125
        v_["ST2"] = 0.125
        v_["ST3"] = 0.125
        v_["ST4"] = 0.125
        v_["ST5"] = 0.125
        v_["ST6"] = 0.125
        v_["ST7"] = 0.125
        v_["ST8"] = 0.125
        v_["MST1"] = 0.0-v_["ST1"]
        v_["MST2"] = 0.0-v_["ST2"]
        v_["MST3"] = 0.0-v_["ST3"]
        v_["MST4"] = 0.0-v_["ST4"]
        v_["MST5"] = 0.0-v_["ST5"]
        v_["MST6"] = 0.0-v_["ST6"]
        v_["MST7"] = 0.0-v_["ST7"]
        v_["MST8"] = 0.0-v_["ST8"]
        v_["TFL1"] = 1.87
        v_["TFH1"] = 4.03
        v_["TF1"] = 0.0+v_["TFL1"]
        v_["TF2"] = 0.0+v_["TFL1"]
        v_["TF3"] = 0.0+v_["TFH1"]
        v_["TF4"] = 0.0+v_["TFH1"]
        v_["TF5"] = 0.0+v_["TFH1"]
        v_["TF6"] = 0.0+v_["TFH1"]
        v_["TF7"] = 0.0+v_["TFH1"]
        v_["TF8"] = 0.0+v_["TFH1"]
        v_["BP1"] = 24.0*v_["ST1"]
        v_["BP2"] = 24.0*v_["ST2"]
        v_["BP3"] = 24.0*v_["ST3"]
        v_["BP4"] = 24.0*v_["ST4"]
        v_["BP5"] = 24.0*v_["ST5"]
        v_["BP6"] = 24.0*v_["ST6"]
        v_["BP7"] = 24.0*v_["ST7"]
        v_["BP8"] = 24.0*v_["ST8"]
        v_["CP1"] = v_["BP1"]*v_["TF1"]
        v_["CP2"] = v_["BP2"]*v_["TF2"]
        v_["CP3"] = v_["BP3"]*v_["TF3"]
        v_["CP4"] = v_["BP4"]*v_["TF4"]
        v_["CP5"] = v_["BP5"]*v_["TF5"]
        v_["CP6"] = v_["BP6"]*v_["TF6"]
        v_["CP7"] = v_["BP7"]*v_["TF7"]
        v_["CP8"] = v_["BP8"]*v_["TF8"]
        v_["NEL"] = 15
        v_["SuVAL"] = 9
        v_["EuVAL"] = 10
        v_["SuPIPE"] = 1
        v_["EuPIPE"] = 8
        v_["MGINV1"] = -3.365170e-3
        v_["MGINV2"] = -2.314284e-2
        v_["MGINV3"] = -6.631203e-3
        v_["MGINV4"] = -1.702093e-3
        v_["MGINV5"] = -1.205983e-2
        v_["MGINV6"] = -9.617776e-4
        v_["MGINV7"] = -1.392046e-5
        v_["MGINV8"] = -4.411625e-3
        v_["MGINV9"] = -2.019250e-3
        v_["MGINV10"] = -2.288437e-3
        v_["NPMP"] = 5
        v_["SuPMP"] = 11
        v_["EuPMP"] = 15
        v_["CONA11"] = -.035520
        v_["CONB11"] = -.054720
        v_["CONC11"] = 99.80
        v_["CONA12"] = -0.07475
        v_["CONB12"] = -9.05
        v_["CONC12"] = 110
        v_["CONA13"] = -.042420
        v_["CONB13"] = -.005370
        v_["CONC13"] = 175.29
        v_["CONA14"] = -.040733
        v_["CONB14"] = -.032036
        v_["CONC14"] = 139.6
        v_["CONA15"] = -.167495
        v_["CONB15"] = -.0019
        v_["CONC15"] = 139.6
        v_["NND"] = 13
        v_["D1"] = 0.0
        v_["D4"] = -33.0
        v_["D2"] = 0.0
        v_["D3"] = -55.0
        v_["D5"] = 0.0
        v_["D6"] = 0.0
        v_["D7"] = 0.0
        v_["D8"] = -25.0
        v_["D9"] = 0.0
        v_["D10"] = -17.0
        v_["D11"] = 0.0
        v_["D12"] = 0.0
        v_["D13"] = 0.0
        v_["SuRES"] = 1
        v_["EuRES"] = 4
        v_["EuRES+1"] = 1+v_["EuRES"]
        v_["HGT1"] = 5.77
        v_["HGT2"] = 3.00
        v_["HGT3"] = 131.08
        v_["HGT4"] = 44.0
        v_["MXHGT1"] = 9.60
        v_["MXHGT2"] = 7.89
        v_["MXHGT3"] = 138.76
        v_["MXHGT4"] = 53.34
        v_["XAR1"] = 1.599
        v_["XAR2"] = 4.6421
        v_["XAR3"] = 30.2307
        v_["XAR4"] = 5.3938
        v_["MXAR1"] = 0.0-v_["XAR1"]
        v_["MXAR2"] = 0.0-v_["XAR2"]
        v_["MXAR3"] = 0.0-v_["XAR3"]
        v_["MXAR4"] = 0.0-v_["XAR4"]
        v_["RXAR1"] = 1.0/v_["XAR1"]
        v_["MRXAR1"] = 0.0-v_["RXAR1"]
        v_["RXAR2"] = 1.0/v_["XAR2"]
        v_["MRXAR2"] = 0.0-v_["RXAR2"]
        v_["RXAR3"] = 1.0/v_["XAR3"]
        v_["MRXAR3"] = 0.0-v_["RXAR3"]
        v_["RXAR4"] = 1.0/v_["XAR4"]
        v_["MRXAR4"] = 0.0-v_["RXAR4"]
        v_["HTXAR1"] = v_["XAR1"]*v_["HGT1"]
        v_["HTXAR2"] = v_["XAR2"]*v_["HGT2"]
        v_["HTXAR3"] = v_["XAR3"]*v_["HGT3"]
        v_["HTXAR4"] = v_["XAR4"]*v_["HGT4"]
        v_["STHGT1"] = 8.5
        v_["STHGT2"] = 6.0
        v_["STHGT3"] = 135.6
        v_["STHGT4"] = 48.5
        v_["STVOL1"] = v_["STHGT1"]*v_["XAR1"]
        v_["STVOL2"] = v_["STHGT2"]*v_["XAR2"]
        v_["STVOL3"] = v_["STHGT3"]*v_["XAR3"]
        v_["STVOL4"] = v_["STHGT4"]*v_["XAR4"]
        v_["MSTVOL1"] = 0.0-v_["STVOL1"]
        v_["MSTVOL2"] = 0.0-v_["STVOL2"]
        v_["MSTVOL3"] = 0.0-v_["STVOL3"]
        v_["MSTVOL4"] = 0.0-v_["STVOL4"]
        v_["WMN1"] = 3.764
        v_["WMN2"] = 11.35
        v_["WMN3"] = 156.648
        v_["WMN4"] = 45.929
        v_["WMX1"] = 5.646
        v_["WMX2"] = 22.133
        v_["WMX3"] = 223.489
        v_["WMX4"] = 61.876
        v_["H1"] = 8.99
        v_["H2"] = 52.84
        v_["H3"] = 138.31
        v_["H4"] = 5.67
        v_["W1"] = 4.728
        v_["W2"] = 15.601
        v_["W3"] = 190.648
        v_["W4"] = 55.00
        v_["SuTW"] = 1
        v_["EuTW"] = 2
        v_["BCTW1"] = 28.34
        v_["BCTW2"] = 18.86
        v_["BCTW3"] = 1.0
        v_["BCTW4"] = 1.0
        v_["BCTW5"] = 28.34
        v_["BCTW6"] = 18.86
        v_["BCTW7"] = 1.0
        v_["BCTW8"] = 1.0
        v_["CTW1,1"] = v_["BCTW1"]*v_["ST1"]
        v_["CTW1,2"] = v_["BCTW1"]*v_["ST2"]
        v_["CTW1,3"] = v_["BCTW1"]*v_["ST3"]
        v_["CTW1,4"] = v_["BCTW1"]*v_["ST4"]
        v_["CTW1,5"] = v_["BCTW1"]*v_["ST5"]
        v_["CTW1,6"] = v_["BCTW1"]*v_["ST6"]
        v_["CTW1,7"] = v_["BCTW1"]*v_["ST7"]
        v_["CTW1,8"] = v_["BCTW1"]*v_["ST8"]
        v_["CTW2,1"] = v_["BCTW2"]*v_["ST1"]
        v_["CTW2,2"] = v_["BCTW2"]*v_["ST2"]
        v_["CTW2,3"] = v_["BCTW2"]*v_["ST3"]
        v_["CTW2,4"] = v_["BCTW2"]*v_["ST4"]
        v_["CTW2,5"] = v_["BCTW2"]*v_["ST5"]
        v_["CTW2,6"] = v_["BCTW2"]*v_["ST6"]
        v_["CTW2,7"] = v_["BCTW2"]*v_["ST7"]
        v_["CTW2,8"] = v_["BCTW2"]*v_["ST8"]
        v_["WS1"] = 94.67
        v_["WS2"] = 30.87
        v_["QSMX1"] = 500.0
        v_["QSMX2"] = 40.0
        v_["SCuQ"] = 1.0
        v_["SCuH"] = 1.0
        v_["SCuV"] = 1.0
        v_["PIPEuSC"] = 1.0
        v_["VALuSC"] = 1.0
        v_["NDuSC"] = 1.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["ONE"]):Int64(v_["NSTEP"])
            for J = Int64(v_["ONE"]):Int64(v_["NEL"])
                iv,ix_,_ = s2mpj_ii("q"*string(J)*","*string(I),ix_)
                arrset(pb.xnames,iv,"q"*string(J)*","*string(I))
            end
            for J = Int64(v_["ONE"]):Int64(v_["NND"])
                iv,ix_,_ = s2mpj_ii("h"*string(J)*","*string(I),ix_)
                arrset(pb.xnames,iv,"h"*string(J)*","*string(I))
            end
            for J = Int64(v_["SuRES"]):Int64(v_["EuRES"])
                iv,ix_,_ = s2mpj_ii("qr"*string(J)*","*string(I),ix_)
                arrset(pb.xnames,iv,"qr"*string(J)*","*string(I))
            end
            for J = Int64(v_["SuTW"]):Int64(v_["EuTW"])
                iv,ix_,_ = s2mpj_ii("qs"*string(J)*","*string(I),ix_)
                arrset(pb.xnames,iv,"qs"*string(J)*","*string(I))
            end
            for J = Int64(v_["SuPMP"]):Int64(v_["EuPMP"])
                iv,ix_,_ = s2mpj_ii("u"*string(J)*","*string(I),ix_)
                arrset(pb.xnames,iv,"u"*string(J)*","*string(I))
            end
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["ONE"]):Int64(v_["NSTEP"])
            ig,ig_,_ = s2mpj_ii("ND"*string(Int64(v_["1"]))*","*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"ND"*string(Int64(v_["1"]))*","*string(I))
            iv = ix_["qr"*string(Int64(v_["1"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["qs"*string(Int64(v_["1"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("ND"*string(Int64(v_["1"]))*","*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"ND"*string(Int64(v_["1"]))*","*string(I))
            iv = ix_["q"*string(Int64(v_["13"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            ig,ig_,_ = s2mpj_ii("ND"*string(Int64(v_["2"]))*","*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"ND"*string(Int64(v_["2"]))*","*string(I))
            iv = ix_["qr"*string(Int64(v_["2"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["qs"*string(Int64(v_["2"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("ND"*string(Int64(v_["2"]))*","*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"ND"*string(Int64(v_["2"]))*","*string(I))
            iv = ix_["q"*string(Int64(v_["14"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["q"*string(Int64(v_["15"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            ig,ig_,_ = s2mpj_ii("ND"*string(Int64(v_["3"]))*","*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"ND"*string(Int64(v_["3"]))*","*string(I))
            iv = ix_["qr"*string(Int64(v_["3"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["q"*string(Int64(v_["5"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("ND"*string(Int64(v_["3"]))*","*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"ND"*string(Int64(v_["3"]))*","*string(I))
            iv = ix_["q"*string(Int64(v_["6"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            ig,ig_,_ = s2mpj_ii("ND"*string(Int64(v_["4"]))*","*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"ND"*string(Int64(v_["4"]))*","*string(I))
            iv = ix_["qr"*string(Int64(v_["4"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["q"*string(Int64(v_["12"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            ig,ig_,_ = s2mpj_ii("ND"*string(Int64(v_["4"]))*","*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"ND"*string(Int64(v_["4"]))*","*string(I))
            iv = ix_["q"*string(Int64(v_["8"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["q"*string(Int64(v_["9"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("ND"*string(Int64(v_["5"]))*","*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"ND"*string(Int64(v_["5"]))*","*string(I))
            iv = ix_["q"*string(Int64(v_["13"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["q"*string(Int64(v_["1"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            ig,ig_,_ = s2mpj_ii("ND"*string(Int64(v_["6"]))*","*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"ND"*string(Int64(v_["6"]))*","*string(I))
            iv = ix_["q"*string(Int64(v_["2"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["q"*string(Int64(v_["3"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            ig,ig_,_ = s2mpj_ii("ND"*string(Int64(v_["6"]))*","*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"ND"*string(Int64(v_["6"]))*","*string(I))
            iv = ix_["q"*string(Int64(v_["12"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("ND"*string(Int64(v_["7"]))*","*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"ND"*string(Int64(v_["7"]))*","*string(I))
            iv = ix_["q"*string(Int64(v_["2"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["q"*string(Int64(v_["7"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("ND"*string(Int64(v_["7"]))*","*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"ND"*string(Int64(v_["7"]))*","*string(I))
            iv = ix_["q"*string(Int64(v_["11"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            ig,ig_,_ = s2mpj_ii("ND"*string(Int64(v_["8"]))*","*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"ND"*string(Int64(v_["8"]))*","*string(I))
            iv = ix_["q"*string(Int64(v_["3"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["q"*string(Int64(v_["7"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            ig,ig_,_ = s2mpj_ii("ND"*string(Int64(v_["9"]))*","*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"ND"*string(Int64(v_["9"]))*","*string(I))
            iv = ix_["q"*string(Int64(v_["4"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["q"*string(Int64(v_["5"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            ig,ig_,_ = s2mpj_ii("ND"*string(Int64(v_["9"]))*","*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"ND"*string(Int64(v_["9"]))*","*string(I))
            iv = ix_["q"*string(Int64(v_["11"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("ND"*string(Int64(v_["10"]))*","*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"ND"*string(Int64(v_["10"]))*","*string(I))
            iv = ix_["q"*string(Int64(v_["4"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["q"*string(Int64(v_["6"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("ND"*string(Int64(v_["11"]))*","*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"ND"*string(Int64(v_["11"]))*","*string(I))
            iv = ix_["q"*string(Int64(v_["10"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["q"*string(Int64(v_["14"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("ND"*string(Int64(v_["11"]))*","*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"ND"*string(Int64(v_["11"]))*","*string(I))
            iv = ix_["q"*string(Int64(v_["15"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("ND"*string(Int64(v_["12"]))*","*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"ND"*string(Int64(v_["12"]))*","*string(I))
            iv = ix_["q"*string(Int64(v_["8"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["q"*string(Int64(v_["10"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("ND"*string(Int64(v_["13"]))*","*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"ND"*string(Int64(v_["13"]))*","*string(I))
            iv = ix_["q"*string(Int64(v_["9"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["q"*string(Int64(v_["1"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
        end
        for I = Int64(v_["1"]):Int64(v_["NSTEP"])
            ig,ig_,_ = s2mpj_ii("EL"*string(Int64(v_["1"]))*","*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"EL"*string(Int64(v_["1"]))*","*string(I))
            iv = ix_["h"*string(Int64(v_["13"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["h"*string(Int64(v_["5"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            ig,ig_,_ = s2mpj_ii("EL"*string(Int64(v_["2"]))*","*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"EL"*string(Int64(v_["2"]))*","*string(I))
            iv = ix_["h"*string(Int64(v_["7"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["h"*string(Int64(v_["6"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            ig,ig_,_ = s2mpj_ii("EL"*string(Int64(v_["3"]))*","*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"EL"*string(Int64(v_["3"]))*","*string(I))
            iv = ix_["h"*string(Int64(v_["8"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["h"*string(Int64(v_["6"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            ig,ig_,_ = s2mpj_ii("EL"*string(Int64(v_["4"]))*","*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"EL"*string(Int64(v_["4"]))*","*string(I))
            iv = ix_["h"*string(Int64(v_["10"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["h"*string(Int64(v_["9"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            ig,ig_,_ = s2mpj_ii("EL"*string(Int64(v_["5"]))*","*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"EL"*string(Int64(v_["5"]))*","*string(I))
            iv = ix_["h"*string(Int64(v_["3"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["h"*string(Int64(v_["9"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            ig,ig_,_ = s2mpj_ii("EL"*string(Int64(v_["6"]))*","*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"EL"*string(Int64(v_["6"]))*","*string(I))
            iv = ix_["h"*string(Int64(v_["10"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["h"*string(Int64(v_["3"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            ig,ig_,_ = s2mpj_ii("EL"*string(Int64(v_["7"]))*","*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"EL"*string(Int64(v_["7"]))*","*string(I))
            iv = ix_["h"*string(Int64(v_["7"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["h"*string(Int64(v_["8"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            ig,ig_,_ = s2mpj_ii("EL"*string(Int64(v_["8"]))*","*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"EL"*string(Int64(v_["8"]))*","*string(I))
            iv = ix_["h"*string(Int64(v_["4"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["h"*string(Int64(v_["12"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            ig,ig_,_ = s2mpj_ii("EL"*string(Int64(v_["9"]))*","*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"EL"*string(Int64(v_["9"]))*","*string(I))
            iv = ix_["h"*string(Int64(v_["4"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["h"*string(Int64(v_["13"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            ig,ig_,_ = s2mpj_ii("EL"*string(Int64(v_["10"]))*","*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"EL"*string(Int64(v_["10"]))*","*string(I))
            iv = ix_["h"*string(Int64(v_["12"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["h"*string(Int64(v_["11"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            for J = Int64(v_["SuPMP"]):Int64(v_["EuPMP"])
                ig,ig_,_ = s2mpj_ii("EL"*string(J)*","*string(I),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"EL"*string(J)*","*string(I))
            end
        end
        for I = Int64(v_["ONE"]):Int64(v_["NSTEP"])
            for J = Int64(v_["SuTW"]):Int64(v_["EuTW"])
                ig,ig_,_ = s2mpj_ii("TWC"*string(J)*","*string(I),ig_)
                arrset(gtype,ig,"<>")
                iv = ix_["qs"*string(J)*","*string(I)]
                pbm.A[ig,iv] += Float64(v_["CTW"*string(J)*","*string(I)])
            end
            for J = Int64(v_["SuPMP"]):Int64(v_["EuPMP"])
                ig,ig_,_ = s2mpj_ii("PC"*string(J)*","*string(I),ig_)
                arrset(gtype,ig,"<>")
            end
        end
        for J = Int64(v_["SuRES"]):Int64(v_["EuRES"])
            ig,ig_,_ = s2mpj_ii("RD"*string(J)*","*string(Int64(v_["1"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"RD"*string(J)*","*string(Int64(v_["1"])))
            iv = ix_["h"*string(J)*","*string(Int64(v_["1"]))]
            pbm.A[ig,iv] += Float64(v_["MXAR"*string(J)])
            ig,ig_,_ = s2mpj_ii("RD"*string(J)*","*string(Int64(v_["1"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"RD"*string(J)*","*string(Int64(v_["1"])))
            iv = ix_["qr"*string(J)*","*string(Int64(v_["1"]))]
            pbm.A[ig,iv] += Float64(v_["MST"*string(Int64(v_["1"]))])
        end
        for I = Int64(v_["2"]):Int64(v_["NSTEP"])
            for J = Int64(v_["SuRES"]):Int64(v_["EuRES"])
                v_["A"] = -1+I
                ig,ig_,_ = s2mpj_ii("RD"*string(J)*","*string(I),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"RD"*string(J)*","*string(I))
                iv = ix_["h"*string(J)*","*string(I)]
                pbm.A[ig,iv] += Float64(v_["MXAR"*string(J)])
                iv = ix_["h"*string(J)*","*string(Int64(v_["A"]))]
                pbm.A[ig,iv] += Float64(v_["XAR"*string(J)])
                iv = ix_["qr"*string(J)*","*string(I)]
                pbm.A[ig,iv] += Float64(v_["MST"*string(I)])
            end
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
        for I = Int64(v_["ONE"]):Int64(v_["NSTEP"])
            for J = Int64(v_["ONE"]):Int64(v_["NND"])
                pbm.gconst[ig_["ND"*string(J)*","*string(I)]] = Float64(v_["D"*string(J)])
            end
        end
        for J = Int64(v_["SuRES"]):Int64(v_["EuRES"])
            pbm.gconst[ig_["RD"*string(J)*","*string(Int64(v_["1"]))]]  = (
                  Float64(v_["MSTVOL"*string(J)]))
        end
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        for I = Int64(v_["ONE"]):Int64(v_["NSTEP"])
            for J = Int64(v_["SuRES"]):Int64(v_["EuRES"])
                pb.xlower[ix_["h"*string(J)*","*string(I)]] = v_["HGT"*string(J)]
                pb.xupper[ix_["h"*string(J)*","*string(I)]] = v_["MXHGT"*string(J)]
                pb.xlower[ix_["qr"*string(J)*","*string(I)]] = -Inf
                pb.xupper[ix_["qr"*string(J)*","*string(I)]] = +Inf
            end
            for J = Int64(v_["EuRES+1"]):Int64(v_["NND"])
                pb.xlower[ix_["h"*string(J)*","*string(I)]] = 0.0
            end
            for J = Int64(v_["SuTW"]):Int64(v_["EuTW"])
                pb.xlower[ix_["qs"*string(J)*","*string(I)]] = 0.0
                pb.xupper[ix_["qs"*string(J)*","*string(I)]] = v_["QSMX"*string(J)]
            end
            for J = Int64(v_["ONE"]):Int64(v_["NEL"])
                pb.xlower[ix_["q"*string(J)*","*string(I)]] = -Inf
                pb.xupper[ix_["q"*string(J)*","*string(I)]] = +Inf
            end
            for J = Int64(v_["SuPMP"]):Int64(v_["EuPMP"])
                pb.xlower[ix_["u"*string(J)*","*string(I)]] = 0.0
                pb.xupper[ix_["u"*string(J)*","*string(I)]] = 7.0
                pb.xlower[ix_["q"*string(J)*","*string(I)]] = 0.0
            end
        end
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(0.0),pb.n)
        pb.y0 = fill(Float64(0.0),pb.m)
        for I = Int64(v_["1"]):Int64(v_["NSTEP"])
            for J = Int64(v_["EuRES+1"]):Int64(v_["NND"])
                pb.x0[ix_["h"*string(J)*","*string(I)]] = Float64(0.0)
            end
            for J = Int64(v_["ONE"]):Int64(v_["NEL"])
                pb.x0[ix_["q"*string(J)*","*string(I)]] = Float64(20.0)
            end
            for J = Int64(v_["SuPMP"]):Int64(v_["EuPMP"])
                pb.x0[ix_["u"*string(J)*","*string(I)]] = Float64(3.0)
            end
            pb.x0[ix_["qs"*string(Int64(v_["1"]))*","*string(I)]] = Float64(25.0)
            pb.x0[ix_["qs"*string(Int64(v_["2"]))*","*string(I)]] = Float64(50.0)
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eXSQ", iet_)
        loaset(elftv,it,1,"X")
        it,iet_,_ = s2mpj_ii( "eXZ", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Z")
        it,iet_,_ = s2mpj_ii( "eXPOW", iet_)
        loaset(elftv,it,1,"X")
        it,iet_,_ = s2mpj_ii( "eXMYZSQ", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        loaset(elftv,it,3,"Z")
        it,iet_,_ = s2mpj_ii( "ePMP1", iet_)
        loaset(elftv,it,1,"Q")
        it,iet_,_ = s2mpj_ii( "ePMP2", iet_)
        loaset(elftv,it,1,"Q")
        it,iet_,_ = s2mpj_ii( "ePMP3", iet_)
        loaset(elftv,it,1,"Q")
        it,iet_,_ = s2mpj_ii( "ePMP4", iet_)
        loaset(elftv,it,1,"Q")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["NSTEP"])
            for J = Int64(v_["SuPIPE"]):Int64(v_["EuPIPE"])
                ename = "EA"*string(J)*","*string(I)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eXPOW")
                arrset(ielftype,ie,iet_["eXPOW"])
                vname = "q"*string(J)*","*string(I)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
                posev = findfirst(x->x=="X",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
            for J = Int64(v_["SuVAL"]):Int64(v_["EuVAL"])
                ename = "EA"*string(J)*","*string(I)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eXPOW")
                arrset(ielftype,ie,iet_["eXPOW"])
                vname = "q"*string(J)*","*string(I)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
                posev = findfirst(x->x=="X",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
            for J = Int64(v_["SuPMP"]):Int64(v_["EuPMP"])
                ename = "EA"*string(J)*","*string(I)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eXSQ")
                arrset(ielftype,ie,iet_["eXSQ"])
                vname = "q"*string(J)*","*string(I)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
                posev = findfirst(x->x=="X",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                ename = "EB"*string(J)*","*string(I)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eXZ")
                arrset(ielftype,ie,iet_["eXZ"])
                vname = "q"*string(J)*","*string(I)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
                posev = findfirst(x->x=="X",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "u"*string(J)*","*string(I)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
                posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                ename = "EC"*string(J)*","*string(I)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eXSQ")
                arrset(ielftype,ie,iet_["eXSQ"])
                vname = "u"*string(J)*","*string(I)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
                posev = findfirst(x->x=="X",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                ename = "EH"*string(J)*","*string(I)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eXMYZSQ")
                arrset(ielftype,ie,iet_["eXMYZSQ"])
                vname = "u"*string(J)*","*string(I)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
                posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
            ename = "EH"*string(Int64(v_["11"]))*","*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "h"*string(Int64(v_["7"]))*","*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "EH"*string(Int64(v_["11"]))*","*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "h"*string(Int64(v_["9"]))*","*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
            posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "EH"*string(Int64(v_["12"]))*","*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "h"*string(Int64(v_["4"]))*","*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "EH"*string(Int64(v_["12"]))*","*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "h"*string(Int64(v_["6"]))*","*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
            posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "EH"*string(Int64(v_["13"]))*","*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "h"*string(Int64(v_["1"]))*","*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "EH"*string(Int64(v_["13"]))*","*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "h"*string(Int64(v_["5"]))*","*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
            posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "EH"*string(Int64(v_["14"]))*","*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "h"*string(Int64(v_["2"]))*","*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "EH"*string(Int64(v_["14"]))*","*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "h"*string(Int64(v_["11"]))*","*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
            posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "EH"*string(Int64(v_["15"]))*","*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "h"*string(Int64(v_["2"]))*","*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "EH"*string(Int64(v_["15"]))*","*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "h"*string(Int64(v_["11"]))*","*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
            posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        for I = Int64(v_["1"]):Int64(v_["NSTEP"])
            ename = "PPW"*string(Int64(v_["11"]))*","*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"ePMP4")
            arrset(ielftype,ie,iet_["ePMP4"])
            ename = "PPW"*string(Int64(v_["12"]))*","*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"ePMP3")
            arrset(ielftype,ie,iet_["ePMP3"])
            ename = "PPW"*string(Int64(v_["13"]))*","*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"ePMP1")
            arrset(ielftype,ie,iet_["ePMP1"])
            ename = "PPW"*string(Int64(v_["14"]))*","*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"ePMP2")
            arrset(ielftype,ie,iet_["ePMP2"])
            ename = "PPW"*string(Int64(v_["15"]))*","*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"ePMP2")
            arrset(ielftype,ie,iet_["ePMP2"])
            for J = Int64(v_["SuPMP"]):Int64(v_["EuPMP"])
                ename = "PPW"*string(J)*","*string(I)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                vname = "q"*string(J)*","*string(I)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
                posev = findfirst(x->x=="Q",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["NSTEP"])
            for J = Int64(v_["SuPIPE"]):Int64(v_["EuPIPE"])
                ig = ig_["EL"*string(J)*","*string(I)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["EA"*string(J)*","*string(I)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["MGINV"*string(J)]))
            end
            for J = Int64(v_["SuVAL"]):Int64(v_["EuVAL"])
                ig = ig_["EL"*string(J)*","*string(I)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["EA"*string(J)*","*string(I)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["MGINV"*string(J)]))
            end
            for J = Int64(v_["SuPMP"]):Int64(v_["EuPMP"])
                ig = ig_["EL"*string(J)*","*string(I)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["EA"*string(J)*","*string(I)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["CONA"*string(J)]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["EB"*string(J)*","*string(I)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["CONB"*string(J)]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["EC"*string(J)*","*string(I)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["CONC"*string(J)]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["EH"*string(J)*","*string(I)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(1.0))
            end
        end
        for J = Int64(v_["SuPMP"]):Int64(v_["EuPMP"])
            ig = ig_["PC"*string(J)*","*string(Int64(v_["2"]))]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["PPW"*string(J)*","*string(Int64(v_["2"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(100.0))
            ig = ig_["PC"*string(J)*","*string(Int64(v_["3"]))]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["PPW"*string(J)*","*string(Int64(v_["3"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(100.0))
            ig = ig_["PC"*string(J)*","*string(Int64(v_["4"]))]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["PPW"*string(J)*","*string(Int64(v_["4"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(100.0))
        end
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
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-SOR2-AY-312-256"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2]
        pb.conderlvl  = pbm.conderlvl;
        return pb, pbm


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "eXSQ"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]*EV_[1]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[1]+EV_[1]
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 2.0
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eXZ"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]*EV_[2]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[2]
            g_[2] = EV_[1]
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = 1.0
                H_[2,1] = H_[1,2]
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eXMYZSQ"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,2,3)
        IV_ =  zeros(Float64,2)
        U_[1,1] = U_[1,1]+1
        U_[1,2] = U_[1,2]-1
        U_[2,3] = U_[2,3]+1
        IV_[1] = dot(U_[1,:],EV_)
        IV_[2] = dot(U_[2,:],EV_)
        f_   = IV_[1]*IV_[2]^2
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = IV_[2]^2
            g_[2] = 2.0*IV_[1]*IV_[2]
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = 2.0*IV_[2]
                H_[2,1] = H_[1,2]
                H_[2,2] = 2.0*IV_[1]
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

    elseif action == "eXPOW"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        MODX = abs(EV_[1])
        ISNEG = EV_[1]<0.0
        ISPOS = EV_[1]>=0.0
        if ISNEG
            SGN = -1.0
        end
        if ISPOS
            SGN = +1.0
        end
        f_   = SGN*MODX^1.852
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 1.852*MODX^0.852
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = SGN*1.577904*MODX^(-0.148)
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "ePMP1"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = 0.074*EV_[1]*EV_[1]+3.062*EV_[1]+50.357
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 0.148*EV_[1]+3.062
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 0.148
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "ePMP2"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = 0.747*EV_[1]*EV_[1]-10.287*EV_[1]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 1.494*EV_[1]-10.287
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 1.494
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "ePMP3"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = 0.034*EV_[1]*EV_[1]+0.220*EV_[1]+6.685
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 0.068*EV_[1]+0.220
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 0.068
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "ePMP4"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = 0.079*EV_[1]*EV_[1]-2.761*EV_[1]+35.014
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 0.158*EV_[1]-2.761
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 0.158
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

