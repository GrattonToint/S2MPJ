function CmRELOAD(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : CmRELOAD
#    *********
# 
#    Source: Nuclear Reactor Core Reload Pattern Optimization
#    A.J. Quist et.al., draft paper, September 1997.
#    (2nd data set implemented here)
#    SIF input: S. Leyffer, November 1997
# 
#    classification = "C-LOR2-MN-342-284"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "CmRELOAD"

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
        v_["N"] = 14
        v_["M"] = 3
        v_["L"] = 4
        v_["T"] = 6
        v_["D11"] = 1
        v_["D12"] = 7
        v_["D21"] = 11
        v_["D22"] = 14
        v_["KFRESH"] = 1.25
        v_["FLIM"] = 1.8
        v_["KEFFuINI"] = 0.956145
        v_["ALPHA"] = 6E-6
        v_["CONSPW"] = 364.0
        v_["CYTIME"] = 350.0
        v_["TT"] = Float64(v_["T"])
        v_["T-1"] = -1.0+v_["TT"]
        v_["DELTAT"] = v_["CYTIME"]/v_["T-1"]
        v_["ACC"] = v_["ALPHA"]*v_["CONSPW"]
        v_["ACC"] = v_["ACC"]*v_["DELTAT"]
        v_["-ACC"] = -1.0*v_["ACC"]
        for I = Int64(v_["1"]):Int64(v_["N"])
            v_["V"*string(I)] = 1.0
        end
        v_["V"*string(Int64(v_["D11"]))] = 0.5
        v_["V"*string(Int64(v_["D12"]))] = 0.5
        v_["V"*string(Int64(v_["D21"]))] = 0.5
        v_["V"*string(Int64(v_["D22"]))] = 0.5
        v_["NROW1"] = 2
        v_["ROWu"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))] = 1
        v_["ROWu"*string(Int64(v_["1"]))*","*string(Int64(v_["2"]))] = 2
        v_["G"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))] = 0.705
        v_["G"*string(Int64(v_["1"]))*","*string(Int64(v_["2"]))] = 0.25
        v_["NROW2"] = 5
        v_["ROWu"*string(Int64(v_["2"]))*","*string(Int64(v_["1"]))] = 1
        v_["ROWu"*string(Int64(v_["2"]))*","*string(Int64(v_["2"]))] = 2
        v_["ROWu"*string(Int64(v_["2"]))*","*string(Int64(v_["3"]))] = 3
        v_["ROWu"*string(Int64(v_["2"]))*","*string(Int64(v_["4"]))] = 7
        v_["ROWu"*string(Int64(v_["2"]))*","*string(Int64(v_["5"]))] = 8
        v_["G"*string(Int64(v_["2"]))*","*string(Int64(v_["1"]))] = 0.125
        v_["G"*string(Int64(v_["2"]))*","*string(Int64(v_["2"]))] = 0.625
        v_["G"*string(Int64(v_["2"]))*","*string(Int64(v_["3"]))] = 0.125
        v_["G"*string(Int64(v_["2"]))*","*string(Int64(v_["7"]))] = 0.08
        v_["G"*string(Int64(v_["2"]))*","*string(Int64(v_["8"]))] = 0.045
        v_["NROW3"] = 6
        v_["ROWu"*string(Int64(v_["3"]))*","*string(Int64(v_["1"]))] = 2
        v_["ROWu"*string(Int64(v_["3"]))*","*string(Int64(v_["2"]))] = 3
        v_["ROWu"*string(Int64(v_["3"]))*","*string(Int64(v_["3"]))] = 4
        v_["ROWu"*string(Int64(v_["3"]))*","*string(Int64(v_["4"]))] = 7
        v_["ROWu"*string(Int64(v_["3"]))*","*string(Int64(v_["5"]))] = 8
        v_["ROWu"*string(Int64(v_["3"]))*","*string(Int64(v_["6"]))] = 9
        v_["G"*string(Int64(v_["3"]))*","*string(Int64(v_["2"]))] = 0.125
        v_["G"*string(Int64(v_["3"]))*","*string(Int64(v_["3"]))] = 0.58
        v_["G"*string(Int64(v_["3"]))*","*string(Int64(v_["4"]))] = 0.125
        v_["G"*string(Int64(v_["3"]))*","*string(Int64(v_["7"]))] = 0.045
        v_["G"*string(Int64(v_["3"]))*","*string(Int64(v_["8"]))] = 0.08
        v_["G"*string(Int64(v_["3"]))*","*string(Int64(v_["9"]))] = 0.045
        v_["NROW4"] = 6
        v_["ROWu"*string(Int64(v_["4"]))*","*string(Int64(v_["1"]))] = 3
        v_["ROWu"*string(Int64(v_["4"]))*","*string(Int64(v_["2"]))] = 4
        v_["ROWu"*string(Int64(v_["4"]))*","*string(Int64(v_["3"]))] = 5
        v_["ROWu"*string(Int64(v_["4"]))*","*string(Int64(v_["4"]))] = 8
        v_["ROWu"*string(Int64(v_["4"]))*","*string(Int64(v_["5"]))] = 9
        v_["ROWu"*string(Int64(v_["4"]))*","*string(Int64(v_["6"]))] = 10
        v_["G"*string(Int64(v_["4"]))*","*string(Int64(v_["3"]))] = 0.125
        v_["G"*string(Int64(v_["4"]))*","*string(Int64(v_["4"]))] = 0.58
        v_["G"*string(Int64(v_["4"]))*","*string(Int64(v_["5"]))] = 0.125
        v_["G"*string(Int64(v_["4"]))*","*string(Int64(v_["8"]))] = 0.045
        v_["G"*string(Int64(v_["4"]))*","*string(Int64(v_["9"]))] = 0.08
        v_["G"*string(Int64(v_["4"]))*","*string(Int64(v_["10"]))] = 0.045
        v_["NROW5"] = 5
        v_["ROWu"*string(Int64(v_["5"]))*","*string(Int64(v_["1"]))] = 4
        v_["ROWu"*string(Int64(v_["5"]))*","*string(Int64(v_["2"]))] = 5
        v_["ROWu"*string(Int64(v_["5"]))*","*string(Int64(v_["3"]))] = 6
        v_["ROWu"*string(Int64(v_["5"]))*","*string(Int64(v_["4"]))] = 9
        v_["ROWu"*string(Int64(v_["5"]))*","*string(Int64(v_["5"]))] = 10
        v_["G"*string(Int64(v_["5"]))*","*string(Int64(v_["4"]))] = 0.125
        v_["G"*string(Int64(v_["5"]))*","*string(Int64(v_["5"]))] = 0.58
        v_["G"*string(Int64(v_["5"]))*","*string(Int64(v_["6"]))] = 0.125
        v_["G"*string(Int64(v_["5"]))*","*string(Int64(v_["9"]))] = 0.045
        v_["G"*string(Int64(v_["5"]))*","*string(Int64(v_["10"]))] = 0.08
        v_["NROW6"] = 3
        v_["ROWu"*string(Int64(v_["6"]))*","*string(Int64(v_["1"]))] = 5
        v_["ROWu"*string(Int64(v_["6"]))*","*string(Int64(v_["2"]))] = 6
        v_["ROWu"*string(Int64(v_["6"]))*","*string(Int64(v_["3"]))] = 10
        v_["G"*string(Int64(v_["6"]))*","*string(Int64(v_["5"]))] = 0.125
        v_["G"*string(Int64(v_["6"]))*","*string(Int64(v_["6"]))] = 0.58
        v_["G"*string(Int64(v_["6"]))*","*string(Int64(v_["10"]))] = 0.045
        v_["NROW7"] = 5
        v_["ROWu"*string(Int64(v_["7"]))*","*string(Int64(v_["1"]))] = 1
        v_["ROWu"*string(Int64(v_["7"]))*","*string(Int64(v_["2"]))] = 2
        v_["ROWu"*string(Int64(v_["7"]))*","*string(Int64(v_["3"]))] = 7
        v_["ROWu"*string(Int64(v_["7"]))*","*string(Int64(v_["4"]))] = 8
        v_["ROWu"*string(Int64(v_["7"]))*","*string(Int64(v_["5"]))] = 11
        v_["G"*string(Int64(v_["7"]))*","*string(Int64(v_["1"]))] = 0.045
        v_["G"*string(Int64(v_["7"]))*","*string(Int64(v_["2"]))] = 0.16
        v_["G"*string(Int64(v_["7"]))*","*string(Int64(v_["7"]))] = 0.5
        v_["G"*string(Int64(v_["7"]))*","*string(Int64(v_["8"]))] = 0.16
        v_["G"*string(Int64(v_["7"]))*","*string(Int64(v_["11"]))] = 0.045
        v_["NROW8"] = 8
        v_["ROWu"*string(Int64(v_["8"]))*","*string(Int64(v_["1"]))] = 2
        v_["ROWu"*string(Int64(v_["8"]))*","*string(Int64(v_["2"]))] = 3
        v_["ROWu"*string(Int64(v_["8"]))*","*string(Int64(v_["3"]))] = 4
        v_["ROWu"*string(Int64(v_["8"]))*","*string(Int64(v_["4"]))] = 7
        v_["ROWu"*string(Int64(v_["8"]))*","*string(Int64(v_["5"]))] = 8
        v_["ROWu"*string(Int64(v_["8"]))*","*string(Int64(v_["6"]))] = 9
        v_["ROWu"*string(Int64(v_["8"]))*","*string(Int64(v_["7"]))] = 11
        v_["ROWu"*string(Int64(v_["8"]))*","*string(Int64(v_["8"]))] = 12
        v_["G"*string(Int64(v_["8"]))*","*string(Int64(v_["2"]))] = 0.045
        v_["G"*string(Int64(v_["8"]))*","*string(Int64(v_["3"]))] = 0.08
        v_["G"*string(Int64(v_["8"]))*","*string(Int64(v_["4"]))] = 0.045
        v_["G"*string(Int64(v_["8"]))*","*string(Int64(v_["7"]))] = 0.08
        v_["G"*string(Int64(v_["8"]))*","*string(Int64(v_["8"]))] = 0.545
        v_["G"*string(Int64(v_["8"]))*","*string(Int64(v_["9"]))] = 0.08
        v_["G"*string(Int64(v_["8"]))*","*string(Int64(v_["11"]))] = 0.08
        v_["G"*string(Int64(v_["8"]))*","*string(Int64(v_["12"]))] = 0.045
        v_["NROW9"] = 9
        v_["ROWu"*string(Int64(v_["9"]))*","*string(Int64(v_["1"]))] = 3
        v_["ROWu"*string(Int64(v_["9"]))*","*string(Int64(v_["2"]))] = 4
        v_["ROWu"*string(Int64(v_["9"]))*","*string(Int64(v_["3"]))] = 5
        v_["ROWu"*string(Int64(v_["9"]))*","*string(Int64(v_["4"]))] = 8
        v_["ROWu"*string(Int64(v_["9"]))*","*string(Int64(v_["5"]))] = 9
        v_["ROWu"*string(Int64(v_["9"]))*","*string(Int64(v_["6"]))] = 10
        v_["ROWu"*string(Int64(v_["9"]))*","*string(Int64(v_["7"]))] = 11
        v_["ROWu"*string(Int64(v_["9"]))*","*string(Int64(v_["8"]))] = 12
        v_["ROWu"*string(Int64(v_["9"]))*","*string(Int64(v_["9"]))] = 13
        v_["G"*string(Int64(v_["9"]))*","*string(Int64(v_["3"]))] = 0.045
        v_["G"*string(Int64(v_["9"]))*","*string(Int64(v_["4"]))] = 0.08
        v_["G"*string(Int64(v_["9"]))*","*string(Int64(v_["5"]))] = 0.045
        v_["G"*string(Int64(v_["9"]))*","*string(Int64(v_["8"]))] = 0.08
        v_["G"*string(Int64(v_["9"]))*","*string(Int64(v_["9"]))] = 0.5
        v_["G"*string(Int64(v_["9"]))*","*string(Int64(v_["10"]))] = 0.08
        v_["G"*string(Int64(v_["9"]))*","*string(Int64(v_["11"]))] = 0.045
        v_["G"*string(Int64(v_["9"]))*","*string(Int64(v_["12"]))] = 0.08
        v_["G"*string(Int64(v_["9"]))*","*string(Int64(v_["13"]))] = 0.045
        v_["NROW10"] = 7
        v_["ROWu"*string(Int64(v_["10"]))*","*string(Int64(v_["1"]))] = 4
        v_["ROWu"*string(Int64(v_["10"]))*","*string(Int64(v_["2"]))] = 5
        v_["ROWu"*string(Int64(v_["10"]))*","*string(Int64(v_["3"]))] = 6
        v_["ROWu"*string(Int64(v_["10"]))*","*string(Int64(v_["4"]))] = 9
        v_["ROWu"*string(Int64(v_["10"]))*","*string(Int64(v_["5"]))] = 10
        v_["ROWu"*string(Int64(v_["10"]))*","*string(Int64(v_["6"]))] = 12
        v_["ROWu"*string(Int64(v_["10"]))*","*string(Int64(v_["7"]))] = 13
        v_["G"*string(Int64(v_["10"]))*","*string(Int64(v_["4"]))] = 0.045
        v_["G"*string(Int64(v_["10"]))*","*string(Int64(v_["5"]))] = 0.08
        v_["G"*string(Int64(v_["10"]))*","*string(Int64(v_["6"]))] = 0.045
        v_["G"*string(Int64(v_["10"]))*","*string(Int64(v_["9"]))] = 0.08
        v_["G"*string(Int64(v_["10"]))*","*string(Int64(v_["10"]))] = 0.5
        v_["G"*string(Int64(v_["10"]))*","*string(Int64(v_["12"]))] = 0.045
        v_["G"*string(Int64(v_["10"]))*","*string(Int64(v_["13"]))] = 0.08
        v_["NROW11"] = 6
        v_["ROWu"*string(Int64(v_["11"]))*","*string(Int64(v_["1"]))] = 7
        v_["ROWu"*string(Int64(v_["11"]))*","*string(Int64(v_["2"]))] = 8
        v_["ROWu"*string(Int64(v_["11"]))*","*string(Int64(v_["3"]))] = 9
        v_["ROWu"*string(Int64(v_["11"]))*","*string(Int64(v_["4"]))] = 11
        v_["ROWu"*string(Int64(v_["11"]))*","*string(Int64(v_["5"]))] = 12
        v_["ROWu"*string(Int64(v_["11"]))*","*string(Int64(v_["6"]))] = 14
        v_["G"*string(Int64(v_["11"]))*","*string(Int64(v_["7"]))] = 0.045
        v_["G"*string(Int64(v_["11"]))*","*string(Int64(v_["8"]))] = 0.125
        v_["G"*string(Int64(v_["11"]))*","*string(Int64(v_["9"]))] = 0.045
        v_["G"*string(Int64(v_["11"]))*","*string(Int64(v_["11"]))] = 0.5
        v_["G"*string(Int64(v_["11"]))*","*string(Int64(v_["12"]))] = 0.16
        v_["G"*string(Int64(v_["11"]))*","*string(Int64(v_["14"]))] = 0.045
        v_["NROW12"] = 7
        v_["ROWu"*string(Int64(v_["12"]))*","*string(Int64(v_["1"]))] = 8
        v_["ROWu"*string(Int64(v_["12"]))*","*string(Int64(v_["2"]))] = 9
        v_["ROWu"*string(Int64(v_["12"]))*","*string(Int64(v_["3"]))] = 10
        v_["ROWu"*string(Int64(v_["12"]))*","*string(Int64(v_["4"]))] = 11
        v_["ROWu"*string(Int64(v_["12"]))*","*string(Int64(v_["5"]))] = 12
        v_["ROWu"*string(Int64(v_["12"]))*","*string(Int64(v_["6"]))] = 13
        v_["ROWu"*string(Int64(v_["12"]))*","*string(Int64(v_["7"]))] = 14
        v_["G"*string(Int64(v_["12"]))*","*string(Int64(v_["8"]))] = 0.045
        v_["G"*string(Int64(v_["12"]))*","*string(Int64(v_["9"]))] = 0.08
        v_["G"*string(Int64(v_["12"]))*","*string(Int64(v_["10"]))] = 0.045
        v_["G"*string(Int64(v_["12"]))*","*string(Int64(v_["11"]))] = 0.08
        v_["G"*string(Int64(v_["12"]))*","*string(Int64(v_["12"]))] = 0.545
        v_["G"*string(Int64(v_["12"]))*","*string(Int64(v_["13"]))] = 0.08
        v_["G"*string(Int64(v_["12"]))*","*string(Int64(v_["14"]))] = 0.08
        v_["NROW13"] = 5
        v_["ROWu"*string(Int64(v_["13"]))*","*string(Int64(v_["1"]))] = 9
        v_["ROWu"*string(Int64(v_["13"]))*","*string(Int64(v_["2"]))] = 10
        v_["ROWu"*string(Int64(v_["13"]))*","*string(Int64(v_["3"]))] = 12
        v_["ROWu"*string(Int64(v_["13"]))*","*string(Int64(v_["4"]))] = 13
        v_["ROWu"*string(Int64(v_["13"]))*","*string(Int64(v_["5"]))] = 14
        v_["G"*string(Int64(v_["13"]))*","*string(Int64(v_["9"]))] = 0.045
        v_["G"*string(Int64(v_["13"]))*","*string(Int64(v_["10"]))] = 0.08
        v_["G"*string(Int64(v_["13"]))*","*string(Int64(v_["12"]))] = 0.08
        v_["G"*string(Int64(v_["13"]))*","*string(Int64(v_["13"]))] = 0.5
        v_["G"*string(Int64(v_["13"]))*","*string(Int64(v_["14"]))] = 0.045
        v_["NROW14"] = 3
        v_["ROWu"*string(Int64(v_["14"]))*","*string(Int64(v_["1"]))] = 11
        v_["ROWu"*string(Int64(v_["14"]))*","*string(Int64(v_["2"]))] = 12
        v_["ROWu"*string(Int64(v_["14"]))*","*string(Int64(v_["3"]))] = 14
        v_["G"*string(Int64(v_["14"]))*","*string(Int64(v_["11"]))] = 0.045
        v_["G"*string(Int64(v_["14"]))*","*string(Int64(v_["12"]))] = 0.125
        v_["G"*string(Int64(v_["14"]))*","*string(Int64(v_["14"]))] = 0.5
        v_["T-1"] = -1+v_["T"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            for K = Int64(v_["1"]):Int64(v_["L"])
                for J = Int64(v_["1"]):Int64(v_["M"])
                    iv,ix_,_ = s2mpj_ii("X"*string(I)*","*string(K)*","*string(J),ix_)
                    arrset(pb.xnames,iv,"X"*string(I)*","*string(K)*","*string(J))
                end
            end
        end
        for S = Int64(v_["1"]):Int64(v_["T"])
            for I = Int64(v_["1"]):Int64(v_["N"])
                iv,ix_,_ = s2mpj_ii("KINF"*string(I)*","*string(S),ix_)
                arrset(pb.xnames,iv,"KINF"*string(I)*","*string(S))
                iv,ix_,_ = s2mpj_ii("PHI"*string(I)*","*string(S),ix_)
                arrset(pb.xnames,iv,"PHI"*string(I)*","*string(S))
            end
            iv,ix_,_ = s2mpj_ii("KEFF"*string(S),ix_)
            arrset(pb.xnames,iv,"KEFF"*string(S))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["KEFF"*string(Int64(v_["T"]))]
        pbm.A[ig,iv] += Float64(-1.0)
        for K = Int64(v_["1"]):Int64(v_["L"])
            for J = Int64(v_["1"]):Int64(v_["M"])
                for I = Int64(v_["1"]):Int64(v_["N"])
                    ig,ig_,_ = s2mpj_ii("SUMI"*string(K)*","*string(J),ig_)
                    arrset(gtype,ig,"==")
                    arrset(pb.cnames,ig,"SUMI"*string(K)*","*string(J))
                    iv = ix_["X"*string(I)*","*string(K)*","*string(J)]
                    pbm.A[ig,iv] += Float64(v_["V"*string(I)])
                end
            end
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            for K = Int64(v_["1"]):Int64(v_["L"])
                for J = Int64(v_["1"]):Int64(v_["M"])
                    ig,ig_,_ = s2mpj_ii("SUMLM"*string(I),ig_)
                    arrset(gtype,ig,"==")
                    arrset(pb.cnames,ig,"SUMLM"*string(I))
                    iv = ix_["X"*string(I)*","*string(K)*","*string(J)]
                    pbm.A[ig,iv] += Float64(1.0)
                end
            end
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            ig,ig_,_ = s2mpj_ii("PLAC"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"PLAC"*string(I))
            iv = ix_["KINF"*string(I)*","*string(Int64(v_["1"]))]
            pbm.A[ig,iv] += Float64(-1.0)
            for J = Int64(v_["1"]):Int64(v_["M"])
                ig,ig_,_ = s2mpj_ii("PLAC"*string(I),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"PLAC"*string(I))
                iv = ix_["X"*string(I)*","*string(Int64(v_["1"]))*","*string(J)]
                pbm.A[ig,iv] += Float64(v_["KFRESH"])
            end
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            for S = Int64(v_["1"]):Int64(v_["T"])
                ig,ig_,_ = s2mpj_ii("KERN"*string(I)*","*string(S),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"KERN"*string(I)*","*string(S))
            end
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            for S = Int64(v_["1"]):Int64(v_["T-1"])
                v_["R"] = 1+S
                ig,ig_,_ = s2mpj_ii("KINFF"*string(I)*","*string(S),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"KINFF"*string(I)*","*string(S))
                iv = ix_["KINF"*string(I)*","*string(Int64(v_["R"]))]
                pbm.A[ig,iv] += Float64(-1.0)
                iv = ix_["KINF"*string(I)*","*string(S)]
                pbm.A[ig,iv] += Float64(1.0)
            end
        end
        for S = Int64(v_["1"]):Int64(v_["T"])
            ig,ig_,_ = s2mpj_ii("CPOW"*string(S),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"CPOW"*string(S))
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            for S = Int64(v_["1"]):Int64(v_["T"])
                ig,ig_,_ = s2mpj_ii("PEAK"*string(I)*","*string(S),ig_)
                arrset(gtype,ig,"<=")
                arrset(pb.cnames,ig,"PEAK"*string(I)*","*string(S))
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
        for K = Int64(v_["1"]):Int64(v_["L"])
            for J = Int64(v_["1"]):Int64(v_["M"])
                pbm.gconst[ig_["SUMI"*string(K)*","*string(J)]] = Float64(1.0)
            end
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            pbm.gconst[ig_["SUMLM"*string(I)]] = Float64(1.0)
        end
        for S = Int64(v_["1"]):Int64(v_["T"])
            pbm.gconst[ig_["CPOW"*string(S)]] = Float64(1.0)
        end
        v_["TEMP"] = 0.0
        for I = Int64(v_["1"]):Int64(v_["N"])
            v_["TEMP"] = v_["TEMP"]+v_["V"*string(I)]
        end
        v_["TEMP"] = v_["FLIM"]/v_["TEMP"]
        for I = Int64(v_["1"]):Int64(v_["N"])
            for S = Int64(v_["1"]):Int64(v_["T"])
                pbm.gconst[ig_["PEAK"*string(I)*","*string(S)]] = Float64(v_["TEMP"])
            end
        end
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        for I = Int64(v_["1"]):Int64(v_["N"])
            for K = Int64(v_["1"]):Int64(v_["L"])
                for J = Int64(v_["1"]):Int64(v_["M"])
                    pb.xupper[ix_["X"*string(I)*","*string(K)*","*string(J)]] = 1.0
                end
            end
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            for S = Int64(v_["1"]):Int64(v_["T"])
                pb.xupper[ix_["KINF"*string(I)*","*string(S)]] = v_["KFRESH"]
            end
        end
        v_["LOuKEFF"] = 0.9
        v_["UPuKEFF"] = 1.5
        v_["TEMP"] = Float64(v_["T"])
        v_["TEMP"] = -0.015*v_["TEMP"]
        v_["LOuKEFF"] = v_["LOuKEFF"]+v_["TEMP"]
        v_["UPuKEFF"] = v_["UPuKEFF"]+v_["TEMP"]
        for S = Int64(v_["1"]):Int64(v_["T"])
            pb.xlower[ix_["KEFF"*string(S)]] = v_["LOuKEFF"]
            pb.xupper[ix_["KEFF"*string(S)]] = v_["UPuKEFF"]
        end
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        v_["R14"] = 14.0
        v_["TEMP"] = v_["KFRESH"]/v_["R14"]
        for S = Int64(v_["1"]):Int64(v_["T"])
            if haskey(ix_,"KEFF"*string(S))
                pb.x0[ix_["KEFF"*string(S)]] = Float64(v_["KEFFuINI"])
            else
                pb.y0[findfirst(x->x==ig_["KEFF"*string(S)],pbm.congrps)]  = (
                      Float64(v_["KEFFuINI"]))
            end
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            for S = Int64(v_["1"]):Int64(v_["T"])
                if haskey(ix_,"KINF"*string(I)*","*string(S))
                    pb.x0[ix_["KINF"*string(I)*","*string(S)]] = Float64(v_["KFRESH"])
                else
                    pb.y0[findfirst(x->x==ig_["KINF"*string(I)*","*string(S)],pbm.congrps)]  = (
                          Float64(v_["KFRESH"]))
                end
                if haskey(ix_,"PHI"*string(I)*","*string(S))
                    pb.x0[ix_["PHI"*string(I)*","*string(S)]] = Float64(v_["TEMP"])
                else
                    pb.y0[findfirst(x->x==ig_["PHI"*string(I)*","*string(S)],pbm.congrps)]  = (
                          Float64(v_["TEMP"]))
                end
            end
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            for K = Int64(v_["1"]):Int64(v_["L"])
                for J = Int64(v_["1"]):Int64(v_["M"])
                    if haskey(ix_,"X"*string(I)*","*string(K)*","*string(J))
                        pb.x0[ix_["X"*string(I)*","*string(K)*","*string(J)]] = Float64(0.5)
                    else
                        pb.y0[findfirst(x->x==ig_["X"*string(I)*","*string(K)*","*string(J)],pbm.congrps)] = Float64(0.5)
                    end
                end
            end
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "en2PROD", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        it,iet_,_ = s2mpj_ii( "en3PROD", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        loaset(elftv,it,3,"V3")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        v_["K"] = 2
        v_["K1"] = -1+v_["K"]
        for J = Int64(v_["1"]):Int64(v_["M"])
            for I = Int64(v_["1"]):Int64(v_["N"])
                for II = Int64(v_["1"]):Int64(v_["N"])
                    ename = "Au"*string(I)*","*string(II)*","*string(J)
                    ie,ie_,_  = s2mpj_ii(ename,ie_)
                    arrset(pbm.elftype,ie,"en3PROD")
                    arrset(ielftype,ie,iet_["en3PROD"])
                    vname = "X"*string(I)*","*string(Int64(v_["K"]))*","*string(J)
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                    vname = "X"*string(II)*","*string(Int64(v_["K1"]))*","*string(J)
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                    vname = "KINF"*string(II)*","*string(Int64(v_["T"]))
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                end
            end
        end
        v_["K"] = 3
        v_["K1"] = -1+v_["K"]
        for J = Int64(v_["1"]):Int64(v_["M"])
            for I = Int64(v_["1"]):Int64(v_["N"])
                for II = Int64(v_["1"]):Int64(v_["N"])
                    ename = "Bu"*string(I)*","*string(II)*","*string(J)
                    ie,ie_,_  = s2mpj_ii(ename,ie_)
                    arrset(pbm.elftype,ie,"en3PROD")
                    arrset(ielftype,ie,iet_["en3PROD"])
                    vname = "X"*string(I)*","*string(Int64(v_["K"]))*","*string(J)
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                    vname = "X"*string(II)*","*string(Int64(v_["K1"]))*","*string(J)
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                    vname = "KINF"*string(II)*","*string(Int64(v_["T"]))
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                end
            end
        end
        v_["K"] = 4
        v_["K1"] = -1+v_["K"]
        for J = Int64(v_["1"]):Int64(v_["M"])
            for I = Int64(v_["1"]):Int64(v_["N"])
                for II = Int64(v_["1"]):Int64(v_["N"])
                    ename = "Cu"*string(I)*","*string(II)*","*string(J)
                    ie,ie_,_  = s2mpj_ii(ename,ie_)
                    arrset(pbm.elftype,ie,"en3PROD")
                    arrset(ielftype,ie,iet_["en3PROD"])
                    vname = "X"*string(I)*","*string(Int64(v_["K"]))*","*string(J)
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                    vname = "X"*string(II)*","*string(Int64(v_["K1"]))*","*string(J)
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                    vname = "KINF"*string(II)*","*string(Int64(v_["T"]))
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                end
            end
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            for S = Int64(v_["1"]):Int64(v_["T"])
                ename = "KTP"*string(I)*","*string(S)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"en2PROD")
                arrset(ielftype,ie,iet_["en2PROD"])
                vname = "KEFF"*string(S)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "PHI"*string(I)*","*string(S)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            for S = Int64(v_["1"]):Int64(v_["T"])
                ename = "P"*string(I)*","*string(S)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"en2PROD")
                arrset(ielftype,ie,iet_["en2PROD"])
                vname = "KINF"*string(I)*","*string(S)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "PHI"*string(I)*","*string(S)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            for J = Int64(v_["1"]):Int64(v_["M"])
                for II = Int64(v_["1"]):Int64(v_["N"])
                    ig = ig_["PLAC"*string(I)]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["Au"*string(I)*","*string(II)*","*string(J)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["V"*string(II)]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["Bu"*string(I)*","*string(II)*","*string(J)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["V"*string(II)]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["Cu"*string(I)*","*string(II)*","*string(J)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["V"*string(II)]))
                end
            end
        end
        for S = Int64(v_["1"]):Int64(v_["T"])
            for I = Int64(v_["1"]):Int64(v_["N"])
                ig = ig_["KERN"*string(I)*","*string(S)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["KTP"*string(I)*","*string(S)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(-1.0))
                v_["TEMP"] = v_["NROW"*string(I)]
                v_["NuROW"] = trunc(Int,v_["TEMP"])
                for II = Int64(v_["1"]):Int64(v_["NuROW"])
                    v_["TEMP"] = v_["ROWu"*string(I)*","*string(II)]
                    v_["III"] = trunc(Int,v_["TEMP"])
                    ig = ig_["KERN"*string(I)*","*string(S)]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["P"*string(Int64(v_["III"]))*","*string(S)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["G"*string(I)*","*string(Int64(v_["III"]))]))
                end
            end
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            for S = Int64(v_["1"]):Int64(v_["T-1"])
                ig = ig_["KINFF"*string(I)*","*string(S)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["P"*string(I)*","*string(S)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-ACC"]))
            end
        end
        for S = Int64(v_["1"]):Int64(v_["T"])
            for I = Int64(v_["1"]):Int64(v_["N"])
                ig = ig_["CPOW"*string(S)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["P"*string(I)*","*string(S)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["V"*string(I)]))
            end
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            for S = Int64(v_["1"]):Int64(v_["T"])
                ig = ig_["PEAK"*string(I)*","*string(S)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["P"*string(I)*","*string(S)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,1.)
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
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
        pb.pbclass = "C-LOR2-MN-342-284"
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

    elseif action == "en2PROD"

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

    elseif action == "en3PROD"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]*EV_[2]*EV_[3]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[2]*EV_[3]
            g_[2] = EV_[1]*EV_[3]
            g_[3] = EV_[1]*EV_[2]
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,2] = EV_[3]
                H_[2,1] = H_[1,2]
                H_[1,3] = EV_[2]
                H_[3,1] = H_[1,3]
                H_[2,3] = EV_[1]
                H_[3,2] = H_[2,3]
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

