function LAKES(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    A problem of water resource management in Canada, which may be 
#    formulated as
#    Min  SUM   SUM  (T(i,j)- R(i,j))^2 + (O(i,j)-R(N+i,j))^2)
#        i=1,N j=1,5 
#    subject to
#    T(i+1,1)-T(i,1)+O(i,1)        =  G(i,1)
#    T(i+1,2)-T(i,2)-O(i,1)+O(i,2) =  G(i,2)
#    T(i+1,3)-T(i,3)-O(i,2)+O(i,3) =  G(i,3)
#    T(i+1,4)-T(i,4)-O(i,3)+O(i,4) =  G(i,4)
#    T(i+1,5)-T(i,5)-O(i,4)+O(i,5) =  G(i,5) 
#    i=1,N and T(N+1,j) = T(1,j)  for j=1,5
#    O(i,2)-a*((T(i,2)/480.8+T(i,3)/4.6)/2-543.4)^2 * 
#    (T(i,2)/480.8-T(i,3)/4.6)^.5=0
#    O(i,3)-b*((T(i,3)/4.6-543.4)^2*(T(i,3)/4.6-T(i,4)/105.15)^0.5) = 0
#    O(i,4)-c*(T(i,4)/105.15-550.11)^2.2 = 0
#    where T(i,j) and O(i,j) are variables, R(i,j) are given and
#    a=.0841168  b=.1280849 and c=0.2605.
#    Extra variables 
#    
#    v(i,2) = T(i,2) / 961.6 + T(i,3) / 9.2 - 543.4
#    w(i,2) = T(i,2) / 480.8 - T(i,3) / 4.6
#    v(i,3) = T(i,3) / 4.6 - 543.4
#    w(i,3) = T(i,3) / 4.6 - T(i,4) / 105.15
#    v(i,4) = T(i,4) / 105.15 - 550.11
#    are introduced so that the nonlinear constraints may be rewritten as
#    O(i,2)-a*v(i,2)^2 * w(i,2)^0.5 = 0 ; w(i,2) > 0
#    O(i,3)-b*v(i,3)^2 * w(i,3)^0.5 = 0 ; w(i,3) > 0
#    O(i,4)-c*v(i,4)^2.2 = 0 ; v(i,4) > 0
#    Source:
#    S Jafar Sadjadi
#    Dept. of Systems Design Engineering
#    University of Waterloo
#    Ontario, N2L 3G1 Canada
# 
#    SIF input: Nick Gould and Jafar Sadjadi, November 1995
# 
#    classification = "C-QOR2-RN-90-78"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "LAKES"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["N"] = 6
        v_["N-1"] = -1+v_["N"]
        v_["N+1"] = 1+v_["N"]
        v_["NN"] = 2*v_["N"]
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
        v_["S1s1"] = 202761.072
        v_["S1s2"] = 277791.816
        v_["S1s3"] = 2636.996
        v_["S1s4"] = 59987.0235
        v_["S1s5"] = 19490.4
        v_["S2s1"] = 202703.646
        v_["S2s2"] = 277849.512
        v_["S2s3"] = 2638.1
        v_["S2s4"] = 59998.59
        v_["S2s5"] = 19555.2
        v_["S3s1"] = 202720.536
        v_["S3s2"] = 277955.288
        v_["S3s3"] = 2639.894
        v_["S3s4"] = 60046.959
        v_["S3s5"] = 19597.6
        v_["S4s1"] = 202808.364
        v_["S4s2"] = 278104.336
        v_["S4s3"] = 2640.906
        v_["S4s4"] = 60074.298
        v_["S4s5"] = 19652.8
        v_["S5s1"] = 202916.46
        v_["S5s2"] = 278224.536
        v_["S5s3"] = 2641.458
        v_["S5s4"] = 60091.122
        v_["S5s5"] = 19708.8
        v_["S6s1"] = 202953.618
        v_["S6s2"] = 278277.424
        v_["S6s3"] = 2641.458
        v_["S6s4"] = 60082.71
        v_["S6s5"] = 19706.4
        v_["O1o1"] = 83.728
        v_["O1o2"] = 174.665
        v_["O1o3"] = 180.539
        v_["O1o4"] = 211.558
        v_["O1o5"] = 232.252
        v_["O2o1"] = 83.789
        v_["O2o2"] = 173.255
        v_["O2o3"] = 179.917
        v_["O2o4"] = 210.585
        v_["O2o5"] = 215.254
        v_["O3o1"] = 82.9160
        v_["O3o2"] = 173.721
        v_["O3o3"] = 182.676
        v_["O3o4"] = 207.838
        v_["O3o5"] = 203.855
        v_["O4o1"] = 80.134
        v_["O4o2"] = 178.654
        v_["O4o3"] = 185.917
        v_["O4o4"] = 206.416
        v_["O4o5"] = 186.308
        v_["O5o1"] = 65.345
        v_["O5o2"] = 188.01
        v_["O5o3"] = 192.568
        v_["O5o4"] = 204.3
        v_["O5o5"] = 201.1
        v_["O6o1"] = 72.005
        v_["O6o2"] = 193.833
        v_["O6o3"] = 196.651
        v_["O6o4"] = 204.25
        v_["O6o5"] = 241.079
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            for J = Int64(v_["1"]):Int64(v_["5"])
                iv,ix_,_ = s2mpj_ii("T"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"T"*string(I)*","*string(J))
                iv,ix_,_ = s2mpj_ii("O"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"O"*string(I)*","*string(J))
            end
            iv,ix_,_ = s2mpj_ii("V"*string(I)*","*string(Int64(v_["2"])),ix_)
            arrset(pb.xnames,iv,"V"*string(I)*","*string(Int64(v_["2"])))
            iv,ix_,_ = s2mpj_ii("W"*string(I)*","*string(Int64(v_["2"])),ix_)
            arrset(pb.xnames,iv,"W"*string(I)*","*string(Int64(v_["2"])))
            iv,ix_,_ = s2mpj_ii("V"*string(I)*","*string(Int64(v_["3"])),ix_)
            arrset(pb.xnames,iv,"V"*string(I)*","*string(Int64(v_["3"])))
            iv,ix_,_ = s2mpj_ii("W"*string(I)*","*string(Int64(v_["3"])),ix_)
            arrset(pb.xnames,iv,"W"*string(I)*","*string(Int64(v_["3"])))
            iv,ix_,_ = s2mpj_ii("V"*string(I)*","*string(Int64(v_["4"])),ix_)
            arrset(pb.xnames,iv,"V"*string(I)*","*string(Int64(v_["4"])))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for i = Int64(v_["1"]):Int64(v_["N"])
            v_["n+i"] = v_["N"]+i
            for j = Int64(v_["1"]):Int64(v_["5"])
                ig,ig_,_ = s2mpj_ii("R"*string(i)*","*string(j),ig_)
                arrset(gtype,ig,"<>")
                iv = ix_["T"*string(i)*","*string(j)]
                pbm.A[ig,iv] += Float64(1.0)
                ig,ig_,_ = s2mpj_ii("R"*string(Int64(v_["n+i"]))*","*string(j),ig_)
                arrset(gtype,ig,"<>")
                iv = ix_["O"*string(i)*","*string(j)]
                pbm.A[ig,iv] += Float64(1.0)
            end
        end
        for i = Int64(v_["1"]):Int64(v_["N-1"])
            v_["k+1"] = 1+i
            for j = Int64(v_["1"]):Int64(v_["5"])
                ig,ig_,_ = s2mpj_ii("G"*string(i)*","*string(j),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"G"*string(i)*","*string(j))
                iv = ix_["T"*string(Int64(v_["k+1"]))*","*string(j)]
                pbm.A[ig,iv] += Float64(1.0)
                iv = ix_["T"*string(i)*","*string(j)]
                pbm.A[ig,iv] += Float64(-1.0)
                iv = ix_["O"*string(i)*","*string(j)]
                pbm.A[ig,iv] += Float64(1.0)
            end
            for j = Int64(v_["2"]):Int64(v_["5"])
                v_["j-1"] = -1+j
                ig,ig_,_ = s2mpj_ii("G"*string(i)*","*string(j),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"G"*string(i)*","*string(j))
                iv = ix_["O"*string(i)*","*string(Int64(v_["j-1"]))]
                pbm.A[ig,iv] += Float64(-1.0)
            end
        end
        for j = Int64(v_["1"]):Int64(v_["5"])
            ig,ig_,_ = s2mpj_ii("G"*string(Int64(v_["N"]))*","*string(j),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"G"*string(Int64(v_["N"]))*","*string(j))
            iv = ix_["T"*string(Int64(v_["1"]))*","*string(j)]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["T"*string(Int64(v_["N"]))*","*string(j)]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("G"*string(Int64(v_["N"]))*","*string(j),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"G"*string(Int64(v_["N"]))*","*string(j))
            iv = ix_["O"*string(Int64(v_["N"]))*","*string(j)]
            pbm.A[ig,iv] += Float64(1.0)
        end
        for j = Int64(v_["2"]):Int64(v_["5"])
            v_["j-1"] = -1+j
            ig,ig_,_ = s2mpj_ii("G"*string(Int64(v_["N"]))*","*string(j),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"G"*string(Int64(v_["N"]))*","*string(j))
            iv = ix_["O"*string(Int64(v_["N"]))*","*string(Int64(v_["j-1"]))]
            pbm.A[ig,iv] += Float64(-1.0)
        end
        for i = Int64(v_["1"]):Int64(v_["N"])
            ig,ig_,_ = s2mpj_ii("A"*string(i)*","*string(Int64(v_["1"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"A"*string(i)*","*string(Int64(v_["1"])))
            iv = ix_["O"*string(i)*","*string(Int64(v_["2"]))]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("A"*string(i)*","*string(Int64(v_["2"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"A"*string(i)*","*string(Int64(v_["2"])))
            iv = ix_["O"*string(i)*","*string(Int64(v_["3"]))]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("A"*string(i)*","*string(Int64(v_["3"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"A"*string(i)*","*string(Int64(v_["3"])))
            iv = ix_["O"*string(i)*","*string(Int64(v_["4"]))]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("V"*string(i)*","*string(Int64(v_["2"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"V"*string(i)*","*string(Int64(v_["2"])))
            iv = ix_["V"*string(i)*","*string(Int64(v_["2"]))]
            pbm.A[ig,iv] += Float64(-1.0)
            v_["C"] = 961.6
            v_["C"] = 1.0/v_["C"]
            ig,ig_,_ = s2mpj_ii("V"*string(i)*","*string(Int64(v_["2"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"V"*string(i)*","*string(Int64(v_["2"])))
            iv = ix_["T"*string(i)*","*string(Int64(v_["2"]))]
            pbm.A[ig,iv] += Float64(v_["C"])
            v_["C"] = 9.2
            v_["C"] = 1.0/v_["C"]
            ig,ig_,_ = s2mpj_ii("V"*string(i)*","*string(Int64(v_["2"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"V"*string(i)*","*string(Int64(v_["2"])))
            iv = ix_["T"*string(i)*","*string(Int64(v_["3"]))]
            pbm.A[ig,iv] += Float64(v_["C"])
            ig,ig_,_ = s2mpj_ii("W"*string(i)*","*string(Int64(v_["2"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"W"*string(i)*","*string(Int64(v_["2"])))
            iv = ix_["W"*string(i)*","*string(Int64(v_["2"]))]
            pbm.A[ig,iv] += Float64(-1.0)
            v_["C"] = 480.8
            v_["C"] = 1.0/v_["C"]
            ig,ig_,_ = s2mpj_ii("W"*string(i)*","*string(Int64(v_["2"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"W"*string(i)*","*string(Int64(v_["2"])))
            iv = ix_["T"*string(i)*","*string(Int64(v_["2"]))]
            pbm.A[ig,iv] += Float64(v_["C"])
            v_["C"] = -4.6
            v_["C"] = 1.0/v_["C"]
            ig,ig_,_ = s2mpj_ii("W"*string(i)*","*string(Int64(v_["2"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"W"*string(i)*","*string(Int64(v_["2"])))
            iv = ix_["T"*string(i)*","*string(Int64(v_["3"]))]
            pbm.A[ig,iv] += Float64(v_["C"])
            ig,ig_,_ = s2mpj_ii("V"*string(i)*","*string(Int64(v_["3"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"V"*string(i)*","*string(Int64(v_["3"])))
            iv = ix_["V"*string(i)*","*string(Int64(v_["3"]))]
            pbm.A[ig,iv] += Float64(-1.0)
            v_["C"] = 4.6
            v_["C"] = 1.0/v_["C"]
            ig,ig_,_ = s2mpj_ii("V"*string(i)*","*string(Int64(v_["3"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"V"*string(i)*","*string(Int64(v_["3"])))
            iv = ix_["T"*string(i)*","*string(Int64(v_["3"]))]
            pbm.A[ig,iv] += Float64(v_["C"])
            ig,ig_,_ = s2mpj_ii("W"*string(i)*","*string(Int64(v_["3"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"W"*string(i)*","*string(Int64(v_["3"])))
            iv = ix_["W"*string(i)*","*string(Int64(v_["3"]))]
            pbm.A[ig,iv] += Float64(-1.0)
            v_["C"] = 4.6
            v_["C"] = 1.0/v_["C"]
            ig,ig_,_ = s2mpj_ii("W"*string(i)*","*string(Int64(v_["3"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"W"*string(i)*","*string(Int64(v_["3"])))
            iv = ix_["T"*string(i)*","*string(Int64(v_["3"]))]
            pbm.A[ig,iv] += Float64(v_["C"])
            v_["C"] = -105.15
            v_["C"] = 1.0/v_["C"]
            ig,ig_,_ = s2mpj_ii("W"*string(i)*","*string(Int64(v_["3"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"W"*string(i)*","*string(Int64(v_["3"])))
            iv = ix_["T"*string(i)*","*string(Int64(v_["4"]))]
            pbm.A[ig,iv] += Float64(v_["C"])
            ig,ig_,_ = s2mpj_ii("V"*string(i)*","*string(Int64(v_["4"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"V"*string(i)*","*string(Int64(v_["4"])))
            iv = ix_["V"*string(i)*","*string(Int64(v_["4"]))]
            pbm.A[ig,iv] += Float64(-1.0)
            v_["C"] = 105.15
            v_["C"] = 1.0/v_["C"]
            ig,ig_,_ = s2mpj_ii("V"*string(i)*","*string(Int64(v_["4"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"V"*string(i)*","*string(Int64(v_["4"])))
            iv = ix_["T"*string(i)*","*string(Int64(v_["4"]))]
            pbm.A[ig,iv] += Float64(v_["C"])
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
        pbm.gconst[ig_["R"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))]]  = (
              Float64(v_["S1s1"]))
        pbm.gconst[ig_["R"*string(Int64(v_["1"]))*","*string(Int64(v_["2"]))]]  = (
              Float64(v_["S1s2"]))
        pbm.gconst[ig_["R"*string(Int64(v_["1"]))*","*string(Int64(v_["3"]))]]  = (
              Float64(v_["S1s3"]))
        pbm.gconst[ig_["R"*string(Int64(v_["1"]))*","*string(Int64(v_["4"]))]]  = (
              Float64(v_["S1s4"]))
        pbm.gconst[ig_["R"*string(Int64(v_["1"]))*","*string(Int64(v_["5"]))]]  = (
              Float64(v_["S1s5"]))
        pbm.gconst[ig_["R"*string(Int64(v_["2"]))*","*string(Int64(v_["1"]))]]  = (
              Float64(v_["S2s1"]))
        pbm.gconst[ig_["R"*string(Int64(v_["2"]))*","*string(Int64(v_["2"]))]]  = (
              Float64(v_["S2s2"]))
        pbm.gconst[ig_["R"*string(Int64(v_["2"]))*","*string(Int64(v_["3"]))]]  = (
              Float64(v_["S2s3"]))
        pbm.gconst[ig_["R"*string(Int64(v_["2"]))*","*string(Int64(v_["4"]))]]  = (
              Float64(v_["S2s4"]))
        pbm.gconst[ig_["R"*string(Int64(v_["2"]))*","*string(Int64(v_["5"]))]]  = (
              Float64(v_["S2s5"]))
        pbm.gconst[ig_["R"*string(Int64(v_["3"]))*","*string(Int64(v_["1"]))]]  = (
              Float64(v_["S3s1"]))
        pbm.gconst[ig_["R"*string(Int64(v_["3"]))*","*string(Int64(v_["2"]))]]  = (
              Float64(v_["S3s2"]))
        pbm.gconst[ig_["R"*string(Int64(v_["3"]))*","*string(Int64(v_["3"]))]]  = (
              Float64(v_["S3s3"]))
        pbm.gconst[ig_["R"*string(Int64(v_["3"]))*","*string(Int64(v_["4"]))]]  = (
              Float64(v_["S3s4"]))
        pbm.gconst[ig_["R"*string(Int64(v_["3"]))*","*string(Int64(v_["5"]))]]  = (
              Float64(v_["S3s5"]))
        pbm.gconst[ig_["R"*string(Int64(v_["4"]))*","*string(Int64(v_["1"]))]]  = (
              Float64(v_["S4s1"]))
        pbm.gconst[ig_["R"*string(Int64(v_["4"]))*","*string(Int64(v_["2"]))]]  = (
              Float64(v_["S4s2"]))
        pbm.gconst[ig_["R"*string(Int64(v_["4"]))*","*string(Int64(v_["3"]))]]  = (
              Float64(v_["S4s3"]))
        pbm.gconst[ig_["R"*string(Int64(v_["4"]))*","*string(Int64(v_["4"]))]]  = (
              Float64(v_["S4s4"]))
        pbm.gconst[ig_["R"*string(Int64(v_["4"]))*","*string(Int64(v_["5"]))]]  = (
              Float64(v_["S4s5"]))
        pbm.gconst[ig_["R"*string(Int64(v_["5"]))*","*string(Int64(v_["1"]))]]  = (
              Float64(v_["S5s1"]))
        pbm.gconst[ig_["R"*string(Int64(v_["5"]))*","*string(Int64(v_["2"]))]]  = (
              Float64(v_["S5s2"]))
        pbm.gconst[ig_["R"*string(Int64(v_["5"]))*","*string(Int64(v_["3"]))]]  = (
              Float64(v_["S5s3"]))
        pbm.gconst[ig_["R"*string(Int64(v_["5"]))*","*string(Int64(v_["4"]))]]  = (
              Float64(v_["S5s4"]))
        pbm.gconst[ig_["R"*string(Int64(v_["5"]))*","*string(Int64(v_["5"]))]]  = (
              Float64(v_["S5s5"]))
        pbm.gconst[ig_["R"*string(Int64(v_["6"]))*","*string(Int64(v_["1"]))]]  = (
              Float64(v_["S6s1"]))
        pbm.gconst[ig_["R"*string(Int64(v_["6"]))*","*string(Int64(v_["2"]))]]  = (
              Float64(v_["S6s2"]))
        pbm.gconst[ig_["R"*string(Int64(v_["6"]))*","*string(Int64(v_["3"]))]]  = (
              Float64(v_["S6s3"]))
        pbm.gconst[ig_["R"*string(Int64(v_["6"]))*","*string(Int64(v_["4"]))]]  = (
              Float64(v_["S6s4"]))
        pbm.gconst[ig_["R"*string(Int64(v_["6"]))*","*string(Int64(v_["5"]))]]  = (
              Float64(v_["S6s5"]))
        pbm.gconst[ig_["R"*string(Int64(v_["7"]))*","*string(Int64(v_["1"]))]]  = (
              Float64(v_["O1o5"]))
        pbm.gconst[ig_["R"*string(Int64(v_["7"]))*","*string(Int64(v_["2"]))]]  = (
              Float64(v_["O1o2"]))
        pbm.gconst[ig_["R"*string(Int64(v_["7"]))*","*string(Int64(v_["3"]))]]  = (
              Float64(v_["O1o3"]))
        pbm.gconst[ig_["R"*string(Int64(v_["7"]))*","*string(Int64(v_["4"]))]]  = (
              Float64(v_["O1o4"]))
        pbm.gconst[ig_["R"*string(Int64(v_["7"]))*","*string(Int64(v_["5"]))]]  = (
              Float64(v_["O1o5"]))
        pbm.gconst[ig_["R"*string(Int64(v_["8"]))*","*string(Int64(v_["1"]))]]  = (
              Float64(v_["O2o1"]))
        pbm.gconst[ig_["R"*string(Int64(v_["8"]))*","*string(Int64(v_["2"]))]]  = (
              Float64(v_["O2o2"]))
        pbm.gconst[ig_["R"*string(Int64(v_["8"]))*","*string(Int64(v_["3"]))]]  = (
              Float64(v_["O2o3"]))
        pbm.gconst[ig_["R"*string(Int64(v_["8"]))*","*string(Int64(v_["4"]))]]  = (
              Float64(v_["O2o4"]))
        pbm.gconst[ig_["R"*string(Int64(v_["8"]))*","*string(Int64(v_["5"]))]]  = (
              Float64(v_["O2o5"]))
        pbm.gconst[ig_["R"*string(Int64(v_["9"]))*","*string(Int64(v_["1"]))]]  = (
              Float64(v_["O3o1"]))
        pbm.gconst[ig_["R"*string(Int64(v_["9"]))*","*string(Int64(v_["2"]))]]  = (
              Float64(v_["O3o2"]))
        pbm.gconst[ig_["R"*string(Int64(v_["9"]))*","*string(Int64(v_["3"]))]]  = (
              Float64(v_["O3o3"]))
        pbm.gconst[ig_["R"*string(Int64(v_["9"]))*","*string(Int64(v_["4"]))]]  = (
              Float64(v_["O3o4"]))
        pbm.gconst[ig_["R"*string(Int64(v_["9"]))*","*string(Int64(v_["5"]))]]  = (
              Float64(v_["O3o5"]))
        pbm.gconst[ig_["R"*string(Int64(v_["10"]))*","*string(Int64(v_["1"]))]]  = (
              Float64(v_["O4o1"]))
        pbm.gconst[ig_["R"*string(Int64(v_["10"]))*","*string(Int64(v_["2"]))]]  = (
              Float64(v_["O4o2"]))
        pbm.gconst[ig_["R"*string(Int64(v_["10"]))*","*string(Int64(v_["3"]))]]  = (
              Float64(v_["O4o3"]))
        pbm.gconst[ig_["R"*string(Int64(v_["10"]))*","*string(Int64(v_["4"]))]]  = (
              Float64(v_["O4o4"]))
        pbm.gconst[ig_["R"*string(Int64(v_["10"]))*","*string(Int64(v_["5"]))]]  = (
              Float64(v_["O4o5"]))
        pbm.gconst[ig_["R"*string(Int64(v_["11"]))*","*string(Int64(v_["1"]))]]  = (
              Float64(v_["O5o1"]))
        pbm.gconst[ig_["R"*string(Int64(v_["11"]))*","*string(Int64(v_["2"]))]]  = (
              Float64(v_["O5o2"]))
        pbm.gconst[ig_["R"*string(Int64(v_["11"]))*","*string(Int64(v_["3"]))]]  = (
              Float64(v_["O5o3"]))
        pbm.gconst[ig_["R"*string(Int64(v_["11"]))*","*string(Int64(v_["4"]))]]  = (
              Float64(v_["O5o4"]))
        pbm.gconst[ig_["R"*string(Int64(v_["11"]))*","*string(Int64(v_["5"]))]]  = (
              Float64(v_["O5o5"]))
        pbm.gconst[ig_["R"*string(Int64(v_["12"]))*","*string(Int64(v_["1"]))]]  = (
              Float64(v_["O6o1"]))
        pbm.gconst[ig_["R"*string(Int64(v_["12"]))*","*string(Int64(v_["2"]))]]  = (
              Float64(v_["O6o2"]))
        pbm.gconst[ig_["R"*string(Int64(v_["12"]))*","*string(Int64(v_["3"]))]]  = (
              Float64(v_["O6o3"]))
        pbm.gconst[ig_["R"*string(Int64(v_["12"]))*","*string(Int64(v_["4"]))]]  = (
              Float64(v_["O6o4"]))
        pbm.gconst[ig_["R"*string(Int64(v_["12"]))*","*string(Int64(v_["5"]))]]  = (
              Float64(v_["O6o5"]))
        pbm.gconst[ig_["G"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))]]  = (
              Float64(-22.0))
        pbm.gconst[ig_["G"*string(Int64(v_["1"]))*","*string(Int64(v_["2"]))]]  = (
              Float64(-1.0))
        pbm.gconst[ig_["G"*string(Int64(v_["1"]))*","*string(Int64(v_["3"]))]]  = (
              Float64(3.0))
        pbm.gconst[ig_["G"*string(Int64(v_["1"]))*","*string(Int64(v_["4"]))]]  = (
              Float64(-27.2))
        pbm.gconst[ig_["G"*string(Int64(v_["1"]))*","*string(Int64(v_["5"]))]]  = (
              Float64(51.5))
        pbm.gconst[ig_["G"*string(Int64(v_["2"]))*","*string(Int64(v_["1"]))]]  = (
              Float64(44.0))
        pbm.gconst[ig_["G"*string(Int64(v_["2"]))*","*string(Int64(v_["2"]))]]  = (
              Float64(162.0))
        pbm.gconst[ig_["G"*string(Int64(v_["2"]))*","*string(Int64(v_["3"]))]]  = (
              Float64(8.0))
        pbm.gconst[ig_["G"*string(Int64(v_["2"]))*","*string(Int64(v_["4"]))]]  = (
              Float64(12.5))
        pbm.gconst[ig_["G"*string(Int64(v_["2"]))*","*string(Int64(v_["5"]))]]  = (
              Float64(53.5))
        pbm.gconst[ig_["G"*string(Int64(v_["3"]))*","*string(Int64(v_["1"]))]]  = (
              Float64(-11.0))
        pbm.gconst[ig_["G"*string(Int64(v_["3"]))*","*string(Int64(v_["2"]))]]  = (
              Float64(60.0))
        pbm.gconst[ig_["G"*string(Int64(v_["3"]))*","*string(Int64(v_["3"]))]]  = (
              Float64(10.0))
        pbm.gconst[ig_["G"*string(Int64(v_["3"]))*","*string(Int64(v_["4"]))]]  = (
              Float64(18.0))
        pbm.gconst[ig_["G"*string(Int64(v_["3"]))*","*string(Int64(v_["5"]))]]  = (
              Float64(39.0))
        pbm.gconst[ig_["G"*string(Int64(v_["4"]))*","*string(Int64(v_["1"]))]]  = (
              Float64(124.0))
        pbm.gconst[ig_["G"*string(Int64(v_["4"]))*","*string(Int64(v_["2"]))]]  = (
              Float64(246.0))
        pbm.gconst[ig_["G"*string(Int64(v_["4"]))*","*string(Int64(v_["3"]))]]  = (
              Float64(6.0))
        pbm.gconst[ig_["G"*string(Int64(v_["4"]))*","*string(Int64(v_["4"]))]]  = (
              Float64(9.7))
        pbm.gconst[ig_["G"*string(Int64(v_["4"]))*","*string(Int64(v_["5"]))]]  = (
              Float64(17.2))
        pbm.gconst[ig_["G"*string(Int64(v_["5"]))*","*string(Int64(v_["1"]))]]  = (
              Float64(127.0))
        pbm.gconst[ig_["G"*string(Int64(v_["5"]))*","*string(Int64(v_["2"]))]]  = (
              Float64(175.0))
        pbm.gconst[ig_["G"*string(Int64(v_["5"]))*","*string(Int64(v_["3"]))]]  = (
              Float64(3.0))
        pbm.gconst[ig_["G"*string(Int64(v_["5"]))*","*string(Int64(v_["4"]))]]  = (
              Float64(10.0))
        pbm.gconst[ig_["G"*string(Int64(v_["5"]))*","*string(Int64(v_["5"]))]]  = (
              Float64(30.2))
        pbm.gconst[ig_["G"*string(Int64(v_["6"]))*","*string(Int64(v_["1"]))]]  = (
              Float64(78.0))
        pbm.gconst[ig_["G"*string(Int64(v_["6"]))*","*string(Int64(v_["2"]))]]  = (
              Float64(156.0))
        pbm.gconst[ig_["G"*string(Int64(v_["6"]))*","*string(Int64(v_["3"]))]]  = (
              Float64(3.0))
        pbm.gconst[ig_["G"*string(Int64(v_["6"]))*","*string(Int64(v_["4"]))]]  = (
              Float64(14.0))
        pbm.gconst[ig_["G"*string(Int64(v_["6"]))*","*string(Int64(v_["5"]))]]  = (
              Float64(23.2))
        for i = Int64(v_["1"]):Int64(v_["N"])
            pbm.gconst[ig_["V"*string(i)*","*string(Int64(v_["2"]))]] = Float64(543.4)
            pbm.gconst[ig_["W"*string(i)*","*string(Int64(v_["2"]))]] = Float64(0.0)
            pbm.gconst[ig_["V"*string(i)*","*string(Int64(v_["3"]))]] = Float64(543.4)
            pbm.gconst[ig_["W"*string(i)*","*string(Int64(v_["3"]))]] = Float64(0.0)
            pbm.gconst[ig_["V"*string(i)*","*string(Int64(v_["4"]))]] = Float64(550.11)
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        for i = Int64(v_["1"]):Int64(v_["N"])
            pb.xlower[ix_["W"*string(i)*","*string(Int64(v_["2"]))]] = 0.0001
            pb.xlower[ix_["W"*string(i)*","*string(Int64(v_["3"]))]] = 0.0001
            pb.xlower[ix_["V"*string(i)*","*string(Int64(v_["4"]))]] = 0.0001
        end
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(1.0),pb.n)
        pb.y0 = fill(Float64(1.0),pb.m)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "en1VAR", iet_)
        loaset(elftv,it,1,"V")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"P")
        it,iet_,_ = s2mpj_ii( "en2VAR", iet_)
        loaset(elftv,it,1,"V")
        loaset(elftv,it,2,"W")
        loaset(elftp,it,1,"P")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for i = Int64(v_["1"]):Int64(v_["N"])
            ename = "B"*string(i)*","*string(Int64(v_["1"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"en2VAR")
            arrset(ielftype,ie,iet_["en2VAR"])
            ename = "B"*string(i)*","*string(Int64(v_["1"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "V"*string(i)*","*string(Int64(v_["2"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(i)*","*string(Int64(v_["1"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "W"*string(i)*","*string(Int64(v_["2"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="W",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(i)*","*string(Int64(v_["1"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            posep = findfirst(x->x=="P",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(0.0841168))
            ename = "B"*string(i)*","*string(Int64(v_["2"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"en2VAR")
            arrset(ielftype,ie,iet_["en2VAR"])
            ename = "B"*string(i)*","*string(Int64(v_["2"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "V"*string(i)*","*string(Int64(v_["3"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(i)*","*string(Int64(v_["2"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "W"*string(i)*","*string(Int64(v_["3"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="W",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(i)*","*string(Int64(v_["2"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            posep = findfirst(x->x=="P",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(0.1280849))
            ename = "B"*string(i)*","*string(Int64(v_["3"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"en1VAR")
            arrset(ielftype,ie,iet_["en1VAR"])
            ename = "B"*string(i)*","*string(Int64(v_["3"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "V"*string(i)*","*string(Int64(v_["4"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(i)*","*string(Int64(v_["3"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            posep = findfirst(x->x=="P",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(0.2605))
        end
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gL2",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for i = Int64(v_["1"]):Int64(v_["N"])
            ig = ig_["A"*string(i)*","*string(Int64(v_["1"]))]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["B"*string(i)*","*string(Int64(v_["1"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
            ig = ig_["A"*string(i)*","*string(Int64(v_["2"]))]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["B"*string(i)*","*string(Int64(v_["2"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
            ig = ig_["A"*string(i)*","*string(Int64(v_["3"]))]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["B"*string(i)*","*string(Int64(v_["3"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        for i = Int64(v_["1"]):Int64(v_["N"])
            v_["n+i"] = v_["N"]+i
            for j = Int64(v_["1"]):Int64(v_["5"])
                ig = ig_["R"*string(i)*","*string(j)]
                arrset(pbm.grftype,ig,"gL2")
                ig = ig_["R"*string(Int64(v_["n+i"]))*","*string(j)]
                arrset(pbm.grftype,ig,"gL2")
            end
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
        pb.pbclass = "C-QOR2-RN-90-78"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2]
        pb.conderlvl  = pbm.conderlvl;
        return pb, pbm


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "en1VAR"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = pbm.elpar[iel_][1]*EV_[1]^2.2
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = pbm.elpar[iel_][1]*2.2*EV_[1]^1.2
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = pbm.elpar[iel_][1]*2.64*EV_[1]^0.2
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "en2VAR"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = pbm.elpar[iel_][1]*EV_[1]^2*EV_[2]^0.5
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = pbm.elpar[iel_][1]*2.0*EV_[1]*EV_[2]^0.5
            g_[2] = pbm.elpar[iel_][1]*0.5*EV_[1]^2/EV_[2]^0.5
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = pbm.elpar[iel_][1]*2.0*EV_[2]^0.5
                H_[1,2] = pbm.elpar[iel_][1]*EV_[1]/EV_[2]^0.5
                H_[2,1] = H_[1,2]
                H_[2,2] = -pbm.elpar[iel_][1]*0.25*EV_[1]^2/EV_[2]^1.5
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    elseif action == "gL2"

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

