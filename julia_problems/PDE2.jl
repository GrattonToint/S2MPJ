function PDE2(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : PDE2
#    *********
# 
#    The pde_2, _20 & _200.mod AMPL models from Hans Mittelmann 
#    (mittelmann@asu.edu)
#    See: http://plato.asu.edu/ftp/barrier/
# 
#    SIF input: Nick Gould, April 25th 2012
#               correction by S. Gratton & Ph. Toint, May 2024
# 
#    classification = "C-LLR2-AN-V-V"
# 
#    the x-y discretization 
# 
#       Alternative values for the SIF file parameters:
# IE N                   3              $-PARAMETER
# IE N                   299            $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "PDE2"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        if nargin<1
            v_["N"] = Int64(6);  #  SIF file default value
        else
            v_["N"] = Int64(args[1]);
        end
# IE N                   2999           $-PARAMETER     pde_2.mod value
# IE N                   2099           $-PARAMETER     pde_20.mod value
# IE N                   1299           $-PARAMETER     pde_200.mod value
        v_["0"] = 0
        v_["1"] = 1
        v_["ONE"] = 1.0
        v_["N1"] = 1+v_["N"]
        v_["RN1"] = Float64(v_["N1"])
        v_["A"] = 0.01
        v_["G"] = 20.0
        v_["H"] = v_["ONE"]/v_["RN1"]
        v_["-H"] = -1.0*v_["H"]
        v_["H2"] = v_["H"]*v_["H"]
        v_["GH2"] = v_["G"]*v_["H2"]
        v_["AH"] = v_["A"]*v_["H"]
        v_["SQRTAH"] = sqrt(v_["AH"])
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["0"]):Int64(v_["N1"])
            for J = Int64(v_["0"]):Int64(v_["N1"])
                iv,ix_,_ = s2mpj_ii("X"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"X"*string(I)*","*string(J))
            end
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            for J = Int64(v_["0"]):Int64(v_["N1"])
                iv,ix_,_ = s2mpj_ii("T"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"T"*string(I)*","*string(J))
            end
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            for J = Int64(v_["0"]):Int64(v_["N1"])
                ig,ig_,_ = s2mpj_ii("OBJ",ig_)
                arrset(gtype,ig,"<>")
                iv = ix_["T"*string(I)*","*string(J)]
                pbm.A[ig,iv] += Float64(1.0)
            end
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            v_["I+"] = 1+I
            v_["I-"] = -1+I
            for J = Int64(v_["1"]):Int64(v_["N"])
                v_["J+"] = 1+J
                v_["J-"] = -1+J
                ig,ig_,_ = s2mpj_ii("P"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"P"*string(I)*","*string(J))
                iv = ix_["X"*string(I)*","*string(J)]
                pbm.A[ig,iv] += Float64(4.0)
                iv = ix_["X"*string(I)*","*string(Int64(v_["J+"]))]
                pbm.A[ig,iv] += Float64(-1.0)
                iv = ix_["X"*string(I)*","*string(Int64(v_["J-"]))]
                pbm.A[ig,iv] += Float64(-1.0)
                iv = ix_["X"*string(Int64(v_["I+"]))*","*string(J)]
                pbm.A[ig,iv] += Float64(-1.0)
                iv = ix_["X"*string(Int64(v_["I-"]))*","*string(J)]
                pbm.A[ig,iv] += Float64(-1.0)
            end
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            for J = Int64(v_["1"]):Int64(v_["N"])
                ig,ig_,_ = s2mpj_ii("A"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,">=")
                arrset(pb.cnames,ig,"A"*string(I)*","*string(J))
                iv = ix_["T"*string(I)*","*string(J)]
                pbm.A[ig,iv] += Float64(1.0)
                iv = ix_["X"*string(I)*","*string(J)]
                pbm.A[ig,iv] += Float64(v_["H"])
                ig,ig_,_ = s2mpj_ii("B"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"<=")
                arrset(pb.cnames,ig,"B"*string(I)*","*string(J))
                iv = ix_["T"*string(I)*","*string(J)]
                pbm.A[ig,iv] += Float64(-1.0)
                iv = ix_["X"*string(I)*","*string(J)]
                pbm.A[ig,iv] += Float64(v_["H"])
            end
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            ig,ig_,_ = s2mpj_ii("C"*string(I)*","*string(Int64(v_["0"])),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"C"*string(I)*","*string(Int64(v_["0"])))
            iv = ix_["T"*string(I)*","*string(Int64(v_["0"]))]
            pbm.A[ig,iv] += Float64(1.0)
            ig,ig_,_ = s2mpj_ii("C"*string(I)*","*string(Int64(v_["0"])),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"C"*string(I)*","*string(Int64(v_["0"])))
            iv = ix_["X"*string(I)*","*string(Int64(v_["0"]))]
            pbm.A[ig,iv] += Float64(v_["SQRTAH"])
            ig,ig_,_ = s2mpj_ii("D"*string(I)*","*string(Int64(v_["0"])),ig_)
            arrset(gtype,ig,"<=")
            arrset(pb.cnames,ig,"D"*string(I)*","*string(Int64(v_["0"])))
            iv = ix_["T"*string(I)*","*string(Int64(v_["0"]))]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("D"*string(I)*","*string(Int64(v_["0"])),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"D"*string(I)*","*string(Int64(v_["0"])))
            iv = ix_["X"*string(I)*","*string(Int64(v_["0"]))]
            pbm.A[ig,iv] += Float64(v_["SQRTAH"])
            ig,ig_,_ = s2mpj_ii("C"*string(I)*","*string(Int64(v_["N1"])),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"C"*string(I)*","*string(Int64(v_["N1"])))
            iv = ix_["T"*string(I)*","*string(Int64(v_["0"]))]
            pbm.A[ig,iv] += Float64(1.0)
            ig,ig_,_ = s2mpj_ii("C"*string(I)*","*string(Int64(v_["N1"])),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"C"*string(I)*","*string(Int64(v_["N1"])))
            iv = ix_["X"*string(I)*","*string(Int64(v_["N1"]))]
            pbm.A[ig,iv] += Float64(v_["SQRTAH"])
            ig,ig_,_ = s2mpj_ii("D"*string(I)*","*string(Int64(v_["N1"])),ig_)
            arrset(gtype,ig,"<=")
            arrset(pb.cnames,ig,"D"*string(I)*","*string(Int64(v_["N1"])))
            iv = ix_["T"*string(I)*","*string(Int64(v_["0"]))]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("D"*string(I)*","*string(Int64(v_["N1"])),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"D"*string(I)*","*string(Int64(v_["N1"])))
            iv = ix_["X"*string(I)*","*string(Int64(v_["N1"]))]
            pbm.A[ig,iv] += Float64(v_["SQRTAH"])
            ig,ig_,_ = s2mpj_ii("E"*string(Int64(v_["0"]))*","*string(I),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"E"*string(Int64(v_["0"]))*","*string(I))
            iv = ix_["T"*string(I)*","*string(Int64(v_["0"]))]
            pbm.A[ig,iv] += Float64(1.0)
            ig,ig_,_ = s2mpj_ii("E"*string(Int64(v_["0"]))*","*string(I),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"E"*string(Int64(v_["0"]))*","*string(I))
            iv = ix_["X"*string(Int64(v_["0"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(v_["SQRTAH"])
            ig,ig_,_ = s2mpj_ii("F"*string(Int64(v_["0"]))*","*string(I),ig_)
            arrset(gtype,ig,"<=")
            arrset(pb.cnames,ig,"F"*string(Int64(v_["0"]))*","*string(I))
            iv = ix_["T"*string(I)*","*string(Int64(v_["0"]))]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("F"*string(Int64(v_["0"]))*","*string(I),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"F"*string(Int64(v_["0"]))*","*string(I))
            iv = ix_["X"*string(Int64(v_["0"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(v_["SQRTAH"])
            ig,ig_,_ = s2mpj_ii("E"*string(Int64(v_["N1"]))*","*string(I),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"E"*string(Int64(v_["N1"]))*","*string(I))
            iv = ix_["T"*string(I)*","*string(Int64(v_["0"]))]
            pbm.A[ig,iv] += Float64(1.0)
            ig,ig_,_ = s2mpj_ii("E"*string(Int64(v_["N1"]))*","*string(I),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"E"*string(Int64(v_["N1"]))*","*string(I))
            iv = ix_["X"*string(Int64(v_["N1"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(v_["SQRTAH"])
            ig,ig_,_ = s2mpj_ii("F"*string(Int64(v_["N1"]))*","*string(I),ig_)
            arrset(gtype,ig,"<=")
            arrset(pb.cnames,ig,"F"*string(Int64(v_["N1"]))*","*string(I))
            iv = ix_["T"*string(I)*","*string(Int64(v_["0"]))]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("F"*string(Int64(v_["N1"]))*","*string(I),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"F"*string(Int64(v_["N1"]))*","*string(I))
            iv = ix_["X"*string(Int64(v_["N1"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(v_["SQRTAH"])
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
        for I = Int64(v_["1"]):Int64(v_["N"])
            for J = Int64(v_["1"]):Int64(v_["N"])
                pbm.gconst[ig_["P"*string(I)*","*string(J)]] = Float64(v_["GH2"])
            end
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            v_["RI"] = Float64(I)
            v_["IH"] = v_["RI"]*v_["H"]
            v_["IH-1"] = -1.0+v_["IH"]
            v_["P"] = v_["RI"]*v_["IH-1"]
            v_["P"] = 5.0*v_["P"]
            for J = Int64(v_["1"]):Int64(v_["N"])
                v_["RJ"] = Float64(J)
                v_["JH"] = v_["RJ"]*v_["H"]
                v_["JH-1"] = -1.0+v_["JH"]
                v_["YD"] = v_["RJ"]*v_["JH-1"]
                v_["YD"] = v_["YD"]*v_["P"]
                v_["YD"] = 3.0+v_["YD"]
                v_["YD"] = v_["YD"]*v_["H"]
                v_["-YD"] = -1.0*v_["YD"]
                pbm.gconst[ig_["A"*string(I)*","*string(J)]] = Float64(v_["YD"])
                pbm.gconst[ig_["B"*string(I)*","*string(J)]] = Float64(v_["-YD"])
            end
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = fill(0.0,pb.n)
        pb.xupper = fill(3.5,pb.n)
        for I = Int64(v_["1"]):Int64(v_["N"])
            pb.xupper[ix_["X"*string(I)*","*string(Int64(v_["0"]))]] = 10.0
            pb.xupper[ix_["X"*string(I)*","*string(Int64(v_["N1"]))]] = 10.0
            pb.xupper[ix_["X"*string(Int64(v_["0"]))*","*string(I)]] = 10.0
            pb.xupper[ix_["X"*string(Int64(v_["N1"]))*","*string(I)]] = 10.0
        end
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
        pb.lincons   = collect(1:length(pbm.congrps))
        pb.pbclass = "C-LLR2-AN-V-V"
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

