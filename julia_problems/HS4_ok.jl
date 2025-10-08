function HS4(action,args...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem: HS4
#    *********
# 
#    Source: problem 4 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: A.R. Conn March 1990
# 
#    classification = "OBR2-AN-2-0"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name         = "HS4" #SG string avec "

    if action == "setup"
        pbm      = PBM(name)
        pb       = PB(name)
        pb.sifpbname = "HS4"  #SG string avec "
        nargin   = length(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}(); #SG parenthèse de trop
        ix_ = Dict{String,Int}(); # SG syntaxe
        ig_ = Dict{String,Int}(); # SG syntaxe
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        iv,ix_,_ = s2x_ii("X1",ix_)
        update_list(pb.xnames,iv,"X1")
        iv,ix_,_ = s2x_ii("X2",ix_)
        update_list(pb.xnames,iv,"X2")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%% %SG commentaire
        ig,ig_,_ = s2x_ii("G1",ig_)
        update_list(pbm.grnames,ig,"G1")
        update_list(gtype,ig,"<>")
        iv = ix_["X1"]
        pbm.A[ig,iv] = 1.0+pbm.A[ig,iv]
        ig,ig_,_ = s2x_ii("G1",ig_)
        update_list(pbm.grnames,ig,"G1")
        update_list(gtype,ig,"<>")
        update_list(pbm.gscale,ig,3.0)
        ig,ig_,_ = s2x_ii("G2",ig_)
        update_list(pbm.grnames,ig,"G2")
        update_list(gtype,ig,"<>")
        iv = ix_["X2"]
        pbm.A[ig,iv] = 1.0+pbm.A[ig,iv]
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        pbm.gconst[ig_["G1"]] = -1.0
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower[ix_["X1"]] = 1.0 # SG parenthèse de trop
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.x0[ix_["X1"]] = 1.125
        pb.x0[ix_["X2"]] = 0.125
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%  #SG commentaires
        iet_ = Dict{String,Int}()
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = Dict{String,Int}()
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2x_ii("CUBE",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = [Int[] for _ in 1:ngrp ]
        ig = ig_["G1"]
        update_list(pbm.grftype,ig,"CUBE")
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 2.66
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        Asave = pbm.A
        pbm.A = Asave[1:ngrp, 1:pb.n]
        #%%%%% RETURN VALUES FROM THE SETUP ACTIONS %%%%%%%
        pb.pbclass = "OBR2-AN-2-0"
        return pb, pbm

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    elseif action == "CUBE"

        GVAR_   = args[1]
        igr_    = args[2]
        nargout = args[3]
        f_= GVAR_ * GVAR_ * GVAR_
        if nargout>1
            g_ = 3.0 * GVAR_ * GVAR_
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_ = 6.0 * GVAR_
            end
        end # SG en manquant
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    #%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    elseif action in  ["fx","fgx","fgHx","cx","cJx","cJHx","cIx","cIJx","cIJHx","cIJxv","fHxv","cJxv","Lxy","Lgxy","LgHxy","LIxy","LIgxy","LIgHxy","LHxyv","LIHxyv"]

        pbm = args[1]  # SG pass pb
        if pbm.name == pbm.name  # SG changed to  sifpbname ????
            pbm.has_globs = [0,0]
            return s2x_eval(action,args...)
        else
            println("ERROR: please run "*pb.sifpbname*" with action = setup")
            return (args[end],) # SG ntuple(1->args[end],undef)
        end

    else
        println("ERROR: unknown action " * action * " requested from " * pb.name * ".jl")
        return (args[end],) # SG  ntuple(1->args[end],undef)
    end

end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

