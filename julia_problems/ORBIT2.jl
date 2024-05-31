function ORBIT2(action,args...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#        A reformulation of the discretized optimal control problem
#        ORBIT
#        Consider minimising the time needed for a spacecraft to
#        move from one circular orbit around the earth to another.
# 
#        This problem can be put in the following form:
#        let (x1,x2,x3) be the position and (x4,x5,x6) the velocity
#        let (u1,u2,u3) be the driving force vector
#        let q be the time required, then our problem becomes:
#        MINIMISE        q
#        with the equations of motion:
#                        dx1/dt = hv*q*x4
#                        dx2/dt = hv*q*x5
#                        dx3/dt = hv*q*x6
#                        dx4/dt = hv*q*(hg*x1/r^3-hf*u(1)/m)     (1)
#                        dx5/dt = hv*q*(hg*x1/r^3-hf*u(2)/m)             
#                        dx6/dt = hv*q*(hg*x1/r^3-hf*u(3)/m)
#        
#        with    m = m0-hm*q*t   ( - the mass variation)
#                        
#                r = sqrt( x1^2 + x2^2 + x3^2 )  -  dist. from the center
#                                                   of the earth
#        't' is a rescaled time varying between 0 and 1, and
#        'hv,hf,hg,hm' are scaling constants.
#        the driving force is bounded by 
#                        u1^2 + u2^2 + u3^2 <= 1                 (2)
#        (the rather arbitrary no. '1' representing the max. power
#        of the spacecraft).
#        We choose the initial conditions:
#                x1 = x2 = 0   ,   x3 = 1   -  initial position
#                x5 = x6 = 0   ,   x4 = Vorb    -   corresponding orbital
#                                                   speed
#        and the final conditions are:
#        x1^2 + x2^2 + x3^2 = Rf^2  -   final orbit radius
#        x4^2 + x5^2 + x6^2 = Vforb^2 -  corresponding orbital
#                                                 speed
#        x1*x4  + x2*x5 + x3*x6 = 0  -  direction must be parallel
#                                        to the velocity
#        we have chosen the constants hg,hf,hv,hm so that the x,u vectors
#        are of order one. These correspond to an initial orbit at 150 km
#        above the earth's surface, and a final orbit at 250 km above the 
#        earth's surface.
#        The reduction to an NLP is done in that same way as for CAR.SIF
#        The time taken should be:
#                        q = 315 secs
# 
#    SIF input: Thomas Felici, Universite de Nancy, France, 
#               October 1993, with modifications to compute derivatives
#               by Nick Gould
#      SAVEs removed December 3rd 2014
# 
#    classification = "LOR1-RN-V-V"
# 
# *******************************************
# 
#        M = Number of time nodes.
# 
#        Change this for different resolution
# 
# IE M                              3   $-PARAMETER n= 25, m= 18
# IE M                             10   $-PARAMETER n= 88, m= 67
# IE M                             30   $-PARAMETER n=268, m=207 original value
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "ORBIT2"

    if action == "setup"
        pbm          = PBM(name)
        pb           = PB(name)
        pb.sifpbname = "ORBIT2"
        nargin       = length(args)
        pbm.call     = eval( Meta.parse( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        if nargin<1
            v_["M"] = 300;  #  SIF file default value
        else
            v_["M"] = args[1];
        end
        v_["NEQ"] = 6
        v_["NCON"] = 1
        v_["NLIM"] = 3
        v_["NLINEQ"] = 0
        v_["NLINCON"] = 0
        v_["NLINLIM"] = 0
        v_["NTEQ"] = v_["NEQ"]+v_["NLINEQ"]
        v_["NTCON"] = v_["NCON"]+v_["NLINCON"]
        v_["NTLIM"] = v_["NLIM"]+v_["NLINLIM"]
        v_["NX"] = 6
        v_["NU"] = 3
        v_["NQ"] = 1
        v_["1"] = 1
        v_["2"] = 2
        v_["3"] = 3
        v_["4"] = 4
        v_["5"] = 5
        v_["6"] = 6
        v_["M-1"] = v_["M"]-v_["1"]
        v_["RM"] = v_["M-1"]
        for I = Int64(v_["1"]):Int64(v_["M"])
            v_["I-1"] = I-v_["1"]
            v_["RI1"] = v_["I-1"]
            v_["T"*string(I)] = v_["RI1"]/v_["RM"]
        end
        v_["RT"] = 6371.0e0
        v_["RT*RT"] = v_["RT"]*v_["RT"]
        v_["MU"] = 9.81e-3*v_["RT*RT"]
        v_["ONE"] = 1.0e+0
        v_["PI/4"] = atan(v_["ONE"])
        v_["PI"] = 4.0e+0*v_["PI/4"]
        v_["R0"] = 1.5e+2+v_["RT"]
        v_["R0*R0"] = v_["R0"]*v_["R0"]
        v_["MU/R0"] = v_["MU"]/v_["R0"]
        v_["VORB"] = sqrt(v_["MU/R0"])
        v_["RF"] = 2.5e+2+v_["RT"]
        v_["M0"] = 3.0e+3
        v_["QT"] = 3.333e+0
        v_["VG"] = 3.0e+0
        v_["TS"] = 1.0e+0
        v_["HG"] = v_["VORB"]*v_["R0*R0"]
        v_["HG"] = v_["MU"]/v_["HG"]
        v_["HG"] = v_["TS"]*v_["HG"]
        v_["HV"] = v_["VORB"]/v_["R0"]
        v_["HV"] = v_["TS"]*v_["HV"]
        v_["HF"] = v_["VORB"]*v_["M0"]
        v_["HF"] = v_["VG"]/v_["HF"]
        v_["HF"] = v_["QT"]*v_["HF"]
        v_["HF"] = v_["TS"]*v_["HF"]
        v_["HM"] = v_["QT"]/v_["M0"]
        v_["HM"] = v_["TS"]*v_["HM"]
        v_["VF"] = v_["MU"]/v_["RF"]
        v_["VF"] = sqrt(v_["VF"])
        v_["VF"] = v_["VF"]/v_["VORB"]
        v_["VFVF"] = v_["VF"]*v_["VF"]
        v_["RF"] = v_["RF"]/v_["R0"]
        v_["RFRF"] = v_["RF"]*v_["RF"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        xscale  = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["M"])
            for J = Int64(v_["1"]):Int64(v_["NX"])
                iv,ix_,_ = s2x_ii("X"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"X"*string(I)*","*string(J))
            end
        end
        for I = Int64(v_["1"]):Int64(v_["M-1"])
            for J = Int64(v_["1"]):Int64(v_["NU"])
                iv,ix_,_ = s2x_ii("U"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"U"*string(I)*","*string(J))
            end
        end
        for J = Int64(v_["1"]):Int64(v_["NQ"])
            iv,ix_,_ = s2x_ii("Q"*string(J),ix_)
            arrset(pb.xnames,iv,"Q"*string(J))
            arrset(xscale,iv,1.0e+2)
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2x_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["Q"*string(Int64(v_["1"]))]
        pbm.A[ig,iv] += 1.0
        for I = Int64(v_["1"]):Int64(v_["M-1"])
            for J = Int64(v_["1"]):Int64(v_["NTEQ"])
                ig,ig_,_ = s2x_ii("K"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"K"*string(I)*","*string(J))
            end
        end
        for J = Int64(v_["1"]):Int64(v_["NTLIM"])
            ig,ig_,_ = s2x_ii("L"*string(J),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"L"*string(J))
        end
        for I = Int64(v_["1"]):Int64(v_["M"])
            for J = Int64(v_["1"]):Int64(v_["NTCON"])
                ig,ig_,_ = s2x_ii("G"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"<=")
                arrset(pb.cnames,ig,"G"*string(I)*","*string(J))
                arrset(pbm.gscale,ig,1.0e+2)
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
        pbm.congrps = findall(x->x!="<>",gtype)
        pb.nob = ngrp-pb.m
        pbm.objgrps = findall(x->x=="<>",gtype)
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        pbm.gconst[ig_["L"*string(Int64(v_["1"]))]] = v_["RFRF"]
        pbm.gconst[ig_["L"*string(Int64(v_["3"]))]] = v_["VFVF"]
        for I = Int64(v_["1"]):Int64(v_["M"])
            pbm.gconst[ig_["G"*string(I)*","*string(Int64(v_["1"]))]] = 1.0e+0
        end
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower[ix_["Q"*string(Int64(v_["1"]))]] = .10000E+03
        pb.xupper[ix_["Q"*string(Int64(v_["1"]))]] = .40000E+26
        pb.xlower[ix_["X"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))]]  = (
              .00000E+00)
        pb.xupper[ix_["X"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))]]  = (
              .00000E+00)
        pb.xlower[ix_["X"*string(Int64(v_["1"]))*","*string(Int64(v_["2"]))]]  = (
              .00000E+00)
        pb.xupper[ix_["X"*string(Int64(v_["1"]))*","*string(Int64(v_["2"]))]]  = (
              .00000E+00)
        pb.xlower[ix_["X"*string(Int64(v_["1"]))*","*string(Int64(v_["3"]))]]  = (
              .10000E+01)
        pb.xupper[ix_["X"*string(Int64(v_["1"]))*","*string(Int64(v_["3"]))]]  = (
              .10000E+01)
        pb.xlower[ix_["X"*string(Int64(v_["1"]))*","*string(Int64(v_["4"]))]]  = (
              .10000E+01)
        pb.xupper[ix_["X"*string(Int64(v_["1"]))*","*string(Int64(v_["4"]))]]  = (
              .10000E+01)
        pb.xlower[ix_["X"*string(Int64(v_["1"]))*","*string(Int64(v_["5"]))]]  = (
              .00000E+00)
        pb.xupper[ix_["X"*string(Int64(v_["1"]))*","*string(Int64(v_["5"]))]]  = (
              .00000E+00)
        pb.xlower[ix_["X"*string(Int64(v_["1"]))*","*string(Int64(v_["6"]))]]  = (
              .00000E+00)
        pb.xupper[ix_["X"*string(Int64(v_["1"]))*","*string(Int64(v_["6"]))]]  = (
              .00000E+00)
        for I = Int64(v_["2"]):Int64(v_["M"])
            pb.xlower[ix_["X"*string(I)*","*string(Int64(v_["1"]))]] = -.10000E+26
            pb.xupper[ix_["X"*string(I)*","*string(Int64(v_["1"]))]] = .10000E+26
            pb.xlower[ix_["X"*string(I)*","*string(Int64(v_["2"]))]] = -.10000E+26
            pb.xupper[ix_["X"*string(I)*","*string(Int64(v_["2"]))]] = .10000E+26
            pb.xlower[ix_["X"*string(I)*","*string(Int64(v_["3"]))]] = -.10000E+26
            pb.xupper[ix_["X"*string(I)*","*string(Int64(v_["3"]))]] = .10000E+26
            pb.xlower[ix_["X"*string(I)*","*string(Int64(v_["4"]))]] = -.10000E+26
            pb.xupper[ix_["X"*string(I)*","*string(Int64(v_["4"]))]] = .10000E+26
            pb.xlower[ix_["X"*string(I)*","*string(Int64(v_["5"]))]] = -.10000E+26
            pb.xupper[ix_["X"*string(I)*","*string(Int64(v_["5"]))]] = .10000E+26
            pb.xlower[ix_["X"*string(I)*","*string(Int64(v_["6"]))]] = -.10000E+26
            pb.xupper[ix_["X"*string(I)*","*string(Int64(v_["6"]))]] = .10000E+26
        end
        for I = Int64(v_["1"]):Int64(v_["M-1"])
            pb.xlower[ix_["U"*string(I)*","*string(Int64(v_["1"]))]] = -.10000E+01
            pb.xupper[ix_["U"*string(I)*","*string(Int64(v_["1"]))]] = .10000E+01
            pb.xlower[ix_["U"*string(I)*","*string(Int64(v_["2"]))]] = -.10000E+01
            pb.xupper[ix_["U"*string(I)*","*string(Int64(v_["2"]))]] = .10000E+01
            pb.xlower[ix_["U"*string(I)*","*string(Int64(v_["3"]))]] = -.10000E+01
            pb.xupper[ix_["U"*string(I)*","*string(Int64(v_["3"]))]] = .10000E+01
        end
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        if haskey(ix_,"Q"*string(Int64(v_["1"])))
            pb.x0[ix_["Q"*string(Int64(v_["1"]))]] = .10000E+03
        else
            pb.y0[findfirst(x->x==ig_["Q"*string(Int64(v_["1"]))],pbm.congrps)]pb.y0[findfirst(x->x==ig_[.10000E+03],pbm.congrps)]
        end
        for I = Int64(v_["1"]):Int64(v_["M"])
            if haskey(ix_,"X"*string(I)*","*string(Int64(v_["1"])))
                pb.x0[ix_["X"*string(I)*","*string(Int64(v_["1"]))]] = .00000E+00
            else
                pb.y0[findfirst(x->x==ig_["X"*string(I)*","*string(Int64(v_["1"]))],pbm.congrps)]pb.y0[findfirst(x->x==ig_[.00000E+00],pbm.congrps)]
            end
            if haskey(ix_,"X"*string(I)*","*string(Int64(v_["2"])))
                pb.x0[ix_["X"*string(I)*","*string(Int64(v_["2"]))]] = .00000E+00
            else
                pb.y0[findfirst(x->x==ig_["X"*string(I)*","*string(Int64(v_["2"]))],pbm.congrps)]pb.y0[findfirst(x->x==ig_[.00000E+00],pbm.congrps)]
            end
            if haskey(ix_,"X"*string(I)*","*string(Int64(v_["3"])))
                pb.x0[ix_["X"*string(I)*","*string(Int64(v_["3"]))]] = .10000E+01
            else
                pb.y0[findfirst(x->x==ig_["X"*string(I)*","*string(Int64(v_["3"]))],pbm.congrps)]pb.y0[findfirst(x->x==ig_[.10000E+01],pbm.congrps)]
            end
            if haskey(ix_,"X"*string(I)*","*string(Int64(v_["4"])))
                pb.x0[ix_["X"*string(I)*","*string(Int64(v_["4"]))]] = .10000E+01
            else
                pb.y0[findfirst(x->x==ig_["X"*string(I)*","*string(Int64(v_["4"]))],pbm.congrps)]pb.y0[findfirst(x->x==ig_[.10000E+01],pbm.congrps)]
            end
            if haskey(ix_,"X"*string(I)*","*string(Int64(v_["5"])))
                pb.x0[ix_["X"*string(I)*","*string(Int64(v_["5"]))]] = .00000E+00
            else
                pb.y0[findfirst(x->x==ig_["X"*string(I)*","*string(Int64(v_["5"]))],pbm.congrps)]pb.y0[findfirst(x->x==ig_[.00000E+00],pbm.congrps)]
            end
            if haskey(ix_,"X"*string(I)*","*string(Int64(v_["6"])))
                pb.x0[ix_["X"*string(I)*","*string(Int64(v_["6"]))]] = .00000E+00
            else
                pb.y0[findfirst(x->x==ig_["X"*string(I)*","*string(Int64(v_["6"]))],pbm.congrps)]pb.y0[findfirst(x->x==ig_[.00000E+00],pbm.congrps)]
            end
        end
        for I = Int64(v_["1"]):Int64(v_["M-1"])
            if haskey(ix_,"U"*string(I)*","*string(Int64(v_["1"])))
                pb.x0[ix_["U"*string(I)*","*string(Int64(v_["1"]))]] = .10000E+01
            else
                pb.y0[findfirst(x->x==ig_["U"*string(I)*","*string(Int64(v_["1"]))],pbm.congrps)]pb.y0[findfirst(x->x==ig_[.10000E+01],pbm.congrps)]
            end
            if haskey(ix_,"U"*string(I)*","*string(Int64(v_["2"])))
                pb.x0[ix_["U"*string(I)*","*string(Int64(v_["2"]))]] = .10000E+01
            else
                pb.y0[findfirst(x->x==ig_["U"*string(I)*","*string(Int64(v_["2"]))],pbm.congrps)]pb.y0[findfirst(x->x==ig_[.10000E+01],pbm.congrps)]
            end
            if haskey(ix_,"U"*string(I)*","*string(Int64(v_["3"])))
                pb.x0[ix_["U"*string(I)*","*string(Int64(v_["3"]))]] = .10000E+01
            else
                pb.y0[findfirst(x->x==ig_["U"*string(I)*","*string(Int64(v_["3"]))],pbm.congrps)]pb.y0[findfirst(x->x==ig_[.10000E+01],pbm.congrps)]
            end
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2x_ii( "SQR", iet_)
        loaset(elftv,it,1,"X")
        it,iet_,_ = s2x_ii( "PROD", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        it,iet_,_ = s2x_ii( "Ktyp", iet_)
        loaset(elftv,it,1,"XI1")
        loaset(elftv,it,2,"XF1")
        loaset(elftv,it,3,"XI2")
        loaset(elftv,it,4,"XF2")
        loaset(elftv,it,5,"XI3")
        loaset(elftv,it,6,"XF3")
        loaset(elftv,it,7,"XI4")
        loaset(elftv,it,8,"XF4")
        loaset(elftv,it,9,"XI5")
        loaset(elftv,it,10,"XF5")
        loaset(elftv,it,11,"XI6")
        loaset(elftv,it,12,"XF6")
        loaset(elftv,it,13,"UV1")
        loaset(elftv,it,14,"UV2")
        loaset(elftv,it,15,"UV3")
        loaset(elftv,it,16,"QV1")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"TN")
        loaset(elftp,it,2,"TN1")
        loaset(elftp,it,3,"ID")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["M-1"])
            v_["S"] = 1+I
            for J = Int64(v_["1"]):Int64(v_["NTEQ"])
                v_["ReJ"] = J
                ename = "Ke"*string(I)*","*string(J)
                ie,ie_,_  = s2x_ii(ename,ie_)
                arrset(pbm.elftype,ie,"Ktyp")
                arrset(ielftype, ie, iet_["Ktyp"])
                posep = findfirst(x->x=="TN",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,v_["T"*string(I)])
                posep = findfirst(x->x=="TN1",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,v_["T"*string(Int64(v_["S"]))])
                posep = findfirst(x->x=="ID",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,v_["ReJ"])
                vname = "X"*string(I)*","*string(Int64(v_["1"]))
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="XI1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(Int64(v_["S"]))*","*string(Int64(v_["1"]))
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="XF1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(I)*","*string(Int64(v_["2"]))
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="XI2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(Int64(v_["S"]))*","*string(Int64(v_["2"]))
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="XF2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(I)*","*string(Int64(v_["3"]))
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="XI3",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(Int64(v_["S"]))*","*string(Int64(v_["3"]))
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="XF3",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(I)*","*string(Int64(v_["4"]))
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="XI4",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(Int64(v_["S"]))*","*string(Int64(v_["4"]))
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="XF4",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(I)*","*string(Int64(v_["5"]))
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="XI5",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(Int64(v_["S"]))*","*string(Int64(v_["5"]))
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="XF5",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(I)*","*string(Int64(v_["6"]))
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="XI6",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(Int64(v_["S"]))*","*string(Int64(v_["6"]))
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="XF6",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "U"*string(I)*","*string(Int64(v_["1"]))
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="UV1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "U"*string(I)*","*string(Int64(v_["2"]))
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="UV2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "U"*string(I)*","*string(Int64(v_["3"]))
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="UV3",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "Q"*string(Int64(v_["1"]))
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="QV1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
        end
        for J = Int64(v_["1"]):Int64(v_["NTCON"])
            for I = Int64(v_["1"]):Int64(v_["M-1"])
                ename = "Ge"*string(I)*","*string(J)*","*string(Int64(v_["1"]))
                ie,ie_,_  = s2x_ii(ename,ie_)
                arrset(pbm.elftype,ie,"SQR")
                arrset(ielftype, ie, iet_["SQR"])
                ename = "Ge"*string(I)*","*string(J)*","*string(Int64(v_["1"]))
                ie,ie_,_  = s2x_ii(ename,ie_)
                vname = "U"*string(I)*","*string(Int64(v_["1"]))
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="X",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                ename = "Ge"*string(I)*","*string(J)*","*string(Int64(v_["2"]))
                ie,ie_,_  = s2x_ii(ename,ie_)
                arrset(pbm.elftype,ie,"SQR")
                arrset(ielftype, ie, iet_["SQR"])
                ename = "Ge"*string(I)*","*string(J)*","*string(Int64(v_["2"]))
                ie,ie_,_  = s2x_ii(ename,ie_)
                vname = "U"*string(I)*","*string(Int64(v_["2"]))
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="X",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                ename = "Ge"*string(I)*","*string(J)*","*string(Int64(v_["3"]))
                ie,ie_,_  = s2x_ii(ename,ie_)
                arrset(pbm.elftype,ie,"SQR")
                arrset(ielftype, ie, iet_["SQR"])
                ename = "Ge"*string(I)*","*string(J)*","*string(Int64(v_["3"]))
                ie,ie_,_  = s2x_ii(ename,ie_)
                vname = "U"*string(I)*","*string(Int64(v_["3"]))
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="X",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
        end
        for J = Int64(v_["1"]):Int64(v_["NTCON"])
            ename = "Ge"*string(Int64(v_["M"]))*","*string(J)*","*string(Int64(v_["1"]))
            ie,ie_,_  = s2x_ii(ename,ie_)
            arrset(pbm.elftype,ie,"SQR")
            arrset(ielftype, ie, iet_["SQR"])
            ename = "Ge"*string(Int64(v_["M"]))*","*string(J)*","*string(Int64(v_["1"]))
            ie,ie_,_  = s2x_ii(ename,ie_)
            vname = "U"*string(Int64(v_["M-1"]))*","*string(Int64(v_["1"]))
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "Ge"*string(Int64(v_["M"]))*","*string(J)*","*string(Int64(v_["2"]))
            ie,ie_,_  = s2x_ii(ename,ie_)
            arrset(pbm.elftype,ie,"SQR")
            arrset(ielftype, ie, iet_["SQR"])
            ename = "Ge"*string(Int64(v_["M"]))*","*string(J)*","*string(Int64(v_["2"]))
            ie,ie_,_  = s2x_ii(ename,ie_)
            vname = "U"*string(Int64(v_["M-1"]))*","*string(Int64(v_["2"]))
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "Ge"*string(Int64(v_["M"]))*","*string(J)*","*string(Int64(v_["3"]))
            ie,ie_,_  = s2x_ii(ename,ie_)
            arrset(pbm.elftype,ie,"SQR")
            arrset(ielftype, ie, iet_["SQR"])
            ename = "Ge"*string(Int64(v_["M"]))*","*string(J)*","*string(Int64(v_["3"]))
            ie,ie_,_  = s2x_ii(ename,ie_)
            vname = "U"*string(Int64(v_["M-1"]))*","*string(Int64(v_["3"]))
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        ename = "Le"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))
        ie,ie_,_  = s2x_ii(ename,ie_)
        arrset(pbm.elftype,ie,"SQR")
        arrset(ielftype, ie, iet_["SQR"])
        ename = "Le"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))
        ie,ie_,_  = s2x_ii(ename,ie_)
        vname = "X"*string(Int64(v_["M"]))*","*string(Int64(v_["1"]))
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "Le"*string(Int64(v_["1"]))*","*string(Int64(v_["2"]))
        ie,ie_,_  = s2x_ii(ename,ie_)
        arrset(pbm.elftype,ie,"SQR")
        arrset(ielftype, ie, iet_["SQR"])
        ename = "Le"*string(Int64(v_["1"]))*","*string(Int64(v_["2"]))
        ie,ie_,_  = s2x_ii(ename,ie_)
        vname = "X"*string(Int64(v_["M"]))*","*string(Int64(v_["2"]))
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "Le"*string(Int64(v_["1"]))*","*string(Int64(v_["3"]))
        ie,ie_,_  = s2x_ii(ename,ie_)
        arrset(pbm.elftype,ie,"SQR")
        arrset(ielftype, ie, iet_["SQR"])
        ename = "Le"*string(Int64(v_["1"]))*","*string(Int64(v_["3"]))
        ie,ie_,_  = s2x_ii(ename,ie_)
        vname = "X"*string(Int64(v_["M"]))*","*string(Int64(v_["3"]))
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "Le"*string(Int64(v_["2"]))*","*string(Int64(v_["1"]))
        ie,ie_,_  = s2x_ii(ename,ie_)
        arrset(pbm.elftype,ie,"PROD")
        arrset(ielftype, ie, iet_["PROD"])
        ename = "Le"*string(Int64(v_["2"]))*","*string(Int64(v_["1"]))
        ie,ie_,_  = s2x_ii(ename,ie_)
        vname = "X"*string(Int64(v_["M"]))*","*string(Int64(v_["1"]))
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "Le"*string(Int64(v_["2"]))*","*string(Int64(v_["1"]))
        ie,ie_,_  = s2x_ii(ename,ie_)
        vname = "X"*string(Int64(v_["M"]))*","*string(Int64(v_["4"]))
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "Le"*string(Int64(v_["2"]))*","*string(Int64(v_["2"]))
        ie,ie_,_  = s2x_ii(ename,ie_)
        arrset(pbm.elftype,ie,"PROD")
        arrset(ielftype, ie, iet_["PROD"])
        ename = "Le"*string(Int64(v_["2"]))*","*string(Int64(v_["2"]))
        ie,ie_,_  = s2x_ii(ename,ie_)
        vname = "X"*string(Int64(v_["M"]))*","*string(Int64(v_["2"]))
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "Le"*string(Int64(v_["2"]))*","*string(Int64(v_["2"]))
        ie,ie_,_  = s2x_ii(ename,ie_)
        vname = "X"*string(Int64(v_["M"]))*","*string(Int64(v_["5"]))
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "Le"*string(Int64(v_["2"]))*","*string(Int64(v_["3"]))
        ie,ie_,_  = s2x_ii(ename,ie_)
        arrset(pbm.elftype,ie,"PROD")
        arrset(ielftype, ie, iet_["PROD"])
        ename = "Le"*string(Int64(v_["2"]))*","*string(Int64(v_["3"]))
        ie,ie_,_  = s2x_ii(ename,ie_)
        vname = "X"*string(Int64(v_["M"]))*","*string(Int64(v_["3"]))
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "Le"*string(Int64(v_["2"]))*","*string(Int64(v_["3"]))
        ie,ie_,_  = s2x_ii(ename,ie_)
        vname = "X"*string(Int64(v_["M"]))*","*string(Int64(v_["6"]))
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "Le"*string(Int64(v_["3"]))*","*string(Int64(v_["1"]))
        ie,ie_,_  = s2x_ii(ename,ie_)
        arrset(pbm.elftype,ie,"SQR")
        arrset(ielftype, ie, iet_["SQR"])
        ename = "Le"*string(Int64(v_["3"]))*","*string(Int64(v_["1"]))
        ie,ie_,_  = s2x_ii(ename,ie_)
        vname = "X"*string(Int64(v_["M"]))*","*string(Int64(v_["4"]))
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "Le"*string(Int64(v_["3"]))*","*string(Int64(v_["2"]))
        ie,ie_,_  = s2x_ii(ename,ie_)
        arrset(pbm.elftype,ie,"SQR")
        arrset(ielftype, ie, iet_["SQR"])
        ename = "Le"*string(Int64(v_["3"]))*","*string(Int64(v_["2"]))
        ie,ie_,_  = s2x_ii(ename,ie_)
        vname = "X"*string(Int64(v_["M"]))*","*string(Int64(v_["5"]))
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "Le"*string(Int64(v_["3"]))*","*string(Int64(v_["3"]))
        ie,ie_,_  = s2x_ii(ename,ie_)
        arrset(pbm.elftype,ie,"SQR")
        arrset(ielftype, ie, iet_["SQR"])
        ename = "Le"*string(Int64(v_["3"]))*","*string(Int64(v_["3"]))
        ie,ie_,_  = s2x_ii(ename,ie_)
        vname = "X"*string(Int64(v_["M"]))*","*string(Int64(v_["6"]))
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["M-1"])
            for J = Int64(v_["1"]):Int64(v_["NTEQ"])
                ig = ig_["K"*string(I)*","*string(J)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["Ke"*string(I)*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,1.)
            end
        end
        for I = Int64(v_["1"]):Int64(v_["M"])
            for J = Int64(v_["1"]):Int64(v_["NTCON"])
                ig = ig_["G"*string(I)*","*string(J)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["Ge"*string(I)*","*string(J)*","*string(Int64(v_["1"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,1.)
                posel = posel+1
                loaset(pbm.grelt,ig,posel,ie_["Ge"*string(I)*","*string(J)*","*string(Int64(v_["2"]))])
                loaset(pbm.grelw,ig,posel,1.)
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["Ge"*string(I)*","*string(J)*","*string(Int64(v_["3"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,1.)
            end
        end
        for I = Int64(v_["1"]):Int64(v_["NTLIM"])
            ig = ig_["L"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["Le"*string(I)*","*string(Int64(v_["1"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["Le"*string(I)*","*string(Int64(v_["2"]))])
            loaset(pbm.grelw,ig,posel,1.)
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["Le"*string(I)*","*string(Int64(v_["3"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.cupper[1:pb.nle] = zeros(Float64,pb.nle)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        #%%%%%%%%%%%%%%% VARIABLES' SCALING %%%%%%%%%%%%%%%
        sA2 = size(pbm.A,2);
        lxs = length(xscale);
        for j = 1:min(sA2,pb.n,length(xscale))
            if xscale[j] != 0.0 && xscale[j] != 1.0
                for i in findall(x->x!=0,pbm.A[:,j])
                      pbm.A[i,j] = pbm.A[i,j]/xscale[j]
                end
            end
        end
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTIONS %%%%%%%
        lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "LOR1-RN-V-V"
        return pb, pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%
