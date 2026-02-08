from s2mpjlib import *
class  BENNETT5(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : BENNETT5
#    *********
# 
#    NIST Data fitting problem BENNETT5 given as an inconsistent set of
#    nonlinear equations.
# 
#    Fit: y = b1 * (b2+x)**(-1/b3) + e
# 
#    Source:  Problem from the NIST nonlinear regression test set
#      http://www.itl.nist.gov/div898/strd/nls/nls_main.shtml
# 
#    Reference:	Bennett, L., L. Swartzendruber, H. Brown, NIST (1994).
#      Superconductivity Magnetization Modeling.
# 
#    SIF input: Nick Gould and Tyrone Rees, Oct 2015
# 
#    classification = "C-CNOR2-MN-3-154"
# 
#    Number of data values
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'BENNETT5'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['M'] = 154
        v_['N'] = 3
        v_['1'] = 1
        v_['X1'] = 7.447168E0
        v_['X2'] = 8.102586E0
        v_['X3'] = 8.452547E0
        v_['X4'] = 8.711278E0
        v_['X5'] = 8.916774E0
        v_['X6'] = 9.087155E0
        v_['X7'] = 9.232590E0
        v_['X8'] = 9.359535E0
        v_['X9'] = 9.472166E0
        v_['X10'] = 9.573384E0
        v_['X11'] = 9.665293E0
        v_['X12'] = 9.749461E0
        v_['X13'] = 9.827092E0
        v_['X14'] = 9.899128E0
        v_['X15'] = 9.966321E0
        v_['X16'] = 10.029280E0
        v_['X17'] = 10.088510E0
        v_['X18'] = 10.144430E0
        v_['X19'] = 10.197380E0
        v_['X20'] = 10.247670E0
        v_['X21'] = 10.295560E0
        v_['X22'] = 10.341250E0
        v_['X23'] = 10.384950E0
        v_['X24'] = 10.426820E0
        v_['X25'] = 10.467000E0
        v_['X26'] = 10.505640E0
        v_['X27'] = 10.542830E0
        v_['X28'] = 10.578690E0
        v_['X29'] = 10.613310E0
        v_['X30'] = 10.646780E0
        v_['X31'] = 10.679150E0
        v_['X32'] = 10.710520E0
        v_['X33'] = 10.740920E0
        v_['X34'] = 10.770440E0
        v_['X35'] = 10.799100E0
        v_['X36'] = 10.826970E0
        v_['X37'] = 10.854080E0
        v_['X38'] = 10.880470E0
        v_['X39'] = 10.906190E0
        v_['X40'] = 10.931260E0
        v_['X41'] = 10.955720E0
        v_['X42'] = 10.979590E0
        v_['X43'] = 11.002910E0
        v_['X44'] = 11.025700E0
        v_['X45'] = 11.047980E0
        v_['X46'] = 11.069770E0
        v_['X47'] = 11.091100E0
        v_['X48'] = 11.111980E0
        v_['X49'] = 11.132440E0
        v_['X50'] = 11.152480E0
        v_['X51'] = 11.172130E0
        v_['X52'] = 11.191410E0
        v_['X53'] = 11.210310E0
        v_['X54'] = 11.228870E0
        v_['X55'] = 11.247090E0
        v_['X56'] = 11.264980E0
        v_['X57'] = 11.282560E0
        v_['X58'] = 11.299840E0
        v_['X59'] = 11.316820E0
        v_['X60'] = 11.333520E0
        v_['X61'] = 11.349940E0
        v_['X62'] = 11.366100E0
        v_['X63'] = 11.382000E0
        v_['X64'] = 11.397660E0
        v_['X65'] = 11.413070E0
        v_['X66'] = 11.428240E0
        v_['X67'] = 11.443200E0
        v_['X68'] = 11.457930E0
        v_['X69'] = 11.472440E0
        v_['X70'] = 11.486750E0
        v_['X71'] = 11.500860E0
        v_['X72'] = 11.514770E0
        v_['X73'] = 11.528490E0
        v_['X74'] = 11.542020E0
        v_['X75'] = 11.555380E0
        v_['X76'] = 11.568550E0
        v_['X77'] = 11.581560E0
        v_['X78'] = 11.594420E0
        v_['X79'] = 11.607121E0
        v_['X80'] = 11.619640E0
        v_['X81'] = 11.632000E0
        v_['X82'] = 11.644210E0
        v_['X83'] = 11.656280E0
        v_['X84'] = 11.668200E0
        v_['X85'] = 11.679980E0
        v_['X86'] = 11.691620E0
        v_['X87'] = 11.703130E0
        v_['X88'] = 11.714510E0
        v_['X89'] = 11.725760E0
        v_['X90'] = 11.736880E0
        v_['X91'] = 11.747890E0
        v_['X92'] = 11.758780E0
        v_['X93'] = 11.769550E0
        v_['X94'] = 11.780200E0
        v_['X95'] = 11.790730E0
        v_['X96'] = 11.801160E0
        v_['X97'] = 11.811480E0
        v_['X98'] = 11.821700E0
        v_['X99'] = 11.831810E0
        v_['X100'] = 11.841820E0
        v_['X101'] = 11.851730E0
        v_['X102'] = 11.861550E0
        v_['X103'] = 11.871270E0
        v_['X104'] = 11.880890E0
        v_['X105'] = 11.890420E0
        v_['X106'] = 11.899870E0
        v_['X107'] = 11.909220E0
        v_['X108'] = 11.918490E0
        v_['X109'] = 11.927680E0
        v_['X110'] = 11.936780E0
        v_['X111'] = 11.945790E0
        v_['X112'] = 11.954730E0
        v_['X113'] = 11.963590E0
        v_['X114'] = 11.972370E0
        v_['X115'] = 11.981070E0
        v_['X116'] = 11.989700E0
        v_['X117'] = 11.998260E0
        v_['X118'] = 12.006740E0
        v_['X119'] = 12.015150E0
        v_['X120'] = 12.023490E0
        v_['X121'] = 12.031760E0
        v_['X122'] = 12.039970E0
        v_['X123'] = 12.048100E0
        v_['X124'] = 12.056170E0
        v_['X125'] = 12.064180E0
        v_['X126'] = 12.072120E0
        v_['X127'] = 12.080010E0
        v_['X128'] = 12.087820E0
        v_['X129'] = 12.095580E0
        v_['X130'] = 12.103280E0
        v_['X131'] = 12.110920E0
        v_['X132'] = 12.118500E0
        v_['X133'] = 12.126030E0
        v_['X134'] = 12.133500E0
        v_['X135'] = 12.140910E0
        v_['X136'] = 12.148270E0
        v_['X137'] = 12.155570E0
        v_['X138'] = 12.162830E0
        v_['X139'] = 12.170030E0
        v_['X140'] = 12.177170E0
        v_['X141'] = 12.184270E0
        v_['X142'] = 12.191320E0
        v_['X143'] = 12.198320E0
        v_['X144'] = 12.205270E0
        v_['X145'] = 12.212170E0
        v_['X146'] = 12.219030E0
        v_['X147'] = 12.225840E0
        v_['X148'] = 12.232600E0
        v_['X149'] = 12.239320E0
        v_['X150'] = 12.245990E0
        v_['X151'] = 12.252620E0
        v_['X152'] = 12.259200E0
        v_['X153'] = 12.265750E0
        v_['X154'] = 12.272240E0
        v_['Y1'] = -34.834702E0
        v_['Y2'] = -34.393200E0
        v_['Y3'] = -34.152901E0
        v_['Y4'] = -33.979099E0
        v_['Y5'] = -33.845901E0
        v_['Y6'] = -33.732899E0
        v_['Y7'] = -33.640301E0
        v_['Y8'] = -33.559200E0
        v_['Y9'] = -33.486801E0
        v_['Y10'] = -33.423100E0
        v_['Y11'] = -33.365101E0
        v_['Y12'] = -33.313000E0
        v_['Y13'] = -33.260899E0
        v_['Y14'] = -33.217400E0
        v_['Y15'] = -33.176899E0
        v_['Y16'] = -33.139198E0
        v_['Y17'] = -33.101601E0
        v_['Y18'] = -33.066799E0
        v_['Y19'] = -33.035000E0
        v_['Y20'] = -33.003101E0
        v_['Y21'] = -32.971298E0
        v_['Y22'] = -32.942299E0
        v_['Y23'] = -32.916302E0
        v_['Y24'] = -32.890202E0
        v_['Y25'] = -32.864101E0
        v_['Y26'] = -32.841000E0
        v_['Y27'] = -32.817799E0
        v_['Y28'] = -32.797501E0
        v_['Y29'] = -32.774300E0
        v_['Y30'] = -32.757000E0
        v_['Y31'] = -32.733799E0
        v_['Y32'] = -32.716400E0
        v_['Y33'] = -32.699100E0
        v_['Y34'] = -32.678799E0
        v_['Y35'] = -32.661400E0
        v_['Y36'] = -32.644001E0
        v_['Y37'] = -32.626701E0
        v_['Y38'] = -32.612202E0
        v_['Y39'] = -32.597698E0
        v_['Y40'] = -32.583199E0
        v_['Y41'] = -32.568699E0
        v_['Y42'] = -32.554298E0
        v_['Y43'] = -32.539799E0
        v_['Y44'] = -32.525299E0
        v_['Y45'] = -32.510799E0
        v_['Y46'] = -32.499199E0
        v_['Y47'] = -32.487598E0
        v_['Y48'] = -32.473202E0
        v_['Y49'] = -32.461601E0
        v_['Y50'] = -32.435501E0
        v_['Y51'] = -32.435501E0
        v_['Y52'] = -32.426800E0
        v_['Y53'] = -32.412300E0
        v_['Y54'] = -32.400799E0
        v_['Y55'] = -32.392101E0
        v_['Y56'] = -32.380501E0
        v_['Y57'] = -32.366001E0
        v_['Y58'] = -32.357300E0
        v_['Y59'] = -32.348598E0
        v_['Y60'] = -32.339901E0
        v_['Y61'] = -32.328400E0
        v_['Y62'] = -32.319698E0
        v_['Y63'] = -32.311001E0
        v_['Y64'] = -32.299400E0
        v_['Y65'] = -32.290699E0
        v_['Y66'] = -32.282001E0
        v_['Y67'] = -32.273300E0
        v_['Y68'] = -32.264599E0
        v_['Y69'] = -32.256001E0
        v_['Y70'] = -32.247299E0
        v_['Y71'] = -32.238602E0
        v_['Y72'] = -32.229900E0
        v_['Y73'] = -32.224098E0
        v_['Y74'] = -32.215401E0
        v_['Y75'] = -32.203800E0
        v_['Y76'] = -32.198002E0
        v_['Y77'] = -32.189400E0
        v_['Y78'] = -32.183601E0
        v_['Y79'] = -32.174900E0
        v_['Y80'] = -32.169102E0
        v_['Y81'] = -32.163300E0
        v_['Y82'] = -32.154598E0
        v_['Y83'] = -32.145901E0
        v_['Y84'] = -32.140099E0
        v_['Y85'] = -32.131401E0
        v_['Y86'] = -32.125599E0
        v_['Y87'] = -32.119801E0
        v_['Y88'] = -32.111198E0
        v_['Y89'] = -32.105400E0
        v_['Y90'] = -32.096699E0
        v_['Y91'] = -32.090900E0
        v_['Y92'] = -32.088001E0
        v_['Y93'] = -32.079300E0
        v_['Y94'] = -32.073502E0
        v_['Y95'] = -32.067699E0
        v_['Y96'] = -32.061901E0
        v_['Y97'] = -32.056099E0
        v_['Y98'] = -32.050301E0
        v_['Y99'] = -32.044498E0
        v_['Y100'] = -32.038799E0
        v_['Y101'] = -32.033001E0
        v_['Y102'] = -32.027199E0
        v_['Y103'] = -32.024300E0
        v_['Y104'] = -32.018501E0
        v_['Y105'] = -32.012699E0
        v_['Y106'] = -32.004002E0
        v_['Y107'] = -32.001099E0
        v_['Y108'] = -31.995300E0
        v_['Y109'] = -31.989500E0
        v_['Y110'] = -31.983700E0
        v_['Y111'] = -31.977900E0
        v_['Y112'] = -31.972099E0
        v_['Y113'] = -31.969299E0
        v_['Y114'] = -31.963501E0
        v_['Y115'] = -31.957701E0
        v_['Y116'] = -31.951900E0
        v_['Y117'] = -31.946100E0
        v_['Y118'] = -31.940300E0
        v_['Y119'] = -31.937401E0
        v_['Y120'] = -31.931601E0
        v_['Y121'] = -31.925800E0
        v_['Y122'] = -31.922899E0
        v_['Y123'] = -31.917101E0
        v_['Y124'] = -31.911301E0
        v_['Y125'] = -31.908400E0
        v_['Y126'] = -31.902599E0
        v_['Y127'] = -31.896900E0
        v_['Y128'] = -31.893999E0
        v_['Y129'] = -31.888201E0
        v_['Y130'] = -31.885300E0
        v_['Y131'] = -31.882401E0
        v_['Y132'] = -31.876600E0
        v_['Y133'] = -31.873699E0
        v_['Y134'] = -31.867901E0
        v_['Y135'] = -31.862101E0
        v_['Y136'] = -31.859200E0
        v_['Y137'] = -31.856300E0
        v_['Y138'] = -31.850500E0
        v_['Y139'] = -31.844700E0
        v_['Y140'] = -31.841801E0
        v_['Y141'] = -31.838900E0
        v_['Y142'] = -31.833099E0
        v_['Y143'] = -31.830200E0
        v_['Y144'] = -31.827299E0
        v_['Y145'] = -31.821600E0
        v_['Y146'] = -31.818701E0
        v_['Y147'] = -31.812901E0
        v_['Y148'] = -31.809999E0
        v_['Y149'] = -31.807100E0
        v_['Y150'] = -31.801300E0
        v_['Y151'] = -31.798401E0
        v_['Y152'] = -31.795500E0
        v_['Y153'] = -31.789700E0
        v_['Y154'] = -31.786800E0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('B'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'B'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['1']),int(v_['M'])+1):
            [ig,ig_,_] = s2mpj_ii('F'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'F'+str(I))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        legrps = np.where(gtype=='<=')[0]
        eqgrps = np.where(gtype=='==')[0]
        gegrps = np.where(gtype=='>=')[0]
        self.nle = len(legrps)
        self.neq = len(eqgrps)
        self.nge = len(gegrps)
        self.m   = self.nle+self.neq+self.nge
        self.congrps = np.concatenate((legrps,eqgrps,gegrps))
        self.cnames = cnames[self.congrps]
        self.nob = ngrp-self.m
        self.objgrps = np.where(gtype=='<>')[0]
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        for I in range(int(v_['1']),int(v_['M'])+1):
            self.gconst = arrset(self.gconst,ig_['F'+str(I)],float(v_['Y'+str(I)]))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        if('B1' in ix_):
            self.x0[ix_['B1']] = float(-2000.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['B1']),float(-2000.0)))
        if('B2' in ix_):
            self.x0[ix_['B2']] = float(50.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['B2']),float(50.0)))
        if('B3' in ix_):
            self.x0[ix_['B3']] = float(0.8)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['B3']),float(0.8)))
        pass
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eE15', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftp = []
        elftp = loaset(elftp,it,0,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for I in range(int(v_['1']),int(v_['M'])+1):
            ename = 'E'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eE15')
            ielftype = arrset(ielftype,ie,iet_["eE15"])
            vname = 'B1'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'B2'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'B3'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V3')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='X')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['X'+str(I)]))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['M'])+1):
            ig = ig_['F'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
        self.objlower = 0.0
#    Solution
# LO SOLTN               
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-CNOR2-MN-3-154"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eE15(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        V3INV = 1.0/EV_[2,0]
        V2PX = EV_[1,0]+self.elpar[iel_][0]
        V2PXL = np.log(V2PX)
        V2PXP = V2PX**V3INV
        V2PXP1 = V2PX**(V3INV+1.0)
        V2PXP2 = V2PX**(V3INV+2.0)
        f_   = EV_[0,0]/V2PXP
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 1.0/V2PXP
            g_[1] = -EV_[0,0]/(EV_[2,0]*V2PXP1)
            g_[2] = EV_[0,0]*V2PXL/(V2PXP*EV_[2,0]**2)
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = -1.0/(EV_[2,0]*V2PXP1)
                H_[1,0] = H_[0,1]
                H_[0,2] = V2PXL/(V2PXP*EV_[2,0]**2)
                H_[2,0] = H_[0,2]
                H_[1,1] = EV_[0,0]*(1.0/EV_[2,0]+1.0)/(EV_[2,0]*V2PXP2)
                H_[1,2]  = (
                      EV_[0,0]/(V2PX*V2PXP*EV_[2,0]**2)-EV_[0,0]*V2PXL/(V2PXP1*EV_[2,0]**3))
                H_[2,1] = H_[1,2]
                H_[2,2] = (EV_[0,0]*V2PXL**2/(V2PXP*EV_[2,0]**4)-2.0*EV_[0,0]*V2PXL/
                     (V2PXP*EV_[2,0]**3))
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

