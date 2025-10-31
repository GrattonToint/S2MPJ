from s2mpjlib import *
class  NELSON(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : NELSON
#    *********
# 
#    NIST Data fitting problem NELSON given as an inconsistent set of
#    nonlinear equations.
# 
#    Fit: log[y] = b1 - b2*x1 * exp[-b3*x2] + e
# 
#    Source:  Problem from the NIST nonlinear regression test set
#      http://www.itl.nist.gov/div898/strd/nls/nls_main.shtml
# 
#    Reference: Nelson, W. (1981).  
#      Analysis of Performance-Degradation Data.  
#      IEEE Transactions on Reliability. Vol. 2, R-30, No. 2, pp. 149-155.
# 
#    SIF input: Nick Gould and Tyrone Rees, Oct 2015
# 
#    classification = "C-CNOR2-MN-3-128"
# 
#    Number of data values
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'NELSON'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['M'] = 128
        v_['N'] = 3
        v_['1'] = 1
        v_['X11'] = 1.0
        v_['X12'] = 1.0
        v_['X13'] = 1.0
        v_['X14'] = 1.0
        v_['X15'] = 1.0
        v_['X16'] = 1.0
        v_['X17'] = 1.0
        v_['X18'] = 1.0
        v_['X19'] = 1.0
        v_['X110'] = 1.0
        v_['X111'] = 1.0
        v_['X112'] = 1.0
        v_['X113'] = 1.0
        v_['X114'] = 1.0
        v_['X115'] = 1.0
        v_['X116'] = 1.0
        v_['X117'] = 2.0
        v_['X118'] = 2.0
        v_['X119'] = 2.0
        v_['X120'] = 2.0
        v_['X121'] = 2.0
        v_['X122'] = 2.0
        v_['X123'] = 2.0
        v_['X124'] = 2.0
        v_['X125'] = 2.0
        v_['X126'] = 2.0
        v_['X127'] = 2.0
        v_['X128'] = 2.0
        v_['X129'] = 2.0
        v_['X130'] = 2.0
        v_['X131'] = 2.0
        v_['X132'] = 2.0
        v_['X133'] = 4.0
        v_['X134'] = 4.0
        v_['X135'] = 4.0
        v_['X136'] = 4.0
        v_['X137'] = 4.0
        v_['X138'] = 4.0
        v_['X139'] = 4.0
        v_['X140'] = 4.0
        v_['X141'] = 4.0
        v_['X142'] = 4.0
        v_['X143'] = 4.0
        v_['X144'] = 4.0
        v_['X145'] = 4.0
        v_['X146'] = 4.0
        v_['X147'] = 4.0
        v_['X148'] = 4.0
        v_['X149'] = 8.0
        v_['X150'] = 8.0
        v_['X151'] = 8.0
        v_['X152'] = 8.0
        v_['X153'] = 8.0
        v_['X154'] = 8.0
        v_['X155'] = 8.0
        v_['X156'] = 8.0
        v_['X157'] = 8.0
        v_['X158'] = 8.0
        v_['X159'] = 8.0
        v_['X160'] = 8.0
        v_['X161'] = 8.0
        v_['X162'] = 8.0
        v_['X163'] = 8.0
        v_['X164'] = 8.0
        v_['X165'] = 16.0
        v_['X166'] = 16.0
        v_['X167'] = 16.0
        v_['X168'] = 16.0
        v_['X169'] = 16.0
        v_['X170'] = 16.0
        v_['X171'] = 16.0
        v_['X172'] = 16.0
        v_['X173'] = 16.0
        v_['X174'] = 16.0
        v_['X175'] = 16.0
        v_['X176'] = 16.0
        v_['X177'] = 16.0
        v_['X178'] = 16.0
        v_['X179'] = 16.0
        v_['X180'] = 16.0
        v_['X181'] = 32.0
        v_['X182'] = 32.0
        v_['X183'] = 32.0
        v_['X184'] = 32.0
        v_['X185'] = 32.0
        v_['X186'] = 32.0
        v_['X187'] = 32.0
        v_['X188'] = 32.0
        v_['X189'] = 32.0
        v_['X190'] = 32.0
        v_['X191'] = 32.0
        v_['X192'] = 32.0
        v_['X193'] = 32.0
        v_['X194'] = 32.0
        v_['X195'] = 32.0
        v_['X196'] = 32.0
        v_['X197'] = 48.0
        v_['X198'] = 48.0
        v_['X199'] = 48.0
        v_['X1100'] = 48.0
        v_['X1101'] = 48.0
        v_['X1102'] = 48.0
        v_['X1103'] = 48.0
        v_['X1104'] = 48.0
        v_['X1105'] = 48.0
        v_['X1106'] = 48.0
        v_['X1107'] = 48.0
        v_['X1108'] = 48.0
        v_['X1109'] = 48.0
        v_['X1110'] = 48.0
        v_['X1111'] = 48.0
        v_['X1112'] = 48.0
        v_['X1113'] = 64.0
        v_['X1114'] = 64.0
        v_['X1115'] = 64.0
        v_['X1116'] = 64.0
        v_['X1117'] = 64.0
        v_['X1118'] = 64.0
        v_['X1119'] = 64.0
        v_['X1120'] = 64.0
        v_['X1121'] = 64.0
        v_['X1122'] = 64.0
        v_['X1123'] = 64.0
        v_['X1124'] = 64.0
        v_['X1125'] = 64.0
        v_['X1126'] = 64.0
        v_['X1127'] = 64.0
        v_['X1128'] = 64.0
        v_['X21'] = 180.0
        v_['X22'] = 180.0
        v_['X23'] = 180.0
        v_['X24'] = 180.0
        v_['X25'] = 225.0
        v_['X26'] = 225.0
        v_['X27'] = 225.0
        v_['X28'] = 225.0
        v_['X29'] = 250.0
        v_['X210'] = 250.0
        v_['X211'] = 250.0
        v_['X212'] = 250.0
        v_['X213'] = 275.0
        v_['X214'] = 275.0
        v_['X215'] = 275.0
        v_['X216'] = 275.0
        v_['X217'] = 180.0
        v_['X218'] = 180.0
        v_['X219'] = 180.0
        v_['X220'] = 180.0
        v_['X221'] = 225.0
        v_['X222'] = 225.0
        v_['X223'] = 225.0
        v_['X224'] = 225.0
        v_['X225'] = 250.0
        v_['X226'] = 250.0
        v_['X227'] = 250.0
        v_['X228'] = 250.0
        v_['X229'] = 275.0
        v_['X230'] = 275.0
        v_['X231'] = 275.0
        v_['X232'] = 275.0
        v_['X233'] = 180.0
        v_['X234'] = 180.0
        v_['X235'] = 180.0
        v_['X236'] = 180.0
        v_['X237'] = 225.0
        v_['X238'] = 225.0
        v_['X239'] = 225.0
        v_['X240'] = 225.0
        v_['X241'] = 250.0
        v_['X242'] = 250.0
        v_['X243'] = 250.0
        v_['X244'] = 250.0
        v_['X245'] = 275.0
        v_['X246'] = 275.0
        v_['X247'] = 275.0
        v_['X248'] = 275.0
        v_['X249'] = 180.0
        v_['X250'] = 180.0
        v_['X251'] = 180.0
        v_['X252'] = 180.0
        v_['X253'] = 225.0
        v_['X254'] = 225.0
        v_['X255'] = 225.0
        v_['X256'] = 225.0
        v_['X257'] = 250.0
        v_['X258'] = 250.0
        v_['X259'] = 250.0
        v_['X260'] = 250.0
        v_['X261'] = 275.0
        v_['X262'] = 275.0
        v_['X263'] = 275.0
        v_['X264'] = 275.0
        v_['X265'] = 180.0
        v_['X266'] = 180.0
        v_['X267'] = 180.0
        v_['X268'] = 180.0
        v_['X269'] = 225.0
        v_['X270'] = 225.0
        v_['X271'] = 225.0
        v_['X272'] = 225.0
        v_['X273'] = 250.0
        v_['X274'] = 250.0
        v_['X275'] = 250.0
        v_['X276'] = 250.0
        v_['X277'] = 275.0
        v_['X278'] = 275.0
        v_['X279'] = 275.0
        v_['X280'] = 275.0
        v_['X281'] = 180.0
        v_['X282'] = 180.0
        v_['X283'] = 180.0
        v_['X284'] = 180.0
        v_['X285'] = 225.0
        v_['X286'] = 225.0
        v_['X287'] = 225.0
        v_['X288'] = 225.0
        v_['X289'] = 250.0
        v_['X290'] = 250.0
        v_['X291'] = 250.0
        v_['X292'] = 250.0
        v_['X293'] = 275.0
        v_['X294'] = 275.0
        v_['X295'] = 275.0
        v_['X296'] = 275.0
        v_['X297'] = 180.0
        v_['X298'] = 180.0
        v_['X299'] = 180.0
        v_['X2100'] = 180.0
        v_['X2101'] = 225.0
        v_['X2102'] = 225.0
        v_['X2103'] = 225.0
        v_['X2104'] = 225.0
        v_['X2105'] = 250.0
        v_['X2106'] = 250.0
        v_['X2107'] = 250.0
        v_['X2108'] = 250.0
        v_['X2109'] = 275.0
        v_['X2110'] = 275.0
        v_['X2111'] = 275.0
        v_['X2112'] = 275.0
        v_['X2113'] = 180.0
        v_['X2114'] = 180.0
        v_['X2115'] = 180.0
        v_['X2116'] = 180.0
        v_['X2117'] = 225.0
        v_['X2118'] = 225.0
        v_['X2119'] = 225.0
        v_['X2120'] = 225.0
        v_['X2121'] = 250.0
        v_['X2122'] = 250.0
        v_['X2123'] = 250.0
        v_['X2124'] = 250.0
        v_['X2125'] = 275.0
        v_['X2126'] = 275.0
        v_['X2127'] = 275.0
        v_['X2128'] = 275.0
        v_['Y1'] = 15.00
        v_['Y2'] = 17.00
        v_['Y3'] = 15.50
        v_['Y4'] = 16.50
        v_['Y5'] = 15.50
        v_['Y6'] = 15.00
        v_['Y7'] = 16.00
        v_['Y8'] = 14.50
        v_['Y9'] = 15.00
        v_['Y10'] = 14.50
        v_['Y11'] = 12.50
        v_['Y12'] = 11.00
        v_['Y13'] = 14.00
        v_['Y14'] = 13.00
        v_['Y15'] = 14.00
        v_['Y16'] = 11.50
        v_['Y17'] = 14.00
        v_['Y18'] = 16.00
        v_['Y19'] = 13.00
        v_['Y20'] = 13.50
        v_['Y21'] = 13.00
        v_['Y22'] = 13.50
        v_['Y23'] = 12.50
        v_['Y24'] = 12.50
        v_['Y25'] = 12.50
        v_['Y26'] = 12.00
        v_['Y27'] = 11.50
        v_['Y28'] = 12.00
        v_['Y29'] = 13.00
        v_['Y30'] = 11.50
        v_['Y31'] = 13.00
        v_['Y32'] = 12.50
        v_['Y33'] = 13.50
        v_['Y34'] = 17.50
        v_['Y35'] = 17.50
        v_['Y36'] = 13.50
        v_['Y37'] = 12.50
        v_['Y38'] = 12.50
        v_['Y39'] = 15.00
        v_['Y40'] = 13.00
        v_['Y41'] = 12.00
        v_['Y42'] = 13.00
        v_['Y43'] = 12.00
        v_['Y44'] = 13.50
        v_['Y45'] = 10.00
        v_['Y46'] = 11.50
        v_['Y47'] = 11.00
        v_['Y48'] = 9.50
        v_['Y49'] = 15.00
        v_['Y50'] = 15.00
        v_['Y51'] = 15.50
        v_['Y52'] = 16.00
        v_['Y53'] = 13.00
        v_['Y54'] = 10.50
        v_['Y55'] = 13.50
        v_['Y56'] = 14.00
        v_['Y57'] = 12.50
        v_['Y58'] = 12.00
        v_['Y59'] = 11.50
        v_['Y60'] = 11.50
        v_['Y61'] = 6.50
        v_['Y62'] = 5.50
        v_['Y63'] = 6.00
        v_['Y64'] = 6.00
        v_['Y65'] = 18.50
        v_['Y66'] = 17.00
        v_['Y67'] = 15.30
        v_['Y68'] = 16.00
        v_['Y69'] = 13.00
        v_['Y70'] = 14.00
        v_['Y71'] = 12.50
        v_['Y72'] = 11.00
        v_['Y73'] = 12.00
        v_['Y74'] = 12.00
        v_['Y75'] = 11.50
        v_['Y76'] = 12.00
        v_['Y77'] = 6.00
        v_['Y78'] = 6.00
        v_['Y79'] = 5.00
        v_['Y80'] = 5.50
        v_['Y81'] = 12.50
        v_['Y82'] = 13.00
        v_['Y83'] = 16.00
        v_['Y84'] = 12.00
        v_['Y85'] = 11.00
        v_['Y86'] = 9.50
        v_['Y87'] = 11.00
        v_['Y88'] = 11.00
        v_['Y89'] = 11.00
        v_['Y90'] = 10.00
        v_['Y91'] = 10.50
        v_['Y92'] = 10.50
        v_['Y93'] = 2.70
        v_['Y94'] = 2.70
        v_['Y95'] = 2.50
        v_['Y96'] = 2.40
        v_['Y97'] = 13.00
        v_['Y98'] = 13.50
        v_['Y99'] = 16.50
        v_['Y100'] = 13.60
        v_['Y101'] = 11.50
        v_['Y102'] = 10.50
        v_['Y103'] = 13.50
        v_['Y104'] = 12.00
        v_['Y105'] = 7.00
        v_['Y106'] = 6.90
        v_['Y107'] = 8.80
        v_['Y108'] = 7.90
        v_['Y109'] = 1.20
        v_['Y110'] = 1.50
        v_['Y111'] = 1.00
        v_['Y112'] = 1.50
        v_['Y113'] = 13.00
        v_['Y114'] = 12.50
        v_['Y115'] = 16.50
        v_['Y116'] = 16.00
        v_['Y117'] = 11.00
        v_['Y118'] = 11.50
        v_['Y119'] = 10.50
        v_['Y120'] = 10.00
        v_['Y121'] = 7.27
        v_['Y122'] = 7.50
        v_['Y123'] = 6.70
        v_['Y124'] = 7.60
        v_['Y125'] = 1.50
        v_['Y126'] = 1.00
        v_['Y127'] = 1.20
        v_['Y128'] = 1.20
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
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['B1']])
            valA = np.append(valA,float(1.0))
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
            v_['LOGY'] = np.log(v_['Y'+str(I)])
            self.gconst = arrset(self.gconst,ig_['F'+str(I)],float(v_['LOGY']))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        if('B1' in ix_):
            self.x0[ix_['B1']] = float(2.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['B1']),float(2.0)))
        if('B2' in ix_):
            self.x0[ix_['B2']] = float(0.0001)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['B2']),float(0.0001)))
        if('B3' in ix_):
            self.x0[ix_['B3']] = float(-0.01)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['B3']),float(-0.01)))
        pass
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eE6', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftp = []
        elftp = loaset(elftp,it,0,'X1')
        elftp = loaset(elftp,it,1,'X2')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for I in range(int(v_['1']),int(v_['M'])+1):
            ename = 'E'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eE6')
            ielftype = arrset(ielftype,ie,iet_["eE6"])
            vname = 'B2'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'B3'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='X1')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['X1'+str(I)]))
            posep = np.where(elftp[ielftype[ie]]=='X2')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['X2'+str(I)]))
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
            self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
        self.objlower = 0.0
#    Solution
# LO SOLTN               
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-CNOR2-MN-3-128"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eE6(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        E = np.exp(-EV_[1]*self.elpar[iel_][1])
        X1E = self.elpar[iel_][0]*E
        V1X1E = EV_[0]*self.elpar[iel_][0]*E
        f_   = V1X1E
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = X1E
            g_[1] = -V1X1E*self.elpar[iel_][1]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = -X1E*self.elpar[iel_][1]
                H_[1,0] = H_[0,1]
                H_[1,1] = V1X1E*self.elpar[iel_][1]**2
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

