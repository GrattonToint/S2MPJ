from s2mpjlib import *
class  ENSOLS(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : ENSOLS
#    *********
# 
#    NIST Data fitting problem ENSO.
# 
#    Fit: y = b1 + b2*cos( 2*pi*x/12 ) + b3*sin( 2*pi*x/12 ) 
#                      + b5*cos( 2*pi*x/b4 ) + b6*sin( 2*pi*x/b4 )
#                      + b8*cos( 2*pi*x/b7 ) + b9*sin( 2*pi*x/b7 ) + e
# 
#    Source:  Problem from the NIST nonlinear regression test set
#      http://www.itl.nist.gov/div898/strd/nls/nls_main.shtml
# 
#    Reference: Kahaner, D., C. Moler, and S. Nash, (1989). 
#     Numerical Methods and Software.  
#     Englewood Cliffs, NJ: Prentice Hall, pp. 441-445.
# 
#    SIF input: Nick Gould and Tyrone Rees, Oct 2015
# 
#    classification = "C-CSUR2-MN-9-0"
# 
#    Number of data values
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'ENSOLS'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['M'] = 168
        v_['N'] = 9
        v_['1'] = 1
        v_['TWELVE'] = 12.0
        v_['PI/4'] = np.arctan(1.0)
        v_['2PI'] = 8.0*v_['PI/4']
        v_['2PIBY12'] = v_['2PI']/v_['TWELVE']
        v_['X1'] = 1.0
        v_['X2'] = 2.0
        v_['X3'] = 3.0
        v_['X4'] = 4.0
        v_['X5'] = 5.0
        v_['X6'] = 6.0
        v_['X7'] = 7.0
        v_['X8'] = 8.0
        v_['X9'] = 9.0
        v_['X10'] = 10.0
        v_['X11'] = 11.0
        v_['X12'] = 12.0
        v_['X13'] = 13.0
        v_['X14'] = 14.0
        v_['X15'] = 15.0
        v_['X16'] = 16.0
        v_['X17'] = 17.0
        v_['X18'] = 18.0
        v_['X19'] = 19.0
        v_['X20'] = 20.0
        v_['X21'] = 21.0
        v_['X22'] = 22.0
        v_['X23'] = 23.0
        v_['X24'] = 24.0
        v_['X25'] = 25.0
        v_['X26'] = 26.0
        v_['X27'] = 27.0
        v_['X28'] = 28.0
        v_['X29'] = 29.0
        v_['X30'] = 30.0
        v_['X31'] = 31.0
        v_['X32'] = 32.0
        v_['X33'] = 33.0
        v_['X34'] = 34.0
        v_['X35'] = 35.0
        v_['X36'] = 36.0
        v_['X37'] = 37.0
        v_['X38'] = 38.0
        v_['X39'] = 39.0
        v_['X40'] = 40.0
        v_['X41'] = 41.0
        v_['X42'] = 42.0
        v_['X43'] = 43.0
        v_['X44'] = 44.0
        v_['X45'] = 45.0
        v_['X46'] = 46.0
        v_['X47'] = 47.0
        v_['X48'] = 48.0
        v_['X49'] = 49.0
        v_['X50'] = 50.0
        v_['X51'] = 51.0
        v_['X52'] = 52.0
        v_['X53'] = 53.0
        v_['X54'] = 54.0
        v_['X55'] = 55.0
        v_['X56'] = 56.0
        v_['X57'] = 57.0
        v_['X58'] = 58.0
        v_['X59'] = 59.0
        v_['X60'] = 60.0
        v_['X61'] = 61.0
        v_['X62'] = 62.0
        v_['X63'] = 63.0
        v_['X64'] = 64.0
        v_['X65'] = 65.0
        v_['X66'] = 66.0
        v_['X67'] = 67.0
        v_['X68'] = 68.0
        v_['X69'] = 69.0
        v_['X70'] = 70.0
        v_['X71'] = 71.0
        v_['X72'] = 72.0
        v_['X73'] = 73.0
        v_['X74'] = 74.0
        v_['X75'] = 75.0
        v_['X76'] = 76.0
        v_['X77'] = 77.0
        v_['X78'] = 78.0
        v_['X79'] = 79.0
        v_['X80'] = 80.0
        v_['X81'] = 81.0
        v_['X82'] = 82.0
        v_['X83'] = 83.0
        v_['X84'] = 84.0
        v_['X85'] = 85.0
        v_['X86'] = 86.0
        v_['X87'] = 87.0
        v_['X88'] = 88.0
        v_['X89'] = 89.0
        v_['X90'] = 90.0
        v_['X91'] = 91.0
        v_['X92'] = 92.0
        v_['X93'] = 93.0
        v_['X94'] = 94.0
        v_['X95'] = 95.0
        v_['X96'] = 96.0
        v_['X97'] = 97.0
        v_['X98'] = 98.0
        v_['X99'] = 99.0
        v_['X100'] = 100.0
        v_['X101'] = 101.0
        v_['X102'] = 102.0
        v_['X103'] = 103.0
        v_['X104'] = 104.0
        v_['X105'] = 105.0
        v_['X106'] = 106.0
        v_['X107'] = 107.0
        v_['X108'] = 108.0
        v_['X109'] = 109.0
        v_['X110'] = 110.0
        v_['X111'] = 111.0
        v_['X112'] = 112.0
        v_['X113'] = 113.0
        v_['X114'] = 114.0
        v_['X115'] = 115.0
        v_['X116'] = 116.0
        v_['X117'] = 117.0
        v_['X118'] = 118.0
        v_['X119'] = 119.0
        v_['X120'] = 120.0
        v_['X121'] = 121.0
        v_['X122'] = 122.0
        v_['X123'] = 123.0
        v_['X124'] = 124.0
        v_['X125'] = 125.0
        v_['X126'] = 126.0
        v_['X127'] = 127.0
        v_['X128'] = 128.0
        v_['X129'] = 129.0
        v_['X130'] = 130.0
        v_['X131'] = 131.0
        v_['X132'] = 132.0
        v_['X133'] = 133.0
        v_['X134'] = 134.0
        v_['X135'] = 135.0
        v_['X136'] = 136.0
        v_['X137'] = 137.0
        v_['X138'] = 138.0
        v_['X139'] = 139.0
        v_['X140'] = 140.0
        v_['X141'] = 141.0
        v_['X142'] = 142.0
        v_['X143'] = 143.0
        v_['X144'] = 144.0
        v_['X145'] = 145.0
        v_['X146'] = 146.0
        v_['X147'] = 147.0
        v_['X148'] = 148.0
        v_['X149'] = 149.0
        v_['X150'] = 150.0
        v_['X151'] = 151.0
        v_['X152'] = 152.0
        v_['X153'] = 153.0
        v_['X154'] = 154.0
        v_['X155'] = 155.0
        v_['X156'] = 156.0
        v_['X157'] = 157.0
        v_['X158'] = 158.0
        v_['X159'] = 159.0
        v_['X160'] = 160.0
        v_['X161'] = 161.0
        v_['X162'] = 162.0
        v_['X163'] = 163.0
        v_['X164'] = 164.0
        v_['X165'] = 165.0
        v_['X166'] = 166.0
        v_['X167'] = 167.0
        v_['X168'] = 168.0
        v_['Y1'] = 12.9
        v_['Y2'] = 11.3
        v_['Y3'] = 10.6
        v_['Y4'] = 11.2
        v_['Y5'] = 10.9
        v_['Y6'] = 7.5
        v_['Y7'] = 7.7
        v_['Y8'] = 11.7
        v_['Y9'] = 12.9
        v_['Y10'] = 14.3
        v_['Y11'] = 10.9
        v_['Y12'] = 13.7
        v_['Y13'] = 17.1
        v_['Y14'] = 14.0
        v_['Y15'] = 15.3
        v_['Y16'] = 8.5
        v_['Y17'] = 5.7
        v_['Y18'] = 5.5
        v_['Y19'] = 7.6
        v_['Y20'] = 8.6
        v_['Y21'] = 7.3
        v_['Y22'] = 7.6
        v_['Y23'] = 12.7
        v_['Y24'] = 11.0
        v_['Y25'] = 12.7
        v_['Y26'] = 12.9
        v_['Y27'] = 13.0
        v_['Y28'] = 10.9
        v_['Y29'] = 10.4
        v_['Y30'] = 10.2
        v_['Y31'] = 8.0
        v_['Y32'] = 10.9
        v_['Y33'] = 13.6
        v_['Y34'] = 10.5
        v_['Y35'] = 9.2
        v_['Y36'] = 12.4
        v_['Y37'] = 12.7
        v_['Y38'] = 13.3
        v_['Y39'] = 10.1
        v_['Y40'] = 7.8
        v_['Y41'] = 4.8
        v_['Y42'] = 3.0
        v_['Y43'] = 2.5
        v_['Y44'] = 6.3
        v_['Y45'] = 9.7
        v_['Y46'] = 11.6
        v_['Y47'] = 8.6
        v_['Y48'] = 12.4
        v_['Y49'] = 10.5
        v_['Y50'] = 13.3
        v_['Y51'] = 10.4
        v_['Y52'] = 8.1
        v_['Y53'] = 3.7
        v_['Y54'] = 10.7
        v_['Y55'] = 5.1
        v_['Y56'] = 10.4
        v_['Y57'] = 10.9
        v_['Y58'] = 11.7
        v_['Y59'] = 11.4
        v_['Y60'] = 13.7
        v_['Y61'] = 14.1
        v_['Y62'] = 14.0
        v_['Y63'] = 12.5
        v_['Y64'] = 6.3
        v_['Y65'] = 9.6
        v_['Y66'] = 11.7
        v_['Y67'] = 5.0
        v_['Y68'] = 10.8
        v_['Y69'] = 12.7
        v_['Y70'] = 10.8
        v_['Y71'] = 11.8
        v_['Y72'] = 12.6
        v_['Y73'] = 15.7
        v_['Y74'] = 12.6
        v_['Y75'] = 14.8
        v_['Y76'] = 7.8
        v_['Y77'] = 7.1
        v_['Y78'] = 11.2
        v_['Y79'] = 8.1
        v_['Y80'] = 6.4
        v_['Y81'] = 5.2
        v_['Y82'] = 12.0
        v_['Y83'] = 10.2
        v_['Y84'] = 12.7
        v_['Y85'] = 10.2
        v_['Y86'] = 14.7
        v_['Y87'] = 12.2
        v_['Y88'] = 7.1
        v_['Y89'] = 5.7
        v_['Y90'] = 6.7
        v_['Y91'] = 3.9
        v_['Y92'] = 8.5
        v_['Y93'] = 8.3
        v_['Y94'] = 10.8
        v_['Y95'] = 16.7
        v_['Y96'] = 12.6
        v_['Y97'] = 12.5
        v_['Y98'] = 12.5
        v_['Y99'] = 9.8
        v_['Y100'] = 7.2
        v_['Y101'] = 4.1
        v_['Y102'] = 10.6
        v_['Y103'] = 10.1
        v_['Y104'] = 10.1
        v_['Y105'] = 11.9
        v_['Y106'] = 13.6
        v_['Y107'] = 16.3
        v_['Y108'] = 17.6
        v_['Y109'] = 15.5
        v_['Y110'] = 16.0
        v_['Y111'] = 15.2
        v_['Y112'] = 11.2
        v_['Y113'] = 14.3
        v_['Y114'] = 14.5
        v_['Y115'] = 8.5
        v_['Y116'] = 12.0
        v_['Y117'] = 12.7
        v_['Y118'] = 11.3
        v_['Y119'] = 14.5
        v_['Y120'] = 15.1
        v_['Y121'] = 10.4
        v_['Y122'] = 11.5
        v_['Y123'] = 13.4
        v_['Y124'] = 7.5
        v_['Y125'] = 0.6
        v_['Y126'] = 0.3
        v_['Y127'] = 5.5
        v_['Y128'] = 5.0
        v_['Y129'] = 4.6
        v_['Y130'] = 8.2
        v_['Y131'] = 9.9
        v_['Y132'] = 9.2
        v_['Y133'] = 12.5
        v_['Y134'] = 10.9
        v_['Y135'] = 9.9
        v_['Y136'] = 8.9
        v_['Y137'] = 7.6
        v_['Y138'] = 9.5
        v_['Y139'] = 8.4
        v_['Y140'] = 10.7
        v_['Y141'] = 13.6
        v_['Y142'] = 13.7
        v_['Y143'] = 13.7
        v_['Y144'] = 16.5
        v_['Y145'] = 16.8
        v_['Y146'] = 17.1
        v_['Y147'] = 15.4
        v_['Y148'] = 9.5
        v_['Y149'] = 6.1
        v_['Y150'] = 10.1
        v_['Y151'] = 9.3
        v_['Y152'] = 5.3
        v_['Y153'] = 11.2
        v_['Y154'] = 16.6
        v_['Y155'] = 15.6
        v_['Y156'] = 12.0
        v_['Y157'] = 11.5
        v_['Y158'] = 8.6
        v_['Y159'] = 13.8
        v_['Y160'] = 8.7
        v_['Y161'] = 8.6
        v_['Y162'] = 8.6
        v_['Y163'] = 8.7
        v_['Y164'] = 12.8
        v_['Y165'] = 13.2
        v_['Y166'] = 14.0
        v_['Y167'] = 13.4
        v_['Y168'] = 14.8
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
            v_['ARG'] = v_['2PIBY12']*v_['X'+str(I)]
            v_['C'] = np.cos(v_['ARG'])
            v_['S'] = np.sin(v_['ARG'])
            [ig,ig_,_] = s2mpj_ii('F'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['B1']])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['B2']])
            valA = np.append(valA,float(v_['C']))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['B3']])
            valA = np.append(valA,float(v_['S']))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        self.objgrps = np.arange(ngrp)
        self.m       = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        for I in range(int(v_['1']),int(v_['M'])+1):
            self.gconst = arrset(self.gconst,ig_['F'+str(I)],float(v_['Y'+str(I)]))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.x0[ix_['B1']] = float(11.0)
        self.x0[ix_['B2']] = float(3.0)
        self.x0[ix_['B3']] = float(0.5)
        self.x0[ix_['B4']] = float(40.0)
        self.x0[ix_['B5']] = float(-0.7)
        self.x0[ix_['B6']] = float(-1.3)
        self.x0[ix_['B7']] = float(25.0)
        self.x0[ix_['B8']] = float(-0.3)
        self.x0[ix_['B9']] = float(1.4)
        pass
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eE8', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftp = []
        elftp = loaset(elftp,it,0,'X')
        [it,iet_,_] = s2mpj_ii( 'eE9', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftp = loaset(elftp,it,0,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for I in range(int(v_['1']),int(v_['M'])+1):
            ename = 'EA'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eE8')
            ielftype = arrset(ielftype,ie,iet_["eE8"])
            vname = 'B4'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'B5'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='X')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['X'+str(I)]))
            ename = 'EB'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eE9')
            ielftype = arrset(ielftype,ie,iet_["eE9"])
            vname = 'B4'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'B6'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='X')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['X'+str(I)]))
            ename = 'EC'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eE8')
            ielftype = arrset(ielftype,ie,iet_["eE8"])
            vname = 'B7'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'B8'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='X')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['X'+str(I)]))
            ename = 'ED'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eE9')
            ielftype = arrset(ielftype,ie,iet_["eE9"])
            vname = 'B7'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'B9'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='X')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['X'+str(I)]))
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gL2',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for ig in range(0,ngrp):
            self.grftype = arrset(self.grftype,ig,'gL2')
        for I in range(int(v_['1']),int(v_['M'])+1):
            ig = ig_['F'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['EA'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,1.)
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['EB'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,1.)
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['EC'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,1.)
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['ED'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
        self.objlower = 0.0
#    Solution
# LO SOLTN               
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass   = "C-CSUR2-MN-9-0"
        self.objderlvl = 2

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def e_globs(self):

        import numpy as np
        self.efpar = np.array([]);
        self.efpar = arrset( self.efpar,0,8.0*np.arctan(1.0e0))
        return pbm

    @staticmethod
    def eE8(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        V12 = EV_[0]*EV_[0]
        V13 = EV_[0]*V12
        V14 = V12*V12
        TPIX = self.efpar[0]*self.elpar[iel_][0]
        TPIXV1 = TPIX/EV_[0]
        C = np.cos(TPIXV1)
        S = np.sin(TPIXV1)
        TPIXS = TPIX*S
        f_   = EV_[1]*C
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = TPIXS*EV_[1]/V12
            g_[1] = C
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = (-self.efpar[0]*self.efpar[0]*EV_[1]*C*self.elpar[iel_][0]**2/
                     V14-2.0*TPIX*EV_[1]*S/V13)
                H_[0,1] = TPIXS/V12
                H_[1,0] = H_[0,1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eE9(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        V12 = EV_[0]*EV_[0]
        V13 = EV_[0]*V12
        V14 = V12*V12
        TPIX = self.efpar[0]*self.elpar[iel_][0]
        TPIXV1 = TPIX/EV_[0]
        C = np.cos(TPIXV1)
        S = np.sin(TPIXV1)
        TPIXC = TPIX*C
        f_   = EV_[1]*S
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = -TPIXC*EV_[1]/V12
            g_[1] = S
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0]  = (
                      2.0*TPIX*EV_[1]*C/V13-self.efpar[0]*self.efpar[0]*EV_[1]*S*self.elpar[iel_][0]**2/V14)
                H_[0,1] = -TPIXC/V12
                H_[1,0] = H_[0,1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gL2(self,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= GVAR_*GVAR_
        if nargout>1:
            g_ = GVAR_+GVAR_
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = 2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

