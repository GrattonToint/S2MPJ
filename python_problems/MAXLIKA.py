from s2mpjlib import *
class  MAXLIKA(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : MAXLIKA
#    *********
# 
#    A variant of Hock and Schittkowski problem 105, where the
#    (inactive) inequality constraint is dropped.
# 
#    Source:
#    Ph. Toint and A. Griewank.
# 
#    SIF input: Ph. Toint, June 1990.
# 
#    classification = "C-COBR2-AY-8-0"
# 
#    Constants
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'MAXLIKA'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['1'] = 1
        v_['235'] = 235
        v_['Y1'] = 95.0
        v_['Y2'] = 105.0
        v_['Y3'] = 110.0
        v_['Y4'] = 110.0
        v_['Y5'] = 110.0
        v_['Y6'] = 110.0
        v_['Y7'] = 115.0
        v_['Y8'] = 115.0
        v_['Y9'] = 115.0
        v_['Y10'] = 115.0
        v_['Y11'] = 120.0
        v_['Y12'] = 120.0
        v_['Y13'] = 120.0
        v_['Y14'] = 120.0
        v_['Y15'] = 120.0
        v_['Y16'] = 120.0
        v_['Y17'] = 120.0
        v_['Y18'] = 120.0
        v_['Y19'] = 120.0
        v_['Y20'] = 120.0
        v_['Y21'] = 120.0
        v_['Y22'] = 120.0
        v_['Y23'] = 120.0
        v_['Y24'] = 120.0
        v_['Y25'] = 120.0
        v_['Y26'] = 125.0
        v_['Y27'] = 125.0
        v_['Y28'] = 125.0
        v_['Y29'] = 125.0
        v_['Y30'] = 125.0
        v_['Y31'] = 125.0
        v_['Y32'] = 125.0
        v_['Y33'] = 125.0
        v_['Y34'] = 125.0
        v_['Y35'] = 125.0
        v_['Y36'] = 125.0
        v_['Y37'] = 125.0
        v_['Y38'] = 125.0
        v_['Y39'] = 125.0
        v_['Y40'] = 125.0
        v_['Y41'] = 130.0
        v_['Y42'] = 130.0
        v_['Y43'] = 130.0
        v_['Y44'] = 130.0
        v_['Y45'] = 130.0
        v_['Y46'] = 130.0
        v_['Y47'] = 130.0
        v_['Y48'] = 130.0
        v_['Y49'] = 130.0
        v_['Y50'] = 130.0
        v_['Y51'] = 130.0
        v_['Y52'] = 130.0
        v_['Y53'] = 130.0
        v_['Y54'] = 130.0
        v_['Y55'] = 130.0
        v_['Y56'] = 135.0
        v_['Y57'] = 135.0
        v_['Y58'] = 135.0
        v_['Y59'] = 135.0
        v_['Y60'] = 135.0
        v_['Y61'] = 135.0
        v_['Y62'] = 135.0
        v_['Y63'] = 135.0
        v_['Y64'] = 135.0
        v_['Y65'] = 135.0
        v_['Y66'] = 135.0
        v_['Y67'] = 135.0
        v_['Y68'] = 135.0
        v_['Y69'] = 140.0
        v_['Y70'] = 140.0
        v_['Y71'] = 140.0
        v_['Y72'] = 140.0
        v_['Y73'] = 140.0
        v_['Y74'] = 140.0
        v_['Y75'] = 140.0
        v_['Y76'] = 140.0
        v_['Y77'] = 140.0
        v_['Y78'] = 140.0
        v_['Y79'] = 140.0
        v_['Y80'] = 140.0
        v_['Y81'] = 140.0
        v_['Y82'] = 140.0
        v_['Y83'] = 140.0
        v_['Y84'] = 140.0
        v_['Y85'] = 140.0
        v_['Y86'] = 140.0
        v_['Y87'] = 140.0
        v_['Y88'] = 140.0
        v_['Y89'] = 140.0
        v_['Y90'] = 145.0
        v_['Y91'] = 145.0
        v_['Y92'] = 145.0
        v_['Y93'] = 145.0
        v_['Y94'] = 145.0
        v_['Y95'] = 145.0
        v_['Y96'] = 145.0
        v_['Y97'] = 145.0
        v_['Y98'] = 145.0
        v_['Y99'] = 145.0
        v_['Y100'] = 145.0
        v_['Y101'] = 145.0
        v_['Y102'] = 150.0
        v_['Y103'] = 150.0
        v_['Y104'] = 150.0
        v_['Y105'] = 150.0
        v_['Y106'] = 150.0
        v_['Y107'] = 150.0
        v_['Y108'] = 150.0
        v_['Y109'] = 150.0
        v_['Y110'] = 150.0
        v_['Y111'] = 150.0
        v_['Y112'] = 150.0
        v_['Y113'] = 150.0
        v_['Y114'] = 150.0
        v_['Y115'] = 150.0
        v_['Y116'] = 150.0
        v_['Y117'] = 150.0
        v_['Y118'] = 150.0
        v_['Y119'] = 155.0
        v_['Y120'] = 155.0
        v_['Y121'] = 155.0
        v_['Y122'] = 155.0
        v_['Y123'] = 160.0
        v_['Y124'] = 160.0
        v_['Y125'] = 160.0
        v_['Y126'] = 160.0
        v_['Y127'] = 160.0
        v_['Y128'] = 160.0
        v_['Y129'] = 160.0
        v_['Y130'] = 160.0
        v_['Y131'] = 160.0
        v_['Y132'] = 160.0
        v_['Y133'] = 160.0
        v_['Y134'] = 160.0
        v_['Y135'] = 160.0
        v_['Y136'] = 160.0
        v_['Y137'] = 160.0
        v_['Y138'] = 160.0
        v_['Y139'] = 160.0
        v_['Y140'] = 160.0
        v_['Y141'] = 160.0
        v_['Y142'] = 160.0
        v_['Y143'] = 165.0
        v_['Y144'] = 165.0
        v_['Y145'] = 165.0
        v_['Y146'] = 165.0
        v_['Y147'] = 165.0
        v_['Y148'] = 165.0
        v_['Y149'] = 165.0
        v_['Y150'] = 165.0
        v_['Y151'] = 170.0
        v_['Y152'] = 170.0
        v_['Y153'] = 170.0
        v_['Y154'] = 170.0
        v_['Y155'] = 170.0
        v_['Y156'] = 170.0
        v_['Y157'] = 170.0
        v_['Y158'] = 170.0
        v_['Y159'] = 170.0
        v_['Y160'] = 170.0
        v_['Y161'] = 170.0
        v_['Y162'] = 170.0
        v_['Y163'] = 170.0
        v_['Y164'] = 170.0
        v_['Y165'] = 170.0
        v_['Y166'] = 170.0
        v_['Y167'] = 170.0
        v_['Y168'] = 175.0
        v_['Y169'] = 175.0
        v_['Y170'] = 175.0
        v_['Y171'] = 175.0
        v_['Y172'] = 175.0
        v_['Y173'] = 175.0
        v_['Y174'] = 175.0
        v_['Y175'] = 175.0
        v_['Y176'] = 180.0
        v_['Y177'] = 180.0
        v_['Y178'] = 180.0
        v_['Y179'] = 180.0
        v_['Y180'] = 180.0
        v_['Y181'] = 180.0
        v_['Y182'] = 185.0
        v_['Y183'] = 185.0
        v_['Y184'] = 185.0
        v_['Y185'] = 185.0
        v_['Y186'] = 185.0
        v_['Y187'] = 185.0
        v_['Y188'] = 190.0
        v_['Y189'] = 190.0
        v_['Y190'] = 190.0
        v_['Y191'] = 190.0
        v_['Y192'] = 190.0
        v_['Y193'] = 190.0
        v_['Y194'] = 190.0
        v_['Y195'] = 195.0
        v_['Y196'] = 195.0
        v_['Y197'] = 195.0
        v_['Y198'] = 195.0
        v_['Y199'] = 200.0
        v_['Y200'] = 200.0
        v_['Y201'] = 200.0
        v_['Y202'] = 205.0
        v_['Y203'] = 205.0
        v_['Y204'] = 205.0
        v_['Y205'] = 210.0
        v_['Y206'] = 210.0
        v_['Y207'] = 210.0
        v_['Y208'] = 210.0
        v_['Y209'] = 210.0
        v_['Y210'] = 210.0
        v_['Y211'] = 210.0
        v_['Y212'] = 210.0
        v_['Y213'] = 215.0
        v_['Y214'] = 220.0
        v_['Y215'] = 220.0
        v_['Y216'] = 220.0
        v_['Y217'] = 220.0
        v_['Y218'] = 220.0
        v_['Y219'] = 220.0
        v_['Y220'] = 230.0
        v_['Y221'] = 230.0
        v_['Y222'] = 230.0
        v_['Y223'] = 230.0
        v_['Y224'] = 230.0
        v_['Y225'] = 235.0
        v_['Y226'] = 240.0
        v_['Y227'] = 240.0
        v_['Y228'] = 240.0
        v_['Y229'] = 240.0
        v_['Y230'] = 240.0
        v_['Y231'] = 240.0
        v_['Y232'] = 240.0
        v_['Y233'] = 245.0
        v_['Y234'] = 250.0
        v_['Y235'] = 250.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        [iv,ix_,_] = s2mpj_ii('X1',ix_)
        self.xnames=arrset(self.xnames,iv,'X1')
        [iv,ix_,_] = s2mpj_ii('X2',ix_)
        self.xnames=arrset(self.xnames,iv,'X2')
        [iv,ix_,_] = s2mpj_ii('X3',ix_)
        self.xnames=arrset(self.xnames,iv,'X3')
        [iv,ix_,_] = s2mpj_ii('X4',ix_)
        self.xnames=arrset(self.xnames,iv,'X4')
        [iv,ix_,_] = s2mpj_ii('X5',ix_)
        self.xnames=arrset(self.xnames,iv,'X5')
        [iv,ix_,_] = s2mpj_ii('X6',ix_)
        self.xnames=arrset(self.xnames,iv,'X6')
        [iv,ix_,_] = s2mpj_ii('X7',ix_)
        self.xnames=arrset(self.xnames,iv,'X7')
        [iv,ix_,_] = s2mpj_ii('X8',ix_)
        self.xnames=arrset(self.xnames,iv,'X8')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['1']),int(v_['235'])+1):
            [ig,ig_,_] = s2mpj_ii('L'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            self.gscale = arrset(self.gscale,ig,float(-1.0))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        self.objgrps = np.arange(ngrp)
        self.m       = 0
        selfnob      = ngrp
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        self.xlower[ix_['X1']] = 0.001
        self.xupper[ix_['X1']] = 0.499
        self.xlower[ix_['X2']] = 0.001
        self.xupper[ix_['X2']] = 0.499
        self.xlower[ix_['X3']] = 100.0
        self.xupper[ix_['X3']] = 180.0
        self.xlower[ix_['X4']] = 130.0
        self.xupper[ix_['X4']] = 210.0
        self.xlower[ix_['X5']] = 170.0
        self.xupper[ix_['X5']] = 240.0
        self.xlower[ix_['X6']] = 5.0
        self.xupper[ix_['X6']] = 25.0
        self.xlower[ix_['X7']] = 5.0
        self.xupper[ix_['X7']] = 25.0
        self.xlower[ix_['X8']] = 5.0
        self.xupper[ix_['X8']] = 25.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.x0[ix_['X1']] = float(0.1)
        self.x0[ix_['X2']] = float(0.2)
        self.x0[ix_['X3']] = float(100.0)
        self.x0[ix_['X4']] = float(125.0)
        self.x0[ix_['X5']] = float(175.0)
        self.x0[ix_['X6']] = float(11.2)
        self.x0[ix_['X7']] = float(13.2)
        self.x0[ix_['X8']] = float(15.8)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eAB', iet_)
        elftv = loaset(elftv,it,0,'U')
        elftv = loaset(elftv,it,1,'V')
        elftv = loaset(elftv,it,2,'W')
        elftp = []
        elftp = loaset(elftp,it,0,'Y')
        [it,iet_,_] = s2mpj_ii( 'eC', iet_)
        elftv = loaset(elftv,it,0,'U')
        elftv = loaset(elftv,it,1,'X')
        elftv = loaset(elftv,it,2,'S')
        elftv = loaset(elftv,it,3,'T')
        elftp = loaset(elftp,it,0,'Y')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for I in range(int(v_['1']),int(v_['235'])+1):
            ename = 'A'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eAB')
            ielftype = arrset(ielftype,ie,iet_["eAB"])
            vname = 'X1'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='U')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X6'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X3'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='W')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='Y')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['Y'+str(I)]))
            ename = 'B'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eAB')
            ielftype = arrset(ielftype,ie,iet_["eAB"])
            vname = 'X2'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='U')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X7'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X4'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='W')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='Y')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['Y'+str(I)]))
            ename = 'C'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eC')
            ielftype = arrset(ielftype,ie,iet_["eC"])
            vname = 'X2'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='U')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X1'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X8'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='S')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X5'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='T')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='Y')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['Y'+str(I)]))
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gLN',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['235'])+1):
            ig = ig_['L'+str(I)]
            self.grftype = arrset(self.grftype,ig,'gLN')
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['A'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,1.)
            posel = posel+1
            self.grelt = loaset(self.grelt,ig,posel,ie_['B'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel, 1.)
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['C'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass   = "C-COBR2-AY-8-0"
        self.objderlvl = 2


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eAB(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        YMW = self.elpar[iel_][0]-EV_[2,0]
        YMWSQ = YMW*YMW
        VSQ = EV_[1,0]*EV_[1,0]
        VCB = VSQ*EV_[1,0]
        A = -YMWSQ/(2.0*VSQ)
        DADV = YMWSQ/VCB
        DADW = YMW/VSQ
        D2ADV2 = -3.0*YMWSQ/(VSQ*VSQ)
        D2ADVW = -2.0*YMW/VCB
        D2ADW2 = -1.0/VSQ
        E = np.exp(A)
        DEDV = E*DADV
        DEDW = E*DADW
        B = EV_[0,0]*E
        DBDV = B*DADV
        DBDW = B*DADW
        D2BDV2 = DBDV*DADV+B*D2ADV2
        D2BDVW = DBDW*DADV+B*D2ADVW
        D2BDW2 = DBDW*DADW+B*D2ADW2
        f_   = B/EV_[1,0]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = E/EV_[1,0]
            g_[1] = (DBDV-B/EV_[1,0])/EV_[1,0]
            g_[2] = DBDW/EV_[1,0]
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = (DEDV-E/EV_[1,0])/EV_[1,0]
                H_[1,0] = H_[0,1]
                H_[0,2] = DEDW/EV_[1,0]
                H_[2,0] = H_[0,2]
                H_[1,1] = (D2BDV2-DBDV/EV_[1,0]+B/VSQ)/EV_[1,0]-(DBDV-B/EV_[1,0])/VSQ
                H_[1,2] = (D2BDVW-DBDW/EV_[1,0])/EV_[1,0]
                H_[2,1] = H_[1,2]
                H_[2,2] = D2BDW2/EV_[1,0]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eC(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((3,4))
        IV_ = np.zeros(3)
        U_[0,0] = U_[0,0]+1
        U_[0,1] = U_[0,1]+1
        U_[1,2] = U_[1,2]+1
        U_[2,3] = U_[2,3]+1
        IV_[0] = to_scalar(U_[0:1,:].dot(EV_))
        IV_[1] = to_scalar(U_[1:2,:].dot(EV_))
        IV_[2] = to_scalar(U_[2:3,:].dot(EV_))
        YMW = self.elpar[iel_][0]-IV_[2]
        YMWSQ = YMW*YMW
        VSQ = IV_[1]*IV_[1]
        VCB = VSQ*IV_[1]
        A = -YMWSQ/(2.0*VSQ)
        DADV = YMWSQ/VCB
        DADW = YMW/VSQ
        D2ADV2 = -3.0*YMWSQ/(VSQ*VSQ)
        D2ADVW = -2.0*YMW/VCB
        D2ADW2 = -1.0/VSQ
        E = np.exp(A)
        DEDV = E*DADV
        DEDW = E*DADW
        B = (1.0-IV_[0])*E
        DBDZ = -E
        DBDV = B*DADV
        DBDW = B*DADW
        D2BDVZ = -DEDV
        D2BDWZ = -DEDW
        D2BDV2 = DBDV*DADV+B*D2ADV2
        D2BDVW = DBDW*DADV+B*D2ADVW
        D2BDW2 = DBDW*DADW+B*D2ADW2
        f_   = B/IV_[1]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = DBDZ/IV_[1]
            g_[1] = (DBDV-B/IV_[1])/IV_[1]
            g_[2] = DBDW/IV_[1]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = (DEDV-E/IV_[1])/IV_[1]
                H_[1,0] = H_[0,1]
                H_[0,2] = D2BDWZ/IV_[1]
                H_[2,0] = H_[0,2]
                H_[1,1] = (D2BDV2-DBDV/IV_[1]+B/VSQ)/IV_[1]-(DBDV-B/IV_[1])/VSQ
                H_[1,2] = (D2BDVW-DBDW/IV_[1])/IV_[1]
                H_[2,1] = H_[1,2]
                H_[2,2] = D2BDW2/IV_[1]
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gLN(self,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= np.log(GVAR_*0.39894228)
        if nargout>1:
            g_ = 1.0/GVAR_
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = -1.0/GVAR_**2
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

