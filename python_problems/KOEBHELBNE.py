from s2mpjlib import *
class  KOEBHELBNE(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : KOEBHELBNE
#    *********
# 
#    An exponential fitting problem arising in the study of the
#    Koebl-Helbling conjecture on energy/time budgets and travel mode.
#    This is the raw data (KHb4).
# 
#    Source:
#    J. P. Hubert and Ph. L . Toint, Summer 2005.
#    SIF input: Ph. Toint, June 2005.
#    Bound-constrained nonlinear equations version: Nick Gould, June 2019.
# 
#    classification = "C-CNOR2-RN-3-156"
# 
#    Useful constants
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 21 VI 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'KOEBHELBNE'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['1'] = 1
        v_['M'] = 156
        v_['X1'] = 0.031181789
        v_['X2'] = 0.046772684
        v_['X3'] = 0.062363579
        v_['X4'] = 0.093545369
        v_['X5'] = 0.109136264
        v_['X6'] = 0.124727159
        v_['X7'] = 0.140318054
        v_['X8'] = 0.155908949
        v_['X9'] = 0.171499844
        v_['X10'] = 0.187090739
        v_['X11'] = 0.202681633
        v_['X12'] = 0.218272528
        v_['X13'] = 0.233863423
        v_['X14'] = 0.249454318
        v_['X15'] = 0.265045213
        v_['X16'] = 0.280636108
        v_['X17'] = 0.296227003
        v_['X18'] = 0.311817898
        v_['X19'] = 0.327408793
        v_['X20'] = 0.358590583
        v_['X21'] = 0.374181478
        v_['X22'] = 0.389772372
        v_['X23'] = 0.405363267
        v_['X24'] = 0.420954162
        v_['X25'] = 0.436545057
        v_['X26'] = 0.452135952
        v_['X27'] = 0.467726847
        v_['X28'] = 0.483317742
        v_['X29'] = 0.498908637
        v_['X30'] = 0.514499532
        v_['X31'] = 0.530090427
        v_['X32'] = 0.545681322
        v_['X33'] = 0.561272217
        v_['X34'] = 0.576863111
        v_['X35'] = 0.592454006
        v_['X36'] = 0.608044901
        v_['X37'] = 0.623635796
        v_['X38'] = 0.639226691
        v_['X39'] = 0.654817586
        v_['X40'] = 0.670408481
        v_['X41'] = 0.685999376
        v_['X42'] = 0.701590271
        v_['X43'] = 0.717181166
        v_['X44'] = 0.732772061
        v_['X45'] = 0.748362956
        v_['X46'] = 0.763953851
        v_['X47'] = 0.779544745
        v_['X48'] = 0.795135640
        v_['X49'] = 0.810726535
        v_['X50'] = 0.826317430
        v_['X51'] = 0.841908325
        v_['X52'] = 0.857499220
        v_['X53'] = 0.873090115
        v_['X54'] = 0.888681010
        v_['X55'] = 0.904271905
        v_['X56'] = 0.919862800
        v_['X57'] = 0.935453695
        v_['X58'] = 0.966635484
        v_['X59'] = 0.997817274
        v_['X60'] = 1.013408169
        v_['X61'] = 1.028999064
        v_['X62'] = 1.044589959
        v_['X63'] = 1.060180854
        v_['X64'] = 1.075771749
        v_['X65'] = 1.091362644
        v_['X66'] = 1.106953539
        v_['X67'] = 1.122544434
        v_['X68'] = 1.138135329
        v_['X69'] = 1.153726223
        v_['X70'] = 1.169317118
        v_['X71'] = 1.184908013
        v_['X72'] = 1.200498908
        v_['X73'] = 1.216089803
        v_['X74'] = 1.247271593
        v_['X75'] = 1.262862488
        v_['X76'] = 1.278453383
        v_['X77'] = 1.294044278
        v_['X78'] = 1.325226068
        v_['X79'] = 1.340816962
        v_['X80'] = 1.356407857
        v_['X81'] = 1.371998752
        v_['X82'] = 1.387589647
        v_['X83'] = 1.403180542
        v_['X84'] = 1.418771437
        v_['X85'] = 1.434362332
        v_['X86'] = 1.481135017
        v_['X87'] = 1.496725912
        v_['X88'] = 1.527907701
        v_['X89'] = 1.543498596
        v_['X90'] = 1.559089491
        v_['X91'] = 1.574680386
        v_['X92'] = 1.590271281
        v_['X93'] = 1.621453071
        v_['X94'] = 1.637043966
        v_['X95'] = 1.652634861
        v_['X96'] = 1.668225756
        v_['X97'] = 1.683816651
        v_['X98'] = 1.699407546
        v_['X99'] = 1.714998440
        v_['X100'] = 1.730589335
        v_['X101'] = 1.746180230
        v_['X102'] = 1.761771125
        v_['X103'] = 1.792952915
        v_['X104'] = 1.808543810
        v_['X105'] = 1.824134705
        v_['X106'] = 1.839725600
        v_['X107'] = 1.855316495
        v_['X108'] = 1.870907390
        v_['X109'] = 1.886498285
        v_['X110'] = 1.902089179
        v_['X111'] = 1.933270969
        v_['X112'] = 1.948861864
        v_['X113'] = 1.964452759
        v_['X114'] = 2.026816339
        v_['X115'] = 2.057998129
        v_['X116'] = 2.089179918
        v_['X117'] = 2.151543498
        v_['X118'] = 2.182725288
        v_['X119'] = 2.198316183
        v_['X120'] = 2.260679763
        v_['X121'] = 2.323043342
        v_['X122'] = 2.338634237
        v_['X123'] = 2.369816027
        v_['X124'] = 2.416588712
        v_['X125'] = 2.432179607
        v_['X126'] = 2.447770502
        v_['X127'] = 2.494543186
        v_['X128'] = 2.541315871
        v_['X129'] = 2.572497661
        v_['X130'] = 2.650452136
        v_['X131'] = 2.728406610
        v_['X132'] = 2.743997505
        v_['X133'] = 2.775179295
        v_['X134'] = 2.806361085
        v_['X135'] = 2.884315559
        v_['X136'] = 2.962270034
        v_['X137'] = 3.040224508
        v_['X138'] = 3.086997193
        v_['X139'] = 3.196133458
        v_['X140'] = 3.274087932
        v_['X141'] = 3.539133146
        v_['X142'] = 3.585905831
        v_['X143'] = 3.897723729
        v_['X144'] = 4.006859993
        v_['X145'] = 4.209541627
        v_['X146'] = 4.287496102
        v_['X147'] = 4.443405051
        v_['X148'] = 5.144995322
        v_['X149'] = 5.378858746
        v_['X150'] = 5.534767695
        v_['X151'] = 5.706267539
        v_['X152'] = 5.924540068
        v_['X153'] = 6.532584970
        v_['X154'] = 6.548175865
        v_['X155'] = 6.906766448
        v_['X156'] = 10.28999064
        v_['Y1'] = 1.395674579
        v_['Y2'] = 0.038471426
        v_['Y3'] = 0.319318940
        v_['Y4'] = 0.545752949
        v_['Y5'] = 0.088398192
        v_['Y6'] = 0.162971961
        v_['Y7'] = 0.515829574
        v_['Y8'] = 3.065817290
        v_['Y9'] = 0.061781114
        v_['Y10'] = 0.422348682
        v_['Y11'] = 0.348789741
        v_['Y12'] = 1.254137154
        v_['Y13'] = 0.806640615
        v_['Y14'] = 0.297642996
        v_['Y15'] = 0.636546644
        v_['Y16'] = 0.731486009
        v_['Y17'] = 0.178811464
        v_['Y18'] = 5.849019044
        v_['Y19'] = 0.211425133
        v_['Y20'] = 0.099010516
        v_['Y21'] = 0.513818504
        v_['Y22'] = 0.659774345
        v_['Y23'] = 0.439481954
        v_['Y24'] = 0.287092436
        v_['Y25'] = 0.427519928
        v_['Y26'] = 0.988931826
        v_['Y27'] = 4.567207169
        v_['Y28'] = 0.143825432
        v_['Y29'] = 0.230813434
        v_['Y30'] = 0.653476950
        v_['Y31'] = 0.349046454
        v_['Y32'] = 0.658389297
        v_['Y33'] = 0.243390186
        v_['Y34'] = 0.198643958
        v_['Y35'] = 0.113079629
        v_['Y36'] = 0.027324724
        v_['Y37'] = 3.309046781
        v_['Y38'] = 0.114929092
        v_['Y39'] = 0.605280127
        v_['Y40'] = 0.274520239
        v_['Y41'] = 0.318091491
        v_['Y42'] = 0.766436148
        v_['Y43'] = 0.869738272
        v_['Y44'] = 0.139373492
        v_['Y45'] = 0.084902595
        v_['Y46'] = 0.181622918
        v_['Y47'] = 1.796726648
        v_['Y48'] = 0.243741458
        v_['Y49'] = 0.293308390
        v_['Y50'] = 0.042250701
        v_['Y51'] = 0.326540975
        v_['Y52'] = 0.615312705
        v_['Y53'] = 0.208960214
        v_['Y54'] = 0.172436454
        v_['Y55'] = 0.296983631
        v_['Y56'] = 0.309713792
        v_['Y57'] = 2.724919333
        v_['Y58'] = 0.329865673
        v_['Y59'] = 0.118646056
        v_['Y60'] = 0.938635931
        v_['Y61'] = 0.193912348
        v_['Y62'] = 0.046748918
        v_['Y63'] = 0.304652956
        v_['Y64'] = 0.128677541
        v_['Y65'] = 1.451844993
        v_['Y66'] = 0.358594758
        v_['Y67'] = 0.089530900
        v_['Y68'] = 0.061870937
        v_['Y69'] = 0.739914176
        v_['Y70'] = 2.155312843
        v_['Y71'] = 0.174792785
        v_['Y72'] = 0.230422808
        v_['Y73'] = 0.169723386
        v_['Y74'] = 1.275355242
        v_['Y75'] = 0.071289699
        v_['Y76'] = 0.061203920
        v_['Y77'] = 0.107195817
        v_['Y78'] = 0.510687662
        v_['Y79'] = 0.266987743
        v_['Y80'] = 0.163949073
        v_['Y81'] = 0.216590731
        v_['Y82'] = 0.116983161
        v_['Y83'] = 0.480761554
        v_['Y84'] = 0.110349980
        v_['Y85'] = 0.117447394
        v_['Y86'] = 0.543860123
        v_['Y87'] = 0.088398192
        v_['Y88'] = 0.018754444
        v_['Y89'] = 0.088398192
        v_['Y90'] = 0.70860539
        v_['Y91'] = 0.059785530
        v_['Y92'] = 0.209830196
        v_['Y93'] = 0.033521364
        v_['Y94'] = 0.486638260
        v_['Y95'] = 0.038248419
        v_['Y96'] = 0.095664501
        v_['Y97'] = 0.054377430
        v_['Y98'] = 0.039112388
        v_['Y99'] = 0.396786332
        v_['Y100'] = 0.435306404
        v_['Y101'] = 0.071289699
        v_['Y102'] = 0.215081426
        v_['Y103'] = 0.127195748
        v_['Y104'] = 0.088398192
        v_['Y105'] = 0.091835487
        v_['Y106'] = 0.437799746
        v_['Y107'] = 0.032532956
        v_['Y108'] = 1.029175100
        v_['Y109'] = 0.067042729
        v_['Y110'] = 0.053346206
        v_['Y111'] = 0.239056856
        v_['Y112'] = 0.141051692
        v_['Y113'] = 0.097071959
        v_['Y114'] = 0.200019167
        v_['Y115'] = 0.072932371
        v_['Y116'] = 0.088398192
        v_['Y117'] = 0.036090863
        v_['Y118'] = 0.324539379
        v_['Y119'] = 0.160828616
        v_['Y120'] = 0.173438346
        v_['Y121'] = 0.164108859
        v_['Y122'] = 0.254930430
        v_['Y123'] = 0.019600012
        v_['Y124'] = 0.050492665
        v_['Y125'] = 0.050950521
        v_['Y126'] = 0.072472146
        v_['Y127'] = 0.269688058
        v_['Y128'] = 0.030918706
        v_['Y129'] = 0.219303253
        v_['Y130'] = 0.135149297
        v_['Y131'] = 0.430161394
        v_['Y132'] = 0.061781114
        v_['Y133'] = 0.180567097
        v_['Y134'] = 0.196062616
        v_['Y135'] = 0.059251516
        v_['Y136'] = 0.081376936
        v_['Y137'] = 0.165941196
        v_['Y138'] = 0.067042729
        v_['Y139'] = 0.115821302
        v_['Y140'] = 0.091364330
        v_['Y141'] = 0.022955501
        v_['Y142'] = 0.038751824
        v_['Y143'] = 0.055819323
        v_['Y144'] = 0.404387698
        v_['Y145'] = 0.077026844
        v_['Y146'] = 0.050492665
        v_['Y147'] = 0.036020536
        v_['Y148'] = 0.130978302
        v_['Y149'] = 0.031490070
        v_['Y150'] = 0.031490070
        v_['Y151'] = 0.147735523
        v_['Y152'] = 0.075884662
        v_['Y153'] = 0.145315971
        v_['Y154'] = 0.837166509
        v_['Y155'] = 0.035501280
        v_['Y156'] = 0.039112388
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        [iv,ix_,_] = s2mpj_ii('N',ix_)
        self.xnames=arrset(self.xnames,iv,'N')
        [iv,ix_,_] = s2mpj_ii('A',ix_)
        self.xnames=arrset(self.xnames,iv,'A')
        [iv,ix_,_] = s2mpj_ii('B',ix_)
        self.xnames=arrset(self.xnames,iv,'B')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['1']),int(v_['M'])+1):
            [ig,ig_,_] = s2mpj_ii('O'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'O'+str(I))
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
            self.gconst = arrset(self.gconst,ig_['O'+str(I)],float(v_['Y'+str(I)]))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        self.xlower[ix_['A']] = -float('Inf')
        self.xupper[ix_['A']] = +float('Inf')
        self.xlower[ix_['B']] = 0.001
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        if('N' in ix_):
            self.x0[ix_['N']] = float(10.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['N']),float(10.0)))
        if('A' in ix_):
            self.x0[ix_['A']] = float(0.01)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['A']),float(0.01)))
        if('B' in ix_):
            self.x0[ix_['B']] = float(0.01)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['B']),float(0.01)))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eKHE', iet_)
        elftv = loaset(elftv,it,0,'VN')
        elftv = loaset(elftv,it,1,'VA')
        elftv = loaset(elftv,it,2,'VB')
        elftp = []
        elftp = loaset(elftp,it,0,'XX')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for I in range(int(v_['1']),int(v_['M'])+1):
            ename = 'E'+str(I)
            [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
            if newelt:
                self.elftype = arrset(self.elftype,ie,'eKHE')
                ielftype = arrset(ielftype,ie,iet_['eKHE'])
            vname = 'N'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='VN')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'A'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='VA')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'B'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='VB')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='XX')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['X'+str(I)]))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['M'])+1):
            ig = ig_['O'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
# LO                       0.0
#    Solution
# LO SOLTN                 77.516347286
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-CNOR2-RN-3-156"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eKHE(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        T = np.exp(-EV_[1]/self.elpar[iel_][0]-self.elpar[iel_][0]/EV_[2])
        M1OX = -1.0/self.elpar[iel_][0]
        XOB2 = self.elpar[iel_][0]/(EV_[2]*EV_[2])
        M2XOB3 = -2.0*self.elpar[iel_][0]/EV_[2]**3
        f_   = EV_[0]*T
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = T
            g_[1] = EV_[0]*T*M1OX
            g_[2] = EV_[0]*T*XOB2
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = T*M1OX
                H_[1,0] = H_[0,1]
                H_[0,2] = T*XOB2
                H_[2,0] = H_[0,2]
                H_[1,1] = EV_[0]*T*M1OX*M1OX
                H_[1,2] = EV_[0]*T*M1OX*XOB2
                H_[2,1] = H_[1,2]
                H_[2,2] = EV_[0]*T*(XOB2*XOB2+M2XOB3)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

