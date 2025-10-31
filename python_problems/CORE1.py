from s2mpjlib import *
class  CORE1(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : CORE1
#    *********
# 
#    A problem from the exploitation of a gas transmission network
#    with consideration of the head loses in the pipes. The aim is
#    to satisfy the demand at several points in the network at a
#    minimal pressure, pumping the gas from a number of different
#    entry points.
# 
#    Sources:
#    D. De Wolf, "Optimisation de reseaux de transport de gas avec
#                 consideration des pertes de charge dans les gazoducs",
#                Ph. D. dissertation, CORE, Belgium, 1992, and
#    D. De Wolf, O. Janssens de Bisthoven and Y. Smeers,
#                "The simplex algorithm extended to piecewise linearly
#                 constrained problems II; an application to the gas
#                 transmission problem", CORE discussion paper 9103, 1991.
# 
# 
#    SDIF input: E. Loute and D. De Wolf, September 1992.
# 
#    classification = "C-CLQI2-RN-65-59"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'CORE1'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        [ig,ig_,_] = s2mpj_ii('COST',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('NODE0001',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'NODE0001')
        [ig,ig_,_] = s2mpj_ii('NODE0002',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'NODE0002')
        [ig,ig_,_] = s2mpj_ii('NODE0003',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'NODE0003')
        [ig,ig_,_] = s2mpj_ii('NODE0004',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'NODE0004')
        [ig,ig_,_] = s2mpj_ii('NODE0005',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'NODE0005')
        [ig,ig_,_] = s2mpj_ii('NODE0006',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'NODE0006')
        [ig,ig_,_] = s2mpj_ii('NODE0007',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'NODE0007')
        [ig,ig_,_] = s2mpj_ii('NODE0008',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'NODE0008')
        [ig,ig_,_] = s2mpj_ii('NODE0009',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'NODE0009')
        [ig,ig_,_] = s2mpj_ii('NODE0010',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'NODE0010')
        [ig,ig_,_] = s2mpj_ii('NODE0011',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'NODE0011')
        [ig,ig_,_] = s2mpj_ii('NODE0012',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'NODE0012')
        [ig,ig_,_] = s2mpj_ii('NODE0013',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'NODE0013')
        [ig,ig_,_] = s2mpj_ii('NODE0014',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'NODE0014')
        [ig,ig_,_] = s2mpj_ii('NODE0015',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'NODE0015')
        [ig,ig_,_] = s2mpj_ii('NODE0016',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'NODE0016')
        [ig,ig_,_] = s2mpj_ii('NODE0017',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'NODE0017')
        [ig,ig_,_] = s2mpj_ii('NODE0018',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'NODE0018')
        [ig,ig_,_] = s2mpj_ii('NODE0019',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'NODE0019')
        [ig,ig_,_] = s2mpj_ii('NODE0020',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'NODE0020')
        [ig,ig_,_] = s2mpj_ii('ARC00001',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'ARC00001')
        [ig,ig_,_] = s2mpj_ii('ARC00002',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'ARC00002')
        [ig,ig_,_] = s2mpj_ii('ARC00003',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'ARC00003')
        [ig,ig_,_] = s2mpj_ii('ARC00004',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'ARC00004')
        [ig,ig_,_] = s2mpj_ii('ARC00005',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'ARC00005')
        [ig,ig_,_] = s2mpj_ii('ARC00006',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'ARC00006')
        [ig,ig_,_] = s2mpj_ii('ARC00007',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'ARC00007')
        [ig,ig_,_] = s2mpj_ii('ARC00008',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'ARC00008')
        [ig,ig_,_] = s2mpj_ii('ARC00009',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'ARC00009')
        [ig,ig_,_] = s2mpj_ii('ARC00010',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'ARC00010')
        [ig,ig_,_] = s2mpj_ii('ARC00011',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'ARC00011')
        [ig,ig_,_] = s2mpj_ii('ARC00012',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'ARC00012')
        [ig,ig_,_] = s2mpj_ii('ARC00013',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'ARC00013')
        [ig,ig_,_] = s2mpj_ii('ARC00014',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'ARC00014')
        [ig,ig_,_] = s2mpj_ii('ARC00015',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'ARC00015')
        [ig,ig_,_] = s2mpj_ii('ARC00016',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'ARC00016')
        [ig,ig_,_] = s2mpj_ii('ARC00017',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'ARC00017')
        [ig,ig_,_] = s2mpj_ii('ARC00018',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'ARC00018')
        [ig,ig_,_] = s2mpj_ii('ARC00019',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'ARC00019')
        [ig,ig_,_] = s2mpj_ii('ARC00020',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'ARC00020')
        [ig,ig_,_] = s2mpj_ii('ARC00021',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'ARC00021')
        [ig,ig_,_] = s2mpj_ii('ARC00022',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'ARC00022')
        [ig,ig_,_] = s2mpj_ii('ARC00023',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'ARC00023')
        [ig,ig_,_] = s2mpj_ii('ARC00024',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'ARC00024')
        [ig,ig_,_] = s2mpj_ii('REGIO001',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'REGIO001')
        [ig,ig_,_] = s2mpj_ii('REGIO002',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'REGIO002')
        [ig,ig_,_] = s2mpj_ii('REGIO003',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'REGIO003')
        [ig,ig_,_] = s2mpj_ii('REGIO004',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'REGIO004')
        [ig,ig_,_] = s2mpj_ii('REGIO005',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'REGIO005')
        [ig,ig_,_] = s2mpj_ii('REGIO006',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'REGIO006')
        [ig,ig_,_] = s2mpj_ii('REGIO007',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'REGIO007')
        [ig,ig_,_] = s2mpj_ii('REGIO008',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'REGIO008')
        [ig,ig_,_] = s2mpj_ii('REGIO009',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'REGIO009')
        [ig,ig_,_] = s2mpj_ii('PROD0001',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'PROD0001')
        [ig,ig_,_] = s2mpj_ii('PROD0002',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'PROD0002')
        [ig,ig_,_] = s2mpj_ii('PROD0003',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'PROD0003')
        [ig,ig_,_] = s2mpj_ii('PROD0004',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'PROD0004')
        [ig,ig_,_] = s2mpj_ii('PROD0005',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'PROD0005')
        [ig,ig_,_] = s2mpj_ii('PROD0006',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'PROD0006')
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        ngrp   = len(ig_)
        [iv,ix_,_] = s2mpj_ii('FLOW0001',ix_)
        self.xnames=arrset(self.xnames,iv,'FLOW0001')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0001']])
        valA = np.append(valA,float(1.00000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0002']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('FLOW0002',ix_)
        self.xnames=arrset(self.xnames,iv,'FLOW0002')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0001']])
        valA = np.append(valA,float(1.00000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0002']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('FLOW0003',ix_)
        self.xnames=arrset(self.xnames,iv,'FLOW0003')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0002']])
        valA = np.append(valA,float(1.00000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0003']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('FLOW0004',ix_)
        self.xnames=arrset(self.xnames,iv,'FLOW0004')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0002']])
        valA = np.append(valA,float(1.00000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0003']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('FLOW0005',ix_)
        self.xnames=arrset(self.xnames,iv,'FLOW0005')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0003']])
        valA = np.append(valA,float(1.00000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0004']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('FLOW0006',ix_)
        self.xnames=arrset(self.xnames,iv,'FLOW0006')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0005']])
        valA = np.append(valA,float(1.00000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0006']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('FLOW0007',ix_)
        self.xnames=arrset(self.xnames,iv,'FLOW0007')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0006']])
        valA = np.append(valA,float(1.00000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0007']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('FLOW0008',ix_)
        self.xnames=arrset(self.xnames,iv,'FLOW0008')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0007']])
        valA = np.append(valA,float(1.00000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0004']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('FLOW0009',ix_)
        self.xnames=arrset(self.xnames,iv,'FLOW0009')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0004']])
        valA = np.append(valA,float(1.00000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0014']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('FLOW0010',ix_)
        self.xnames=arrset(self.xnames,iv,'FLOW0010')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0008']])
        valA = np.append(valA,float(1.00000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0009']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('FLOW0011',ix_)
        self.xnames=arrset(self.xnames,iv,'FLOW0011')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0008']])
        valA = np.append(valA,float(1.00000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0009']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('FLOW0012',ix_)
        self.xnames=arrset(self.xnames,iv,'FLOW0012')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0009']])
        valA = np.append(valA,float(1.00000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0010']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('FLOW0013',ix_)
        self.xnames=arrset(self.xnames,iv,'FLOW0013')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0009']])
        valA = np.append(valA,float(1.00000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0010']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('FLOW0014',ix_)
        self.xnames=arrset(self.xnames,iv,'FLOW0014')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0010']])
        valA = np.append(valA,float(1.00000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0011']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('FLOW0015',ix_)
        self.xnames=arrset(self.xnames,iv,'FLOW0015')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0010']])
        valA = np.append(valA,float(1.00000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0011']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('FLOW0016',ix_)
        self.xnames=arrset(self.xnames,iv,'FLOW0016')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0011']])
        valA = np.append(valA,float(1.00000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0012']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('FLOW0017',ix_)
        self.xnames=arrset(self.xnames,iv,'FLOW0017')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0012']])
        valA = np.append(valA,float(1.00000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0013']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('FLOW0018',ix_)
        self.xnames=arrset(self.xnames,iv,'FLOW0018')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0013']])
        valA = np.append(valA,float(1.00000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0014']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('FLOW0019',ix_)
        self.xnames=arrset(self.xnames,iv,'FLOW0019')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0014']])
        valA = np.append(valA,float(1.00000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0015']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('FLOW0020',ix_)
        self.xnames=arrset(self.xnames,iv,'FLOW0020')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0015']])
        valA = np.append(valA,float(1.00000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0016']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('FLOW0021',ix_)
        self.xnames=arrset(self.xnames,iv,'FLOW0021')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0011']])
        valA = np.append(valA,float(1.00000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0017']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('FLOW0022',ix_)
        self.xnames=arrset(self.xnames,iv,'FLOW0022')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0017']])
        valA = np.append(valA,float(1.00000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0018']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('FLOW0023',ix_)
        self.xnames=arrset(self.xnames,iv,'FLOW0023')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0018']])
        valA = np.append(valA,float(1.00000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0019']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('FLOW0024',ix_)
        self.xnames=arrset(self.xnames,iv,'FLOW0024')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0019']])
        valA = np.append(valA,float(1.00000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0020']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('DEM00001',ix_)
        self.xnames=arrset(self.xnames,iv,'DEM00001')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0003']])
        valA = np.append(valA,float(1.00000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['REGIO001']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('DEM00002',ix_)
        self.xnames=arrset(self.xnames,iv,'DEM00002')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0006']])
        valA = np.append(valA,float(1.00000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['REGIO002']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('DEM00003',ix_)
        self.xnames=arrset(self.xnames,iv,'DEM00003')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0007']])
        valA = np.append(valA,float(1.00000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['REGIO003']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('DEM00004',ix_)
        self.xnames=arrset(self.xnames,iv,'DEM00004')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0010']])
        valA = np.append(valA,float(1.00000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['REGIO004']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('DEM00005',ix_)
        self.xnames=arrset(self.xnames,iv,'DEM00005')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0012']])
        valA = np.append(valA,float(1.00000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['REGIO005']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('DEM00006',ix_)
        self.xnames=arrset(self.xnames,iv,'DEM00006')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0015']])
        valA = np.append(valA,float(1.00000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['REGIO006']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('DEM00007',ix_)
        self.xnames=arrset(self.xnames,iv,'DEM00007')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0016']])
        valA = np.append(valA,float(1.00000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['REGIO007']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('DEM00008',ix_)
        self.xnames=arrset(self.xnames,iv,'DEM00008')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0019']])
        valA = np.append(valA,float(1.00000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['REGIO008']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('DEM00009',ix_)
        self.xnames=arrset(self.xnames,iv,'DEM00009')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0020']])
        valA = np.append(valA,float(1.00000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['REGIO009']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('SUPP0001',ix_)
        self.xnames=arrset(self.xnames,iv,'SUPP0001')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['PROD0001']])
        valA = np.append(valA,float(1.00000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0001']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('SUPP0002',ix_)
        self.xnames=arrset(self.xnames,iv,'SUPP0002')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['PROD0002']])
        valA = np.append(valA,float(1.00000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0002']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('SUPP0003',ix_)
        self.xnames=arrset(self.xnames,iv,'SUPP0003')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['PROD0003']])
        valA = np.append(valA,float(1.00000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0005']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('SUPP0004',ix_)
        self.xnames=arrset(self.xnames,iv,'SUPP0004')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['PROD0004']])
        valA = np.append(valA,float(1.00000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0008']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('SUPP0005',ix_)
        self.xnames=arrset(self.xnames,iv,'SUPP0005')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['PROD0005']])
        valA = np.append(valA,float(1.00000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0013']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('SUPP0006',ix_)
        self.xnames=arrset(self.xnames,iv,'SUPP0006')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['PROD0006']])
        valA = np.append(valA,float(1.00000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['NODE0014']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('PROD0001',ix_)
        self.xnames=arrset(self.xnames,iv,'PROD0001')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['COST']])
        valA = np.append(valA,float(2.28000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['PROD0001']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('PROD0002',ix_)
        self.xnames=arrset(self.xnames,iv,'PROD0002')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['COST']])
        valA = np.append(valA,float(2.28000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['PROD0002']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('PROD0003',ix_)
        self.xnames=arrset(self.xnames,iv,'PROD0003')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['COST']])
        valA = np.append(valA,float(2.28000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['PROD0003']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('PROD0004',ix_)
        self.xnames=arrset(self.xnames,iv,'PROD0004')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['COST']])
        valA = np.append(valA,float(1.68000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['PROD0004']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('PROD0005',ix_)
        self.xnames=arrset(self.xnames,iv,'PROD0005')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['COST']])
        valA = np.append(valA,float(1.68000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['PROD0005']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('PROD0006',ix_)
        self.xnames=arrset(self.xnames,iv,'PROD0006')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['COST']])
        valA = np.append(valA,float(1.68000E+00))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['PROD0006']])
        valA = np.append(valA,float(-1.00000E+00))
        [iv,ix_,_] = s2mpj_ii('PI000001',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000001')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00001']])
        valA = np.append(valA,float(-9.07027E+00))
        [iv,ix_,_] = s2mpj_ii('PI000001',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000001')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00002']])
        valA = np.append(valA,float(-9.07027E+00))
        [iv,ix_,_] = s2mpj_ii('PI000002',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000002')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00001']])
        valA = np.append(valA,float(9.07027E+00))
        [iv,ix_,_] = s2mpj_ii('PI000002',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000002')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00002']])
        valA = np.append(valA,float(9.07027E+00))
        [iv,ix_,_] = s2mpj_ii('PI000002',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000002')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00003']])
        valA = np.append(valA,float(-6.04685E+00))
        [iv,ix_,_] = s2mpj_ii('PI000002',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000002')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00004']])
        valA = np.append(valA,float(-6.04685E+00))
        [iv,ix_,_] = s2mpj_ii('PI000003',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000003')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00003']])
        valA = np.append(valA,float(6.04685E+00))
        [iv,ix_,_] = s2mpj_ii('PI000003',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000003')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00004']])
        valA = np.append(valA,float(6.04685E+00))
        [iv,ix_,_] = s2mpj_ii('PI000003',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000003')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00005']])
        valA = np.append(valA,float(-1.39543E+00))
        [iv,ix_,_] = s2mpj_ii('PI000004',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000004')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00005']])
        valA = np.append(valA,float(1.39543E+00))
        [iv,ix_,_] = s2mpj_ii('PI000004',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000004')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00008']])
        valA = np.append(valA,float(2.26895E-01))
        [iv,ix_,_] = s2mpj_ii('PI000004',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000004')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00009']])
        valA = np.append(valA,float(-6.59656E-01))
        [iv,ix_,_] = s2mpj_ii('PI000005',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000005')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00006']])
        valA = np.append(valA,float(-1.00256E-01))
        [iv,ix_,_] = s2mpj_ii('PI000006',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000006')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00006']])
        valA = np.append(valA,float(1.00256E-01))
        [iv,ix_,_] = s2mpj_ii('PI000006',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000006')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00007']])
        valA = np.append(valA,float(-1.48655E-01))
        [iv,ix_,_] = s2mpj_ii('PI000007',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000007')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00007']])
        valA = np.append(valA,float(1.48655E-01))
        [iv,ix_,_] = s2mpj_ii('PI000007',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000007')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00008']])
        valA = np.append(valA,float(-2.26895E-01))
        [iv,ix_,_] = s2mpj_ii('PI000008',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000008')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00010']])
        valA = np.append(valA,float(-7.25622E+00))
        [iv,ix_,_] = s2mpj_ii('PI000008',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000008')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00011']])
        valA = np.append(valA,float(-1.08033E-01))
        [iv,ix_,_] = s2mpj_ii('PI000009',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000009')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00010']])
        valA = np.append(valA,float(7.25622E+00))
        [iv,ix_,_] = s2mpj_ii('PI000009',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000009')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00011']])
        valA = np.append(valA,float(1.08033E-01))
        [iv,ix_,_] = s2mpj_ii('PI000009',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000009')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00012']])
        valA = np.append(valA,float(-1.81405E+00))
        [iv,ix_,_] = s2mpj_ii('PI000009',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000009')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00013']])
        valA = np.append(valA,float(-2.70084E-02))
        [iv,ix_,_] = s2mpj_ii('PI000010',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000010')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00012']])
        valA = np.append(valA,float(1.81405E+00))
        [iv,ix_,_] = s2mpj_ii('PI000010',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000010')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00013']])
        valA = np.append(valA,float(2.70084E-02))
        [iv,ix_,_] = s2mpj_ii('PI000010',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000010')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00014']])
        valA = np.append(valA,float(-1.45124E+00))
        [iv,ix_,_] = s2mpj_ii('PI000010',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000010')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00015']])
        valA = np.append(valA,float(-2.16067E-02))
        [iv,ix_,_] = s2mpj_ii('PI000011',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000011')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00014']])
        valA = np.append(valA,float(1.45124E+00))
        [iv,ix_,_] = s2mpj_ii('PI000011',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000011')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00015']])
        valA = np.append(valA,float(2.16067E-02))
        [iv,ix_,_] = s2mpj_ii('PI000011',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000011')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00016']])
        valA = np.append(valA,float(-8.63836E-01))
        [iv,ix_,_] = s2mpj_ii('PI000011',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000011')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00021']])
        valA = np.append(valA,float(-5.14445E-02))
        [iv,ix_,_] = s2mpj_ii('PI000012',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000012')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00016']])
        valA = np.append(valA,float(8.63836E-01))
        [iv,ix_,_] = s2mpj_ii('PI000012',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000012')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00017']])
        valA = np.append(valA,float(-9.07027E-01))
        [iv,ix_,_] = s2mpj_ii('PI000013',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000013')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00017']])
        valA = np.append(valA,float(9.07027E-01))
        [iv,ix_,_] = s2mpj_ii('PI000013',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000013')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00018']])
        valA = np.append(valA,float(-7.25622E+00))
        [iv,ix_,_] = s2mpj_ii('PI000014',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000014')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00009']])
        valA = np.append(valA,float(6.59656E-01))
        [iv,ix_,_] = s2mpj_ii('PI000014',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000014')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00018']])
        valA = np.append(valA,float(7.25622E+00))
        [iv,ix_,_] = s2mpj_ii('PI000014',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000014')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00019']])
        valA = np.append(valA,float(-3.62811E+00))
        [iv,ix_,_] = s2mpj_ii('PI000015',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000015')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00019']])
        valA = np.append(valA,float(3.62811E+00))
        [iv,ix_,_] = s2mpj_ii('PI000015',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000015')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00020']])
        valA = np.append(valA,float(-1.45124E+00))
        [iv,ix_,_] = s2mpj_ii('PI000016',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000016')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00020']])
        valA = np.append(valA,float(1.45124E+00))
        [iv,ix_,_] = s2mpj_ii('PI000017',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000017')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00021']])
        valA = np.append(valA,float(5.14445E-02))
        [iv,ix_,_] = s2mpj_ii('PI000017',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000017')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00022']])
        valA = np.append(valA,float(-6.41977E-03))
        [iv,ix_,_] = s2mpj_ii('PI000018',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000018')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00022']])
        valA = np.append(valA,float(6.41977E-03))
        [iv,ix_,_] = s2mpj_ii('PI000018',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000018')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00023']])
        valA = np.append(valA,float(-1.70320E-03))
        [iv,ix_,_] = s2mpj_ii('PI000019',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000019')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00023']])
        valA = np.append(valA,float(1.70320E-03))
        [iv,ix_,_] = s2mpj_ii('PI000019',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000019')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00024']])
        valA = np.append(valA,float(-2.78190E-02))
        [iv,ix_,_] = s2mpj_ii('PI000020',ix_)
        self.xnames=arrset(self.xnames,iv,'PI000020')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['ARC00024']])
        valA = np.append(valA,float(2.78190E-02))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
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
        self.gconst = arrset(self.gconst,ig_['REGIO001'],float(-3.91800E+00))
        self.gconst = arrset(self.gconst,ig_['REGIO002'],float(-4.03400E+00))
        self.gconst = arrset(self.gconst,ig_['REGIO003'],float(-5.25600E+00))
        self.gconst = arrset(self.gconst,ig_['REGIO004'],float(-6.36500E+00))
        self.gconst = arrset(self.gconst,ig_['REGIO005'],float(-2.12000E+00))
        self.gconst = arrset(self.gconst,ig_['REGIO006'],float(-6.84800E+00))
        self.gconst = arrset(self.gconst,ig_['REGIO007'],float(-1.56160E+01))
        self.gconst = arrset(self.gconst,ig_['REGIO008'],float(-2.22000E-01))
        self.gconst = arrset(self.gconst,ig_['REGIO009'],float(-1.91900E+00))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        self.xlower[ix_['FLOW0001']] = -2.20120E+02
        self.xupper[ix_['FLOW0001']] = 2.20120E+02
        self.xlower[ix_['FLOW0002']] = -2.20120E+02
        self.xupper[ix_['FLOW0002']] = 2.20120E+02
        self.xlower[ix_['FLOW0003']] = -2.20120E+02
        self.xupper[ix_['FLOW0003']] = 2.20120E+02
        self.xlower[ix_['FLOW0004']] = -2.20120E+02
        self.xupper[ix_['FLOW0004']] = 2.20120E+02
        self.xlower[ix_['FLOW0005']] = -2.20120E+02
        self.xupper[ix_['FLOW0005']] = 2.20120E+02
        self.xlower[ix_['FLOW0006']] = -2.20120E+02
        self.xupper[ix_['FLOW0006']] = 2.20120E+02
        self.xlower[ix_['FLOW0007']] = -2.20120E+02
        self.xupper[ix_['FLOW0007']] = 2.20120E+02
        self.xlower[ix_['FLOW0008']] = -2.20120E+02
        self.xupper[ix_['FLOW0008']] = 2.20120E+02
        self.xlower[ix_['FLOW0009']] = -2.20120E+02
        self.xupper[ix_['FLOW0009']] = 2.20120E+02
        self.xlower[ix_['FLOW0010']] = 0.00000E+00
        self.xupper[ix_['FLOW0010']] = 2.20120E+02
        self.xlower[ix_['FLOW0011']] = 0.00000E+00
        self.xupper[ix_['FLOW0011']] = 2.20120E+02
        self.xlower[ix_['FLOW0012']] = -2.20120E+02
        self.xupper[ix_['FLOW0012']] = 2.20120E+02
        self.xlower[ix_['FLOW0013']] = -2.20120E+02
        self.xupper[ix_['FLOW0013']] = 2.20120E+02
        self.xlower[ix_['FLOW0014']] = -2.20120E+02
        self.xupper[ix_['FLOW0014']] = 2.20120E+02
        self.xlower[ix_['FLOW0015']] = -2.20120E+02
        self.xupper[ix_['FLOW0015']] = 2.20120E+02
        self.xlower[ix_['FLOW0016']] = -2.20120E+02
        self.xupper[ix_['FLOW0016']] = 2.20120E+02
        self.xlower[ix_['FLOW0017']] = -2.20120E+02
        self.xupper[ix_['FLOW0017']] = 2.20120E+02
        self.xlower[ix_['FLOW0018']] = -2.20120E+02
        self.xupper[ix_['FLOW0018']] = 2.20120E+02
        self.xlower[ix_['FLOW0019']] = -2.20120E+02
        self.xupper[ix_['FLOW0019']] = 2.20120E+02
        self.xlower[ix_['FLOW0020']] = -2.20120E+02
        self.xupper[ix_['FLOW0020']] = 2.20120E+02
        self.xlower[ix_['FLOW0021']] = -2.20120E+02
        self.xupper[ix_['FLOW0021']] = 2.20120E+02
        self.xlower[ix_['FLOW0022']] = 0.00000E+00
        self.xupper[ix_['FLOW0022']] = 2.20120E+02
        self.xlower[ix_['FLOW0023']] = -2.20120E+02
        self.xupper[ix_['FLOW0023']] = 2.20120E+02
        self.xlower[ix_['FLOW0024']] = -2.20120E+02
        self.xupper[ix_['FLOW0024']] = 2.20120E+02
        self.xupper[ix_['SUPP0001']] = 1.15940E+01
        self.xupper[ix_['SUPP0002']] = 8.40000E+00
        self.xupper[ix_['SUPP0003']] = 4.80000E+00
        self.xupper[ix_['SUPP0004']] = 2.20120E+01
        self.xupper[ix_['SUPP0005']] = 1.20000E+00
        self.xupper[ix_['SUPP0006']] = 9.60000E-01
        self.xlower[ix_['PROD0001']] = 8.87000E+00
        self.xupper[ix_['PROD0001']] = 1.15940E+01
        self.xlower[ix_['PROD0002']] = 0.00000E+00
        self.xupper[ix_['PROD0002']] = 8.40000E+00
        self.xlower[ix_['PROD0003']] = 0.00000E+00
        self.xupper[ix_['PROD0003']] = 4.80000E+00
        self.xlower[ix_['PROD0004']] = 2.03440E+01
        self.xupper[ix_['PROD0004']] = 2.20120E+01
        self.xlower[ix_['PROD0005']] = 0.00000E+00
        self.xupper[ix_['PROD0005']] = 1.20000E+00
        self.xlower[ix_['PROD0006']] = 0.00000E+00
        self.xupper[ix_['PROD0006']] = 9.60000E-01
        self.xlower[ix_['PI000001']] = 0.00000E+00
        self.xupper[ix_['PI000001']] = 5.92900E+03
        self.xlower[ix_['PI000002']] = 0.00000E+00
        self.xupper[ix_['PI000002']] = 5.92900E+03
        self.xlower[ix_['PI000003']] = 9.00000E+02
        self.xupper[ix_['PI000003']] = 6.40000E+03
        self.xlower[ix_['PI000004']] = 0.00000E+00
        self.xupper[ix_['PI000004']] = 6.40000E+03
        self.xlower[ix_['PI000005']] = 0.00000E+00
        self.xupper[ix_['PI000005']] = 5.92900E+03
        self.xlower[ix_['PI000006']] = 9.00000E+02
        self.xupper[ix_['PI000006']] = 6.40000E+03
        self.xlower[ix_['PI000007']] = 9.00000E+02
        self.xupper[ix_['PI000007']] = 6.40000E+03
        self.xlower[ix_['PI000008']] = 2.50000E+03
        self.xupper[ix_['PI000008']] = 4.38244E+03
        self.xlower[ix_['PI000009']] = 0.00000E+00
        self.xupper[ix_['PI000009']] = 4.38244E+03
        self.xlower[ix_['PI000010']] = 9.00000E+02
        self.xupper[ix_['PI000010']] = 4.38244E+03
        self.xlower[ix_['PI000011']] = 0.00000E+00
        self.xupper[ix_['PI000011']] = 4.38244E+03
        self.xlower[ix_['PI000012']] = 9.00000E+02
        self.xupper[ix_['PI000012']] = 4.38244E+03
        self.xlower[ix_['PI000013']] = 0.00000E+00
        self.xupper[ix_['PI000013']] = 4.38244E+03
        self.xlower[ix_['PI000014']] = 0.00000E+00
        self.xupper[ix_['PI000014']] = 4.38244E+03
        self.xlower[ix_['PI000015']] = 9.00000E+02
        self.xupper[ix_['PI000015']] = 4.38244E+03
        self.xlower[ix_['PI000016']] = 2.50000E+03
        self.xupper[ix_['PI000016']] = 4.38244E+03
        self.xlower[ix_['PI000017']] = 0.00000E+00
        self.xupper[ix_['PI000017']] = 4.38244E+03
        self.xlower[ix_['PI000018']] = 0.00000E+00
        self.xupper[ix_['PI000018']] = 3.96900E+03
        self.xlower[ix_['PI000019']] = 6.25000E+02
        self.xupper[ix_['PI000019']] = 4.38244E+03
        self.xlower[ix_['PI000020']] = 6.25000E+02
        self.xupper[ix_['PI000020']] = 4.38244E+03
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        self.x0[ix_['FLOW0001']] = float(5.79700E+00)
        self.x0[ix_['FLOW0002']] = float(5.79700E+00)
        self.x0[ix_['FLOW0003']] = float(9.99700E+00)
        self.x0[ix_['FLOW0004']] = float(9.99700E+00)
        self.x0[ix_['FLOW0005']] = float(1.60760E+01)
        self.x0[ix_['FLOW0006']] = float(4.80000E+00)
        self.x0[ix_['FLOW0007']] = float(7.66000E-01)
        self.x0[ix_['FLOW0008']] = float(-4.49000E+00)
        self.x0[ix_['FLOW0009']] = float(1.15860E+01)
        self.x0[ix_['FLOW0010']] = float(1.72404E+01)
        self.x0[ix_['FLOW0011']] = float(2.10363E+00)
        self.x0[ix_['FLOW0012']] = float(1.72404E+01)
        self.x0[ix_['FLOW0013']] = float(2.10363E+00)
        self.x0[ix_['FLOW0014']] = float(1.15676E+01)
        self.x0[ix_['FLOW0015']] = float(1.41145E+00)
        self.x0[ix_['FLOW0016']] = float(1.08380E+01)
        self.x0[ix_['FLOW0017']] = float(8.71800E+00)
        self.x0[ix_['FLOW0018']] = float(9.91800E+00)
        self.x0[ix_['FLOW0019']] = float(2.24640E+01)
        self.x0[ix_['FLOW0020']] = float(1.56160E+01)
        self.x0[ix_['FLOW0021']] = float(2.14100E+00)
        self.x0[ix_['FLOW0022']] = float(2.14100E+00)
        self.x0[ix_['FLOW0023']] = float(2.14100E+00)
        self.x0[ix_['FLOW0024']] = float(1.91900E+00)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQR', iet_)
        elftv = loaset(elftv,it,0,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        ename = 'F00001SQ'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQR')
        ielftype = arrset(ielftype,ie,iet_["eSQR"])
        vname = 'FLOW0001'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F00002SQ'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQR')
        ielftype = arrset(ielftype,ie,iet_["eSQR"])
        vname = 'FLOW0002'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F00003SQ'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQR')
        ielftype = arrset(ielftype,ie,iet_["eSQR"])
        vname = 'FLOW0003'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F00004SQ'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQR')
        ielftype = arrset(ielftype,ie,iet_["eSQR"])
        vname = 'FLOW0004'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F00005SQ'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQR')
        ielftype = arrset(ielftype,ie,iet_["eSQR"])
        vname = 'FLOW0005'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F00006SQ'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQR')
        ielftype = arrset(ielftype,ie,iet_["eSQR"])
        vname = 'FLOW0006'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F00007SQ'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQR')
        ielftype = arrset(ielftype,ie,iet_["eSQR"])
        vname = 'FLOW0007'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F00008SQ'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQR')
        ielftype = arrset(ielftype,ie,iet_["eSQR"])
        vname = 'FLOW0008'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F00009SQ'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQR')
        ielftype = arrset(ielftype,ie,iet_["eSQR"])
        vname = 'FLOW0009'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F00010SQ'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQR')
        ielftype = arrset(ielftype,ie,iet_["eSQR"])
        vname = 'FLOW0010'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F00011SQ'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQR')
        ielftype = arrset(ielftype,ie,iet_["eSQR"])
        vname = 'FLOW0011'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F00012SQ'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQR')
        ielftype = arrset(ielftype,ie,iet_["eSQR"])
        vname = 'FLOW0012'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F00013SQ'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQR')
        ielftype = arrset(ielftype,ie,iet_["eSQR"])
        vname = 'FLOW0013'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F00014SQ'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQR')
        ielftype = arrset(ielftype,ie,iet_["eSQR"])
        vname = 'FLOW0014'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F00015SQ'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQR')
        ielftype = arrset(ielftype,ie,iet_["eSQR"])
        vname = 'FLOW0015'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F00016SQ'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQR')
        ielftype = arrset(ielftype,ie,iet_["eSQR"])
        vname = 'FLOW0016'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F00017SQ'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQR')
        ielftype = arrset(ielftype,ie,iet_["eSQR"])
        vname = 'FLOW0017'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F00018SQ'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQR')
        ielftype = arrset(ielftype,ie,iet_["eSQR"])
        vname = 'FLOW0018'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F00019SQ'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQR')
        ielftype = arrset(ielftype,ie,iet_["eSQR"])
        vname = 'FLOW0019'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F00020SQ'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQR')
        ielftype = arrset(ielftype,ie,iet_["eSQR"])
        vname = 'FLOW0020'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F00021SQ'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQR')
        ielftype = arrset(ielftype,ie,iet_["eSQR"])
        vname = 'FLOW0021'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F00022SQ'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQR')
        ielftype = arrset(ielftype,ie,iet_["eSQR"])
        vname = 'FLOW0022'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F00023SQ'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQR')
        ielftype = arrset(ielftype,ie,iet_["eSQR"])
        vname = 'FLOW0023'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F00024SQ'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQR')
        ielftype = arrset(ielftype,ie,iet_["eSQR"])
        vname = 'FLOW0024'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        ig = ig_['ARC00001']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F00001SQ'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['ARC00002']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F00002SQ'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['ARC00003']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F00003SQ'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['ARC00004']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F00004SQ'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['ARC00005']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F00005SQ'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['ARC00006']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F00006SQ'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['ARC00007']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F00007SQ'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['ARC00008']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F00008SQ'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['ARC00009']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F00009SQ'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['ARC00010']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F00010SQ'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['ARC00011']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F00011SQ'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['ARC00012']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F00012SQ'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['ARC00013']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F00013SQ'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['ARC00014']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F00014SQ'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['ARC00015']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F00015SQ'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['ARC00016']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F00016SQ'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['ARC00017']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F00017SQ'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['ARC00018']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F00018SQ'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['ARC00019']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F00019SQ'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['ARC00020']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F00020SQ'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['ARC00021']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F00021SQ'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['ARC00022']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F00022SQ'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['ARC00023']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F00023SQ'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['ARC00024']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F00024SQ'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.cupper[np.arange(self.nle)] = np.zeros((self.nle,1))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.clower[np.arange(self.nle+self.neq,self.m)] = np.zeros((self.nge,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-CLQI2-RN-65-59"
        self.objderlvl = 2
        self.conderlvl = [2]


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eSQR(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        XGE0 = EV_[0]>=0.0e+0
        if XGE0!=0:
            G = 2.0e+0*EV_[0]
        if XGE0==0:
            G = -2.0e+0*EV_[0]
        if XGE0!=0:
            H = 2.0e+0
        if XGE0==0:
            H = -2.0e+0
        f_   = EV_[0]*np.absolute(EV_[0])
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = G
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = H
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

