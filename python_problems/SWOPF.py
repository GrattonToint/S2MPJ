from s2mpjlib import *
class  SWOPF(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    An optimal electrical powerflow system design problem from Switzerland.
# 
#    Source:
#    a contribution to fullfill the LANCELOT academic licence agreement.
# 
#    SIF input: R. Bacher, Dept of Electrical Engineering, ETH Zurich, 
#               November 1994.
# 
#    classification = "LQR2-RN-83-92"
# 
#    Number of nodes       =   7
#    Number of branches    =   7
#    Number of generators  =   3
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'SWOPF'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'SWOPF'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['FIRST'] = 1
        v_['NOBRANCHES'] = 7
        v_['NOSHUNTS'] = 0
        v_['NOTRAFOS'] = 3
        v_['NOBUSSES'] = 7
        v_['NOGEN'] = 3
        v_['NOGENBK'] = 3
        v_['NOAREAS'] = 0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2mpj_ii('VE0001',ix_)
        pb.xnames=arrset(pb.xnames,iv,'VE0001')
        [iv,ix_,_] = s2mpj_ii('VF0001',ix_)
        pb.xnames=arrset(pb.xnames,iv,'VF0001')
        [iv,ix_,_] = s2mpj_ii('V20001',ix_)
        pb.xnames=arrset(pb.xnames,iv,'V20001')
        [iv,ix_,_] = s2mpj_ii('VE0002',ix_)
        pb.xnames=arrset(pb.xnames,iv,'VE0002')
        [iv,ix_,_] = s2mpj_ii('VF0002',ix_)
        pb.xnames=arrset(pb.xnames,iv,'VF0002')
        [iv,ix_,_] = s2mpj_ii('V20002',ix_)
        pb.xnames=arrset(pb.xnames,iv,'V20002')
        [iv,ix_,_] = s2mpj_ii('VE0003',ix_)
        pb.xnames=arrset(pb.xnames,iv,'VE0003')
        [iv,ix_,_] = s2mpj_ii('VF0003',ix_)
        pb.xnames=arrset(pb.xnames,iv,'VF0003')
        [iv,ix_,_] = s2mpj_ii('V20003',ix_)
        pb.xnames=arrset(pb.xnames,iv,'V20003')
        [iv,ix_,_] = s2mpj_ii('VE0004',ix_)
        pb.xnames=arrset(pb.xnames,iv,'VE0004')
        [iv,ix_,_] = s2mpj_ii('VF0004',ix_)
        pb.xnames=arrset(pb.xnames,iv,'VF0004')
        [iv,ix_,_] = s2mpj_ii('V20004',ix_)
        pb.xnames=arrset(pb.xnames,iv,'V20004')
        [iv,ix_,_] = s2mpj_ii('VE0005',ix_)
        pb.xnames=arrset(pb.xnames,iv,'VE0005')
        [iv,ix_,_] = s2mpj_ii('VF0005',ix_)
        pb.xnames=arrset(pb.xnames,iv,'VF0005')
        [iv,ix_,_] = s2mpj_ii('V20005',ix_)
        pb.xnames=arrset(pb.xnames,iv,'V20005')
        [iv,ix_,_] = s2mpj_ii('VE0006',ix_)
        pb.xnames=arrset(pb.xnames,iv,'VE0006')
        [iv,ix_,_] = s2mpj_ii('VF0006',ix_)
        pb.xnames=arrset(pb.xnames,iv,'VF0006')
        [iv,ix_,_] = s2mpj_ii('V20006',ix_)
        pb.xnames=arrset(pb.xnames,iv,'V20006')
        [iv,ix_,_] = s2mpj_ii('VE0007',ix_)
        pb.xnames=arrset(pb.xnames,iv,'VE0007')
        [iv,ix_,_] = s2mpj_ii('VF0007',ix_)
        pb.xnames=arrset(pb.xnames,iv,'VF0007')
        [iv,ix_,_] = s2mpj_ii('V20007',ix_)
        pb.xnames=arrset(pb.xnames,iv,'V20007')
        [iv,ix_,_] = s2mpj_ii('EI0001',ix_)
        pb.xnames=arrset(pb.xnames,iv,'EI0001')
        [iv,ix_,_] = s2mpj_ii('FI0001',ix_)
        pb.xnames=arrset(pb.xnames,iv,'FI0001')
        [iv,ix_,_] = s2mpj_ii('EJ0001',ix_)
        pb.xnames=arrset(pb.xnames,iv,'EJ0001')
        [iv,ix_,_] = s2mpj_ii('FJ0001',ix_)
        pb.xnames=arrset(pb.xnames,iv,'FJ0001')
        [iv,ix_,_] = s2mpj_ii('PI0001',ix_)
        pb.xnames=arrset(pb.xnames,iv,'PI0001')
        [iv,ix_,_] = s2mpj_ii('QI0001',ix_)
        pb.xnames=arrset(pb.xnames,iv,'QI0001')
        [iv,ix_,_] = s2mpj_ii('PJ0001',ix_)
        pb.xnames=arrset(pb.xnames,iv,'PJ0001')
        [iv,ix_,_] = s2mpj_ii('QJ0001',ix_)
        pb.xnames=arrset(pb.xnames,iv,'QJ0001')
        [iv,ix_,_] = s2mpj_ii('EI0002',ix_)
        pb.xnames=arrset(pb.xnames,iv,'EI0002')
        [iv,ix_,_] = s2mpj_ii('FI0002',ix_)
        pb.xnames=arrset(pb.xnames,iv,'FI0002')
        [iv,ix_,_] = s2mpj_ii('EJ0002',ix_)
        pb.xnames=arrset(pb.xnames,iv,'EJ0002')
        [iv,ix_,_] = s2mpj_ii('FJ0002',ix_)
        pb.xnames=arrset(pb.xnames,iv,'FJ0002')
        [iv,ix_,_] = s2mpj_ii('PI0002',ix_)
        pb.xnames=arrset(pb.xnames,iv,'PI0002')
        [iv,ix_,_] = s2mpj_ii('QI0002',ix_)
        pb.xnames=arrset(pb.xnames,iv,'QI0002')
        [iv,ix_,_] = s2mpj_ii('PJ0002',ix_)
        pb.xnames=arrset(pb.xnames,iv,'PJ0002')
        [iv,ix_,_] = s2mpj_ii('QJ0002',ix_)
        pb.xnames=arrset(pb.xnames,iv,'QJ0002')
        [iv,ix_,_] = s2mpj_ii('EI0003',ix_)
        pb.xnames=arrset(pb.xnames,iv,'EI0003')
        [iv,ix_,_] = s2mpj_ii('FI0003',ix_)
        pb.xnames=arrset(pb.xnames,iv,'FI0003')
        [iv,ix_,_] = s2mpj_ii('EJ0003',ix_)
        pb.xnames=arrset(pb.xnames,iv,'EJ0003')
        [iv,ix_,_] = s2mpj_ii('FJ0003',ix_)
        pb.xnames=arrset(pb.xnames,iv,'FJ0003')
        [iv,ix_,_] = s2mpj_ii('PI0003',ix_)
        pb.xnames=arrset(pb.xnames,iv,'PI0003')
        [iv,ix_,_] = s2mpj_ii('QI0003',ix_)
        pb.xnames=arrset(pb.xnames,iv,'QI0003')
        [iv,ix_,_] = s2mpj_ii('PJ0003',ix_)
        pb.xnames=arrset(pb.xnames,iv,'PJ0003')
        [iv,ix_,_] = s2mpj_ii('QJ0003',ix_)
        pb.xnames=arrset(pb.xnames,iv,'QJ0003')
        [iv,ix_,_] = s2mpj_ii('EI0004',ix_)
        pb.xnames=arrset(pb.xnames,iv,'EI0004')
        [iv,ix_,_] = s2mpj_ii('FI0004',ix_)
        pb.xnames=arrset(pb.xnames,iv,'FI0004')
        [iv,ix_,_] = s2mpj_ii('EJ0004',ix_)
        pb.xnames=arrset(pb.xnames,iv,'EJ0004')
        [iv,ix_,_] = s2mpj_ii('FJ0004',ix_)
        pb.xnames=arrset(pb.xnames,iv,'FJ0004')
        [iv,ix_,_] = s2mpj_ii('PI0004',ix_)
        pb.xnames=arrset(pb.xnames,iv,'PI0004')
        [iv,ix_,_] = s2mpj_ii('QI0004',ix_)
        pb.xnames=arrset(pb.xnames,iv,'QI0004')
        [iv,ix_,_] = s2mpj_ii('PJ0004',ix_)
        pb.xnames=arrset(pb.xnames,iv,'PJ0004')
        [iv,ix_,_] = s2mpj_ii('QJ0004',ix_)
        pb.xnames=arrset(pb.xnames,iv,'QJ0004')
        [iv,ix_,_] = s2mpj_ii('EI0005',ix_)
        pb.xnames=arrset(pb.xnames,iv,'EI0005')
        [iv,ix_,_] = s2mpj_ii('FI0005',ix_)
        pb.xnames=arrset(pb.xnames,iv,'FI0005')
        [iv,ix_,_] = s2mpj_ii('EJ0005',ix_)
        pb.xnames=arrset(pb.xnames,iv,'EJ0005')
        [iv,ix_,_] = s2mpj_ii('FJ0005',ix_)
        pb.xnames=arrset(pb.xnames,iv,'FJ0005')
        [iv,ix_,_] = s2mpj_ii('PI0005',ix_)
        pb.xnames=arrset(pb.xnames,iv,'PI0005')
        [iv,ix_,_] = s2mpj_ii('QI0005',ix_)
        pb.xnames=arrset(pb.xnames,iv,'QI0005')
        [iv,ix_,_] = s2mpj_ii('PJ0005',ix_)
        pb.xnames=arrset(pb.xnames,iv,'PJ0005')
        [iv,ix_,_] = s2mpj_ii('QJ0005',ix_)
        pb.xnames=arrset(pb.xnames,iv,'QJ0005')
        [iv,ix_,_] = s2mpj_ii('EI0006',ix_)
        pb.xnames=arrset(pb.xnames,iv,'EI0006')
        [iv,ix_,_] = s2mpj_ii('FI0006',ix_)
        pb.xnames=arrset(pb.xnames,iv,'FI0006')
        [iv,ix_,_] = s2mpj_ii('EJ0006',ix_)
        pb.xnames=arrset(pb.xnames,iv,'EJ0006')
        [iv,ix_,_] = s2mpj_ii('FJ0006',ix_)
        pb.xnames=arrset(pb.xnames,iv,'FJ0006')
        [iv,ix_,_] = s2mpj_ii('PI0006',ix_)
        pb.xnames=arrset(pb.xnames,iv,'PI0006')
        [iv,ix_,_] = s2mpj_ii('QI0006',ix_)
        pb.xnames=arrset(pb.xnames,iv,'QI0006')
        [iv,ix_,_] = s2mpj_ii('PJ0006',ix_)
        pb.xnames=arrset(pb.xnames,iv,'PJ0006')
        [iv,ix_,_] = s2mpj_ii('QJ0006',ix_)
        pb.xnames=arrset(pb.xnames,iv,'QJ0006')
        [iv,ix_,_] = s2mpj_ii('EI0007',ix_)
        pb.xnames=arrset(pb.xnames,iv,'EI0007')
        [iv,ix_,_] = s2mpj_ii('FI0007',ix_)
        pb.xnames=arrset(pb.xnames,iv,'FI0007')
        [iv,ix_,_] = s2mpj_ii('EJ0007',ix_)
        pb.xnames=arrset(pb.xnames,iv,'EJ0007')
        [iv,ix_,_] = s2mpj_ii('FJ0007',ix_)
        pb.xnames=arrset(pb.xnames,iv,'FJ0007')
        [iv,ix_,_] = s2mpj_ii('PI0007',ix_)
        pb.xnames=arrset(pb.xnames,iv,'PI0007')
        [iv,ix_,_] = s2mpj_ii('QI0007',ix_)
        pb.xnames=arrset(pb.xnames,iv,'QI0007')
        [iv,ix_,_] = s2mpj_ii('PJ0007',ix_)
        pb.xnames=arrset(pb.xnames,iv,'PJ0007')
        [iv,ix_,_] = s2mpj_ii('QJ0007',ix_)
        pb.xnames=arrset(pb.xnames,iv,'QJ0007')
        [iv,ix_,_] = s2mpj_ii('PG0001',ix_)
        pb.xnames=arrset(pb.xnames,iv,'PG0001')
        [iv,ix_,_] = s2mpj_ii('PG0002',ix_)
        pb.xnames=arrset(pb.xnames,iv,'PG0002')
        [iv,ix_,_] = s2mpj_ii('PG0003',ix_)
        pb.xnames=arrset(pb.xnames,iv,'PG0003')
        [iv,ix_,_] = s2mpj_ii('QG0001',ix_)
        pb.xnames=arrset(pb.xnames,iv,'QG0001')
        [iv,ix_,_] = s2mpj_ii('QG0002',ix_)
        pb.xnames=arrset(pb.xnames,iv,'QG0002')
        [iv,ix_,_] = s2mpj_ii('QG0003',ix_)
        pb.xnames=arrset(pb.xnames,iv,'QG0003')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('GV20001',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GV20001')
        iv = ix_['V20001']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('SLF0000',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'SLF0000')
        iv = ix_['VF0001']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GV20002',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GV20002')
        iv = ix_['V20002']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GV20003',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GV20003')
        iv = ix_['V20003']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GV20004',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GV20004')
        iv = ix_['V20004']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GV20005',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GV20005')
        iv = ix_['V20005']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GV20006',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GV20006')
        iv = ix_['V20006']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GV20007',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GV20007')
        iv = ix_['V20007']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('LOSS0000',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['PI0001']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        iv = ix_['PJ0001']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GEI0001',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GEI0001')
        iv = ix_['EI0001']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GFI0001',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GFI0001')
        iv = ix_['FI0001']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GEJ0001',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GEJ0001')
        iv = ix_['EJ0001']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GFJ0001',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GFJ0001')
        iv = ix_['FJ0001']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GEI0001',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GEI0001')
        iv = ix_['VE0001']
        pbm.A[ig,iv] = float(-5.299)+pbm.A[ig,iv]
        iv = ix_['VF0001']
        pbm.A[ig,iv] = float(-66.243)+pbm.A[ig,iv]
        iv = ix_['VE0002']
        pbm.A[ig,iv] = float(5.299)+pbm.A[ig,iv]
        iv = ix_['VF0002']
        pbm.A[ig,iv] = float(66.243)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GFI0001',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GFI0001')
        iv = ix_['VE0001']
        pbm.A[ig,iv] = float(66.243)+pbm.A[ig,iv]
        iv = ix_['VF0001']
        pbm.A[ig,iv] = float(-5.299)+pbm.A[ig,iv]
        iv = ix_['VE0002']
        pbm.A[ig,iv] = float(-66.243)+pbm.A[ig,iv]
        iv = ix_['VF0002']
        pbm.A[ig,iv] = float(5.299)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GEJ0001',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GEJ0001')
        iv = ix_['VE0002']
        pbm.A[ig,iv] = float(-5.299)+pbm.A[ig,iv]
        iv = ix_['VF0002']
        pbm.A[ig,iv] = float(-66.243)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GFJ0001',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GFJ0001')
        iv = ix_['VE0002']
        pbm.A[ig,iv] = float(66.243)+pbm.A[ig,iv]
        iv = ix_['VF0002']
        pbm.A[ig,iv] = float(-5.299)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GEJ0001',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GEJ0001')
        iv = ix_['VE0001']
        pbm.A[ig,iv] = float(5.299)+pbm.A[ig,iv]
        iv = ix_['VF0001']
        pbm.A[ig,iv] = float(66.243)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GFJ0001',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GFJ0001')
        iv = ix_['VE0001']
        pbm.A[ig,iv] = float(-66.243)+pbm.A[ig,iv]
        iv = ix_['VF0001']
        pbm.A[ig,iv] = float(5.299)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GPI0001',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GPI0001')
        iv = ix_['PI0001']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GQI0001',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GQI0001')
        iv = ix_['QI0001']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GPJ0001',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GPJ0001')
        iv = ix_['PJ0001']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GQJ0001',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GQJ0001')
        iv = ix_['QJ0001']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GPNI0001',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GPNI0001')
        iv = ix_['PI0001']
        pbm.A[ig,iv] = float(-1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GPNI0002',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GPNI0002')
        iv = ix_['PJ0001']
        pbm.A[ig,iv] = float(-1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GQNI0001',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GQNI0001')
        iv = ix_['QI0001']
        pbm.A[ig,iv] = float(-1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GQNI0002',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GQNI0002')
        iv = ix_['QJ0001']
        pbm.A[ig,iv] = float(-1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GMXI0001',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'GMXI0001')
        [ig,ig_,_] = s2mpj_ii('GMXJ0001',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'GMXJ0001')
        [ig,ig_,_] = s2mpj_ii('LOSS0000',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['PI0002']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        iv = ix_['PJ0002']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GEI0002',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GEI0002')
        iv = ix_['EI0002']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GFI0002',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GFI0002')
        iv = ix_['FI0002']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GEJ0002',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GEJ0002')
        iv = ix_['EJ0002']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GFJ0002',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GFJ0002')
        iv = ix_['FJ0002']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GEI0002',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GEI0002')
        iv = ix_['VE0002']
        pbm.A[ig,iv] = float(-1.175)+pbm.A[ig,iv]
        iv = ix_['VF0002']
        pbm.A[ig,iv] = float(-6.915)+pbm.A[ig,iv]
        iv = ix_['VE0006']
        pbm.A[ig,iv] = float(1.175)+pbm.A[ig,iv]
        iv = ix_['VF0006']
        pbm.A[ig,iv] = float(7.051)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GFI0002',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GFI0002')
        iv = ix_['VE0002']
        pbm.A[ig,iv] = float(6.915)+pbm.A[ig,iv]
        iv = ix_['VF0002']
        pbm.A[ig,iv] = float(-1.175)+pbm.A[ig,iv]
        iv = ix_['VE0006']
        pbm.A[ig,iv] = float(-7.051)+pbm.A[ig,iv]
        iv = ix_['VF0006']
        pbm.A[ig,iv] = float(1.175)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GEJ0002',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GEJ0002')
        iv = ix_['VE0006']
        pbm.A[ig,iv] = float(-1.175)+pbm.A[ig,iv]
        iv = ix_['VF0006']
        pbm.A[ig,iv] = float(-6.915)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GFJ0002',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GFJ0002')
        iv = ix_['VE0006']
        pbm.A[ig,iv] = float(6.915)+pbm.A[ig,iv]
        iv = ix_['VF0006']
        pbm.A[ig,iv] = float(-1.175)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GEJ0002',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GEJ0002')
        iv = ix_['VE0002']
        pbm.A[ig,iv] = float(1.175)+pbm.A[ig,iv]
        iv = ix_['VF0002']
        pbm.A[ig,iv] = float(7.051)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GFJ0002',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GFJ0002')
        iv = ix_['VE0002']
        pbm.A[ig,iv] = float(-7.051)+pbm.A[ig,iv]
        iv = ix_['VF0002']
        pbm.A[ig,iv] = float(1.175)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GPI0002',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GPI0002')
        iv = ix_['PI0002']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GQI0002',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GQI0002')
        iv = ix_['QI0002']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GPJ0002',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GPJ0002')
        iv = ix_['PJ0002']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GQJ0002',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GQJ0002')
        iv = ix_['QJ0002']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GPNI0002',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GPNI0002')
        iv = ix_['PI0002']
        pbm.A[ig,iv] = float(-1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GPNI0006',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GPNI0006')
        iv = ix_['PJ0002']
        pbm.A[ig,iv] = float(-1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GQNI0002',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GQNI0002')
        iv = ix_['QI0002']
        pbm.A[ig,iv] = float(-1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GQNI0006',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GQNI0006')
        iv = ix_['QJ0002']
        pbm.A[ig,iv] = float(-1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GMXI0002',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'GMXI0002')
        [ig,ig_,_] = s2mpj_ii('GMXJ0002',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'GMXJ0002')
        [ig,ig_,_] = s2mpj_ii('LOSS0000',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['PI0003']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        iv = ix_['PJ0003']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GEI0003',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GEI0003')
        iv = ix_['EI0003']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GFI0003',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GFI0003')
        iv = ix_['FI0003']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GEJ0003',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GEJ0003')
        iv = ix_['EJ0003']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GFJ0003',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GFJ0003')
        iv = ix_['FJ0003']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GEI0003',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GEI0003')
        iv = ix_['VE0002']
        pbm.A[ig,iv] = float(-1.726)+pbm.A[ig,iv]
        iv = ix_['VF0002']
        pbm.A[ig,iv] = float(-10.498)+pbm.A[ig,iv]
        iv = ix_['VE0004']
        pbm.A[ig,iv] = float(1.726)+pbm.A[ig,iv]
        iv = ix_['VF0004']
        pbm.A[ig,iv] = float(10.588)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GFI0003',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GFI0003')
        iv = ix_['VE0002']
        pbm.A[ig,iv] = float(10.498)+pbm.A[ig,iv]
        iv = ix_['VF0002']
        pbm.A[ig,iv] = float(-1.726)+pbm.A[ig,iv]
        iv = ix_['VE0004']
        pbm.A[ig,iv] = float(-10.588)+pbm.A[ig,iv]
        iv = ix_['VF0004']
        pbm.A[ig,iv] = float(1.726)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GEJ0003',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GEJ0003')
        iv = ix_['VE0004']
        pbm.A[ig,iv] = float(-1.726)+pbm.A[ig,iv]
        iv = ix_['VF0004']
        pbm.A[ig,iv] = float(-10.498)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GFJ0003',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GFJ0003')
        iv = ix_['VE0004']
        pbm.A[ig,iv] = float(10.498)+pbm.A[ig,iv]
        iv = ix_['VF0004']
        pbm.A[ig,iv] = float(-1.726)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GEJ0003',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GEJ0003')
        iv = ix_['VE0002']
        pbm.A[ig,iv] = float(1.726)+pbm.A[ig,iv]
        iv = ix_['VF0002']
        pbm.A[ig,iv] = float(10.588)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GFJ0003',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GFJ0003')
        iv = ix_['VE0002']
        pbm.A[ig,iv] = float(-10.588)+pbm.A[ig,iv]
        iv = ix_['VF0002']
        pbm.A[ig,iv] = float(1.726)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GPI0003',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GPI0003')
        iv = ix_['PI0003']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GQI0003',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GQI0003')
        iv = ix_['QI0003']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GPJ0003',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GPJ0003')
        iv = ix_['PJ0003']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GQJ0003',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GQJ0003')
        iv = ix_['QJ0003']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GPNI0002',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GPNI0002')
        iv = ix_['PI0003']
        pbm.A[ig,iv] = float(-1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GPNI0004',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GPNI0004')
        iv = ix_['PJ0003']
        pbm.A[ig,iv] = float(-1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GQNI0002',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GQNI0002')
        iv = ix_['QI0003']
        pbm.A[ig,iv] = float(-1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GQNI0004',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GQNI0004')
        iv = ix_['QJ0003']
        pbm.A[ig,iv] = float(-1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GMXI0003',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'GMXI0003')
        [ig,ig_,_] = s2mpj_ii('GMXJ0003',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'GMXJ0003')
        [ig,ig_,_] = s2mpj_ii('LOSS0000',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['PI0004']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        iv = ix_['PJ0004']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GEI0004',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GEI0004')
        iv = ix_['EI0004']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GFI0004',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GFI0004')
        iv = ix_['FI0004']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GEJ0004',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GEJ0004')
        iv = ix_['EJ0004']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GFJ0004',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GFJ0004')
        iv = ix_['FJ0004']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GEI0004',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GEI0004')
        iv = ix_['VE0003']
        pbm.A[ig,iv] = float(-6.897)+pbm.A[ig,iv]
        iv = ix_['VF0003']
        pbm.A[ig,iv] = float(-82.759)+pbm.A[ig,iv]
        iv = ix_['VE0004']
        pbm.A[ig,iv] = float(6.897)+pbm.A[ig,iv]
        iv = ix_['VF0004']
        pbm.A[ig,iv] = float(82.759)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GFI0004',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GFI0004')
        iv = ix_['VE0003']
        pbm.A[ig,iv] = float(82.759)+pbm.A[ig,iv]
        iv = ix_['VF0003']
        pbm.A[ig,iv] = float(-6.897)+pbm.A[ig,iv]
        iv = ix_['VE0004']
        pbm.A[ig,iv] = float(-82.759)+pbm.A[ig,iv]
        iv = ix_['VF0004']
        pbm.A[ig,iv] = float(6.897)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GEJ0004',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GEJ0004')
        iv = ix_['VE0004']
        pbm.A[ig,iv] = float(-6.897)+pbm.A[ig,iv]
        iv = ix_['VF0004']
        pbm.A[ig,iv] = float(-82.759)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GFJ0004',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GFJ0004')
        iv = ix_['VE0004']
        pbm.A[ig,iv] = float(82.759)+pbm.A[ig,iv]
        iv = ix_['VF0004']
        pbm.A[ig,iv] = float(-6.897)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GEJ0004',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GEJ0004')
        iv = ix_['VE0003']
        pbm.A[ig,iv] = float(6.897)+pbm.A[ig,iv]
        iv = ix_['VF0003']
        pbm.A[ig,iv] = float(82.759)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GFJ0004',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GFJ0004')
        iv = ix_['VE0003']
        pbm.A[ig,iv] = float(-82.759)+pbm.A[ig,iv]
        iv = ix_['VF0003']
        pbm.A[ig,iv] = float(6.897)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GPI0004',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GPI0004')
        iv = ix_['PI0004']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GQI0004',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GQI0004')
        iv = ix_['QI0004']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GPJ0004',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GPJ0004')
        iv = ix_['PJ0004']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GQJ0004',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GQJ0004')
        iv = ix_['QJ0004']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GPNI0003',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GPNI0003')
        iv = ix_['PI0004']
        pbm.A[ig,iv] = float(-1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GPNI0004',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GPNI0004')
        iv = ix_['PJ0004']
        pbm.A[ig,iv] = float(-1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GQNI0003',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GQNI0003')
        iv = ix_['QI0004']
        pbm.A[ig,iv] = float(-1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GQNI0004',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GQNI0004')
        iv = ix_['QJ0004']
        pbm.A[ig,iv] = float(-1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GMXI0004',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'GMXI0004')
        [ig,ig_,_] = s2mpj_ii('GMXJ0004',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'GMXJ0004')
        [ig,ig_,_] = s2mpj_ii('LOSS0000',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['PI0005']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        iv = ix_['PJ0005']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GEI0005',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GEI0005')
        iv = ix_['EI0005']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GFI0005',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GFI0005')
        iv = ix_['FI0005']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GEJ0005',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GEJ0005')
        iv = ix_['EJ0005']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GFJ0005',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GFJ0005')
        iv = ix_['FJ0005']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GEI0005',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GEI0005')
        iv = ix_['VE0004']
        pbm.A[ig,iv] = float(-1.175)+pbm.A[ig,iv]
        iv = ix_['VF0004']
        pbm.A[ig,iv] = float(-6.915)+pbm.A[ig,iv]
        iv = ix_['VE0007']
        pbm.A[ig,iv] = float(1.175)+pbm.A[ig,iv]
        iv = ix_['VF0007']
        pbm.A[ig,iv] = float(7.051)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GFI0005',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GFI0005')
        iv = ix_['VE0004']
        pbm.A[ig,iv] = float(6.915)+pbm.A[ig,iv]
        iv = ix_['VF0004']
        pbm.A[ig,iv] = float(-1.175)+pbm.A[ig,iv]
        iv = ix_['VE0007']
        pbm.A[ig,iv] = float(-7.051)+pbm.A[ig,iv]
        iv = ix_['VF0007']
        pbm.A[ig,iv] = float(1.175)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GEJ0005',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GEJ0005')
        iv = ix_['VE0007']
        pbm.A[ig,iv] = float(-1.175)+pbm.A[ig,iv]
        iv = ix_['VF0007']
        pbm.A[ig,iv] = float(-6.915)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GFJ0005',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GFJ0005')
        iv = ix_['VE0007']
        pbm.A[ig,iv] = float(6.915)+pbm.A[ig,iv]
        iv = ix_['VF0007']
        pbm.A[ig,iv] = float(-1.175)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GEJ0005',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GEJ0005')
        iv = ix_['VE0004']
        pbm.A[ig,iv] = float(1.175)+pbm.A[ig,iv]
        iv = ix_['VF0004']
        pbm.A[ig,iv] = float(7.051)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GFJ0005',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GFJ0005')
        iv = ix_['VE0004']
        pbm.A[ig,iv] = float(-7.051)+pbm.A[ig,iv]
        iv = ix_['VF0004']
        pbm.A[ig,iv] = float(1.175)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GPI0005',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GPI0005')
        iv = ix_['PI0005']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GQI0005',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GQI0005')
        iv = ix_['QI0005']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GPJ0005',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GPJ0005')
        iv = ix_['PJ0005']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GQJ0005',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GQJ0005')
        iv = ix_['QJ0005']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GPNI0004',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GPNI0004')
        iv = ix_['PI0005']
        pbm.A[ig,iv] = float(-1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GPNI0007',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GPNI0007')
        iv = ix_['PJ0005']
        pbm.A[ig,iv] = float(-1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GQNI0004',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GQNI0004')
        iv = ix_['QI0005']
        pbm.A[ig,iv] = float(-1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GQNI0007',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GQNI0007')
        iv = ix_['QJ0005']
        pbm.A[ig,iv] = float(-1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GMXI0005',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'GMXI0005')
        [ig,ig_,_] = s2mpj_ii('GMXJ0005',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'GMXJ0005')
        [ig,ig_,_] = s2mpj_ii('LOSS0000',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['PI0006']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        iv = ix_['PJ0006']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GEI0006',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GEI0006')
        iv = ix_['EI0006']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GFI0006',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GFI0006')
        iv = ix_['FI0006']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GEJ0006',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GEJ0006')
        iv = ix_['EJ0006']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GFJ0006',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GFJ0006')
        iv = ix_['FJ0006']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GEI0006',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GEI0006')
        iv = ix_['VE0005']
        pbm.A[ig,iv] = float(-3.448)+pbm.A[ig,iv]
        iv = ix_['VF0005']
        pbm.A[ig,iv] = float(-41.379)+pbm.A[ig,iv]
        iv = ix_['VE0006']
        pbm.A[ig,iv] = float(3.448)+pbm.A[ig,iv]
        iv = ix_['VF0006']
        pbm.A[ig,iv] = float(41.379)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GFI0006',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GFI0006')
        iv = ix_['VE0005']
        pbm.A[ig,iv] = float(41.379)+pbm.A[ig,iv]
        iv = ix_['VF0005']
        pbm.A[ig,iv] = float(-3.448)+pbm.A[ig,iv]
        iv = ix_['VE0006']
        pbm.A[ig,iv] = float(-41.379)+pbm.A[ig,iv]
        iv = ix_['VF0006']
        pbm.A[ig,iv] = float(3.448)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GEJ0006',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GEJ0006')
        iv = ix_['VE0006']
        pbm.A[ig,iv] = float(-3.448)+pbm.A[ig,iv]
        iv = ix_['VF0006']
        pbm.A[ig,iv] = float(-41.379)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GFJ0006',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GFJ0006')
        iv = ix_['VE0006']
        pbm.A[ig,iv] = float(41.379)+pbm.A[ig,iv]
        iv = ix_['VF0006']
        pbm.A[ig,iv] = float(-3.448)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GEJ0006',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GEJ0006')
        iv = ix_['VE0005']
        pbm.A[ig,iv] = float(3.448)+pbm.A[ig,iv]
        iv = ix_['VF0005']
        pbm.A[ig,iv] = float(41.379)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GFJ0006',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GFJ0006')
        iv = ix_['VE0005']
        pbm.A[ig,iv] = float(-41.379)+pbm.A[ig,iv]
        iv = ix_['VF0005']
        pbm.A[ig,iv] = float(3.448)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GPI0006',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GPI0006')
        iv = ix_['PI0006']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GQI0006',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GQI0006')
        iv = ix_['QI0006']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GPJ0006',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GPJ0006')
        iv = ix_['PJ0006']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GQJ0006',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GQJ0006')
        iv = ix_['QJ0006']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GPNI0005',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GPNI0005')
        iv = ix_['PI0006']
        pbm.A[ig,iv] = float(-1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GPNI0006',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GPNI0006')
        iv = ix_['PJ0006']
        pbm.A[ig,iv] = float(-1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GQNI0005',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GQNI0005')
        iv = ix_['QI0006']
        pbm.A[ig,iv] = float(-1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GQNI0006',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GQNI0006')
        iv = ix_['QJ0006']
        pbm.A[ig,iv] = float(-1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GMXI0006',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'GMXI0006')
        [ig,ig_,_] = s2mpj_ii('GMXJ0006',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'GMXJ0006')
        [ig,ig_,_] = s2mpj_ii('LOSS0000',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['PI0007']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        iv = ix_['PJ0007']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GEI0007',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GEI0007')
        iv = ix_['EI0007']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GFI0007',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GFI0007')
        iv = ix_['FI0007']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GEJ0007',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GEJ0007')
        iv = ix_['EJ0007']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GFJ0007',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GFJ0007')
        iv = ix_['FJ0007']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GEI0007',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GEI0007')
        iv = ix_['VE0006']
        pbm.A[ig,iv] = float(-1.726)+pbm.A[ig,iv]
        iv = ix_['VF0006']
        pbm.A[ig,iv] = float(-10.498)+pbm.A[ig,iv]
        iv = ix_['VE0007']
        pbm.A[ig,iv] = float(1.726)+pbm.A[ig,iv]
        iv = ix_['VF0007']
        pbm.A[ig,iv] = float(10.588)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GFI0007',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GFI0007')
        iv = ix_['VE0006']
        pbm.A[ig,iv] = float(10.498)+pbm.A[ig,iv]
        iv = ix_['VF0006']
        pbm.A[ig,iv] = float(-1.726)+pbm.A[ig,iv]
        iv = ix_['VE0007']
        pbm.A[ig,iv] = float(-10.588)+pbm.A[ig,iv]
        iv = ix_['VF0007']
        pbm.A[ig,iv] = float(1.726)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GEJ0007',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GEJ0007')
        iv = ix_['VE0007']
        pbm.A[ig,iv] = float(-1.726)+pbm.A[ig,iv]
        iv = ix_['VF0007']
        pbm.A[ig,iv] = float(-10.498)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GFJ0007',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GFJ0007')
        iv = ix_['VE0007']
        pbm.A[ig,iv] = float(10.498)+pbm.A[ig,iv]
        iv = ix_['VF0007']
        pbm.A[ig,iv] = float(-1.726)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GEJ0007',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GEJ0007')
        iv = ix_['VE0006']
        pbm.A[ig,iv] = float(1.726)+pbm.A[ig,iv]
        iv = ix_['VF0006']
        pbm.A[ig,iv] = float(10.588)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GFJ0007',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GFJ0007')
        iv = ix_['VE0006']
        pbm.A[ig,iv] = float(-10.588)+pbm.A[ig,iv]
        iv = ix_['VF0006']
        pbm.A[ig,iv] = float(1.726)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GPI0007',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GPI0007')
        iv = ix_['PI0007']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GQI0007',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GQI0007')
        iv = ix_['QI0007']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GPJ0007',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GPJ0007')
        iv = ix_['PJ0007']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GQJ0007',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GQJ0007')
        iv = ix_['QJ0007']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GPNI0006',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GPNI0006')
        iv = ix_['PI0007']
        pbm.A[ig,iv] = float(-1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GPNI0007',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GPNI0007')
        iv = ix_['PJ0007']
        pbm.A[ig,iv] = float(-1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GQNI0006',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GQNI0006')
        iv = ix_['QI0007']
        pbm.A[ig,iv] = float(-1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GQNI0007',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GQNI0007')
        iv = ix_['QJ0007']
        pbm.A[ig,iv] = float(-1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GMXI0007',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'GMXI0007')
        [ig,ig_,_] = s2mpj_ii('GMXJ0007',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'GMXJ0007')
        [ig,ig_,_] = s2mpj_ii('GPNI0001',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GPNI0001')
        iv = ix_['PG0001']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GPNI0003',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GPNI0003')
        iv = ix_['PG0002']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GPNI0005',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GPNI0005')
        iv = ix_['PG0003']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GQNI0001',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GQNI0001')
        iv = ix_['QG0001']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GQNI0003',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GQNI0003')
        iv = ix_['QG0002']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GQNI0005',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GQNI0005')
        iv = ix_['QG0003']
        pbm.A[ig,iv] = float(1.000)+pbm.A[ig,iv]
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        legrps = find(gtype,lambda x:x=='<=')
        eqgrps = find(gtype,lambda x:x=='==')
        gegrps = find(gtype,lambda x:x=='>=')
        pb.nle = len(legrps)
        pb.neq = len(eqgrps)
        pb.nge = len(gegrps)
        pb.m   = pb.nle+pb.neq+pb.nge
        pbm.congrps = find(gtype,lambda x:(x=='<=' or x=='==' or x=='>='))
        pb.cnames= cnames[pbm.congrps]
        pb.nob = ngrp-pb.m
        pbm.objgrps = find(gtype,lambda x:x=='<>')
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        pbm.gconst = arrset(pbm.gconst,ig_['GPNI0001'],float(0.000))
        pbm.gconst = arrset(pbm.gconst,ig_['GQNI0001'],float(0.000))
        pbm.gconst = arrset(pbm.gconst,ig_['SLF0000'],float(0.000))
        pbm.gconst = arrset(pbm.gconst,ig_['GPNI0002'],float(2.000))
        pbm.gconst = arrset(pbm.gconst,ig_['GQNI0002'],float(3.000))
        pbm.gconst = arrset(pbm.gconst,ig_['GPNI0003'],float(0.600))
        pbm.gconst = arrset(pbm.gconst,ig_['GQNI0003'],float(0.080))
        pbm.gconst = arrset(pbm.gconst,ig_['GPNI0004'],float(2.000))
        pbm.gconst = arrset(pbm.gconst,ig_['GQNI0004'],float(0.200))
        pbm.gconst = arrset(pbm.gconst,ig_['GPNI0005'],float(0.500))
        pbm.gconst = arrset(pbm.gconst,ig_['GQNI0005'],float(0.050))
        pbm.gconst = arrset(pbm.gconst,ig_['GPNI0006'],float(1.000))
        pbm.gconst = arrset(pbm.gconst,ig_['GQNI0006'],float(0.300))
        pbm.gconst = arrset(pbm.gconst,ig_['GPNI0007'],float(2.000))
        pbm.gconst = arrset(pbm.gconst,ig_['GQNI0007'],float(1.000))
        pbm.gconst = arrset(pbm.gconst,ig_['GMXI0001'],float(16.000))
        pbm.gconst = arrset(pbm.gconst,ig_['GMXJ0001'],float(16.000))
        pbm.gconst = arrset(pbm.gconst,ig_['GMXI0002'],float(4.000))
        pbm.gconst = arrset(pbm.gconst,ig_['GMXJ0002'],float(4.000))
        pbm.gconst = arrset(pbm.gconst,ig_['GMXI0003'],float(4.000))
        pbm.gconst = arrset(pbm.gconst,ig_['GMXJ0003'],float(4.000))
        pbm.gconst = arrset(pbm.gconst,ig_['GMXI0004'],float(25.000))
        pbm.gconst = arrset(pbm.gconst,ig_['GMXJ0004'],float(25.000))
        pbm.gconst = arrset(pbm.gconst,ig_['GMXI0005'],float(4.000))
        pbm.gconst = arrset(pbm.gconst,ig_['GMXJ0005'],float(4.000))
        pbm.gconst = arrset(pbm.gconst,ig_['GMXI0006'],float(6.250))
        pbm.gconst = arrset(pbm.gconst,ig_['GMXJ0006'],float(6.250))
        pbm.gconst = arrset(pbm.gconst,ig_['GMXI0007'],float(4.000))
        pbm.gconst = arrset(pbm.gconst,ig_['GMXJ0007'],float(4.000))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        pb.xlower[ix_['V20001']] = 0.810
        pb.xupper[ix_['V20001']] = 1.210
        pb.xlower[ix_['V20002']] = 0.810
        pb.xupper[ix_['V20002']] = 1.210
        pb.xlower[ix_['V20003']] = 0.941
        pb.xupper[ix_['V20003']] = 1.210
        pb.xlower[ix_['V20004']] = 0.941
        pb.xupper[ix_['V20004']] = 1.210
        pb.xlower[ix_['V20005']] = 0.941
        pb.xupper[ix_['V20005']] = 1.210
        pb.xlower[ix_['V20006']] = 0.941
        pb.xupper[ix_['V20006']] = 1.210
        pb.xlower[ix_['V20007']] = 0.941
        pb.xupper[ix_['V20007']] = 1.210
        pb.xlower[ix_['PG0001']] = 0.500
        pb.xupper[ix_['PG0001']] = 10.000
        pb.xlower[ix_['PG0002']] = 0.500
        pb.xupper[ix_['PG0002']] = 10.000
        pb.xlower[ix_['PG0003']] = 0.200
        pb.xupper[ix_['PG0003']] = 4.000
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        pb.x0[ix_['VE0001']] = float(1.000)
        pb.x0[ix_['VF0001']] = float(0.000)
        pb.x0[ix_['V20001']] = float(1.000)
        pb.x0[ix_['VE0002']] = float(1.001)
        pb.x0[ix_['VF0002']] = float(0.000)
        pb.x0[ix_['V20002']] = float(1.002)
        pb.x0[ix_['VE0003']] = float(1.050)
        pb.x0[ix_['VF0003']] = float(0.000)
        pb.x0[ix_['V20003']] = float(1.102)
        pb.x0[ix_['VE0004']] = float(1.001)
        pb.x0[ix_['VF0004']] = float(0.000)
        pb.x0[ix_['V20004']] = float(1.002)
        pb.x0[ix_['VE0005']] = float(1.050)
        pb.x0[ix_['VF0005']] = float(0.000)
        pb.x0[ix_['V20005']] = float(1.102)
        pb.x0[ix_['VE0006']] = float(1.001)
        pb.x0[ix_['VF0006']] = float(0.000)
        pb.x0[ix_['V20006']] = float(1.002)
        pb.x0[ix_['VE0007']] = float(1.001)
        pb.x0[ix_['VF0007']] = float(0.000)
        pb.x0[ix_['V20007']] = float(1.002)
        pb.x0[ix_['EI0001']] = float(-0.005)
        pb.x0[ix_['FI0001']] = float(0.066)
        pb.x0[ix_['EJ0001']] = float(0.005)
        pb.x0[ix_['FJ0001']] = float(-0.066)
        pb.x0[ix_['PI0001']] = float(-0.005)
        pb.x0[ix_['QI0001']] = float(-0.066)
        pb.x0[ix_['PJ0001']] = float(0.005)
        pb.x0[ix_['QJ0001']] = float(0.066)
        pb.x0[ix_['EI0002']] = float(0.000)
        pb.x0[ix_['FI0002']] = float(0.136)
        pb.x0[ix_['EJ0002']] = float(0.000)
        pb.x0[ix_['FJ0002']] = float(0.136)
        pb.x0[ix_['PI0002']] = float(0.000)
        pb.x0[ix_['QI0002']] = float(-0.136)
        pb.x0[ix_['PJ0002']] = float(0.000)
        pb.x0[ix_['QJ0002']] = float(-0.136)
        pb.x0[ix_['EI0003']] = float(0.000)
        pb.x0[ix_['FI0003']] = float(0.091)
        pb.x0[ix_['EJ0003']] = float(0.000)
        pb.x0[ix_['FJ0003']] = float(0.091)
        pb.x0[ix_['PI0003']] = float(0.000)
        pb.x0[ix_['QI0003']] = float(-0.091)
        pb.x0[ix_['PJ0003']] = float(0.000)
        pb.x0[ix_['QJ0003']] = float(-0.091)
        pb.x0[ix_['EI0004']] = float(0.338)
        pb.x0[ix_['FI0004']] = float(-4.055)
        pb.x0[ix_['EJ0004']] = float(-0.338)
        pb.x0[ix_['FJ0004']] = float(4.055)
        pb.x0[ix_['PI0004']] = float(0.355)
        pb.x0[ix_['QI0004']] = float(4.258)
        pb.x0[ix_['PJ0004']] = float(-0.338)
        pb.x0[ix_['QJ0004']] = float(-4.059)
        pb.x0[ix_['EI0005']] = float(0.000)
        pb.x0[ix_['FI0005']] = float(0.136)
        pb.x0[ix_['EJ0005']] = float(0.000)
        pb.x0[ix_['FJ0005']] = float(0.136)
        pb.x0[ix_['PI0005']] = float(0.000)
        pb.x0[ix_['QI0005']] = float(-0.136)
        pb.x0[ix_['PJ0005']] = float(0.000)
        pb.x0[ix_['QJ0005']] = float(-0.136)
        pb.x0[ix_['EI0006']] = float(0.169)
        pb.x0[ix_['FI0006']] = float(-2.028)
        pb.x0[ix_['EJ0006']] = float(-0.169)
        pb.x0[ix_['FJ0006']] = float(2.028)
        pb.x0[ix_['PI0006']] = float(0.177)
        pb.x0[ix_['QI0006']] = float(2.129)
        pb.x0[ix_['PJ0006']] = float(-0.169)
        pb.x0[ix_['QJ0006']] = float(-2.030)
        pb.x0[ix_['EI0007']] = float(0.000)
        pb.x0[ix_['FI0007']] = float(0.091)
        pb.x0[ix_['EJ0007']] = float(0.000)
        pb.x0[ix_['FJ0007']] = float(0.091)
        pb.x0[ix_['PI0007']] = float(0.000)
        pb.x0[ix_['QI0007']] = float(-0.091)
        pb.x0[ix_['PJ0007']] = float(0.000)
        pb.x0[ix_['QJ0007']] = float(-0.091)
        pb.x0[ix_['PG0001']] = float(3.000)
        pb.x0[ix_['PG0002']] = float(5.000)
        pb.x0[ix_['PG0003']] = float(2.000)
        pb.x0[ix_['QG0001']] = float(0.000)
        pb.x0[ix_['QG0002']] = float(0.000)
        pb.x0[ix_['QG0003']] = float(0.000)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eXTIMESY', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        [it,iet_,_] = s2mpj_ii( 'eXSQUARE', iet_)
        elftv = loaset(elftv,it,0,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        ename = 'E20001'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eXSQUARE"])
        vname = 'VE0001'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F20001'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eXSQUARE"])
        vname = 'VF0001'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E20002'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eXSQUARE"])
        vname = 'VE0002'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F20002'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eXSQUARE"])
        vname = 'VF0002'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E20003'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eXSQUARE"])
        vname = 'VE0003'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F20003'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eXSQUARE"])
        vname = 'VF0003'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E20004'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eXSQUARE"])
        vname = 'VE0004'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F20004'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eXSQUARE"])
        vname = 'VF0004'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E20005'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eXSQUARE"])
        vname = 'VE0005'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F20005'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eXSQUARE"])
        vname = 'VF0005'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E20006'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eXSQUARE"])
        vname = 'VE0006'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F20006'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eXSQUARE"])
        vname = 'VF0006'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E20007'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eXSQUARE"])
        vname = 'VE0007'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F20007'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eXSQUARE"])
        vname = 'VF0007'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EIEI0001'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'EI0001'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VE0001'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'FIFI0001'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'FI0001'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VF0001'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EIFI0001'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'EI0001'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VF0001'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'FIEI0001'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'FI0001'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VE0001'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EJEJ0001'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'EJ0001'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VE0002'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'FJFJ0001'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'FJ0001'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VF0002'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EJFJ0001'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'EJ0001'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VF0002'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'FJEJ0001'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'FJ0001'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VE0002'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'PI20001'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eXSQUARE"])
        vname = 'PI0001'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'QI20001'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eXSQUARE"])
        vname = 'QI0001'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'PJ20001'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eXSQUARE"])
        vname = 'PJ0001'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'QJ20001'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eXSQUARE"])
        vname = 'QJ0001'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EIEI0002'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'EI0002'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VE0002'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'FIFI0002'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'FI0002'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VF0002'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EIFI0002'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'EI0002'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VF0002'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'FIEI0002'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'FI0002'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VE0002'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EJEJ0002'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'EJ0002'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VE0006'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'FJFJ0002'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'FJ0002'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VF0006'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EJFJ0002'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'EJ0002'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VF0006'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'FJEJ0002'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'FJ0002'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VE0006'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'PI20002'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eXSQUARE"])
        vname = 'PI0002'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'QI20002'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eXSQUARE"])
        vname = 'QI0002'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'PJ20002'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eXSQUARE"])
        vname = 'PJ0002'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'QJ20002'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eXSQUARE"])
        vname = 'QJ0002'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EIEI0003'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'EI0003'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VE0002'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'FIFI0003'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'FI0003'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VF0002'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EIFI0003'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'EI0003'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VF0002'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'FIEI0003'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'FI0003'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VE0002'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EJEJ0003'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'EJ0003'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VE0004'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'FJFJ0003'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'FJ0003'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VF0004'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EJFJ0003'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'EJ0003'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VF0004'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'FJEJ0003'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'FJ0003'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VE0004'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'PI20003'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eXSQUARE"])
        vname = 'PI0003'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'QI20003'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eXSQUARE"])
        vname = 'QI0003'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'PJ20003'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eXSQUARE"])
        vname = 'PJ0003'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'QJ20003'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eXSQUARE"])
        vname = 'QJ0003'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EIEI0004'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'EI0004'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VE0003'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'FIFI0004'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'FI0004'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VF0003'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EIFI0004'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'EI0004'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VF0003'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'FIEI0004'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'FI0004'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VE0003'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EJEJ0004'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'EJ0004'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VE0004'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'FJFJ0004'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'FJ0004'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VF0004'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EJFJ0004'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'EJ0004'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VF0004'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'FJEJ0004'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'FJ0004'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VE0004'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'PI20004'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eXSQUARE"])
        vname = 'PI0004'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'QI20004'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eXSQUARE"])
        vname = 'QI0004'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'PJ20004'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eXSQUARE"])
        vname = 'PJ0004'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'QJ20004'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eXSQUARE"])
        vname = 'QJ0004'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EIEI0005'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'EI0005'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VE0004'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'FIFI0005'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'FI0005'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VF0004'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EIFI0005'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'EI0005'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VF0004'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'FIEI0005'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'FI0005'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VE0004'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EJEJ0005'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'EJ0005'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VE0007'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'FJFJ0005'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'FJ0005'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VF0007'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EJFJ0005'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'EJ0005'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VF0007'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'FJEJ0005'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'FJ0005'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VE0007'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'PI20005'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eXSQUARE"])
        vname = 'PI0005'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'QI20005'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eXSQUARE"])
        vname = 'QI0005'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'PJ20005'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eXSQUARE"])
        vname = 'PJ0005'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'QJ20005'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eXSQUARE"])
        vname = 'QJ0005'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EIEI0006'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'EI0006'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VE0005'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'FIFI0006'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'FI0006'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VF0005'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EIFI0006'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'EI0006'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VF0005'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'FIEI0006'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'FI0006'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VE0005'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EJEJ0006'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'EJ0006'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VE0006'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'FJFJ0006'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'FJ0006'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VF0006'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EJFJ0006'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'EJ0006'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VF0006'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'FJEJ0006'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'FJ0006'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VE0006'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'PI20006'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eXSQUARE"])
        vname = 'PI0006'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'QI20006'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eXSQUARE"])
        vname = 'QI0006'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'PJ20006'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eXSQUARE"])
        vname = 'PJ0006'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'QJ20006'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eXSQUARE"])
        vname = 'QJ0006'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EIEI0007'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'EI0007'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VE0006'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'FIFI0007'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'FI0007'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VF0006'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EIFI0007'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'EI0007'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VF0006'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'FIEI0007'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'FI0007'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VE0006'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EJEJ0007'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'EJ0007'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VE0007'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'FJFJ0007'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'FJ0007'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VF0007'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EJFJ0007'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'EJ0007'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VF0007'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'FJEJ0007'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXTIMESY')
        ielftype = arrset(ielftype, ie, iet_["eXTIMESY"])
        vname = 'FJ0007'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'VE0007'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'PI20007'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eXSQUARE"])
        vname = 'PI0007'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'QI20007'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eXSQUARE"])
        vname = 'QI0007'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'PJ20007'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eXSQUARE"])
        vname = 'PJ0007'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'QJ20007'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eXSQUARE"])
        vname = 'QJ0007'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        ig = ig_['GV20001']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E20001'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F20001'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        ig = ig_['GV20002']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E20002'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F20002'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        ig = ig_['GV20003']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E20003'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F20003'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        ig = ig_['GV20004']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E20004'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F20004'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        ig = ig_['GV20005']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E20005'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F20005'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        ig = ig_['GV20006']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E20006'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F20006'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        ig = ig_['GV20007']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E20007'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F20007'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        ig = ig_['GPI0001']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EIEI0001'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['FIFI0001'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        ig = ig_['GQI0001']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EIFI0001'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['FIEI0001'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.000))
        ig = ig_['GPJ0001']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EJEJ0001'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['FJFJ0001'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        ig = ig_['GQJ0001']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EJFJ0001'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['FJEJ0001'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.000))
        ig = ig_['GMXI0001']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['PI20001'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['QI20001'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.000))
        ig = ig_['GMXJ0001']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['PJ20001'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['QJ20001'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.000))
        ig = ig_['GPI0002']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EIEI0002'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['FIFI0002'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        ig = ig_['GQI0002']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EIFI0002'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['FIEI0002'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.000))
        ig = ig_['GPJ0002']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EJEJ0002'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['FJFJ0002'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        ig = ig_['GQJ0002']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EJFJ0002'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['FJEJ0002'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.000))
        ig = ig_['GMXI0002']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['PI20002'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['QI20002'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.000))
        ig = ig_['GMXJ0002']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['PJ20002'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['QJ20002'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.000))
        ig = ig_['GPI0003']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EIEI0003'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['FIFI0003'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        ig = ig_['GQI0003']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EIFI0003'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['FIEI0003'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.000))
        ig = ig_['GPJ0003']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EJEJ0003'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['FJFJ0003'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        ig = ig_['GQJ0003']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EJFJ0003'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['FJEJ0003'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.000))
        ig = ig_['GMXI0003']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['PI20003'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['QI20003'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.000))
        ig = ig_['GMXJ0003']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['PJ20003'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['QJ20003'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.000))
        ig = ig_['GPI0004']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EIEI0004'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['FIFI0004'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        ig = ig_['GQI0004']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EIFI0004'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['FIEI0004'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.000))
        ig = ig_['GPJ0004']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EJEJ0004'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['FJFJ0004'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        ig = ig_['GQJ0004']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EJFJ0004'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['FJEJ0004'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.000))
        ig = ig_['GMXI0004']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['PI20004'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['QI20004'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.000))
        ig = ig_['GMXJ0004']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['PJ20004'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['QJ20004'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.000))
        ig = ig_['GPI0005']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EIEI0005'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['FIFI0005'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        ig = ig_['GQI0005']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EIFI0005'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['FIEI0005'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.000))
        ig = ig_['GPJ0005']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EJEJ0005'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['FJFJ0005'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        ig = ig_['GQJ0005']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EJFJ0005'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['FJEJ0005'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.000))
        ig = ig_['GMXI0005']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['PI20005'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['QI20005'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.000))
        ig = ig_['GMXJ0005']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['PJ20005'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['QJ20005'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.000))
        ig = ig_['GPI0006']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EIEI0006'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['FIFI0006'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        ig = ig_['GQI0006']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EIFI0006'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['FIEI0006'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.000))
        ig = ig_['GPJ0006']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EJEJ0006'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['FJFJ0006'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        ig = ig_['GQJ0006']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EJFJ0006'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['FJEJ0006'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.000))
        ig = ig_['GMXI0006']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['PI20006'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['QI20006'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.000))
        ig = ig_['GMXJ0006']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['PJ20006'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['QJ20006'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.000))
        ig = ig_['GPI0007']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EIEI0007'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['FIFI0007'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        ig = ig_['GQI0007']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EIFI0007'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['FIEI0007'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.000))
        ig = ig_['GPJ0007']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EJEJ0007'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['FJFJ0007'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        ig = ig_['GQJ0007']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EJFJ0007'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['FJEJ0007'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.000))
        ig = ig_['GMXI0007']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['PI20007'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['QI20007'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.000))
        ig = ig_['GMXJ0007']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['PJ20007'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.000))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['QJ20007'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.000))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solutions
# LO SOLTN               6.78605619D-2
# LO SOLTN               6.78422730D-2
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.cupper[np.arange(pb.nle)] = np.zeros((pb.nle,1))
        pb.clower[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.cupper[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "LQR2-RN-83-92"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eXTIMESY(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[1]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]
            g_[1] = EV_[0]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = 1.0
                H_[1,0] = H_[0,1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eXSQUARE(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[0]+EV_[0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

