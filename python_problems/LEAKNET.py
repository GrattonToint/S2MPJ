from s2mpjlib import *
class  LEAKNET(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : LEAKNET
#    *********
# 
#    The British Gas leaknet problem.
# 
#    The problem is to minimize the gas leakage in a natural gas network
#    by adjusting the gauge pressures (the P variables), the pipe flows
#    (the Q variables) and the source flows (the S variables).  There are a
#    set of nonlinear constraints corresponding to each pipe (the PIP
#    constraints); These relate the pressures at the start and end of the
#    pipe to the leakage from the pipe. There are also conservation
#    equations (the linear N constraints) at each node (flow in = flow
#    out). Finally, the pressures and source flows are restricted.
# 
#    Source:
#    British Gas, private communication.
# 
#    SIF input: Nick Gould, 25th June 1990.
# 
#    classification = "C-CLOR2-RN-156-153"
# 
#    network data
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'LEAKNET'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['NODES'] = 73
        v_['PIPES'] = 80
        v_['1'] = 1
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('N      1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N      1')
        [ig,ig_,_] = s2mpj_ii('N      2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N      2')
        [ig,ig_,_] = s2mpj_ii('N      3',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N      3')
        [ig,ig_,_] = s2mpj_ii('N      4',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N      4')
        [ig,ig_,_] = s2mpj_ii('N      5',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N      5')
        [ig,ig_,_] = s2mpj_ii('N      6',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N      6')
        [ig,ig_,_] = s2mpj_ii('N      9',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N      9')
        [ig,ig_,_] = s2mpj_ii('N     10',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N     10')
        [ig,ig_,_] = s2mpj_ii('N     12',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N     12')
        [ig,ig_,_] = s2mpj_ii('N     13',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N     13')
        [ig,ig_,_] = s2mpj_ii('N     14',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N     14')
        [ig,ig_,_] = s2mpj_ii('N     15',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N     15')
        [ig,ig_,_] = s2mpj_ii('N     16',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N     16')
        [ig,ig_,_] = s2mpj_ii('N     17',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N     17')
        [ig,ig_,_] = s2mpj_ii('N     18',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N     18')
        [ig,ig_,_] = s2mpj_ii('N     19',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N     19')
        [ig,ig_,_] = s2mpj_ii('N     20',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N     20')
        [ig,ig_,_] = s2mpj_ii('N     21',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N     21')
        [ig,ig_,_] = s2mpj_ii('N     22',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N     22')
        [ig,ig_,_] = s2mpj_ii('N     23',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N     23')
        [ig,ig_,_] = s2mpj_ii('N     26',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N     26')
        [ig,ig_,_] = s2mpj_ii('N     27',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N     27')
        [ig,ig_,_] = s2mpj_ii('N    101',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    101')
        [ig,ig_,_] = s2mpj_ii('N    102',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    102')
        [ig,ig_,_] = s2mpj_ii('N    103',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    103')
        [ig,ig_,_] = s2mpj_ii('N    104',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    104')
        [ig,ig_,_] = s2mpj_ii('N    105',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    105')
        [ig,ig_,_] = s2mpj_ii('N    106',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    106')
        [ig,ig_,_] = s2mpj_ii('N    107',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    107')
        [ig,ig_,_] = s2mpj_ii('N    108',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    108')
        [ig,ig_,_] = s2mpj_ii('N    109',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    109')
        [ig,ig_,_] = s2mpj_ii('N    110',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    110')
        [ig,ig_,_] = s2mpj_ii('N    111',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    111')
        [ig,ig_,_] = s2mpj_ii('N    112',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    112')
        [ig,ig_,_] = s2mpj_ii('N    201',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    201')
        [ig,ig_,_] = s2mpj_ii('N    202',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    202')
        [ig,ig_,_] = s2mpj_ii('N    203',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    203')
        [ig,ig_,_] = s2mpj_ii('N    204',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    204')
        [ig,ig_,_] = s2mpj_ii('N    205',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    205')
        [ig,ig_,_] = s2mpj_ii('N    206',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    206')
        [ig,ig_,_] = s2mpj_ii('N    207',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    207')
        [ig,ig_,_] = s2mpj_ii('N    208',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    208')
        [ig,ig_,_] = s2mpj_ii('N    209',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    209')
        [ig,ig_,_] = s2mpj_ii('N    210',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    210')
        [ig,ig_,_] = s2mpj_ii('N    211',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    211')
        [ig,ig_,_] = s2mpj_ii('N    212',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    212')
        [ig,ig_,_] = s2mpj_ii('N    301',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    301')
        [ig,ig_,_] = s2mpj_ii('N    302',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    302')
        [ig,ig_,_] = s2mpj_ii('N    303',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    303')
        [ig,ig_,_] = s2mpj_ii('N    304',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    304')
        [ig,ig_,_] = s2mpj_ii('N    305',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    305')
        [ig,ig_,_] = s2mpj_ii('N    306',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    306')
        [ig,ig_,_] = s2mpj_ii('N    307',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    307')
        [ig,ig_,_] = s2mpj_ii('N    308',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    308')
        [ig,ig_,_] = s2mpj_ii('N    309',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    309')
        [ig,ig_,_] = s2mpj_ii('N    401',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    401')
        [ig,ig_,_] = s2mpj_ii('N    402',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    402')
        [ig,ig_,_] = s2mpj_ii('N    403',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    403')
        [ig,ig_,_] = s2mpj_ii('N    404',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    404')
        [ig,ig_,_] = s2mpj_ii('N    405',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    405')
        [ig,ig_,_] = s2mpj_ii('N    406',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    406')
        [ig,ig_,_] = s2mpj_ii('N    407',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    407')
        [ig,ig_,_] = s2mpj_ii('N    501',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    501')
        [ig,ig_,_] = s2mpj_ii('N    502',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    502')
        [ig,ig_,_] = s2mpj_ii('N    503',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    503')
        [ig,ig_,_] = s2mpj_ii('N    504',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    504')
        [ig,ig_,_] = s2mpj_ii('N    505',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    505')
        [ig,ig_,_] = s2mpj_ii('N    506',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    506')
        [ig,ig_,_] = s2mpj_ii('N    507',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    507')
        [ig,ig_,_] = s2mpj_ii('N    508',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    508')
        [ig,ig_,_] = s2mpj_ii('N    509',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    509')
        [ig,ig_,_] = s2mpj_ii('N    510',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    510')
        [ig,ig_,_] = s2mpj_ii('N    511',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N    511')
        [ig,ig_,_] = s2mpj_ii('PIP    1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP    1')
        [ig,ig_,_] = s2mpj_ii('PIP    2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP    2')
        [ig,ig_,_] = s2mpj_ii('PIP    3',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP    3')
        [ig,ig_,_] = s2mpj_ii('PIP    4',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP    4')
        [ig,ig_,_] = s2mpj_ii('PIP    5',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP    5')
        [ig,ig_,_] = s2mpj_ii('PIP    6',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP    6')
        [ig,ig_,_] = s2mpj_ii('PIP    7',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP    7')
        [ig,ig_,_] = s2mpj_ii('PIP    8',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP    8')
        [ig,ig_,_] = s2mpj_ii('PIP    9',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP    9')
        [ig,ig_,_] = s2mpj_ii('PIP   10',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   10')
        [ig,ig_,_] = s2mpj_ii('PIP   11',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   11')
        [ig,ig_,_] = s2mpj_ii('PIP   12',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   12')
        [ig,ig_,_] = s2mpj_ii('PIP   13',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   13')
        [ig,ig_,_] = s2mpj_ii('PIP   14',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   14')
        [ig,ig_,_] = s2mpj_ii('PIP   15',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   15')
        [ig,ig_,_] = s2mpj_ii('PIP   16',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   16')
        [ig,ig_,_] = s2mpj_ii('PIP   17',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   17')
        [ig,ig_,_] = s2mpj_ii('PIP   18',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   18')
        [ig,ig_,_] = s2mpj_ii('PIP   19',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   19')
        [ig,ig_,_] = s2mpj_ii('PIP   20',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   20')
        [ig,ig_,_] = s2mpj_ii('PIP   21',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   21')
        [ig,ig_,_] = s2mpj_ii('PIP   22',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   22')
        [ig,ig_,_] = s2mpj_ii('PIP   23',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   23')
        [ig,ig_,_] = s2mpj_ii('PIP   24',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   24')
        [ig,ig_,_] = s2mpj_ii('PIP   25',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   25')
        [ig,ig_,_] = s2mpj_ii('PIP   26',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   26')
        [ig,ig_,_] = s2mpj_ii('PIP   27',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   27')
        [ig,ig_,_] = s2mpj_ii('PIP   28',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   28')
        [ig,ig_,_] = s2mpj_ii('PIP   29',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   29')
        [ig,ig_,_] = s2mpj_ii('PIP   30',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   30')
        [ig,ig_,_] = s2mpj_ii('PIP   31',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   31')
        [ig,ig_,_] = s2mpj_ii('PIP   32',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   32')
        [ig,ig_,_] = s2mpj_ii('PIP   33',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   33')
        [ig,ig_,_] = s2mpj_ii('PIP   34',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   34')
        [ig,ig_,_] = s2mpj_ii('PIP   35',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   35')
        [ig,ig_,_] = s2mpj_ii('PIP   36',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   36')
        [ig,ig_,_] = s2mpj_ii('PIP   37',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   37')
        [ig,ig_,_] = s2mpj_ii('PIP   38',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   38')
        [ig,ig_,_] = s2mpj_ii('PIP   39',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   39')
        [ig,ig_,_] = s2mpj_ii('PIP   40',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   40')
        [ig,ig_,_] = s2mpj_ii('PIP   41',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   41')
        [ig,ig_,_] = s2mpj_ii('PIP   42',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   42')
        [ig,ig_,_] = s2mpj_ii('PIP   43',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   43')
        [ig,ig_,_] = s2mpj_ii('PIP   44',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   44')
        [ig,ig_,_] = s2mpj_ii('PIP   45',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   45')
        [ig,ig_,_] = s2mpj_ii('PIP   46',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   46')
        [ig,ig_,_] = s2mpj_ii('PIP   47',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   47')
        [ig,ig_,_] = s2mpj_ii('PIP   48',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   48')
        [ig,ig_,_] = s2mpj_ii('PIP   49',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   49')
        [ig,ig_,_] = s2mpj_ii('PIP   50',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   50')
        [ig,ig_,_] = s2mpj_ii('PIP   51',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   51')
        [ig,ig_,_] = s2mpj_ii('PIP   52',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   52')
        [ig,ig_,_] = s2mpj_ii('PIP   53',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   53')
        [ig,ig_,_] = s2mpj_ii('PIP   54',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   54')
        [ig,ig_,_] = s2mpj_ii('PIP   55',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   55')
        [ig,ig_,_] = s2mpj_ii('PIP   56',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   56')
        [ig,ig_,_] = s2mpj_ii('PIP   57',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   57')
        [ig,ig_,_] = s2mpj_ii('PIP   58',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   58')
        [ig,ig_,_] = s2mpj_ii('PIP   59',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   59')
        [ig,ig_,_] = s2mpj_ii('PIP   60',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   60')
        [ig,ig_,_] = s2mpj_ii('PIP   61',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   61')
        [ig,ig_,_] = s2mpj_ii('PIP   62',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   62')
        [ig,ig_,_] = s2mpj_ii('PIP   63',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   63')
        [ig,ig_,_] = s2mpj_ii('PIP   64',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   64')
        [ig,ig_,_] = s2mpj_ii('PIP   65',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   65')
        [ig,ig_,_] = s2mpj_ii('PIP   66',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   66')
        [ig,ig_,_] = s2mpj_ii('PIP   67',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   67')
        [ig,ig_,_] = s2mpj_ii('PIP   68',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   68')
        [ig,ig_,_] = s2mpj_ii('PIP   69',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   69')
        [ig,ig_,_] = s2mpj_ii('PIP   70',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   70')
        [ig,ig_,_] = s2mpj_ii('PIP   71',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   71')
        [ig,ig_,_] = s2mpj_ii('PIP   72',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   72')
        [ig,ig_,_] = s2mpj_ii('PIP   73',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   73')
        [ig,ig_,_] = s2mpj_ii('PIP   74',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   74')
        [ig,ig_,_] = s2mpj_ii('PIP   75',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   75')
        [ig,ig_,_] = s2mpj_ii('PIP   76',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   76')
        [ig,ig_,_] = s2mpj_ii('PIP   77',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   77')
        [ig,ig_,_] = s2mpj_ii('PIP   78',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   78')
        [ig,ig_,_] = s2mpj_ii('PIP   79',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   79')
        [ig,ig_,_] = s2mpj_ii('PIP   80',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PIP   80')
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        ngrp   = len(ig_)
        [iv,ix_,_] = s2mpj_ii('P1',ix_)
        self.xnames=arrset(self.xnames,iv,'P1')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.12000e-03))
        [iv,ix_,_] = s2mpj_ii('P2',ix_)
        self.xnames=arrset(self.xnames,iv,'P2')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.12000e-03))
        [iv,ix_,_] = s2mpj_ii('P2',ix_)
        self.xnames=arrset(self.xnames,iv,'P2')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.12000e-03))
        [iv,ix_,_] = s2mpj_ii('P3',ix_)
        self.xnames=arrset(self.xnames,iv,'P3')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.12000e-03))
        [iv,ix_,_] = s2mpj_ii('P3',ix_)
        self.xnames=arrset(self.xnames,iv,'P3')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(8.00000e-04))
        [iv,ix_,_] = s2mpj_ii('P4',ix_)
        self.xnames=arrset(self.xnames,iv,'P4')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(8.00000e-04))
        [iv,ix_,_] = s2mpj_ii('P4',ix_)
        self.xnames=arrset(self.xnames,iv,'P4')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(4.56000e-03))
        [iv,ix_,_] = s2mpj_ii('P5',ix_)
        self.xnames=arrset(self.xnames,iv,'P5')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(4.56000e-03))
        [iv,ix_,_] = s2mpj_ii('P5',ix_)
        self.xnames=arrset(self.xnames,iv,'P5')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(5.60000e-04))
        [iv,ix_,_] = s2mpj_ii('P6',ix_)
        self.xnames=arrset(self.xnames,iv,'P6')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(5.60000e-04))
        [iv,ix_,_] = s2mpj_ii('P5',ix_)
        self.xnames=arrset(self.xnames,iv,'P5')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(3.68000e-03))
        [iv,ix_,_] = s2mpj_ii('P26',ix_)
        self.xnames=arrset(self.xnames,iv,'P26')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(3.68000e-03))
        [iv,ix_,_] = s2mpj_ii('P6',ix_)
        self.xnames=arrset(self.xnames,iv,'P6')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(3.60000e-03))
        [iv,ix_,_] = s2mpj_ii('P9',ix_)
        self.xnames=arrset(self.xnames,iv,'P9')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(3.60000e-03))
        [iv,ix_,_] = s2mpj_ii('P6',ix_)
        self.xnames=arrset(self.xnames,iv,'P6')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(7.84000e-03))
        [iv,ix_,_] = s2mpj_ii('P304',ix_)
        self.xnames=arrset(self.xnames,iv,'P304')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(7.84000e-03))
        [iv,ix_,_] = s2mpj_ii('P9',ix_)
        self.xnames=arrset(self.xnames,iv,'P9')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.60000e-03))
        [iv,ix_,_] = s2mpj_ii('P10',ix_)
        self.xnames=arrset(self.xnames,iv,'P10')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.60000e-03))
        [iv,ix_,_] = s2mpj_ii('P10',ix_)
        self.xnames=arrset(self.xnames,iv,'P10')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(6.08000e-03))
        [iv,ix_,_] = s2mpj_ii('P12',ix_)
        self.xnames=arrset(self.xnames,iv,'P12')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(6.08000e-03))
        [iv,ix_,_] = s2mpj_ii('P10',ix_)
        self.xnames=arrset(self.xnames,iv,'P10')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(4.00000e-03))
        [iv,ix_,_] = s2mpj_ii('P27',ix_)
        self.xnames=arrset(self.xnames,iv,'P27')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(4.00000e-03))
        [iv,ix_,_] = s2mpj_ii('P12',ix_)
        self.xnames=arrset(self.xnames,iv,'P12')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.40000e-04))
        [iv,ix_,_] = s2mpj_ii('P13',ix_)
        self.xnames=arrset(self.xnames,iv,'P13')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.40000e-04))
        [iv,ix_,_] = s2mpj_ii('P13',ix_)
        self.xnames=arrset(self.xnames,iv,'P13')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.86000e-03))
        [iv,ix_,_] = s2mpj_ii('P14',ix_)
        self.xnames=arrset(self.xnames,iv,'P14')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.86000e-03))
        [iv,ix_,_] = s2mpj_ii('P13',ix_)
        self.xnames=arrset(self.xnames,iv,'P13')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.32000e-03))
        [iv,ix_,_] = s2mpj_ii('P19',ix_)
        self.xnames=arrset(self.xnames,iv,'P19')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.32000e-03))
        [iv,ix_,_] = s2mpj_ii('P14',ix_)
        self.xnames=arrset(self.xnames,iv,'P14')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.34000e-03))
        [iv,ix_,_] = s2mpj_ii('P15',ix_)
        self.xnames=arrset(self.xnames,iv,'P15')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.34000e-03))
        [iv,ix_,_] = s2mpj_ii('P16',ix_)
        self.xnames=arrset(self.xnames,iv,'P16')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(9.40000e-04))
        [iv,ix_,_] = s2mpj_ii('P17',ix_)
        self.xnames=arrset(self.xnames,iv,'P17')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(9.40000e-04))
        [iv,ix_,_] = s2mpj_ii('P16',ix_)
        self.xnames=arrset(self.xnames,iv,'P16')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.60000e-04))
        [iv,ix_,_] = s2mpj_ii('P18',ix_)
        self.xnames=arrset(self.xnames,iv,'P18')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.60000e-04))
        [iv,ix_,_] = s2mpj_ii('P16',ix_)
        self.xnames=arrset(self.xnames,iv,'P16')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(4.80000e-03))
        [iv,ix_,_] = s2mpj_ii('P26',ix_)
        self.xnames=arrset(self.xnames,iv,'P26')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(4.80000e-03))
        [iv,ix_,_] = s2mpj_ii('P18',ix_)
        self.xnames=arrset(self.xnames,iv,'P18')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.12000e-03))
        [iv,ix_,_] = s2mpj_ii('P19',ix_)
        self.xnames=arrset(self.xnames,iv,'P19')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.12000e-03))
        [iv,ix_,_] = s2mpj_ii('P19',ix_)
        self.xnames=arrset(self.xnames,iv,'P19')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(5.60000e-04))
        [iv,ix_,_] = s2mpj_ii('P20',ix_)
        self.xnames=arrset(self.xnames,iv,'P20')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(5.60000e-04))
        [iv,ix_,_] = s2mpj_ii('P20',ix_)
        self.xnames=arrset(self.xnames,iv,'P20')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(2.88000e-03))
        [iv,ix_,_] = s2mpj_ii('P21',ix_)
        self.xnames=arrset(self.xnames,iv,'P21')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(2.88000e-03))
        [iv,ix_,_] = s2mpj_ii('P22',ix_)
        self.xnames=arrset(self.xnames,iv,'P22')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.60000e-04))
        [iv,ix_,_] = s2mpj_ii('P404',ix_)
        self.xnames=arrset(self.xnames,iv,'P404')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.60000e-04))
        [iv,ix_,_] = s2mpj_ii('P23',ix_)
        self.xnames=arrset(self.xnames,iv,'P23')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.20000e-03))
        [iv,ix_,_] = s2mpj_ii('P404',ix_)
        self.xnames=arrset(self.xnames,iv,'P404')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.20000e-03))
        [iv,ix_,_] = s2mpj_ii('P27',ix_)
        self.xnames=arrset(self.xnames,iv,'P27')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(5.28000e-03))
        [iv,ix_,_] = s2mpj_ii('P404',ix_)
        self.xnames=arrset(self.xnames,iv,'P404')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(5.28000e-03))
        [iv,ix_,_] = s2mpj_ii('P101',ix_)
        self.xnames=arrset(self.xnames,iv,'P101')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.76000e-03))
        [iv,ix_,_] = s2mpj_ii('P102',ix_)
        self.xnames=arrset(self.xnames,iv,'P102')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.76000e-03))
        [iv,ix_,_] = s2mpj_ii('P102',ix_)
        self.xnames=arrset(self.xnames,iv,'P102')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.28000e-03))
        [iv,ix_,_] = s2mpj_ii('P103',ix_)
        self.xnames=arrset(self.xnames,iv,'P103')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.28000e-03))
        [iv,ix_,_] = s2mpj_ii('P103',ix_)
        self.xnames=arrset(self.xnames,iv,'P103')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.40000e-03))
        [iv,ix_,_] = s2mpj_ii('P104',ix_)
        self.xnames=arrset(self.xnames,iv,'P104')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.40000e-03))
        [iv,ix_,_] = s2mpj_ii('P103',ix_)
        self.xnames=arrset(self.xnames,iv,'P103')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.05000e-03))
        [iv,ix_,_] = s2mpj_ii('P111',ix_)
        self.xnames=arrset(self.xnames,iv,'P111')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.05000e-03))
        [iv,ix_,_] = s2mpj_ii('P104',ix_)
        self.xnames=arrset(self.xnames,iv,'P104')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(2.48000e-03))
        [iv,ix_,_] = s2mpj_ii('P105',ix_)
        self.xnames=arrset(self.xnames,iv,'P105')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(2.48000e-03))
        [iv,ix_,_] = s2mpj_ii('P104',ix_)
        self.xnames=arrset(self.xnames,iv,'P104')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(2.45000e-03))
        [iv,ix_,_] = s2mpj_ii('P110',ix_)
        self.xnames=arrset(self.xnames,iv,'P110')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(2.45000e-03))
        [iv,ix_,_] = s2mpj_ii('P105',ix_)
        self.xnames=arrset(self.xnames,iv,'P105')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.40000e-03))
        [iv,ix_,_] = s2mpj_ii('P106',ix_)
        self.xnames=arrset(self.xnames,iv,'P106')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.40000e-03))
        [iv,ix_,_] = s2mpj_ii('P105',ix_)
        self.xnames=arrset(self.xnames,iv,'P105')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.10000e-03))
        [iv,ix_,_] = s2mpj_ii('P112',ix_)
        self.xnames=arrset(self.xnames,iv,'P112')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.10000e-03))
        [iv,ix_,_] = s2mpj_ii('P106',ix_)
        self.xnames=arrset(self.xnames,iv,'P106')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(9.60000e-04))
        [iv,ix_,_] = s2mpj_ii('P107',ix_)
        self.xnames=arrset(self.xnames,iv,'P107')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(9.60000e-04))
        [iv,ix_,_] = s2mpj_ii('P106',ix_)
        self.xnames=arrset(self.xnames,iv,'P106')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(2.02000e-03))
        [iv,ix_,_] = s2mpj_ii('P109',ix_)
        self.xnames=arrset(self.xnames,iv,'P109')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(2.02000e-03))
        [iv,ix_,_] = s2mpj_ii('P107',ix_)
        self.xnames=arrset(self.xnames,iv,'P107')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.48000e-03))
        [iv,ix_,_] = s2mpj_ii('P201',ix_)
        self.xnames=arrset(self.xnames,iv,'P201')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.48000e-03))
        [iv,ix_,_] = s2mpj_ii('P108',ix_)
        self.xnames=arrset(self.xnames,iv,'P108')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(4.90000e-04))
        [iv,ix_,_] = s2mpj_ii('P109',ix_)
        self.xnames=arrset(self.xnames,iv,'P109')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(4.90000e-04))
        [iv,ix_,_] = s2mpj_ii('P108',ix_)
        self.xnames=arrset(self.xnames,iv,'P108')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.20000e-03))
        [iv,ix_,_] = s2mpj_ii('P210',ix_)
        self.xnames=arrset(self.xnames,iv,'P210')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.20000e-03))
        [iv,ix_,_] = s2mpj_ii('P112',ix_)
        self.xnames=arrset(self.xnames,iv,'P112')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(4.90000e-04))
        [iv,ix_,_] = s2mpj_ii('P509',ix_)
        self.xnames=arrset(self.xnames,iv,'P509')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(4.90000e-04))
        [iv,ix_,_] = s2mpj_ii('P201',ix_)
        self.xnames=arrset(self.xnames,iv,'P201')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.40000e-03))
        [iv,ix_,_] = s2mpj_ii('P202',ix_)
        self.xnames=arrset(self.xnames,iv,'P202')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.40000e-03))
        [iv,ix_,_] = s2mpj_ii('P201',ix_)
        self.xnames=arrset(self.xnames,iv,'P201')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.96000e-03))
        [iv,ix_,_] = s2mpj_ii('P510',ix_)
        self.xnames=arrset(self.xnames,iv,'P510')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.96000e-03))
        [iv,ix_,_] = s2mpj_ii('P202',ix_)
        self.xnames=arrset(self.xnames,iv,'P202')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(3.16000e-03))
        [iv,ix_,_] = s2mpj_ii('P203',ix_)
        self.xnames=arrset(self.xnames,iv,'P203')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(3.16000e-03))
        [iv,ix_,_] = s2mpj_ii('P202',ix_)
        self.xnames=arrset(self.xnames,iv,'P202')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(4.48000e-03))
        [iv,ix_,_] = s2mpj_ii('P211',ix_)
        self.xnames=arrset(self.xnames,iv,'P211')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(4.48000e-03))
        [iv,ix_,_] = s2mpj_ii('P203',ix_)
        self.xnames=arrset(self.xnames,iv,'P203')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.72000e-03))
        [iv,ix_,_] = s2mpj_ii('P204',ix_)
        self.xnames=arrset(self.xnames,iv,'P204')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.72000e-03))
        [iv,ix_,_] = s2mpj_ii('P203',ix_)
        self.xnames=arrset(self.xnames,iv,'P203')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(2.96000e-03))
        [iv,ix_,_] = s2mpj_ii('P502',ix_)
        self.xnames=arrset(self.xnames,iv,'P502')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(2.96000e-03))
        [iv,ix_,_] = s2mpj_ii('P204',ix_)
        self.xnames=arrset(self.xnames,iv,'P204')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(6.00000e-04))
        [iv,ix_,_] = s2mpj_ii('P205',ix_)
        self.xnames=arrset(self.xnames,iv,'P205')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(6.00000e-04))
        [iv,ix_,_] = s2mpj_ii('P204',ix_)
        self.xnames=arrset(self.xnames,iv,'P204')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(8.50000e-04))
        [iv,ix_,_] = s2mpj_ii('P208',ix_)
        self.xnames=arrset(self.xnames,iv,'P208')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(8.50000e-04))
        [iv,ix_,_] = s2mpj_ii('P205',ix_)
        self.xnames=arrset(self.xnames,iv,'P205')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(5.20000e-04))
        [iv,ix_,_] = s2mpj_ii('P206',ix_)
        self.xnames=arrset(self.xnames,iv,'P206')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(5.20000e-04))
        [iv,ix_,_] = s2mpj_ii('P205',ix_)
        self.xnames=arrset(self.xnames,iv,'P205')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(6.70000e-04))
        [iv,ix_,_] = s2mpj_ii('P207',ix_)
        self.xnames=arrset(self.xnames,iv,'P207')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(6.70000e-04))
        [iv,ix_,_] = s2mpj_ii('P206',ix_)
        self.xnames=arrset(self.xnames,iv,'P206')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(2.88000e-03))
        [iv,ix_,_] = s2mpj_ii('P301',ix_)
        self.xnames=arrset(self.xnames,iv,'P301')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(2.88000e-03))
        [iv,ix_,_] = s2mpj_ii('P208',ix_)
        self.xnames=arrset(self.xnames,iv,'P208')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(7.00000e-04))
        [iv,ix_,_] = s2mpj_ii('P209',ix_)
        self.xnames=arrset(self.xnames,iv,'P209')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(7.00000e-04))
        [iv,ix_,_] = s2mpj_ii('P208',ix_)
        self.xnames=arrset(self.xnames,iv,'P208')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(5.48000e-03))
        [iv,ix_,_] = s2mpj_ii('P210',ix_)
        self.xnames=arrset(self.xnames,iv,'P210')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(5.48000e-03))
        [iv,ix_,_] = s2mpj_ii('P210',ix_)
        self.xnames=arrset(self.xnames,iv,'P210')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.88000e-03))
        [iv,ix_,_] = s2mpj_ii('P211',ix_)
        self.xnames=arrset(self.xnames,iv,'P211')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.88000e-03))
        [iv,ix_,_] = s2mpj_ii('P211',ix_)
        self.xnames=arrset(self.xnames,iv,'P211')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(6.40000e-04))
        [iv,ix_,_] = s2mpj_ii('P212',ix_)
        self.xnames=arrset(self.xnames,iv,'P212')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(6.40000e-04))
        [iv,ix_,_] = s2mpj_ii('P301',ix_)
        self.xnames=arrset(self.xnames,iv,'P301')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(4.68000e-03))
        [iv,ix_,_] = s2mpj_ii('P302',ix_)
        self.xnames=arrset(self.xnames,iv,'P302')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(4.68000e-03))
        [iv,ix_,_] = s2mpj_ii('P301',ix_)
        self.xnames=arrset(self.xnames,iv,'P301')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(3.80000e-03))
        [iv,ix_,_] = s2mpj_ii('P304',ix_)
        self.xnames=arrset(self.xnames,iv,'P304')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(3.80000e-03))
        [iv,ix_,_] = s2mpj_ii('P302',ix_)
        self.xnames=arrset(self.xnames,iv,'P302')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(2.40000e-03))
        [iv,ix_,_] = s2mpj_ii('P303',ix_)
        self.xnames=arrset(self.xnames,iv,'P303')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(2.40000e-03))
        [iv,ix_,_] = s2mpj_ii('P302',ix_)
        self.xnames=arrset(self.xnames,iv,'P302')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(5.32000e-03))
        [iv,ix_,_] = s2mpj_ii('P305',ix_)
        self.xnames=arrset(self.xnames,iv,'P305')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(5.32000e-03))
        [iv,ix_,_] = s2mpj_ii('P303',ix_)
        self.xnames=arrset(self.xnames,iv,'P303')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.08000e-03))
        [iv,ix_,_] = s2mpj_ii('P401',ix_)
        self.xnames=arrset(self.xnames,iv,'P401')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.08000e-03))
        [iv,ix_,_] = s2mpj_ii('P305',ix_)
        self.xnames=arrset(self.xnames,iv,'P305')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.68000e-03))
        [iv,ix_,_] = s2mpj_ii('P306',ix_)
        self.xnames=arrset(self.xnames,iv,'P306')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.68000e-03))
        [iv,ix_,_] = s2mpj_ii('P305',ix_)
        self.xnames=arrset(self.xnames,iv,'P305')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(7.90000e-04))
        [iv,ix_,_] = s2mpj_ii('P309',ix_)
        self.xnames=arrset(self.xnames,iv,'P309')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(7.90000e-04))
        [iv,ix_,_] = s2mpj_ii('P306',ix_)
        self.xnames=arrset(self.xnames,iv,'P306')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(4.70000e-04))
        [iv,ix_,_] = s2mpj_ii('P307',ix_)
        self.xnames=arrset(self.xnames,iv,'P307')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(4.70000e-04))
        [iv,ix_,_] = s2mpj_ii('P306',ix_)
        self.xnames=arrset(self.xnames,iv,'P306')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(5.20000e-04))
        [iv,ix_,_] = s2mpj_ii('P308',ix_)
        self.xnames=arrset(self.xnames,iv,'P308')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(5.20000e-04))
        [iv,ix_,_] = s2mpj_ii('P307',ix_)
        self.xnames=arrset(self.xnames,iv,'P307')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(3.50000e-04))
        [iv,ix_,_] = s2mpj_ii('P503',ix_)
        self.xnames=arrset(self.xnames,iv,'P503')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(3.50000e-04))
        [iv,ix_,_] = s2mpj_ii('P401',ix_)
        self.xnames=arrset(self.xnames,iv,'P401')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(7.20000e-04))
        [iv,ix_,_] = s2mpj_ii('P402',ix_)
        self.xnames=arrset(self.xnames,iv,'P402')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(7.20000e-04))
        [iv,ix_,_] = s2mpj_ii('P401',ix_)
        self.xnames=arrset(self.xnames,iv,'P401')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(5.44000e-03))
        [iv,ix_,_] = s2mpj_ii('P403',ix_)
        self.xnames=arrset(self.xnames,iv,'P403')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(5.44000e-03))
        [iv,ix_,_] = s2mpj_ii('P403',ix_)
        self.xnames=arrset(self.xnames,iv,'P403')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(3.40000e-03))
        [iv,ix_,_] = s2mpj_ii('P404',ix_)
        self.xnames=arrset(self.xnames,iv,'P404')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(3.40000e-03))
        [iv,ix_,_] = s2mpj_ii('P403',ix_)
        self.xnames=arrset(self.xnames,iv,'P403')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(5.80000e-03))
        [iv,ix_,_] = s2mpj_ii('P405',ix_)
        self.xnames=arrset(self.xnames,iv,'P405')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(5.80000e-03))
        [iv,ix_,_] = s2mpj_ii('P405',ix_)
        self.xnames=arrset(self.xnames,iv,'P405')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(7.90000e-04))
        [iv,ix_,_] = s2mpj_ii('P406',ix_)
        self.xnames=arrset(self.xnames,iv,'P406')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(7.90000e-04))
        [iv,ix_,_] = s2mpj_ii('P405',ix_)
        self.xnames=arrset(self.xnames,iv,'P405')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(3.88000e-03))
        [iv,ix_,_] = s2mpj_ii('P407',ix_)
        self.xnames=arrset(self.xnames,iv,'P407')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(3.88000e-03))
        [iv,ix_,_] = s2mpj_ii('P407',ix_)
        self.xnames=arrset(self.xnames,iv,'P407')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(6.40000e-04))
        [iv,ix_,_] = s2mpj_ii('P501',ix_)
        self.xnames=arrset(self.xnames,iv,'P501')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(6.40000e-04))
        [iv,ix_,_] = s2mpj_ii('P501',ix_)
        self.xnames=arrset(self.xnames,iv,'P501')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.10000e-03))
        [iv,ix_,_] = s2mpj_ii('P502',ix_)
        self.xnames=arrset(self.xnames,iv,'P502')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(1.10000e-03))
        [iv,ix_,_] = s2mpj_ii('P501',ix_)
        self.xnames=arrset(self.xnames,iv,'P501')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(4.20000e-04))
        [iv,ix_,_] = s2mpj_ii('P505',ix_)
        self.xnames=arrset(self.xnames,iv,'P505')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(4.20000e-04))
        [iv,ix_,_] = s2mpj_ii('P502',ix_)
        self.xnames=arrset(self.xnames,iv,'P502')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(3.70000e-04))
        [iv,ix_,_] = s2mpj_ii('P503',ix_)
        self.xnames=arrset(self.xnames,iv,'P503')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(3.70000e-04))
        [iv,ix_,_] = s2mpj_ii('P503',ix_)
        self.xnames=arrset(self.xnames,iv,'P503')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(3.50000e-04))
        [iv,ix_,_] = s2mpj_ii('P504',ix_)
        self.xnames=arrset(self.xnames,iv,'P504')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(3.50000e-04))
        [iv,ix_,_] = s2mpj_ii('P505',ix_)
        self.xnames=arrset(self.xnames,iv,'P505')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(6.20000e-04))
        [iv,ix_,_] = s2mpj_ii('P506',ix_)
        self.xnames=arrset(self.xnames,iv,'P506')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(6.20000e-04))
        [iv,ix_,_] = s2mpj_ii('P506',ix_)
        self.xnames=arrset(self.xnames,iv,'P506')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(7.20000e-04))
        [iv,ix_,_] = s2mpj_ii('P507',ix_)
        self.xnames=arrset(self.xnames,iv,'P507')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(7.20000e-04))
        [iv,ix_,_] = s2mpj_ii('P506',ix_)
        self.xnames=arrset(self.xnames,iv,'P506')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(9.20000e-04))
        [iv,ix_,_] = s2mpj_ii('P508',ix_)
        self.xnames=arrset(self.xnames,iv,'P508')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(9.20000e-04))
        [iv,ix_,_] = s2mpj_ii('P508',ix_)
        self.xnames=arrset(self.xnames,iv,'P508')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(9.00000e-04))
        [iv,ix_,_] = s2mpj_ii('P509',ix_)
        self.xnames=arrset(self.xnames,iv,'P509')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(9.00000e-04))
        [iv,ix_,_] = s2mpj_ii('P508',ix_)
        self.xnames=arrset(self.xnames,iv,'P508')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(3.88000e-03))
        [iv,ix_,_] = s2mpj_ii('P510',ix_)
        self.xnames=arrset(self.xnames,iv,'P510')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(3.88000e-03))
        [iv,ix_,_] = s2mpj_ii('P510',ix_)
        self.xnames=arrset(self.xnames,iv,'P510')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(9.90000e-04))
        [iv,ix_,_] = s2mpj_ii('P511',ix_)
        self.xnames=arrset(self.xnames,iv,'P511')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(9.90000e-04))
        [iv,ix_,_] = s2mpj_ii('Q1',ix_)
        self.xnames=arrset(self.xnames,iv,'Q1')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N      1']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N      2']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q2',ix_)
        self.xnames=arrset(self.xnames,iv,'Q2')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N      2']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N      3']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q3',ix_)
        self.xnames=arrset(self.xnames,iv,'Q3')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N      3']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N      4']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q4',ix_)
        self.xnames=arrset(self.xnames,iv,'Q4')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N      4']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N      5']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q5',ix_)
        self.xnames=arrset(self.xnames,iv,'Q5')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N      5']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N      6']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q6',ix_)
        self.xnames=arrset(self.xnames,iv,'Q6')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N      5']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N     26']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q7',ix_)
        self.xnames=arrset(self.xnames,iv,'Q7')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N      6']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N      9']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q8',ix_)
        self.xnames=arrset(self.xnames,iv,'Q8')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N      6']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    304']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q9',ix_)
        self.xnames=arrset(self.xnames,iv,'Q9')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N      9']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N     10']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q10',ix_)
        self.xnames=arrset(self.xnames,iv,'Q10')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N     10']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N     12']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q11',ix_)
        self.xnames=arrset(self.xnames,iv,'Q11')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N     10']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N     27']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q12',ix_)
        self.xnames=arrset(self.xnames,iv,'Q12')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N     12']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N     13']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q13',ix_)
        self.xnames=arrset(self.xnames,iv,'Q13')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N     13']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N     14']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q14',ix_)
        self.xnames=arrset(self.xnames,iv,'Q14')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N     13']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N     19']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q15',ix_)
        self.xnames=arrset(self.xnames,iv,'Q15')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N     14']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N     15']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q16',ix_)
        self.xnames=arrset(self.xnames,iv,'Q16')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N     16']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N     17']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q17',ix_)
        self.xnames=arrset(self.xnames,iv,'Q17')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N     16']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N     18']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q18',ix_)
        self.xnames=arrset(self.xnames,iv,'Q18')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N     16']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N     26']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q19',ix_)
        self.xnames=arrset(self.xnames,iv,'Q19')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N     18']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N     19']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q20',ix_)
        self.xnames=arrset(self.xnames,iv,'Q20')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N     19']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N     20']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q21',ix_)
        self.xnames=arrset(self.xnames,iv,'Q21')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N     20']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N     21']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q22',ix_)
        self.xnames=arrset(self.xnames,iv,'Q22')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N     22']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    404']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q23',ix_)
        self.xnames=arrset(self.xnames,iv,'Q23')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N     23']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    404']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q24',ix_)
        self.xnames=arrset(self.xnames,iv,'Q24')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N     27']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    404']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q25',ix_)
        self.xnames=arrset(self.xnames,iv,'Q25')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    101']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    102']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q26',ix_)
        self.xnames=arrset(self.xnames,iv,'Q26')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    102']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    103']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q27',ix_)
        self.xnames=arrset(self.xnames,iv,'Q27')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    103']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    104']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q28',ix_)
        self.xnames=arrset(self.xnames,iv,'Q28')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    103']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    111']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q29',ix_)
        self.xnames=arrset(self.xnames,iv,'Q29')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    104']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    105']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q30',ix_)
        self.xnames=arrset(self.xnames,iv,'Q30')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    104']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    110']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q31',ix_)
        self.xnames=arrset(self.xnames,iv,'Q31')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    105']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    106']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q32',ix_)
        self.xnames=arrset(self.xnames,iv,'Q32')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    105']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    112']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q33',ix_)
        self.xnames=arrset(self.xnames,iv,'Q33')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    106']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    107']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q34',ix_)
        self.xnames=arrset(self.xnames,iv,'Q34')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    106']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    109']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q35',ix_)
        self.xnames=arrset(self.xnames,iv,'Q35')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    107']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    201']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q36',ix_)
        self.xnames=arrset(self.xnames,iv,'Q36')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    108']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    109']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q37',ix_)
        self.xnames=arrset(self.xnames,iv,'Q37')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    108']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    210']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q38',ix_)
        self.xnames=arrset(self.xnames,iv,'Q38')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    112']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    509']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q39',ix_)
        self.xnames=arrset(self.xnames,iv,'Q39')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    201']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    202']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q40',ix_)
        self.xnames=arrset(self.xnames,iv,'Q40')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    201']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    510']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q41',ix_)
        self.xnames=arrset(self.xnames,iv,'Q41')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    202']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    203']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q42',ix_)
        self.xnames=arrset(self.xnames,iv,'Q42')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    202']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    211']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q43',ix_)
        self.xnames=arrset(self.xnames,iv,'Q43')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    203']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    204']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q44',ix_)
        self.xnames=arrset(self.xnames,iv,'Q44')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    203']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    502']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q45',ix_)
        self.xnames=arrset(self.xnames,iv,'Q45')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    204']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    205']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q46',ix_)
        self.xnames=arrset(self.xnames,iv,'Q46')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    204']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    208']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q47',ix_)
        self.xnames=arrset(self.xnames,iv,'Q47')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    205']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    206']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q48',ix_)
        self.xnames=arrset(self.xnames,iv,'Q48')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    205']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    207']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q49',ix_)
        self.xnames=arrset(self.xnames,iv,'Q49')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    206']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    301']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q50',ix_)
        self.xnames=arrset(self.xnames,iv,'Q50')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    208']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    209']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q51',ix_)
        self.xnames=arrset(self.xnames,iv,'Q51')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    208']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    210']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q52',ix_)
        self.xnames=arrset(self.xnames,iv,'Q52')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    210']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    211']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q53',ix_)
        self.xnames=arrset(self.xnames,iv,'Q53')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    211']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    212']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q54',ix_)
        self.xnames=arrset(self.xnames,iv,'Q54')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    301']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    302']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q55',ix_)
        self.xnames=arrset(self.xnames,iv,'Q55')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    301']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    304']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q56',ix_)
        self.xnames=arrset(self.xnames,iv,'Q56')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    302']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    303']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q57',ix_)
        self.xnames=arrset(self.xnames,iv,'Q57')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    302']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    305']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q58',ix_)
        self.xnames=arrset(self.xnames,iv,'Q58')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    303']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    401']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q59',ix_)
        self.xnames=arrset(self.xnames,iv,'Q59')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    305']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    306']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q60',ix_)
        self.xnames=arrset(self.xnames,iv,'Q60')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    305']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    309']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q61',ix_)
        self.xnames=arrset(self.xnames,iv,'Q61')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    306']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    307']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q62',ix_)
        self.xnames=arrset(self.xnames,iv,'Q62')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    306']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    308']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q63',ix_)
        self.xnames=arrset(self.xnames,iv,'Q63')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    307']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    503']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q64',ix_)
        self.xnames=arrset(self.xnames,iv,'Q64')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    401']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    402']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q65',ix_)
        self.xnames=arrset(self.xnames,iv,'Q65')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    401']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    403']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q66',ix_)
        self.xnames=arrset(self.xnames,iv,'Q66')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    403']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    404']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q67',ix_)
        self.xnames=arrset(self.xnames,iv,'Q67')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    403']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    405']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q68',ix_)
        self.xnames=arrset(self.xnames,iv,'Q68')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    405']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    406']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q69',ix_)
        self.xnames=arrset(self.xnames,iv,'Q69')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    405']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    407']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q70',ix_)
        self.xnames=arrset(self.xnames,iv,'Q70')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    407']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    501']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q71',ix_)
        self.xnames=arrset(self.xnames,iv,'Q71')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    501']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    502']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q72',ix_)
        self.xnames=arrset(self.xnames,iv,'Q72')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    501']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    505']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q73',ix_)
        self.xnames=arrset(self.xnames,iv,'Q73')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    502']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    503']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q74',ix_)
        self.xnames=arrset(self.xnames,iv,'Q74')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    503']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    504']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q75',ix_)
        self.xnames=arrset(self.xnames,iv,'Q75')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    505']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    506']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q76',ix_)
        self.xnames=arrset(self.xnames,iv,'Q76')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    506']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    507']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q77',ix_)
        self.xnames=arrset(self.xnames,iv,'Q77')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    506']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    508']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q78',ix_)
        self.xnames=arrset(self.xnames,iv,'Q78')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    508']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    509']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q79',ix_)
        self.xnames=arrset(self.xnames,iv,'Q79')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    508']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    510']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('Q80',ix_)
        self.xnames=arrset(self.xnames,iv,'Q80')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    510']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    511']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('S1',ix_)
        self.xnames=arrset(self.xnames,iv,'S1')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N      1']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('S23',ix_)
        self.xnames=arrset(self.xnames,iv,'S23')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N     23']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('S101',ix_)
        self.xnames=arrset(self.xnames,iv,'S101')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N    101']])
        valA = np.append(valA,float(1.0))
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
        self.gconst = arrset(self.gconst,ig_['N      2'],float(7.3070))
        self.gconst = arrset(self.gconst,ig_['N      3'],float(2.4360))
        self.gconst = arrset(self.gconst,ig_['N      4'],float(3.6530))
        self.gconst = arrset(self.gconst,ig_['N      5'],float(4.8710))
        self.gconst = arrset(self.gconst,ig_['N      6'],float(7.3070))
        self.gconst = arrset(self.gconst,ig_['N      9'],float(12.178))
        self.gconst = arrset(self.gconst,ig_['N     10'],float(14.613))
        self.gconst = arrset(self.gconst,ig_['N     12'],float(6.0890))
        self.gconst = arrset(self.gconst,ig_['N     13'],float(9.7420))
        self.gconst = arrset(self.gconst,ig_['N     14'],float(13.395))
        self.gconst = arrset(self.gconst,ig_['N     15'],float(28.008))
        self.gconst = arrset(self.gconst,ig_['N     16'],float(4.8710))
        self.gconst = arrset(self.gconst,ig_['N     17'],float(19.484))
        self.gconst = arrset(self.gconst,ig_['N     18'],float(9.7420))
        self.gconst = arrset(self.gconst,ig_['N     19'],float(6.0890))
        self.gconst = arrset(self.gconst,ig_['N     20'],float(6.0890))
        self.gconst = arrset(self.gconst,ig_['N     21'],float(26.971))
        self.gconst = arrset(self.gconst,ig_['N     22'],float(14.613))
        self.gconst = arrset(self.gconst,ig_['N     26'],float(4.8710))
        self.gconst = arrset(self.gconst,ig_['N     27'],float(8.5240))
        self.gconst = arrset(self.gconst,ig_['N    102'],float(7.0000))
        self.gconst = arrset(self.gconst,ig_['N    103'],float(35.000))
        self.gconst = arrset(self.gconst,ig_['N    104'],float(62.000))
        self.gconst = arrset(self.gconst,ig_['N    105'],float(41.000))
        self.gconst = arrset(self.gconst,ig_['N    106'],float(44.000))
        self.gconst = arrset(self.gconst,ig_['N    107'],float(12.000))
        self.gconst = arrset(self.gconst,ig_['N    108'],float(28.000))
        self.gconst = arrset(self.gconst,ig_['N    109'],float(53.000))
        self.gconst = arrset(self.gconst,ig_['N    110'],float(56.000))
        self.gconst = arrset(self.gconst,ig_['N    111'],float(21.000))
        self.gconst = arrset(self.gconst,ig_['N    112'],float(28.000))
        self.gconst = arrset(self.gconst,ig_['N    201'],float(21.000))
        self.gconst = arrset(self.gconst,ig_['N    202'],float(41.000))
        self.gconst = arrset(self.gconst,ig_['N    203'],float(39.000))
        self.gconst = arrset(self.gconst,ig_['N    204'],float(42.000))
        self.gconst = arrset(self.gconst,ig_['N    205'],float(30.000))
        self.gconst = arrset(self.gconst,ig_['N    206'],float(26.000))
        self.gconst = arrset(self.gconst,ig_['N    207'],float(16.000))
        self.gconst = arrset(self.gconst,ig_['N    208'],float(44.000))
        self.gconst = arrset(self.gconst,ig_['N    209'],float(21.000))
        self.gconst = arrset(self.gconst,ig_['N    210'],float(55.000))
        self.gconst = arrset(self.gconst,ig_['N    211'],float(35.000))
        self.gconst = arrset(self.gconst,ig_['N    212'],float(19.000))
        self.gconst = arrset(self.gconst,ig_['N    301'],float(60.000))
        self.gconst = arrset(self.gconst,ig_['N    302'],float(78.000))
        self.gconst = arrset(self.gconst,ig_['N    303'],float(25.000))
        self.gconst = arrset(self.gconst,ig_['N    304'],float(15.831))
        self.gconst = arrset(self.gconst,ig_['N    305'],float(60.000))
        self.gconst = arrset(self.gconst,ig_['N    306'],float(35.000))
        self.gconst = arrset(self.gconst,ig_['N    307'],float(19.000))
        self.gconst = arrset(self.gconst,ig_['N    308'],float(21.000))
        self.gconst = arrset(self.gconst,ig_['N    309'],float(21.000))
        self.gconst = arrset(self.gconst,ig_['N    401'],float(53.000))
        self.gconst = arrset(self.gconst,ig_['N    402'],float(32.000))
        self.gconst = arrset(self.gconst,ig_['N    403'],float(94.000))
        self.gconst = arrset(self.gconst,ig_['N    404'],float(7.3070))
        self.gconst = arrset(self.gconst,ig_['N    405'],float(88.000))
        self.gconst = arrset(self.gconst,ig_['N    406'],float(21.000))
        self.gconst = arrset(self.gconst,ig_['N    407'],float(37.000))
        self.gconst = arrset(self.gconst,ig_['N    501'],float(35.000))
        self.gconst = arrset(self.gconst,ig_['N    502'],float(32.000))
        self.gconst = arrset(self.gconst,ig_['N    503'],float(14.000))
        self.gconst = arrset(self.gconst,ig_['N    504'],float(7.0000))
        self.gconst = arrset(self.gconst,ig_['N    505'],float(18.000))
        self.gconst = arrset(self.gconst,ig_['N    506'],float(30.000))
        self.gconst = arrset(self.gconst,ig_['N    507'],float(14.000))
        self.gconst = arrset(self.gconst,ig_['N    508'],float(46.000))
        self.gconst = arrset(self.gconst,ig_['N    509'],float(30.000))
        self.gconst = arrset(self.gconst,ig_['N    510'],float(34.000))
        self.gconst = arrset(self.gconst,ig_['N    511'],float(23.000))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        self.xlower[ix_['P1']] = 21.999956
        self.xlower[ix_['P2']] = 21.999956
        self.xlower[ix_['P3']] = 21.999956
        self.xlower[ix_['P4']] = 21.999956
        self.xlower[ix_['P5']] = 21.999956
        self.xlower[ix_['P6']] = 21.999956
        self.xlower[ix_['P9']] = 21.999956
        self.xlower[ix_['P10']] = 21.999956
        self.xlower[ix_['P12']] = 21.999956
        self.xlower[ix_['P13']] = 21.999956
        self.xlower[ix_['P14']] = 21.999956
        self.xlower[ix_['P15']] = 21.999956
        self.xlower[ix_['P16']] = 21.999956
        self.xlower[ix_['P17']] = 21.999956
        self.xlower[ix_['P18']] = 21.999956
        self.xlower[ix_['P19']] = 21.999956
        self.xlower[ix_['P20']] = 21.999956
        self.xlower[ix_['P21']] = 21.999956
        self.xlower[ix_['P22']] = 21.999956
        self.xlower[ix_['P23']] = 21.999956
        self.xlower[ix_['P26']] = 21.999956
        self.xlower[ix_['P27']] = 21.999956
        self.xlower[ix_['P101']] = 21.999956
        self.xlower[ix_['P102']] = 21.999956
        self.xlower[ix_['P103']] = 21.999956
        self.xlower[ix_['P104']] = 21.999956
        self.xlower[ix_['P105']] = 21.999956
        self.xlower[ix_['P106']] = 21.999956
        self.xlower[ix_['P107']] = 21.999956
        self.xlower[ix_['P108']] = 21.999956
        self.xlower[ix_['P109']] = 21.999956
        self.xlower[ix_['P110']] = 21.999956
        self.xlower[ix_['P111']] = 21.999956
        self.xlower[ix_['P112']] = 21.999956
        self.xlower[ix_['P201']] = 21.999956
        self.xlower[ix_['P202']] = 21.999956
        self.xlower[ix_['P203']] = 21.999956
        self.xlower[ix_['P204']] = 21.999956
        self.xlower[ix_['P205']] = 21.999956
        self.xlower[ix_['P206']] = 21.999956
        self.xlower[ix_['P207']] = 21.999956
        self.xlower[ix_['P208']] = 21.999956
        self.xlower[ix_['P209']] = 21.999956
        self.xlower[ix_['P210']] = 21.999956
        self.xlower[ix_['P211']] = 21.999956
        self.xlower[ix_['P212']] = 21.999956
        self.xlower[ix_['P301']] = 21.999956
        self.xlower[ix_['P302']] = 21.999956
        self.xlower[ix_['P303']] = 21.999956
        self.xlower[ix_['P304']] = 21.999956
        self.xlower[ix_['P305']] = 21.999956
        self.xlower[ix_['P306']] = 21.999956
        self.xlower[ix_['P307']] = 21.999956
        self.xlower[ix_['P308']] = 21.999956
        self.xlower[ix_['P309']] = 21.999956
        self.xlower[ix_['P401']] = 21.999956
        self.xlower[ix_['P402']] = 21.999956
        self.xlower[ix_['P403']] = 21.999956
        self.xlower[ix_['P404']] = 21.999956
        self.xlower[ix_['P405']] = 21.999956
        self.xlower[ix_['P406']] = 21.999956
        self.xlower[ix_['P407']] = 21.999956
        self.xlower[ix_['P501']] = 21.999956
        self.xlower[ix_['P502']] = 21.999956
        self.xlower[ix_['P503']] = 21.999956
        self.xlower[ix_['P504']] = 21.999956
        self.xlower[ix_['P505']] = 21.999956
        self.xlower[ix_['P506']] = 21.999956
        self.xlower[ix_['P507']] = 21.999956
        self.xlower[ix_['P508']] = 21.999956
        self.xlower[ix_['P509']] = 21.999956
        self.xlower[ix_['P510']] = 21.999956
        self.xlower[ix_['P511']] = 21.999956
        self.xupper[ix_['P1']] = 50.000000
        self.xupper[ix_['P23']] = 50.000000
        self.xupper[ix_['P101']] = 50.000000
        self.xlower[ix_['S1']] = 0.00000E+00
        self.xupper[ix_['S1']] = 0.10000E+08
        self.xlower[ix_['S23']] = 0.00000E+00
        self.xupper[ix_['S23']] = 0.10000E+08
        self.xlower[ix_['S101']] = 0.00000E+00
        self.xupper[ix_['S101']] = 0.10000E+08
        for I in range(int(v_['1']),int(v_['PIPES'])+1):
            self.xlower[ix_['Q'+str(I)]] = -float('Inf')
            self.xupper[ix_['Q'+str(I)]] = +float('Inf')
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        if('P1' in ix_):
            self.x0[ix_['P1']] = float(49.9999542)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P1']),float(49.9999542)))
        if('P2' in ix_):
            self.x0[ix_['P2']] = float(49.9235382)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P2']),float(49.9235382)))
        if('P3' in ix_):
            self.x0[ix_['P3']] = float(49.8521347)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P3']),float(49.8521347)))
        if('P4' in ix_):
            self.x0[ix_['P4']] = float(49.8023033)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P4']),float(49.8023033)))
        if('P5' in ix_):
            self.x0[ix_['P5']] = float(49.5280037)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P5']),float(49.5280037)))
        if('P6' in ix_):
            self.x0[ix_['P6']] = float(49.4430084)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P6']),float(49.4430084)))
        if('P9' in ix_):
            self.x0[ix_['P9']] = float(49.4192848)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P9']),float(49.4192848)))
        if('P10' in ix_):
            self.x0[ix_['P10']] = float(48.9165802)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P10']),float(48.9165802)))
        if('P12' in ix_):
            self.x0[ix_['P12']] = float(48.4542847)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P12']),float(48.4542847)))
        if('P13' in ix_):
            self.x0[ix_['P13']] = float(48.0059395)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P13']),float(48.0059395)))
        if('P14' in ix_):
            self.x0[ix_['P14']] = float(45.8674431)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P14']),float(45.8674431)))
        if('P15' in ix_):
            self.x0[ix_['P15']] = float(45.0843582)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P15']),float(45.0843582)))
        if('P16' in ix_):
            self.x0[ix_['P16']] = float(49.2248535)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P16']),float(49.2248535)))
        if('P17' in ix_):
            self.x0[ix_['P17']] = float(47.6257820)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P17']),float(47.6257820)))
        if('P18' in ix_):
            self.x0[ix_['P18']] = float(48.7873573)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P18']),float(48.7873573)))
        if('P19' in ix_):
            self.x0[ix_['P19']] = float(47.4473228)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P19']),float(47.4473228)))
        if('P20' in ix_):
            self.x0[ix_['P20']] = float(47.0119705)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P20']),float(47.0119705)))
        if('P21' in ix_):
            self.x0[ix_['P21']] = float(45.4362640)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P21']),float(45.4362640)))
        if('P22' in ix_):
            self.x0[ix_['P22']] = float(49.6542473)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P22']),float(49.6542473)))
        if('P23' in ix_):
            self.x0[ix_['P23']] = float(49.9999542)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P23']),float(49.9999542)))
        if('P26' in ix_):
            self.x0[ix_['P26']] = float(49.3843575)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P26']),float(49.3843575)))
        if('P27' in ix_):
            self.x0[ix_['P27']] = float(49.2706299)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P27']),float(49.2706299)))
        if('P101' in ix_):
            self.x0[ix_['P101']] = float(49.9999542)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P101']),float(49.9999542)))
        if('P102' in ix_):
            self.x0[ix_['P102']] = float(49.7797737)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P102']),float(49.7797737)))
        if('P103' in ix_):
            self.x0[ix_['P103']] = float(49.6217003)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P103']),float(49.6217003)))
        if('P104' in ix_):
            self.x0[ix_['P104']] = float(49.3556252)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P104']),float(49.3556252)))
        if('P105' in ix_):
            self.x0[ix_['P105']] = float(48.9891777)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P105']),float(48.9891777)))
        if('P106' in ix_):
            self.x0[ix_['P106']] = float(48.8152504)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P106']),float(48.8152504)))
        if('P107' in ix_):
            self.x0[ix_['P107']] = float(48.7247696)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P107']),float(48.7247696)))
        if('P108' in ix_):
            self.x0[ix_['P108']] = float(44.6206322)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P108']),float(44.6206322)))
        if('P109' in ix_):
            self.x0[ix_['P109']] = float(44.6906090)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P109']),float(44.6906090)))
        if('P110' in ix_):
            self.x0[ix_['P110']] = float(44.5836792)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P110']),float(44.5836792)))
        if('P111' in ix_):
            self.x0[ix_['P111']] = float(47.5885887)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P111']),float(47.5885887)))
        if('P112' in ix_):
            self.x0[ix_['P112']] = float(44.7669029)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P112']),float(44.7669029)))
        if('P201' in ix_):
            self.x0[ix_['P201']] = float(48.5901833)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P201']),float(48.5901833)))
        if('P202' in ix_):
            self.x0[ix_['P202']] = float(48.2493629)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P202']),float(48.2493629)))
        if('P203' in ix_):
            self.x0[ix_['P203']] = float(47.5757141)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P203']),float(47.5757141)))
        if('P204' in ix_):
            self.x0[ix_['P204']] = float(46.9474792)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P204']),float(46.9474792)))
        if('P205' in ix_):
            self.x0[ix_['P205']] = float(47.2141495)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P205']),float(47.2141495)))
        if('P206' in ix_):
            self.x0[ix_['P206']] = float(48.1905937)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P206']),float(48.1905937)))
        if('P207' in ix_):
            self.x0[ix_['P207']] = float(46.4013824)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P207']),float(46.4013824)))
        if('P208' in ix_):
            self.x0[ix_['P208']] = float(47.0579872)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P208']),float(47.0579872)))
        if('P209' in ix_):
            self.x0[ix_['P209']] = float(45.6997147)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P209']),float(45.6997147)))
        if('P210' in ix_):
            self.x0[ix_['P210']] = float(47.4423180)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P210']),float(47.4423180)))
        if('P211' in ix_):
            self.x0[ix_['P211']] = float(47.8950729)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P211']),float(47.8950729)))
        if('P212' in ix_):
            self.x0[ix_['P212']] = float(46.8515167)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P212']),float(46.8515167)))
        if('P301' in ix_):
            self.x0[ix_['P301']] = float(48.4113693)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P301']),float(48.4113693)))
        if('P302' in ix_):
            self.x0[ix_['P302']] = float(48.6494293)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P302']),float(48.6494293)))
        if('P303' in ix_):
            self.x0[ix_['P303']] = float(49.0014572)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P303']),float(49.0014572)))
        if('P304' in ix_):
            self.x0[ix_['P304']] = float(48.6789932)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P304']),float(48.6789932)))
        if('P305' in ix_):
            self.x0[ix_['P305']] = float(47.3707924)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P305']),float(47.3707924)))
        if('P306' in ix_):
            self.x0[ix_['P306']] = float(47.2600479)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P306']),float(47.2600479)))
        if('P307' in ix_):
            self.x0[ix_['P307']] = float(46.5539703)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P307']),float(46.5539703)))
        if('P308' in ix_):
            self.x0[ix_['P308']] = float(46.2514153)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P308']),float(46.2514153)))
        if('P309' in ix_):
            self.x0[ix_['P309']] = float(45.8382378)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P309']),float(45.8382378)))
        if('P401' in ix_):
            self.x0[ix_['P401']] = float(49.1842041)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P401']),float(49.1842041)))
        if('P402' in ix_):
            self.x0[ix_['P402']] = float(46.2833633)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P402']),float(46.2833633)))
        if('P403' in ix_):
            self.x0[ix_['P403']] = float(49.4383583)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P403']),float(49.4383583)))
        if('P404' in ix_):
            self.x0[ix_['P404']] = float(49.8198280)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P404']),float(49.8198280)))
        if('P405' in ix_):
            self.x0[ix_['P405']] = float(47.8698006)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P405']),float(47.8698006)))
        if('P406' in ix_):
            self.x0[ix_['P406']] = float(46.3379631)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P406']),float(46.3379631)))
        if('P407' in ix_):
            self.x0[ix_['P407']] = float(47.7081528)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P407']),float(47.7081528)))
        if('P501' in ix_):
            self.x0[ix_['P501']] = float(46.5787659)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P501']),float(46.5787659)))
        if('P502' in ix_):
            self.x0[ix_['P502']] = float(47.3301430)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P502']),float(47.3301430)))
        if('P503' in ix_):
            self.x0[ix_['P503']] = float(46.5575447)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P503']),float(46.5575447)))
        if('P504' in ix_):
            self.x0[ix_['P504']] = float(46.4538345)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P504']),float(46.4538345)))
        if('P505' in ix_):
            self.x0[ix_['P505']] = float(46.0914383)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P505']),float(46.0914383)))
        if('P506' in ix_):
            self.x0[ix_['P506']] = float(46.1212387)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P506']),float(46.1212387)))
        if('P507' in ix_):
            self.x0[ix_['P507']] = float(45.4262505)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P507']),float(45.4262505)))
        if('P508' in ix_):
            self.x0[ix_['P508']] = float(47.4120369)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P508']),float(47.4120369)))
        if('P509' in ix_):
            self.x0[ix_['P509']] = float(44.7292328)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P509']),float(44.7292328)))
        if('P510' in ix_):
            self.x0[ix_['P510']] = float(47.9998589)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P510']),float(47.9998589)))
        if('P511' in ix_):
            self.x0[ix_['P511']] = float(45.7520485)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P511']),float(45.7520485)))
        if('Q1' in ix_):
            self.x0[ix_['Q1']] = float(0.192685E+03)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q1']),float(0.192685E+03)))
        if('Q2' in ix_):
            self.x0[ix_['Q2']] = float(0.185378E+03)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q2']),float(0.185378E+03)))
        if('Q3' in ix_):
            self.x0[ix_['Q3']] = float(0.182942E+03)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q3']),float(0.182942E+03)))
        if('Q4' in ix_):
            self.x0[ix_['Q4']] = float(0.179289E+03)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q4']),float(0.179289E+03)))
        if('Q5' in ix_):
            self.x0[ix_['Q5']] = float(0.119560E+03)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q5']),float(0.119560E+03)))
        if('Q6' in ix_):
            self.x0[ix_['Q6']] = float(0.548582E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q6']),float(0.548582E+02)))
        if('Q7' in ix_):
            self.x0[ix_['Q7']] = float(0.194484E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q7']),float(0.194484E+02)))
        if('Q8' in ix_):
            self.x0[ix_['Q8']] = float(0.928046E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q8']),float(0.928046E+02)))
        if('Q9' in ix_):
            self.x0[ix_['Q9']] = float(0.727037E+01)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q9']),float(0.727037E+01)))
        if('Q10' in ix_):
            self.x0[ix_['Q10']] = float(0.804928E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q10']),float(0.804928E+02)))
        if('Q11' in ix_):
            self.x0[ix_['Q11']] = float(-.878354E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q11']),float(-.878354E+02)))
        if('Q12' in ix_):
            self.x0[ix_['Q12']] = float(0.744038E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q12']),float(0.744038E+02)))
        if('Q13' in ix_):
            self.x0[ix_['Q13']] = float(0.414030E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q13']),float(0.414030E+02)))
        if('Q14' in ix_):
            self.x0[ix_['Q14']] = float(0.232588E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q14']),float(0.232588E+02)))
        if('Q15' in ix_):
            self.x0[ix_['Q15']] = float(0.280080E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q15']),float(0.280080E+02)))
        if('Q16' in ix_):
            self.x0[ix_['Q16']] = float(0.194840E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q16']),float(0.194840E+02)))
        if('Q17' in ix_):
            self.x0[ix_['Q17']] = float(0.256322E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q17']),float(0.256322E+02)))
        if('Q18' in ix_):
            self.x0[ix_['Q18']] = float(-.499872E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q18']),float(-.499872E+02)))
        if('Q19' in ix_):
            self.x0[ix_['Q19']] = float(0.158902E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q19']),float(0.158902E+02)))
        if('Q20' in ix_):
            self.x0[ix_['Q20']] = float(0.330600E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q20']),float(0.330600E+02)))
        if('Q21' in ix_):
            self.x0[ix_['Q21']] = float(0.269710E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q21']),float(0.269710E+02)))
        if('Q22' in ix_):
            self.x0[ix_['Q22']] = float(-.146130E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q22']),float(-.146130E+02)))
        if('Q23' in ix_):
            self.x0[ix_['Q23']] = float(0.785198E+03)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q23']),float(0.785198E+03)))
        if('Q24' in ix_):
            self.x0[ix_['Q24']] = float(-.963594E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q24']),float(-.963594E+02)))
        if('Q25' in ix_):
            self.x0[ix_['Q25']] = float(0.959108E+03)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q25']),float(0.959108E+03)))
        if('Q26' in ix_):
            self.x0[ix_['Q26']] = float(0.952108E+03)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q26']),float(0.952108E+03)))
        if('Q27' in ix_):
            self.x0[ix_['Q27']] = float(0.896108E+03)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q27']),float(0.896108E+03)))
        if('Q28' in ix_):
            self.x0[ix_['Q28']] = float(0.210000E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q28']),float(0.210000E+02)))
        if('Q29' in ix_):
            self.x0[ix_['Q29']] = float(0.778108E+03)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q29']),float(0.778108E+03)))
        if('Q30' in ix_):
            self.x0[ix_['Q30']] = float(0.560000E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q30']),float(0.560000E+02)))
        if('Q31' in ix_):
            self.x0[ix_['Q31']] = float(0.705999E+03)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q31']),float(0.705999E+03)))
        if('Q32' in ix_):
            self.x0[ix_['Q32']] = float(0.311093E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q32']),float(0.311093E+02)))
        if('Q33' in ix_):
            self.x0[ix_['Q33']] = float(0.604474E+03)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q33']),float(0.604474E+03)))
        if('Q34' in ix_):
            self.x0[ix_['Q34']] = float(0.575246E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q34']),float(0.575246E+02)))
        if('Q35' in ix_):
            self.x0[ix_['Q35']] = float(0.592474E+03)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q35']),float(0.592474E+03)))
        if('Q36' in ix_):
            self.x0[ix_['Q36']] = float(-.452460E+01)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q36']),float(-.452460E+01)))
        if('Q37' in ix_):
            self.x0[ix_['Q37']] = float(-.234754E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q37']),float(-.234754E+02)))
        if('Q38' in ix_):
            self.x0[ix_['Q38']] = float(0.310933E+01)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q38']),float(0.310933E+01)))
        if('Q39' in ix_):
            self.x0[ix_['Q39']] = float(0.395173E+03)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q39']),float(0.395173E+03)))
        if('Q40' in ix_):
            self.x0[ix_['Q40']] = float(0.176301E+03)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q40']),float(0.176301E+03)))
        if('Q41' in ix_):
            self.x0[ix_['Q41']] = float(0.144907E+03)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q41']),float(0.144907E+03)))
        if('Q42' in ix_):
            self.x0[ix_['Q42']] = float(0.209265E+03)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q42']),float(0.209265E+03)))
        if('Q43' in ix_):
            self.x0[ix_['Q43']] = float(0.213451E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q43']),float(0.213451E+02)))
        if('Q44' in ix_):
            self.x0[ix_['Q44']] = float(0.845622E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q44']),float(0.845622E+02)))
        if('Q45' in ix_):
            self.x0[ix_['Q45']] = float(-.886488E+01)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q45']),float(-.886488E+01)))
        if('Q46' in ix_):
            self.x0[ix_['Q46']] = float(-.117900E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q46']),float(-.117900E+02)))
        if('Q47' in ix_):
            self.x0[ix_['Q47']] = float(-.548649E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q47']),float(-.548649E+02)))
        if('Q48' in ix_):
            self.x0[ix_['Q48']] = float(0.160000E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q48']),float(0.160000E+02)))
        if('Q49' in ix_):
            self.x0[ix_['Q49']] = float(-.808649E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q49']),float(-.808649E+02)))
        if('Q50' in ix_):
            self.x0[ix_['Q50']] = float(0.210000E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q50']),float(0.210000E+02)))
        if('Q51' in ix_):
            self.x0[ix_['Q51']] = float(-.767900E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q51']),float(-.767900E+02)))
        if('Q52' in ix_):
            self.x0[ix_['Q52']] = float(-.155265E+03)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q52']),float(-.155265E+03)))
        if('Q53' in ix_):
            self.x0[ix_['Q53']] = float(0.190000E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q53']),float(0.190000E+02)))
        if('Q54' in ix_):
            self.x0[ix_['Q54']] = float(-.638913E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q54']),float(-.638913E+02)))
        if('Q55' in ix_):
            self.x0[ix_['Q55']] = float(-.769736E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q55']),float(-.769736E+02)))
        if('Q56' in ix_):
            self.x0[ix_['Q56']] = float(-.297006E+03)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q56']),float(-.297006E+03)))
        if('Q57' in ix_):
            self.x0[ix_['Q57']] = float(0.155115E+03)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q57']),float(0.155115E+03)))
        if('Q58' in ix_):
            self.x0[ix_['Q58']] = float(-.322006E+03)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q58']),float(-.322006E+03)))
        if('Q59' in ix_):
            self.x0[ix_['Q59']] = float(0.741150E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q59']),float(0.741150E+02)))
        if('Q60' in ix_):
            self.x0[ix_['Q60']] = float(0.210000E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q60']),float(0.210000E+02)))
        if('Q61' in ix_):
            self.x0[ix_['Q61']] = float(0.181150E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q61']),float(0.181150E+02)))
        if('Q62' in ix_):
            self.x0[ix_['Q62']] = float(0.210000E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q62']),float(0.210000E+02)))
        if('Q63' in ix_):
            self.x0[ix_['Q63']] = float(-.884952E+00)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q63']),float(-.884952E+00)))
        if('Q64' in ix_):
            self.x0[ix_['Q64']] = float(0.320000E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q64']),float(0.320000E+02)))
        if('Q65' in ix_):
            self.x0[ix_['Q65']] = float(-.407006E+03)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q65']),float(-.407006E+03)))
        if('Q66' in ix_):
            self.x0[ix_['Q66']] = float(-.666918E+03)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q66']),float(-.666918E+03)))
        if('Q67' in ix_):
            self.x0[ix_['Q67']] = float(0.165912E+03)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q67']),float(0.165912E+03)))
        if('Q68' in ix_):
            self.x0[ix_['Q68']] = float(0.210000E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q68']),float(0.210000E+02)))
        if('Q69' in ix_):
            self.x0[ix_['Q69']] = float(0.569121E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q69']),float(0.569121E+02)))
        if('Q70' in ix_):
            self.x0[ix_['Q70']] = float(0.199121E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q70']),float(0.199121E+02)))
        if('Q71' in ix_):
            self.x0[ix_['Q71']] = float(-.306772E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q71']),float(-.306772E+02)))
        if('Q72' in ix_):
            self.x0[ix_['Q72']] = float(0.155893E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q72']),float(0.155893E+02)))
        if('Q73' in ix_):
            self.x0[ix_['Q73']] = float(0.218850E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q73']),float(0.218850E+02)))
        if('Q74' in ix_):
            self.x0[ix_['Q74']] = float(0.700000E+01)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q74']),float(0.700000E+01)))
        if('Q75' in ix_):
            self.x0[ix_['Q75']] = float(-.241070E+01)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q75']),float(-.241070E+01)))
        if('Q76' in ix_):
            self.x0[ix_['Q76']] = float(0.140000E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q76']),float(0.140000E+02)))
        if('Q77' in ix_):
            self.x0[ix_['Q77']] = float(-.464107E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q77']),float(-.464107E+02)))
        if('Q78' in ix_):
            self.x0[ix_['Q78']] = float(0.268907E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q78']),float(0.268907E+02)))
        if('Q79' in ix_):
            self.x0[ix_['Q79']] = float(-.119301E+03)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q79']),float(-.119301E+03)))
        if('Q80' in ix_):
            self.x0[ix_['Q80']] = float(0.230000E+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Q80']),float(0.230000E+02)))
        if('S1' in ix_):
            self.x0[ix_['S1']] = float(0.192685E+03)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['S1']),float(0.192685E+03)))
        if('S23' in ix_):
            self.x0[ix_['S23']] = float(0.785198E+03)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['S23']),float(0.785198E+03)))
        if('S101' in ix_):
            self.x0[ix_['S101']] = float(0.959108E+03)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['S101']),float(0.959108E+03)))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQUARE', iet_)
        elftv = loaset(elftv,it,0,'P')
        [it,iet_,_] = s2mpj_ii( 'ePANHAN', iet_)
        elftv = loaset(elftv,it,0,'Q')
        elftp = []
        elftp = loaset(elftp,it,0,'A1')
        elftp = loaset(elftp,it,1,'A2')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        ename = 'PSQR1'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR2'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR3'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR4'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR5'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR6'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR9'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P9'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR10'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P10'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR12'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P12'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR13'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P13'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR14'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P14'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR15'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P15'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR16'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P16'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR17'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P17'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR18'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P18'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR19'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P19'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR20'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P20'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR21'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P21'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR22'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P22'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR23'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P23'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR26'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P26'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR27'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P27'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR101'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P101'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR102'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P102'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR103'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P103'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR104'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P104'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR105'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P105'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR106'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P106'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR107'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P107'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR108'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P108'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR109'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P109'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR110'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P110'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR111'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P111'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR112'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P112'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR201'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P201'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR202'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P202'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR203'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P203'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR204'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P204'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR205'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P205'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR206'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P206'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR207'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P207'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR208'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P208'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR209'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P209'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR210'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P210'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR211'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P211'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR212'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P212'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR301'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P301'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR302'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P302'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR303'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P303'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR304'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P304'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR305'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P305'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR306'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P306'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR307'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P307'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR308'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P308'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR309'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P309'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR401'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P401'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR402'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P402'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR403'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P403'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR404'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P404'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR405'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P405'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR406'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P406'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR407'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P407'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR501'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P501'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR502'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P502'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR503'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P503'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR504'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P504'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR505'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P505'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR506'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P506'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR507'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P507'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR508'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P508'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR509'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P509'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR510'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P510'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PSQR511'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        vname = 'P511'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='P')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        for I in range(int(v_['1']),int(v_['PIPES'])+1):
            ename = 'PANH'+str(I)
            [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'ePANHAN')
            ielftype = arrset(ielftype,ie,iet_["ePANHAN"])
            vname = 'Q'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='Q')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PANH1'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-7.07259e-07))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(1.37621e+02))
        ename = 'PANH2'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-7.07259e-07))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(1.37621e+02))
        ename = 'PANH3'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-5.05185e-07))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(1.37621e+02))
        ename = 'PANH4'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-2.87955e-06))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(1.37621e+02))
        ename = 'PANH5'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.97677e-06))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(1.94163e+02))
        ename = 'PANH6'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.29902e-05))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(1.94163e+02))
        ename = 'PANH7'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.27078e-05))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(1.94163e+02))
        ename = 'PANH8'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-2.76748e-05))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(1.94163e+02))
        ename = 'PANH9'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.89567e-03))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(4.92083e+02))
        ename = 'PANH10'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-2.14621e-05))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(1.94163e+02))
        ename = 'PANH11'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.41198e-05))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(1.94163e+02))
        ename = 'PANH12'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-2.75247e-05))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(3.43580e+02))
        ename = 'PANH13'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-3.65685e-04))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(3.43580e+02))
        ename = 'PANH14'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-2.59518e-04))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(3.43580e+02))
        ename = 'PANH15'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-2.63450e-04))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(3.43580e+02))
        ename = 'PANH16'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.11371e-03))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(4.92083e+02))
        ename = 'PANH17'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.89567e-04))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(4.92083e+02))
        ename = 'PANH18'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.69438e-05))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(1.94163e+02))
        ename = 'PANH19'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.32697e-03))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(4.92083e+02))
        ename = 'PANH20'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.10099e-04))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(3.43580e+02))
        ename = 'PANH21'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-5.66222e-04))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(3.43580e+02))
        ename = 'PANH22'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.89567e-04))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(4.92083e+02))
        ename = 'PANH23'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.27360e-07))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(9.63346e+01))
        ename = 'PANH24'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.86381e-05))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(1.94163e+02))
        ename = 'PANH25'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.06357e-07))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(8.60722e+01))
        ename = 'PANH26'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-7.73503e-08))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(8.60722e+01))
        ename = 'PANH27'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.48586e-07))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(9.63346e+01))
        ename = 'PANH28'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.24403e-03))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(4.92083e+02))
        ename = 'PANH29'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-2.63210e-07))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(9.63346e+01))
        ename = 'PANH30'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-4.81682e-04))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(3.43580e+02))
        ename = 'PANH31'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.48586e-07))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(9.63346e+01))
        ename = 'PANH32'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.30327e-03))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(4.92083e+02))
        ename = 'PANH33'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.01888e-07))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(9.63346e+01))
        ename = 'PANH34'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-3.97142e-04))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(3.43580e+02))
        ename = 'PANH35'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.57077e-07))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(9.63346e+01))
        ename = 'PANH36'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-5.80549e-04))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(4.92083e+02))
        ename = 'PANH37'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.42175e-03))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(4.92083e+02))
        ename = 'PANH38'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-5.80549e-04))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(4.92083e+02))
        ename = 'PANH39'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-8.84073e-07))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(1.37621e+02))
        ename = 'PANH40'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-6.91870e-06))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(1.94163e+02))
        ename = 'PANH41'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.11546e-05))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(1.94163e+02))
        ename = 'PANH42'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-2.82903e-06))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(1.37621e+02))
        ename = 'PANH43'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-3.38160e-04))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(3.43580e+02))
        ename = 'PANH44'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.04487e-05))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(1.94163e+02))
        ename = 'PANH45'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-7.10876e-04))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(4.92083e+02))
        ename = 'PANH46'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.67114e-04))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(3.43580e+02))
        ename = 'PANH47'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.02234e-04))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(3.43580e+02))
        ename = 'PANH48'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-7.93812e-04))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(4.92083e+02))
        ename = 'PANH49'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.01663e-05))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(1.94163e+02))
        ename = 'PANH50'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-8.29355e-04))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(4.92083e+02))
        ename = 'PANH51'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.93441e-05))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(1.94163e+02))
        ename = 'PANH52'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-6.63631e-06))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(1.94163e+02))
        ename = 'PANH53'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-7.58268e-04))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(4.92083e+02))
        ename = 'PANH54'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.65202e-05))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(1.94163e+02))
        ename = 'PANH55'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.34138e-05))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(1.94163e+02))
        ename = 'PANH56'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.51555e-06))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(1.37621e+02))
        ename = 'PANH57'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.87793e-05))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(1.94163e+02))
        ename = 'PANH58'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-6.81999e-07))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(1.37621e+02))
        ename = 'PANH59'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-5.93032e-06))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(1.94163e+02))
        ename = 'PANH60'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-9.35987e-04))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(4.92083e+02))
        ename = 'PANH61'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-5.56853e-04))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(4.92083e+02))
        ename = 'PANH62'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-6.16093e-04))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(4.92083e+02))
        ename = 'PANH63'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-4.14678e-04))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(4.92083e+02))
        ename = 'PANH64'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-8.53051e-04))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(4.92083e+02))
        ename = 'PANH65'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-5.77363e-07))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(9.63346e+01))
        ename = 'PANH66'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-3.60852e-07))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(9.63346e+01))
        ename = 'PANH67'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-2.04737e-05))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(1.94163e+02))
        ename = 'PANH68'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-9.35987e-04))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(4.92083e+02))
        ename = 'PANH69'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.36962e-05))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(1.94163e+02))
        ename = 'PANH70'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-7.58268e-04))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(4.92083e+02))
        ename = 'PANH71'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-2.16265e-04))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(3.43580e+02))
        ename = 'PANH72'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-4.97613e-04))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(4.92083e+02))
        ename = 'PANH73'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-4.38374e-04))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(4.92083e+02))
        ename = 'PANH74'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-4.14678e-04))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(4.92083e+02))
        ename = 'PANH75'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-7.34572e-04))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(4.92083e+02))
        ename = 'PANH76'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-8.53051e-04))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(4.92083e+02))
        ename = 'PANH77'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.80876e-04))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(3.43580e+02))
        ename = 'PANH78'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.06631e-03))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(4.92083e+02))
        ename = 'PANH79'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.36962e-05))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(1.94163e+02))
        ename = 'PANH80'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_['eSQUARE'])
        posep = np.where(elftp[ielftype[ie]]=='A1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.17295e-03))
        posep = np.where(elftp[ielftype[ie]]=='A2')[0]
        loaset(self.elpar,ie,posep[0],float(4.92083e+02))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        ig = ig_['PIP    1']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR1'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR2'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH1'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP    2']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR2'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR3'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH2'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP    3']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR3'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR4'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH3'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP    4']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR4'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR5'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH4'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP    5']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR5'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR6'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH5'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP    6']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR5'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR26'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH6'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP    7']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR6'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR9'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH7'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP    8']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR6'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR304'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH8'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP    9']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR9'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR10'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH9'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   10']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR10'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR12'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH10'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   11']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR10'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR27'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH11'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   12']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR12'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR13'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH12'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   13']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR13'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR14'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH13'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   14']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR13'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR19'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH14'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   15']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR14'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR15'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH15'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   16']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR16'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR17'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH16'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   17']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR16'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR18'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH17'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   18']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR16'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR26'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH18'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   19']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR18'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR19'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH19'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   20']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR19'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR20'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH20'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   21']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR20'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR21'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH21'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   22']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR22'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR404'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH22'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   23']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR23'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR404'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH23'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   24']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR27'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR404'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH24'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   25']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR101'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR102'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH25'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   26']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR102'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR103'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH26'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   27']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR103'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR104'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH27'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   28']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR103'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR111'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH28'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   29']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR104'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR105'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH29'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   30']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR104'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR110'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH30'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   31']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR105'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR106'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH31'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   32']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR105'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR112'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH32'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   33']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR106'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR107'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH33'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   34']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR106'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR109'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH34'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   35']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR107'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR201'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH35'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   36']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR108'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR109'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH36'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   37']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR108'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR210'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH37'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   38']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR112'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR509'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH38'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   39']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR201'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR202'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH39'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   40']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR201'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR510'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH40'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   41']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR202'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR203'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH41'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   42']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR202'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR211'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH42'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   43']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR203'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR204'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH43'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   44']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR203'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR502'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH44'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   45']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR204'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR205'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH45'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   46']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR204'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR208'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH46'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   47']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR205'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR206'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH47'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   48']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR205'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR207'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH48'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   49']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR206'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR301'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH49'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   50']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR208'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR209'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH50'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   51']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR208'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR210'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH51'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   52']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR210'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR211'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH52'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   53']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR211'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR212'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH53'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   54']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR301'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR302'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH54'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   55']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR301'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR304'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH55'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   56']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR302'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR303'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH56'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   57']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR302'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR305'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH57'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   58']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR303'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR401'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH58'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   59']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR305'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR306'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH59'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   60']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR305'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR309'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH60'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   61']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR306'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR307'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH61'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   62']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR306'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR308'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH62'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   63']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR307'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR503'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH63'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   64']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR401'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR402'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH64'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   65']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR401'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR403'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH65'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   66']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR403'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR404'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH66'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   67']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR403'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR405'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH67'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   68']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR405'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR406'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH68'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   69']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR405'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR407'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH69'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   70']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR407'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR501'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH70'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   71']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR501'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR502'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH71'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   72']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR501'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR505'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH72'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   73']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR502'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR503'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH73'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   74']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR503'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR504'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH74'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   75']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR505'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR506'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH75'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   76']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR506'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR507'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH76'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   77']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR506'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR508'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH77'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   78']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR508'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR509'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH78'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   79']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR508'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR510'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH79'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['PIP   80']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR510'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQR511'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PANH80'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        self.objlower = 0.0
#    Solution
# LO SOLUTN               8.00028D+00
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
        self.pbclass   = "C-CLOR2-RN-156-153"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def e_globs(self):

        import numpy as np
        self.efpar = np.array([]);
        self.efpar = arrset( self.efpar,0,1.0e-3)
        self.efpar = arrset( self.efpar,1,1.01325)
        self.efpar = arrset( self.efpar,2,3.62e-2)
        self.efpar = arrset( self.efpar,3,3.5657e+0)
        self.efpar = arrset( self.efpar,4,1.47519e+1)
        self.efpar = arrset( self.efpar,5,1.0e-1)
        self.efpar = arrset( self.efpar,6,np.log10(np.exp(1.0e+0)))
        return pbm

    @staticmethod
    def eSQUARE(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = (self.efpar[0]*EV_[0]+self.efpar[1])**2
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0e+0*self.efpar[0]*(self.efpar[0]*EV_[0]+self.efpar[1])
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0e+0*self.efpar[0]*self.efpar[0]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def ePANHAN(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        QGE = EV_[0]>=self.efpar[5]
        QLE = EV_[0]<=-self.efpar[5]
        QELSE = not(QGE or QLE)
        if QELSE!=0:
            QRATIO = EV_[0]/self.efpar[5]
        if QGE!=0:
            H = EV_[0]
        if QLE!=0:
            H = -EV_[0]
        if QELSE!=0:
            H = self.efpar[5]*(3.75e-1+7.5e-1*QRATIO**2-1.25e-1*QRATIO**4)
        if QGE!=0:
            H1 = 1.0e+0
        if QLE!=0:
            H1 = -1.0e+0
        if QELSE!=0:
            H1 = 1.5e+0*QRATIO-5.0e-1*QRATIO**3
        if QGE!=0:
            H2 = 0.0e+0
        if QLE!=0:
            H2 = 0.0e+0
        if QELSE!=0:
            H2 = (1.5e+0-1.5e+0*QRATIO**2)/self.efpar[5]
        ARGLOG = self.elpar[iel_][1]*H
        X = np.log10(ARGLOG)-5.0e+0
        X1 = self.efpar[6]*H1/H
        X2 = self.efpar[6]*(H2/H-(H1/H)**2)
        FROOT = 1.0/((self.efpar[2]*X+self.efpar[3])*X+self.efpar[4])
        DERIV = 2.0*self.efpar[2]*X+self.efpar[3]
        F = FROOT*FROOT
        F1 = -2.0*DERIV*X1*FROOT**3
        F2 = -2.0*FROOT**3*(DERIV*X2+2.0*self.efpar[2]*X1**2-0.75*F1**2/FROOT**5)
        f_   = self.elpar[iel_][0]*F*EV_[0]*H
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = self.elpar[iel_][0]*(F*H+EV_[0]*(F1*H+F*H1))
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = self.elpar[iel_][0]*(2.0*(F1*H+F*H1)+EV_[0]*(F2*H+2.0*F1*H1+F*H2))
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

