function varargout = SPMSRTLS(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : SPMSRTLS
%    *********
% 
%    Liu and Nocedal tridiagonal matrix square root problem.
% 
%    Source:  problem 151 (p. 93) in
%    A.R. Buckley,
%    "Test functions for unconstrained minimization",
%    TR 1989CS-3, Mathematics, statistics and computing centre,
%    Dalhousie University, Halifax (CDN), 1989.
% 
%    This is a least-squares variant of problem SPMSQRT.
% 
%    SIF input: Ph. Toint, Dec 1989.
% 
%    classification = 'C-CSUR2-AN-V-V'
% 
%    M is the dimension of the matrix
%    The number of variables is 3*M-2
% 
%       Alternative values for the SIF file parameters:
% IE M                   10             $-PARAMETER n = 28     original value
% IE M                   34             $-PARAMETER n = 100
% IE M                   167            $-PARAMETER n = 499
% IE M                   334            $-PARAMETER n = 1000
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 21 VI 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'SPMSRTLS';

switch(action)

    case {'setup','setup_redprec'}

        if(strcmp(action,'setup_redprec'))
            if(isfield(pbm,'ndigs'))
                rmfield(pbm,'ndigs');
            end
            pbm.ndigs = max(1,min(15,varargin{end}));
            nargs     = nargin-2;
        else
            nargs = nargin-1;
        end
        pb.name   = name;
        pbm.name  = name;
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = containers.Map('KeyType','char', 'ValueType', 'double');
        ix_ = containers.Map('KeyType','char', 'ValueType', 'double');
        ig_ = containers.Map('KeyType','char', 'ValueType', 'double');
        if(nargs<1)
            v_('M') = 1667;  %  SIF file default value
        else
            v_('M') = varargin{1};
        end
% IE M                   3334           $-PARAMETER n = 10000
        v_('M-1') = -1+v_('M');
        v_('M-2') = -2+v_('M');
        v_('M-3') = -3+v_('M');
        v_('1') = 1;
        v_('2') = 2;
        v_('3') = 3;
        v_('4') = 4;
        v_('B1,1') = sin(1.0);
        v_('B1,2') = sin(4.0);
        v_('K') = 2;
        for I=v_('2'):v_('M-1')
            v_('I-1') = -1+I;
            v_('I+1') = 1+I;
            v_('K') = 1+v_('K');
            v_('RK') = v_('K');
            v_('RKSQR') = v_('RK')*v_('RK');
            v_(['B',int2str(I),',',int2str(round(v_('I-1')))]) = sin(v_('RKSQR'));
            v_('K') = 1+v_('K');
            v_('RK') = v_('K');
            v_('RKSQR') = v_('RK')*v_('RK');
            v_(['B',int2str(I),',',int2str(I)]) = sin(v_('RKSQR'));
            v_('K') = 1+v_('K');
            v_('RK') = v_('K');
            v_('RKSQR') = v_('RK')*v_('RK');
            v_(['B',int2str(I),',',int2str(round(v_('I+1')))]) = sin(v_('RKSQR'));
        end
        v_('K') = 1+v_('K');
        v_('RK') = v_('K');
        v_('RKSQR') = v_('RK')*v_('RK');
        v_(['B',int2str(round(v_('M'))),',',int2str(round(v_('M-1')))]) =...
              sin(v_('RKSQR'));
        v_('K') = 1+v_('K');
        v_('RK') = v_('K');
        v_('RKSQR') = v_('RK')*v_('RK');
        v_(['B',int2str(round(v_('M'))),',',int2str(round(v_('M')))]) =...
              sin(v_('RKSQR'));
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        [iv,ix_] =...
              s2mpjlib('ii',['X',int2str(round(v_('1'))),',',int2str(round(v_('1')))],ix_);
        pb.xnames{iv} = ['X',int2str(round(v_('1'))),',',int2str(round(v_('1')))];
        [iv,ix_] =...
              s2mpjlib('ii',['X',int2str(round(v_('1'))),',',int2str(round(v_('2')))],ix_);
        pb.xnames{iv} = ['X',int2str(round(v_('1'))),',',int2str(round(v_('2')))];
        for I=v_('2'):v_('M-1')
            v_('I-1') = -1+I;
            v_('I+1') = 1+I;
            [iv,ix_] =...
                  s2mpjlib('ii',['X',int2str(I),',',int2str(round(v_('I-1')))],ix_);
            pb.xnames{iv} = ['X',int2str(I),',',int2str(round(v_('I-1')))];
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I),',',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I),',',int2str(I)];
            [iv,ix_] =...
                  s2mpjlib('ii',['X',int2str(I),',',int2str(round(v_('I+1')))],ix_);
            pb.xnames{iv} = ['X',int2str(I),',',int2str(round(v_('I+1')))];
        end
        [iv,ix_] =...
              s2mpjlib('ii',['X',int2str(round(v_('M'))),',',int2str(round(v_('M-1')))],ix_);
        pb.xnames{iv} = ['X',int2str(round(v_('M'))),',',int2str(round(v_('M-1')))];
        [iv,ix_] =...
              s2mpjlib('ii',['X',int2str(round(v_('M'))),',',int2str(round(v_('M')))],ix_);
        pb.xnames{iv} = ['X',int2str(round(v_('M'))),',',int2str(round(v_('M')))];
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        for J=v_('1'):v_('3')
            [ig,ig_] = s2mpjlib('ii',['E',int2str(round(v_('1'))),',',int2str(J)],ig_);
            gtype{ig} = '<>';
        end
        for J=v_('1'):v_('4')
            [ig,ig_] = s2mpjlib('ii',['E',int2str(round(v_('2'))),',',int2str(J)],ig_);
            gtype{ig} = '<>';
        end
        for I=v_('3'):v_('M-2')
            v_('I-2') = -2+I;
            v_('I+2') = 2+I;
            for J=v_('I-2'):v_('I+2')
                [ig,ig_] = s2mpjlib('ii',['E',int2str(I),',',int2str(J)],ig_);
                gtype{ig} = '<>';
            end
        end
        for J=v_('M-3'):v_('M')
            [ig,ig_] =...
                  s2mpjlib('ii',['E',int2str(round(v_('M-1'))),',',int2str(J)],ig_);
            gtype{ig} = '<>';
        end
        for J=v_('M-2'):v_('M')
            [ig,ig_] = s2mpjlib('ii',['E',int2str(round(v_('M'))),',',int2str(J)],ig_);
            gtype{ig} = '<>';
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = ix_.Count;
        ngrp   = ig_.Count;
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        v_('ENTRY') = v_('B1,1')*v_('B1,1');
        v_('PROD') = v_('B1,2')*v_('B2,1');
        v_('ENTRY') = v_('ENTRY')+v_('PROD');
        pbm.gconst(ig_('E1,1')) = v_('ENTRY');
        for I=v_('2'):v_('M-1')
            v_('I-1') = -1+I;
            v_('I+1') = 1+I;
            v_('ENTRY') =...
                  v_(['B',int2str(I),',',int2str(I)])*v_(['B',int2str(I),',',int2str(I)]);
            v_('PROD') =...
                  v_(['B',int2str(round(v_('I-1'))),',',int2str(I)])*v_(['B',int2str(I),',',int2str(round(v_('I-1')))]);
            v_('ENTRY') = v_('ENTRY')+v_('PROD');
            v_('PROD') =...
                  v_(['B',int2str(round(v_('I+1'))),',',int2str(I)])*v_(['B',int2str(I),',',int2str(round(v_('I+1')))]);
            v_('ENTRY') = v_('ENTRY')+v_('PROD');
            pbm.gconst(ig_(['E',int2str(I),',',int2str(I)])) = v_('ENTRY');
        end
        v_('ENTRY') = v_(['B',int2str(round(v_('M'))),',',int2str(round(v_('M')))])*...
             v_(['B',int2str(round(v_('M'))),',',int2str(round(v_('M')))]);
        v_('PROD') = v_(['B',int2str(round(v_('M-1'))),',',int2str(round(v_('M')))])*...
             v_(['B',int2str(round(v_('M'))),',',int2str(round(v_('M-1')))]);
        v_('ENTRY') = v_('ENTRY')+v_('PROD');
        pbm.gconst(ig_(['E',int2str(round(v_('M'))),',',int2str(round(v_('M')))])) =...
              v_('ENTRY');
        for I=v_('1'):v_('M-1')
            v_('I+1') = 1+I;
            v_('ENTRY') = v_(['B',int2str(round(v_('I+1'))),',',int2str(I)])*...
                 v_(['B',int2str(I),',',int2str(I)]);
            v_('PROD') = v_(['B',int2str(round(v_('I+1'))),',',int2str(round(v_('I+1')))])*...
                 v_(['B',int2str(round(v_('I+1'))),',',int2str(I)]);
            v_('ENTRY') = v_('ENTRY')+v_('PROD');
            pbm.gconst(ig_(['E',int2str(round(v_('I+1'))),',',int2str(I)])) =...
                  v_('ENTRY');
        end
        for I=v_('2'):v_('M')
            v_('I-1') = -1+I;
            v_('ENTRY') = v_(['B',int2str(round(v_('I-1'))),',',int2str(I)])*...
                 v_(['B',int2str(I),',',int2str(I)]);
            v_('PROD') = v_(['B',int2str(round(v_('I-1'))),',',int2str(round(v_('I-1')))])*...
                 v_(['B',int2str(round(v_('I-1'))),',',int2str(I)]);
            v_('ENTRY') = v_('ENTRY')+v_('PROD');
            pbm.gconst(ig_(['E',int2str(round(v_('I-1'))),',',int2str(I)])) =...
                  v_('ENTRY');
        end
        for I=v_('2'):v_('M-1')
            v_('I-1') = -1+I;
            v_('I+1') = 1+I;
            v_('ENTRY') = v_(['B',int2str(round(v_('I+1'))),',',int2str(I)])*...
                 v_(['B',int2str(I),',',int2str(round(v_('I-1')))]);
            pbm.gconst(ig_(['E',int2str(round(v_('I+1'))),',',int2str(round(v_('I-1')))])) = v_('ENTRY');
        end
        for I=v_('2'):v_('M-1')
            v_('I-1') = -1+I;
            v_('I+1') = 1+I;
            v_('ENTRY') = v_(['B',int2str(round(v_('I-1'))),',',int2str(I)])*...
                 v_(['B',int2str(I),',',int2str(round(v_('I+1')))]);
            pbm.gconst(ig_(['E',int2str(round(v_('I-1'))),',',int2str(round(v_('I+1')))])) = v_('ENTRY');
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        v_('PROD') = 0.2*v_('B1,1');
        pb.x0(ix_('X1,1'),1) = v_('PROD');
        v_('PROD') = 0.2*v_('B1,2');
        pb.x0(ix_('X1,2'),1) = v_('PROD');
        for I=v_('2'):v_('M-1')
            v_('I-1') = -1+I;
            v_('I+1') = 1+I;
            v_('PROD') = 0.2*v_(['B',int2str(I),',',int2str(round(v_('I-1')))]);
            pb.x0(ix_(['X',int2str(I),',',int2str(round(v_('I-1')))]),1) = v_('PROD');
            v_('PROD') = 0.2*v_(['B',int2str(I),',',int2str(I)]);
            pb.x0(ix_(['X',int2str(I),',',int2str(I)]),1) = v_('PROD');
            v_('PROD') = 0.2*v_(['B',int2str(I),',',int2str(round(v_('I+1')))]);
            pb.x0(ix_(['X',int2str(I),',',int2str(round(v_('I+1')))]),1) = v_('PROD');
        end
        v_('PROD') =...
              0.2*v_(['B',int2str(round(v_('M'))),',',int2str(round(v_('M-1')))]);
        pb.x0(ix_(['X',int2str(round(v_('M'))),',',int2str(round(v_('M-1')))]),1) =...
              v_('PROD');
        v_('PROD') =...
              0.2*v_(['B',int2str(round(v_('M'))),',',int2str(round(v_('M')))]);
        pb.x0(ix_(['X',int2str(round(v_('M'))),',',int2str(round(v_('M')))]),1) =...
              v_('PROD');
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'ePROD2',iet_);
        elftv{it}{1} = 'VI';
        elftv{it}{2} = 'VJ';
        [it,iet_] = s2mpjlib( 'ii', 'eSQ',iet_);
        elftv{it}{1} = 'V';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = ['D',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQ';
        ielftype(ie) = iet_('eSQ');
        ename = ['D',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('1'))),',',int2str(round(v_('1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['G',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD2';
        ielftype(ie) = iet_('ePROD2');
        ename = ['G',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('1'))),',',int2str(round(v_('2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VI',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['G',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('2'))),',',int2str(round(v_('1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VJ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['H',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD2';
        ielftype(ie) = iet_('ePROD2');
        ename = ['H',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('1'))),',',int2str(round(v_('1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VI',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['H',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('1'))),',',int2str(round(v_('2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VJ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['R',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD2';
        ielftype(ie) = iet_('ePROD2');
        ename = ['R',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('1'))),',',int2str(round(v_('2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VI',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['R',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('2'))),',',int2str(round(v_('2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VJ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['S',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD2';
        ielftype(ie) = iet_('ePROD2');
        ename = ['S',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('1'))),',',int2str(round(v_('2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VI',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['S',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('2'))),',',int2str(round(v_('3')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VJ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['B',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD2';
        ielftype(ie) = iet_('ePROD2');
        ename = ['B',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('2'))),',',int2str(round(v_('1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VI',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['B',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('1'))),',',int2str(round(v_('1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VJ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD2';
        ielftype(ie) = iet_('ePROD2');
        ename = ['C',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('2'))),',',int2str(round(v_('2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VI',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('2'))),',',int2str(round(v_('1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VJ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['F',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD2';
        ielftype(ie) = iet_('ePROD2');
        ename = ['F',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('2'))),',',int2str(round(v_('1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VI',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['F',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('1'))),',',int2str(round(v_('2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VJ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['D',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQ';
        ielftype(ie) = iet_('eSQ');
        ename = ['D',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('2'))),',',int2str(round(v_('2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['G',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD2';
        ielftype(ie) = iet_('ePROD2');
        ename = ['G',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('2'))),',',int2str(round(v_('3')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VI',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['G',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('3'))),',',int2str(round(v_('2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VJ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['H',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD2';
        ielftype(ie) = iet_('ePROD2');
        ename = ['H',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('2'))),',',int2str(round(v_('2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VI',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['H',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('2'))),',',int2str(round(v_('3')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VJ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['R',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD2';
        ielftype(ie) = iet_('ePROD2');
        ename = ['R',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('2'))),',',int2str(round(v_('3')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VI',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['R',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('3'))),',',int2str(round(v_('3')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VJ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['S',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD2';
        ielftype(ie) = iet_('ePROD2');
        ename = ['S',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('2'))),',',int2str(round(v_('3')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VI',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['S',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('3'))),',',int2str(round(v_('4')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VJ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        for I=v_('3'):v_('M-2')
            v_('I-2') = -2+I;
            v_('I-1') = -1+I;
            v_('I+1') = 1+I;
            v_('I+2') = 2+I;
            ename = ['A',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'ePROD2';
            ielftype(ie) = iet_('ePROD2');
            vname = ['X',int2str(I),',',int2str(round(v_('I-1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('VI',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('I-1'))),',',int2str(round(v_('I-2')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('VJ',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'ePROD2';
            ielftype(ie) = iet_('ePROD2');
            vname = ['X',int2str(I),',',int2str(round(v_('I-1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('VI',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('I-1'))),',',int2str(round(v_('I-1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('VJ',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['C',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'ePROD2';
            ielftype(ie) = iet_('ePROD2');
            vname = ['X',int2str(I),',',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('VI',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(I),',',int2str(round(v_('I-1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('VJ',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['F',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'ePROD2';
            ielftype(ie) = iet_('ePROD2');
            vname = ['X',int2str(I),',',int2str(round(v_('I-1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('VI',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('I-1'))),',',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('VJ',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['D',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQ';
            ielftype(ie) = iet_('eSQ');
            vname = ['X',int2str(I),',',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['G',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'ePROD2';
            ielftype(ie) = iet_('ePROD2');
            vname = ['X',int2str(I),',',int2str(round(v_('I+1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('VI',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('I+1'))),',',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('VJ',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['H',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'ePROD2';
            ielftype(ie) = iet_('ePROD2');
            vname = ['X',int2str(I),',',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('VI',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(I),',',int2str(round(v_('I+1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('VJ',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['R',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'ePROD2';
            ielftype(ie) = iet_('ePROD2');
            vname = ['X',int2str(I),',',int2str(round(v_('I+1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('VI',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('I+1'))),',',int2str(round(v_('I+1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('VJ',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['S',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'ePROD2';
            ielftype(ie) = iet_('ePROD2');
            vname = ['X',int2str(I),',',int2str(round(v_('I+1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('VI',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('I+1'))),',',int2str(round(v_('I+2')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('VJ',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        ename = ['A',int2str(round(v_('M-1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD2';
        ielftype(ie) = iet_('ePROD2');
        ename = ['A',int2str(round(v_('M-1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('M-1'))),',',int2str(round(v_('M-2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VI',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['A',int2str(round(v_('M-1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('M-2'))),',',int2str(round(v_('M-3')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VJ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['B',int2str(round(v_('M-1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD2';
        ielftype(ie) = iet_('ePROD2');
        ename = ['B',int2str(round(v_('M-1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('M-1'))),',',int2str(round(v_('M-2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VI',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['B',int2str(round(v_('M-1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('M-2'))),',',int2str(round(v_('M-2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VJ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C',int2str(round(v_('M-1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD2';
        ielftype(ie) = iet_('ePROD2');
        ename = ['C',int2str(round(v_('M-1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('M-1'))),',',int2str(round(v_('M-1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VI',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C',int2str(round(v_('M-1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('M-1'))),',',int2str(round(v_('M-2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VJ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['F',int2str(round(v_('M-1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD2';
        ielftype(ie) = iet_('ePROD2');
        ename = ['F',int2str(round(v_('M-1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('M-1'))),',',int2str(round(v_('M-2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VI',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['F',int2str(round(v_('M-1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('M-2'))),',',int2str(round(v_('M-1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VJ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['D',int2str(round(v_('M-1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQ';
        ielftype(ie) = iet_('eSQ');
        ename = ['D',int2str(round(v_('M-1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('M-1'))),',',int2str(round(v_('M-1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['G',int2str(round(v_('M-1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD2';
        ielftype(ie) = iet_('ePROD2');
        ename = ['G',int2str(round(v_('M-1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('M-1'))),',',int2str(round(v_('M')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VI',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['G',int2str(round(v_('M-1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('M'))),',',int2str(round(v_('M-1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VJ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['H',int2str(round(v_('M-1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD2';
        ielftype(ie) = iet_('ePROD2');
        ename = ['H',int2str(round(v_('M-1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('M-1'))),',',int2str(round(v_('M-1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VI',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['H',int2str(round(v_('M-1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('M-1'))),',',int2str(round(v_('M')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VJ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['R',int2str(round(v_('M-1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD2';
        ielftype(ie) = iet_('ePROD2');
        ename = ['R',int2str(round(v_('M-1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('M-1'))),',',int2str(round(v_('M')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VI',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['R',int2str(round(v_('M-1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('M'))),',',int2str(round(v_('M')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VJ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['A',int2str(round(v_('M')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD2';
        ielftype(ie) = iet_('ePROD2');
        ename = ['A',int2str(round(v_('M')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('M'))),',',int2str(round(v_('M-1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VI',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['A',int2str(round(v_('M')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('M-1'))),',',int2str(round(v_('M-2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VJ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['B',int2str(round(v_('M')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD2';
        ielftype(ie) = iet_('ePROD2');
        ename = ['B',int2str(round(v_('M')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('M'))),',',int2str(round(v_('M-1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VI',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['B',int2str(round(v_('M')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('M-1'))),',',int2str(round(v_('M-1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VJ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C',int2str(round(v_('M')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD2';
        ielftype(ie) = iet_('ePROD2');
        ename = ['C',int2str(round(v_('M')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('M'))),',',int2str(round(v_('M')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VI',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C',int2str(round(v_('M')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('M'))),',',int2str(round(v_('M-1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VJ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['F',int2str(round(v_('M')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD2';
        ielftype(ie) = iet_('ePROD2');
        ename = ['F',int2str(round(v_('M')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('M'))),',',int2str(round(v_('M-1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VI',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['F',int2str(round(v_('M')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('M-1'))),',',int2str(round(v_('M')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('VJ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['D',int2str(round(v_('M')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQ';
        ielftype(ie) = iet_('eSQ');
        ename = ['D',int2str(round(v_('M')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('M'))),',',int2str(round(v_('M')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for ig = 1:ngrp
            pbm.grftype{ig} = 'gL2';
        end
        ig = ig_(['E',int2str(round(v_('1'))),',',int2str(round(v_('1')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['D',int2str(round(v_('1')))]);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['G',int2str(round(v_('1')))]);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_(['E',int2str(round(v_('1'))),',',int2str(round(v_('2')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['H',int2str(round(v_('1')))]);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['R',int2str(round(v_('1')))]);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_(['E',int2str(round(v_('1'))),',',int2str(round(v_('3')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['S',int2str(round(v_('1')))]);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_(['E',int2str(round(v_('2'))),',',int2str(round(v_('1')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['B',int2str(round(v_('2')))]);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['C',int2str(round(v_('2')))]);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_(['E',int2str(round(v_('2'))),',',int2str(round(v_('2')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['F',int2str(round(v_('2')))]);
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['D',int2str(round(v_('2')))]);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['G',int2str(round(v_('2')))]);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_(['E',int2str(round(v_('2'))),',',int2str(round(v_('3')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['H',int2str(round(v_('2')))]);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['R',int2str(round(v_('2')))]);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_(['E',int2str(round(v_('2'))),',',int2str(round(v_('4')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['S',int2str(round(v_('2')))]);
        pbm.grelw{ig}(posel) = 1.;
        for I=v_('3'):v_('M-2')
            v_('I-2') = -2+I;
            v_('I-1') = -1+I;
            v_('I+1') = 1+I;
            v_('I+2') = 2+I;
            ig = ig_(['E',int2str(I),',',int2str(round(v_('I-2')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['A',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
            ig = ig_(['E',int2str(I),',',int2str(round(v_('I-1')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['B',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['C',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
            ig = ig_(['E',int2str(I),',',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['F',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['D',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['G',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
            ig = ig_(['E',int2str(I),',',int2str(round(v_('I+1')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['H',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['R',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
            ig = ig_(['E',int2str(I),',',int2str(round(v_('I+2')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['S',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
        end
        ig = ig_(['E',int2str(round(v_('M-1'))),',',int2str(round(v_('M-3')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['A',int2str(round(v_('M-1')))]);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_(['E',int2str(round(v_('M-1'))),',',int2str(round(v_('M-2')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['B',int2str(round(v_('M-1')))]);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['C',int2str(round(v_('M-1')))]);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_(['E',int2str(round(v_('M-1'))),',',int2str(round(v_('M-1')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['F',int2str(round(v_('M-1')))]);
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['D',int2str(round(v_('M-1')))]);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['G',int2str(round(v_('M-1')))]);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_(['E',int2str(round(v_('M-1'))),',',int2str(round(v_('M')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['H',int2str(round(v_('M-1')))]);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['R',int2str(round(v_('M-1')))]);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_(['E',int2str(round(v_('M'))),',',int2str(round(v_('M-2')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['A',int2str(round(v_('M')))]);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_(['E',int2str(round(v_('M'))),',',int2str(round(v_('M-1')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['B',int2str(round(v_('M')))]);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['C',int2str(round(v_('M')))]);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_(['E',int2str(round(v_('M'))),',',int2str(round(v_('M')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['F',int2str(round(v_('M')))]);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['D',int2str(round(v_('M')))]);
        pbm.grelw{ig}(posel) = 1.;
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'C-CSUR2-AN-V-V';
        pbm.objderlvl = 2;
        pb.objderlvl = pbm.objderlvl;
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb, pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm,pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end


    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'ePROD2'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(2);
        if(nargout>1)
            g_(1,1) = EV_(2);
            g_(2,1) = EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = 1.0;
                H_(2,1) = H_(1,2);
                varargout{3} = H_;
            end
        end

    case 'eSQ'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(1);
        if(nargout>1)
            g_(1,1) = EV_(1)+EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0;
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gL2'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = GVAR_*GVAR_;
        if(nargout>1)
            g_ = GVAR_+GVAR_;
            varargout{2} = g_;
            if(nargout>2)
                H_ = 2.0;
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    case {'fx','fgx','fgHx','cx','cJx','cJHx','cIx','cIJx','cIJHx','cIJxv','fHxv',...
          'cJxv','cJtxv','cIJtxv','Lxy','Lgxy','LgHxy','LIxy','LIgxy','LIgHxy',...
          'LHxyv','LIHxyv'}

        if(isfield(pbm,'name')&&strcmp(pbm.name,name))
            pbm.has_globs = [0,0];
            [varargout{1:max(1,nargout)}] = s2mpjlib(action,pbm,varargin{:});
        else
            disp(['ERROR: please run ',name,' with action = setup'])
            [varargout{1:nargout}] = deal(NaN);
        end

    otherwise
        disp([' ERROR: action ',action,' unavailable for problem ',name,'.m'])
    end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

