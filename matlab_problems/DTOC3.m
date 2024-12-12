function varargout = DTOC3(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : DTOC3
%    *********
% 
%    This is a discrete time optimal control (DTOC) problem.  
%    The system has N time periods, 1 control variable and 2 state variables.
% 
%    The problem is convex.
% 
%    Sources: problem 3 in
%    T.F. Coleman and A. Liao,
%    "An Efficient Trust Region Method for Unconstrained Discret-Time Optimal
%    Control Problems",
%    Tech. Report, ctc93tr144,  Advanced Computing Research Institute, 
%    Cornell University, 1992.
% 
%    D.P. Bertsekas,
%    "Projected Newton methods for optimization problems with simple
%    constraints", 
%    SIAM J. Control and Optimization 20, pp. 221-246, 1982.
% 
%    SIF input: Ph. Toint, August 1993
% 
%    classification = 'C-CQLR2-AN-V-V'
% 
%    Problem variants: they are identified by the value of the parameter N.
% 
%    The problem has 3N-1  variables (of which 2 are fixed),
%    and 2(N-1) constraints
% 
%       Alternative values for the SIF file parameters:
% IE N                   10             $-PARAMETER  n=   29,m= 18 original value
% IE N                   50             $-PARAMETER  n=  149,m= 98
% IE N                   100            $-PARAMETER  n=  299,m=198
% IE N                   500            $-PARAMETER  n= 1499,m=998
% IE N                   1000           $-PARAMETER  n= 2999,m=1998
% IE N                   1500           $-PARAMETER  n= 4499,m=2998
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'DTOC3';

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
            v_('N') = 10;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
% IE N                   5000           $-PARAMETER  n=14999,m=9998
        v_('N-1') = -1+v_('N');
        v_('1') = 1;
        v_('2') = 2;
        v_('RN') = v_('N');
        v_('S') = 1.0/v_('RN');
        v_('2/S') = 2.0/v_('S');
        v_('-S') = -1.0*v_('S');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for T=v_('1'):v_('N-1')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(T)],ix_);
            pb.xnames{iv} = ['X',int2str(T)];
        end
        for T=v_('1'):v_('N')
            for I=v_('1'):v_('2')
                [iv,ix_] = s2mpjlib('ii',['Y',int2str(T),',',int2str(I)],ix_);
                pb.xnames{iv} = ['Y',int2str(T),',',int2str(I)];
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        for T=v_('1'):v_('N-1')
            [ig,ig_] = s2mpjlib('ii',['O',int2str(T)],ig_);
            gtype{ig} = '<>';
            pbm.gscale(ig,1) = v_('2/S');
        end
        for T=v_('1'):v_('N-1')
            v_('T+1') = 1+T;
            [ig,ig_] = s2mpjlib('ii',['TT',int2str(T),',',int2str(round(v_('1')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['TT',int2str(T),',',int2str(round(v_('1')))];
            irA(end+1)  = ig;
            icA(end+1)  =...
                  ix_(['Y',int2str(round(v_('T+1'))),',',int2str(round(v_('1')))]);
            valA(end+1) = -1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['Y',int2str(T),',',int2str(round(v_('1')))]);
            valA(end+1) = 1.0;
            [ig,ig_] = s2mpjlib('ii',['TT',int2str(T),',',int2str(round(v_('1')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['TT',int2str(T),',',int2str(round(v_('1')))];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['Y',int2str(T),',',int2str(round(v_('2')))]);
            valA(end+1) = v_('S');
            [ig,ig_] = s2mpjlib('ii',['TT',int2str(T),',',int2str(round(v_('2')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['TT',int2str(T),',',int2str(round(v_('2')))];
            irA(end+1)  = ig;
            icA(end+1)  =...
                  ix_(['Y',int2str(round(v_('T+1'))),',',int2str(round(v_('2')))]);
            valA(end+1) = -1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['Y',int2str(T),',',int2str(round(v_('2')))]);
            valA(end+1) = 1.0;
            [ig,ig_] = s2mpjlib('ii',['TT',int2str(T),',',int2str(round(v_('2')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['TT',int2str(T),',',int2str(round(v_('2')))];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['Y',int2str(T),',',int2str(round(v_('1')))]);
            valA(end+1) = v_('-S');
            [ig,ig_] = s2mpjlib('ii',['TT',int2str(T),',',int2str(round(v_('2')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['TT',int2str(T),',',int2str(round(v_('2')))];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(T)]);
            valA(end+1) = v_('S');
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = ix_.Count;
        ngrp   = ig_.Count;
        legrps = find(strcmp(gtype,'<='));
        eqgrps = find(strcmp(gtype,'=='));
        gegrps = find(strcmp(gtype,'>='));
        pb.nle = length(legrps);
        pb.neq = length(eqgrps);
        pb.nge = length(gegrps);
        pb.m   = pb.nle+pb.neq+pb.nge;
        pbm.congrps = [ legrps, eqgrps, gegrps ];
        [pb.cnames{1:pb.m}] = deal(cnames{pbm.congrps});
        pb.nob = ngrp-pb.m;
        pbm.objgrps = find(strcmp(gtype,'<>'));
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower(ix_(['Y',int2str(round(v_('1'))),',',int2str(round(v_('1')))]),1) = 15.0;
        pb.xupper(ix_(['Y',int2str(round(v_('1'))),',',int2str(round(v_('1')))]),1) = 15.0;
        pb.xlower(ix_(['Y',int2str(round(v_('1'))),',',int2str(round(v_('2')))]),1) = 5.0;
        pb.xupper(ix_(['Y',int2str(round(v_('1'))),',',int2str(round(v_('2')))]),1) = 5.0;
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        pb.x0(ix_(['Y',int2str(round(v_('1'))),',',int2str(round(v_('1')))]),1) =...
              15.0;
        pb.x0(ix_(['Y',int2str(round(v_('1'))),',',int2str(round(v_('2')))]),1) =...
              5.0;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSQ',iet_);
        elftv{it}{1} = 'YY';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for T=v_('2'):v_('N')
            ename = ['Y1SQ',int2str(T)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQ';
            ielftype(ie) = iet_('eSQ');
            vname = ['Y',int2str(T),',',int2str(round(v_('1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('YY',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['Y2SQ',int2str(T)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQ';
            ielftype(ie) = iet_('eSQ');
            vname = ['Y',int2str(T),',',int2str(round(v_('2')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('YY',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        for T=v_('1'):v_('N-1')
            ename = ['XSQ',int2str(T)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQ';
            ielftype(ie) = iet_('eSQ');
            vname = ['X',int2str(T)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('YY',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for T=v_('1'):v_('N-1')
            v_('T+1') = 1+T;
            ig = ig_(['O',int2str(T)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['Y1SQ',int2str(round(v_('T+1')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 2.0;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['Y2SQ',int2str(round(v_('T+1')))]);
            pbm.grelw{ig}(posel) = 1.0;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['XSQ',int2str(T)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 6.0;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
% LO SOLUTION(  10)      224.590381002
% LO SOLUTION(  50)      233.278523083
% LO SOLUTION( 100)      234.286202920
% LO SOLUTION( 500)      235.084407947
% LO SOLUTION(1000)      235.182824435
% LO SOLUTION(5000)      235.154640099
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CQLR2-AN-V-V';
        pbm.objderlvl = 2;
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2];
        pb.conderlvl  = pbm.conderlvl;
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb, pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm,pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end

% **********************
%  SET UP THE FUNCTION *
%  AND RANGE ROUTINES  *
% **********************

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

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

