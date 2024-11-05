function varargout = MINC44(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : MINC44
%    *********
% 
%    Minimize the permanent of a doubly stochastic n dimensional square matrix
%    whose trace is zero.
%    The conjecture is that the minimum is achieved when all non-diagonal
%    entries of the matrix are equal to 1/(n-1).
% 
%    Source: conjecture 44 in
%    H. Minc,
%    "Theory of Permanents 1982-1985",
%    Linear Algebra and Applications, vol. 21, pp. 109-148, 1987.
% 
%    SIF input: Ph. Toint, April 1992.
%               minor correction by Ph. Shott, Jan 1995.
% 
%    classification = 'C-CLQR2-AN-V-V'
% 
%    Size of matrix
% 
%       Alternative values for the SIF file parameters:
% IE N                   2              $-PARAMETER n = 5
% IE N                   3              $-PARAMETER n = 13
% IE N                   4              $-PARAMETER n = 27
% IE N                   5              $-PARAMETER n = 51
% IE N                   6              $-PARAMETER n = 93
% IE N                   7              $-PARAMETER n = 169
% IE N                   8              $-PARAMETER n = 311
% IE N                   9              $-PARAMETER n = 583
% IE N                   10             $-PARAMETER n = 1113
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 17 X 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'MINC44';

switch(action)

    case {'setup','setup_redprec'}

        if(isfield(pbm,'ndigs'))
            rmfield(pbm,'ndigs');
        end
        if(strcmp(action,'setup_redprec'))
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
            v_('N') = 5;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
        v_('0') = 0;
        v_('1') = 1;
        v_('2') = 2;
        v_('N+1') = 1+v_('N');
        v_('N-1') = -1+v_('N');
        v_('2**N') = 1;
        for I=v_('1'):v_('N')
            v_('N-I+1') = v_('N+1')-I;
            v_('I-1') = - 1+I;
            v_('R2**N') = v_('2**N');
            v_(['S',int2str(round(v_('N-I+1')))]) = 0.1+v_('R2**N');
            v_(['T',int2str(round(v_('I-1')))]) = 0.1+v_('R2**N');
            v_('2**N') = v_('2**N')*v_('2');
        end
        v_('2**N-1') = - 1+v_('2**N');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for M=v_('1'):v_('N-1')
            v_('RK1') = v_(['T',int2str(M)]);
            v_('K1') = fix(v_('RK1'));
            v_('K2') = 2*v_('K1');
            v_('K1') = 1+v_('K1');
            v_('K2') = - 1+v_('K2');
            for K=v_('K1'):v_('K2')
                [iv,ix_] = s2mpjlib('ii',['P',int2str(K)],ix_);
                pb.xnames{iv} = ['P',int2str(K)];
            end
        end
        for I=v_('1'):v_('N')
            for J=v_('1'):v_('N')
                [iv,ix_] = s2mpjlib('ii',['A',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['A',int2str(I),',',int2str(J)];
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        iv = ix_(['P',int2str(round(v_('2**N-1')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        for M=v_('1'):v_('N-1')
            v_('RK1') = v_(['T',int2str(M)]);
            v_('K1') = fix(v_('RK1'));
            v_('K2') = 2*v_('K1');
            v_('K1') = 1+v_('K1');
            v_('K2') = - 1+v_('K2');
            for K=v_('K1'):v_('K2')
                [ig,ig_] = s2mpjlib('ii',['PE',int2str(K)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['PE',int2str(K)];
                iv = ix_(['P',int2str(K)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = - 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = - 1.0;
                end
            end
        end
        for I=v_('1'):v_('N')
            for J=v_('1'):v_('N')
                [ig,ig_] = s2mpjlib('ii',['C',int2str(J)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['C',int2str(J)];
                iv = ix_(['A',int2str(I),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
            end
        end
        for I=v_('1'):v_('N-1')
            for J=v_('1'):v_('N')
                [ig,ig_] = s2mpjlib('ii',['R',int2str(I)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['R',int2str(I)];
                iv = ix_(['A',int2str(I),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
            end
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
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        for J=v_('1'):v_('N')
            pbm.gconst(ig_(['C',int2str(J)])) = 1.0;
        end
        for I=v_('1'):v_('N-1')
            pbm.gconst(ig_(['R',int2str(I)])) = 1.0;
        end
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        for I=v_('1'):v_('N')
            for J=v_('1'):v_('N')
                pb.xupper(ix_(['A',int2str(I),',',int2str(J)])) = 1.0;
            end
        end
        for I=v_('1'):v_('N')
            pb.xlower(ix_(['A',int2str(I),',',int2str(I)]),1) = 0.0;
            pb.xupper(ix_(['A',int2str(I),',',int2str(I)]),1) = 0.0;
        end
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'en2PR',iet_);
        elftv{it}{1} = 'A';
        elftv{it}{2} = 'P';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for M=v_('1'):v_('N-1')
            v_('RK1') = v_(['T',int2str(M)]);
            v_('K1') = fix(v_('RK1'));
            v_('K2') = 2*v_('K1');
            v_('K1') = 1+v_('K1');
            v_('K2') = - 1+v_('K2');
            for K=v_('K1'):v_('K2')
                v_('ID') = 0;
                v_('PT') = 1;
                v_('KK') = K;
                for I=v_('1'):v_('N')
                    v_('SI') = v_(['S',int2str(I)]);
                    v_('ISI') = fix(v_('SI'));
                    v_('BI') = fix(v_('KK')/v_('ISI'));
                    v_('ID') = v_('ID')+v_('BI');
                    v_('BISI') = v_('BI')*v_('ISI');
                    v_('KK') = v_('KK')-v_('BISI');
                    v_('RI') = I;
                    v_(['RNZ',int2str(round(v_('PT')))]) = 0.1+v_('RI');
                    v_('PT') = v_('PT')+v_('BI');
                end
                v_('I1') = v_('0');
                v_('I2') = v_('1');
                v_('ID-2') = - 2+v_('ID');
                for I=v_('1'):v_('ID-2')
                    v_('I1') = v_('ID');
                    v_('I2') = v_('0');
                end
                for I=v_('1'):v_('I1')
                    v_('RJ') = v_(['RNZ',int2str(I)]);
                    v_('J') = fix(v_('RJ'));
                    v_('SI') = v_(['S',int2str(round(v_('J')))]);
                    v_('ISI') = fix(v_('SI'));
                    v_('IPP') = K-v_('ISI');
                    ename = ['E',int2str(K),',',int2str(I)];
                    [ie,ie_] = s2mpjlib('ii',ename,ie_);
                    pbm.elftype{ie} = 'en2PR';
                    ielftype(ie) = iet_('en2PR');
                    vname = ['A',int2str(round(v_('ID'))),',',int2str(round(v_('J')))];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('A',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                    vname = ['P',int2str(round(v_('IPP')))];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('P',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                end
                for I=v_('1'):v_('I2')
                    v_('RJ') = v_(['RNZ',int2str(round(v_('1')))]);
                    v_('J') = fix(v_('RJ'));
                    v_('RJJ') = v_(['RNZ',int2str(round(v_('2')))]);
                    v_('JJ') = fix(v_('RJJ'));
                    ename = ['E',int2str(K),',',int2str(round(v_('1')))];
                    [ie,ie_] = s2mpjlib('ii',ename,ie_);
                    pbm.elftype{ie} = 'en2PR';
                    ielftype(ie) = iet_('en2PR');
                    ename = ['E',int2str(K),',',int2str(round(v_('1')))];
                    [ie,ie_] = s2mpjlib('ii',ename,ie_);
                    vname = ['A',int2str(round(v_('2'))),',',int2str(round(v_('J')))];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('A',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                    ename = ['E',int2str(K),',',int2str(round(v_('1')))];
                    [ie,ie_] = s2mpjlib('ii',ename,ie_);
                    vname = ['A',int2str(round(v_('1'))),',',int2str(round(v_('JJ')))];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('P',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                    ename = ['E',int2str(K),',',int2str(round(v_('2')))];
                    [ie,ie_] = s2mpjlib('ii',ename,ie_);
                    pbm.elftype{ie} = 'en2PR';
                    ielftype(ie) = iet_('en2PR');
                    ename = ['E',int2str(K),',',int2str(round(v_('2')))];
                    [ie,ie_] = s2mpjlib('ii',ename,ie_);
                    vname = ['A',int2str(round(v_('2'))),',',int2str(round(v_('JJ')))];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('A',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                    ename = ['E',int2str(K),',',int2str(round(v_('2')))];
                    [ie,ie_] = s2mpjlib('ii',ename,ie_);
                    vname = ['A',int2str(round(v_('1'))),',',int2str(round(v_('J')))];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('P',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                end
                v_('RD') = v_('ID');
                v_(['D',int2str(K)]) = 0.1+v_('RD');
            end
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for M=v_('1'):v_('N-1')
            v_('RK1') = v_(['T',int2str(M)]);
            v_('K1') = fix(v_('RK1'));
            v_('K2') = 2*v_('K1');
            v_('K1') = 1+v_('K1');
            v_('K2') = - 1+v_('K2');
            for K=v_('K1'):v_('K2')
                v_('RD') = v_(['D',int2str(K)]);
                v_('ID') = fix(v_('RD'));
                for I=v_('1'):v_('ID')
                    ig = ig_(['PE',int2str(K)]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) = ie_(['E',int2str(K),',',int2str(I)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = 1.;
                end
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
%    Solution
% LO SOLTN(2)           1.0
% LO SOLTN(3)           0.25
% LO SOLTN(4)           0.11111111
% LO SOLTN(5)           0.04296835
% LO SOLTN(6)           0.01695926
% LO SOLTN(7)           6.62293832D-03
% LO SOLTN(8)           2.57309338D-03
% LO SOLTN(9)           9.94617795D-04
% LO SOLTN(10)          3.83144655D-04
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CLQR2-AN-V-V';
        pb.x0          = zeros(pb.n,1);
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

    case 'en2PR'

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

