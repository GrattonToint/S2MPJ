function varargout = DISCS(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : DISCS
%    *********
% 
%    The problem is to replace the nodes of a planar graph by discs
%    such that adjacent nodes have touching discs and non adjacent nodes
%    disjoint discs.  The smallest disc is constrained to have radius equal to
%    one and is centered at the origin.  The second discs is also
%    constrained to have its Y = 0.0, in order to cancel invariance under
%    rotation.
% 
%    Source:
%    W. Pulleyblank,
%    private communication, 1991.
% 
%    classification = 'C-CLQR2-MY-36-66'
% 
%    SIF input: A.R. Conn and Ph. Toint, April 1991.
% 
%    Number of nodes
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 21 VI 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'DISCS';

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
        v_('NNODES') = 12;
        v_('EPSIL') = 0.0001;
        v_('1') = 1;
        v_('2') = 2;
        for I=v_('1'):v_('NNODES')
            v_('I-1') = -1+I;
            for J=v_('1'):v_('I-1')
                v_(['A',int2str(J),',',int2str(I)]) = 0.0;
            end
        end
        v_('A1,2') = 1.01;
        v_('A1,7') = 1.01;
        v_('A2,3') = 1.01;
        v_('A2,4') = 1.01;
        v_('A3,4') = 1.01;
        v_('A4,5') = 1.01;
        v_('A4,6') = 1.01;
        v_('A5,6') = 1.01;
        v_('A5,11') = 1.01;
        v_('A6,7') = 1.01;
        v_('A7,8') = 1.01;
        v_('A7,9') = 1.01;
        v_('A8,9') = 1.01;
        v_('A8,10') = 1.01;
        v_('A9,10') = 1.01;
        v_('A10,11') = 1.01;
        v_('A10,12') = 1.01;
        v_('A11,12') = 1.01;
        v_('RNODES') = v_('NNODES');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for I=v_('1'):v_('NNODES')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
            [iv,ix_] = s2mpjlib('ii',['Y',int2str(I)],ix_);
            pb.xnames{iv} = ['Y',int2str(I)];
            [iv,ix_] = s2mpjlib('ii',['R',int2str(I)],ix_);
            pb.xnames{iv} = ['R',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        for I=v_('1'):v_('NNODES')
            [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
            gtype{ig} = '<>';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['R',int2str(I)]);
            valA(end+1) = 1.0;
        end
        for I=v_('2'):v_('NNODES')
            v_('I-1') = -1+I;
            for J=v_('1'):v_('I-1')
                v_('RAIJ') = v_(['A',int2str(J),',',int2str(I)]);
                v_('AIJ') = fix(v_('RAIJ'));
                for K=v_('1'):v_('AIJ')
                    [ig,ig_] = s2mpjlib('ii',['S',int2str(J),',',int2str(I)],ig_);
                    gtype{ig}  = '==';
                    cnames{ig} = ['S',int2str(J),',',int2str(I)];
                end
                v_('AIJ-1') = -1+v_('AIJ');
                v_('NAIJ') = -1*v_('AIJ-1');
                for K=v_('1'):v_('NAIJ')
                    [ig,ig_] = s2mpjlib('ii',['S',int2str(J),',',int2str(I)],ig_);
                    gtype{ig}  = '<=';
                    cnames{ig} = ['S',int2str(J),',',int2str(I)];
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
        v_('-EPSIL') = -1.0*v_('EPSIL');
        for I=v_('1'):v_('NNODES')
            v_('I-1') = -1+I;
            for J=v_('1'):v_('I-1')
                v_('RAIJ') = v_(['A',int2str(J),',',int2str(I)]);
                v_('AIJ') = fix(v_('RAIJ'));
                v_('AIJ-1') = -1+v_('AIJ');
                v_('NAIJ') = -1*v_('AIJ-1');
                for K=v_('1'):v_('NAIJ')
                    pbm.gconst(ig_(['S',int2str(J),',',int2str(I)])) = v_('-EPSIL');
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        for I=v_('1'):v_('NNODES')
            pb.xlower(ix_(['R',int2str(I)]),1) = 1.0;
        end
        pb.xlower(ix_(['X',int2str(round(v_('1')))]),1) = 0.0;
        pb.xupper(ix_(['X',int2str(round(v_('1')))]),1) = 0.0;
        pb.xlower(ix_(['Y',int2str(round(v_('1')))]),1) = 0.0;
        pb.xupper(ix_(['Y',int2str(round(v_('1')))]),1) = 0.0;
        pb.xlower(ix_(['Y',int2str(round(v_('2')))]),1) = 0.0;
        pb.xupper(ix_(['Y',int2str(round(v_('2')))]),1) = 0.0;
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        for I=v_('1'):v_('NNODES')
            v_('RI') = I;
            v_('RM') = 0.03*v_('RI');
            v_('RP') = 0.0055555*v_('RI');
            pb.x0(ix_(['R',int2str(I)]),1) = v_('RI');
            pb.x0(ix_(['X',int2str(I)]),1) = v_('RM');
            pb.x0(ix_(['Y',int2str(I)]),1) = v_('RP');
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eIPSQ',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        [it,iet_] = s2mpjlib( 'ii', 'eIMSQ',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('2'):v_('NNODES')
            v_('I-1') = -1+I;
            for J=v_('1'):v_('I-1')
                ename = ['DR',int2str(J),',',int2str(I)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eIPSQ';
                ielftype(ie) = iet_('eIPSQ');
                vname = ['R',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('X',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['R',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('Y',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['DX',int2str(J),',',int2str(I)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eIMSQ';
                ielftype(ie) = iet_('eIMSQ');
                vname = ['X',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('X',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('Y',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['DY',int2str(J),',',int2str(I)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eIMSQ';
                ielftype(ie) = iet_('eIMSQ');
                vname = ['Y',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('X',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['Y',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('Y',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('2'):v_('NNODES')
            v_('I-1') = -1+I;
            for J=v_('1'):v_('I-1')
                ig = ig_(['S',int2str(J),',',int2str(I)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['DR',int2str(J),',',int2str(I)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['DX',int2str(J),',',int2str(I)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = -1.0;
                posel = posel+1;
                pbm.grelt{ig}(posel) = ie_(['DY',int2str(J),',',int2str(I)]);
                pbm.grelw{ig}(posel) = -1.0;
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
% ZL DISCS                              RNODES
%    Solution
% LO SOLTN(12)           20.46122911
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CLQR2-MY-36-66';
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

    case 'eIPSQ'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(1,2);
        U_(1,1) = U_(1,1)+1;
        U_(1,2) = U_(1,2)+1;
        IV_(1) = U_(1,:)*EV_;
        varargout{1} = IV_(1)*IV_(1);
        if(nargout>1)
            g_(1,1) = IV_(1)+IV_(1);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_(1,1) = 2.0;
                varargout{3} = U_.'*H_*U_;
            end
        end

    case 'eIMSQ'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(1,2);
        U_(1,1) = U_(1,1)+1;
        U_(1,2) = U_(1,2)-1;
        IV_(1) = U_(1,:)*EV_;
        varargout{1} = IV_(1)*IV_(1);
        if(nargout>1)
            g_(1,1) = IV_(1)+IV_(1);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_(1,1) = 2.0;
                varargout{3} = U_.'*H_*U_;
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

