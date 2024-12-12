function varargout = HIMMELBFNE(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HIMMELBFNE
%    *********
% 
%    A 4 variables data fitting problems by Himmelblau.
% 
%    Source: problem 32 in
%    D.H. Himmelblau,
%    "Applied nonlinear programming",
%    McGraw-Hill, New-York, 1972.
% 
%    See Buckley#76 (p. 66)
% 
%    SIF input: Ph. Toint, Dec 1989.
%    Nonlinear-equations version of HIMMELBF.SIF, Nick Gould, Jan 2020.
% 
%    classification = 'C-CNOR2-AN-4-7'
% 
%    Problem data
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HIMMELBFNE';

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
        v_('A1') = 0.0;
        v_('A2') = 0.000428;
        v_('A3') = 0.001000;
        v_('A4') = 0.001610;
        v_('A5') = 0.002090;
        v_('A6') = 0.003480;
        v_('A7') = 0.005250;
        v_('B1') = 7.391;
        v_('B2') = 11.18;
        v_('B3') = 16.44;
        v_('B4') = 16.20;
        v_('B5') = 22.20;
        v_('B6') = 24.02;
        v_('B7') = 31.32;
        v_('1') = 1;
        v_('7') = 7;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        [iv,ix_] = s2mpjlib('ii','X1',ix_);
        pb.xnames{iv} = 'X1';
        [iv,ix_] = s2mpjlib('ii','X2',ix_);
        pb.xnames{iv} = 'X2';
        [iv,ix_] = s2mpjlib('ii','X3',ix_);
        pb.xnames{iv} = 'X3';
        [iv,ix_] = s2mpjlib('ii','X4',ix_);
        pb.xnames{iv} = 'X4';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        for I=v_('1'):v_('7')
            [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['G',int2str(I)];
            pbm.gscale(ig,1) = 0.0001;
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
        %%%%%%%%%%%%%%%%%%%  CONSTANTS %%%%%%%%%%%%%%%%%%%
        pbm.gconst = 1.0*ones(ngrp,1);
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,'X1'))
            pb.x0(ix_('X1'),1) = 2.7;
        else
            pb.y0(find(pbm.congrps==ig_('X1')),1) = 2.7;
        end
        if(isKey(ix_,'X2'))
            pb.x0(ix_('X2'),1) = 90.0;
        else
            pb.y0(find(pbm.congrps==ig_('X2')),1) = 90.0;
        end
        if(isKey(ix_,'X3'))
            pb.x0(ix_('X3'),1) = 1500.0;
        else
            pb.y0(find(pbm.congrps==ig_('X3')),1) = 1500.0;
        end
        if(isKey(ix_,'X4'))
            pb.x0(ix_('X4'),1) = 10.0;
        else
            pb.y0(find(pbm.congrps==ig_('X4')),1) = 10.0;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eHF',iet_);
        elftv{it}{1} = 'XA';
        elftv{it}{2} = 'XB';
        elftv{it}{3} = 'XC';
        elftv{it}{4} = 'XD';
        elftp{it}{1} = 'A';
        elftp{it}{2} = 'B';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('1'):v_('7')
            ename = ['E',int2str(I)];
            [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
            if(newelt)
                pbm.elftype{ie} = 'eHF';
                ielftype(ie) = iet_('eHF');
            end
            vname = 'X1';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('XA',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X2';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('XB',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X3';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('XC',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X4';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('XD',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('A',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['A',int2str(I)]);
            [~,posep] = ismember('B',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['B',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('7')
            ig = ig_(['G',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               318.572
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CNOR2-AN-4-7';
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

    case 'eHF'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U = EV_(1)*EV_(1)+pbm.elpar{iel_}(1)*EV_(2)*EV_(2)+pbm.elpar{iel_}(1)*...
             pbm.elpar{iel_}(1)*EV_(3)*EV_(3);
        V = pbm.elpar{iel_}(2)*(1.0+pbm.elpar{iel_}(1)*EV_(4)*EV_(4));
        V2 = V*V;
        AB = pbm.elpar{iel_}(1)*pbm.elpar{iel_}(2);
        A2 = pbm.elpar{iel_}(1)*pbm.elpar{iel_}(1);
        T = -4.0*AB/V2;
        varargout{1} = U/V;
        if(nargout>1)
            g_(1,1) = 2.0*EV_(1)/V;
            g_(2,1) = 2.0*pbm.elpar{iel_}(1)*EV_(2)/V;
            g_(3,1) = 2.0*A2*EV_(3)/V;
            g_(4,1) = -2.0*AB*EV_(4)*U/V2;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(4,4);
                H_(1,1) = 2.0/V;
                H_(1,4) = T*EV_(4)*EV_(1);
                H_(4,1) = H_(1,4);
                H_(2,2) = 2.0*pbm.elpar{iel_}(1)/V;
                H_(2,4) = T*pbm.elpar{iel_}(1)*EV_(4)*EV_(2);
                H_(4,2) = H_(2,4);
                H_(3,3) = 2.0*A2/V;
                H_(3,4) = T*pbm.elpar{iel_}(1)*EV_(4)*EV_(3);
                H_(4,3) = H_(3,4);
                H_(4,4) = -2.0*AB*U/V2+8.0*(AB*EV_(4))^2*U/(V2*V);
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

