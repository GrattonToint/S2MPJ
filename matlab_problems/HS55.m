function varargout = HS55(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HS55
%    *********
% 
%    Source: problem 55 in
%    W. Hock and K. Schittkowski,
%    "Test examples for nonlinear programming codes",
%    Lectures Notes in Economics and Mathematical Systems 187, Springer
%    Verlag, Heidelberg, 1981.
% 
%    SIF input: A.R. Conn, April 1990
% 
%    classification = 'C-COLR2-AN-6-6'
% 
%    some useful parameters, including N, the number of variables.
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 21 VI 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HS55';

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
        v_('N') = 6;
        v_('1') = 1;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2');
        valA(end+1) = 2.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X5');
        valA(end+1) = 4.0;
        [ig,ig_] = s2mpjlib('ii','CON1',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'CON1';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2');
        valA(end+1) = 2.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X5');
        valA(end+1) = 5.0;
        [ig,ig_] = s2mpjlib('ii','CON2',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'CON2';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X3');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','CON3',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'CON3';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X4');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X5');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X6');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','CON4',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'CON4';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X4');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','CON5',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'CON5';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X5');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','CON6',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'CON6';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X3');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X6');
        valA(end+1) = 1.0;
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
        pbm.gconst(ig_('CON1')) = 6.0;
        pbm.gconst(ig_('CON2')) = 3.0;
        pbm.gconst(ig_('CON3')) = 2.0;
        pbm.gconst(ig_('CON4')) = 1.0;
        pbm.gconst(ig_('CON5')) = 2.0;
        pbm.gconst(ig_('CON6')) = 2.0;
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xupper(ix_('X1')) = 1.0;
        pb.xupper(ix_('X4')) = 1.0;
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 0.0*ones(pb.n,1);
        if(isKey(ix_,'X1'))
            pb.x0(ix_('X1'),1) = 1.0;
        else
            pb.y0(find(pbm.congrps==ig_('X1')),1) = 1.0;
        end
        if(isKey(ix_,'X2'))
            pb.x0(ix_('X2'),1) = 2.0;
        else
            pb.y0(find(pbm.congrps==ig_('X2')),1) = 2.0;
        end
        if(isKey(ix_,'X6'))
            pb.x0(ix_('X6'),1) = 2.0;
        else
            pb.y0(find(pbm.congrps==ig_('X6')),1) = 2.0;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eEPROD',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = 'E1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eEPROD';
        ielftype(ie) = iet_('eEPROD');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        ig = ig_('OBJ');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN                6.66666666
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-COLR2-AN-6-6';
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

    case 'eEPROD'

        EV_  = varargin{1};
        iel_ = varargin{2};
        PROD = exp(EV_(1)*EV_(2));
        varargout{1} = PROD;
        if(nargout>1)
            g_(1,1) = EV_(2)*PROD;
            g_(2,1) = EV_(1)*PROD;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = EV_(2)*EV_(2)*PROD;
                H_(1,2) = PROD+EV_(1)*EV_(2)*PROD;
                H_(2,1) = H_(1,2);
                H_(2,2) = EV_(1)*EV_(1)*PROD;
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

