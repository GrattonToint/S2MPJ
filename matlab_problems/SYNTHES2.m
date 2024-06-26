function varargout = SYNTHES2(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : SYNTHES2
%    *********
% 
%    Source: Test problem 2 (Synthesis of processing system) in 
%    M. Duran & I.E. Grossmann,
%    "An outer approximation algorithm for a class of mixed integer nonlinear
%     programs", Mathematical Programming 36, pp. 307-339, 1986.
% 
%    SIF input: S. Leyffer, October 1997
% 
%    classification = 'OOR2-AN-11-14'
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'SYNTHES2';

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
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        v_('1') = 1;
        v_('5') = 5;
        v_('6') = 6;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('6')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        for I=v_('1'):v_('5')
            [iv,ix_] = s2mpjlib('ii',['Y',int2str(I)],ix_);
            pb.xnames{iv} = ['Y',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        iv = ix_('Y1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 5.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 5.0;
        end
        iv = ix_('Y2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 8.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 8.0;
        end
        iv = ix_('Y3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 6.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 6.0;
        end
        iv = ix_('Y4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 10.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 10.0;
        end
        iv = ix_('Y5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 6.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 6.0;
        end
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -10.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -10.0;
        end
        iv = ix_('X2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -15.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -15.0;
        end
        iv = ix_('X3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -15.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -15.0;
        end
        iv = ix_('X4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 15.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 15.0;
        end
        iv = ix_('X5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 5.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 5.0;
        end
        iv = ix_('X6');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -20.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -20.0;
        end
        [ig,ig_] = s2mpjlib('ii','N1',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'N1';
        [ig,ig_] = s2mpjlib('ii','N2',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'N2';
        iv = ix_('Y1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -10.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -10.0;
        end
        [ig,ig_] = s2mpjlib('ii','N3',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'N3';
        iv = ix_('Y2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -10.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -10.0;
        end
        [ig,ig_] = s2mpjlib('ii','L1',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'L1';
        iv = ix_('X3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.25+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.25;
        end
        iv = ix_('Y3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -10.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -10.0;
        end
        [ig,ig_] = s2mpjlib('ii','L2',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'L2';
        iv = ix_('X4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('Y4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -10.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -10.0;
        end
        [ig,ig_] = s2mpjlib('ii','L3',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'L3';
        iv = ix_('X3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -2.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -2.0;
        end
        iv = ix_('X6');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2.0;
        end
        iv = ix_('Y5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -10.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -10.0;
        end
        [ig,ig_] = s2mpjlib('ii','L4',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'L4';
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -2.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -2.0;
        end
        iv = ix_('X4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X6');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2.0;
        end
        [ig,ig_] = s2mpjlib('ii','L5',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'L5';
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.75+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.75;
        end
        iv = ix_('X4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X6');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2.0;
        end
        [ig,ig_] = s2mpjlib('ii','L6',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'L6';
        iv = ix_('X3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X6');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [ig,ig_] = s2mpjlib('ii','L7',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'L7';
        iv = ix_('X3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2.0;
        end
        iv = ix_('X4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X6');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -2.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -2.0;
        end
        [ig,ig_] = s2mpjlib('ii','L8',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'L8';
        iv = ix_('X4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.5+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.5;
        end
        iv = ix_('X5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii','L9',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'L9';
        iv = ix_('X4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.2+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.2;
        end
        iv = ix_('X5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [ig,ig_] = s2mpjlib('ii','L10',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'L10';
        iv = ix_('Y1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('Y2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii','L11',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'L11';
        iv = ix_('Y4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('Y5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
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
        pbm.gconst(ig_('OBJ')) = -140.0;
        pbm.gconst(ig_('N2')) = 1.0;
        pbm.gconst(ig_('N3')) = 1.0;
        pbm.gconst(ig_('L10')) = 1.0;
        pbm.gconst(ig_('L11')) = 1.0;
        %%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%
        grange(legrps,1) = Inf*ones(pb.nle,1);
        grange(ig_('L11')) = -1.0;
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xupper(ix_('X1')) = 2.0;
        pb.xupper(ix_('X2')) = 2.0;
        pb.xupper(ix_('X3')) = 2.0;
        pb.xupper(ix_('X6')) = 3.0;
        for I=v_('1'):v_('5')
            pb.xupper(ix_(['Y',int2str(I)])) = 1.0;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eLOGSUM',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        [it,iet_] = s2mpjlib( 'ii', 'eEXPA',iet_);
        elftv{it}{1} = 'X';
        elftp{it}{1} = 'A';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        ename = 'LOGX4X5';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eLOGSUM';
        ielftype(ie) = iet_('eLOGSUM');
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EXPX1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eEXPA';
        ielftype(ie) = iet_('eEXPA');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('A',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        ename = 'EXPX2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eEXPA';
        ielftype(ie) = iet_('eEXPA');
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('A',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.2;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        ig = ig_('OBJ');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EXPX1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('EXPX2');
        pbm.grelw{ig}(posel) = 1.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('LOGX4X5');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -60.0;
        ig = ig_('N1');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('LOGX4X5');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0;
        ig = ig_('N2');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EXPX1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        ig = ig_('N3');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EXPX2');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = grange(legrps);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'OOR2-AN-11-14';
        pb.x0          = zeros(pb.n,1);
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb,  pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm, pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end
% **********************
%  SET UP THE FUNCTION *
%  AND RANGE ROUTINES  *
% **********************

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eLOGSUM'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = log(EV_(1)+EV_(2)+1.0);
        DX = 1.0/(EV_(1)+EV_(2)+1.0);
        DXDX = -1.0/(EV_(1)+EV_(2)+1.0)^2;
        if(nargout>1)
            g_(1,1) = DX;
            g_(2,1) = DX;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = DXDX;
                H_(1,2) = DXDX;
                H_(2,1) = H_(1,2);
                H_(2,2) = DXDX;
                varargout{3} = H_;
            end
        end

    case 'eEXPA'

        EV_  = varargin{1};
        iel_ = varargin{2};
        EXPXA = exp(EV_(1)/pbm.elpar{iel_}(1));
        varargout{1} = EXPXA;
        if(nargout>1)
            g_(1,1) = EXPXA/pbm.elpar{iel_}(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = EXPXA/pbm.elpar{iel_}(1)/pbm.elpar{iel_}(1);
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    case {'fx','fgx','fgHx','cx','cJx','cJHx','cIx','cIJx','cIJHx','cIJxv','fHxv',...
          'cJxv','Lxy','Lgxy','LgHxy','LIxy','LIgxy','LIgHxy','LHxyv','LIHxyv'}

        if(isfield(pbm,'name')&&strcmp(pbm.name,name))
            pbm.has_globs = [0,0];
            [varargout{1:max(1,nargout)}] = s2mpjlib(action,pbm,varargin{:});
        else
            disp(['ERROR: please run ',name,' with action = setup'])
            [varargout{1:nargout}] = deal(repmat(NaN,1:nargout));
        end

    otherwise
        disp([' ERROR: unknown action ',action,' requested from ',name,'.m'])
    end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

