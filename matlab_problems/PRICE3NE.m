function varargout = PRICE3NE(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : PRICE3NE
%    *********
% 
%    SCIPY global optimization benchmark example Price03
% 
%    Fit: y  =  (10(x_1^2-x_2),sqrt(6)(6.4(x_2-0.5)^2 -x_1)) +  e
% 
%    Source:  Problem from the SCIPY benchmark set
%      https://github.com/scipy/scipy/tree/master/benchmarks/ ...
%              benchmarks/go_benchmark_functions
% 
%    Nonlinear-equation formulation of PRICE3.SIF
% 
%    SIF input: Nick Gould, Jan 2020
% 
%    classification = 'NOR2-MN-2-2'
% 
%    Number of data values
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'PRICE3NE';

switch(action)

    case 'setup'

    pb.name      = 'PRICE3NE';
    pb.sifpbname = 'PRICE3NE';
    pbm.name     = 'PRICE3NE';
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        v_('M') = 2;
        v_('N') = 2;
        v_('1') = 1;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2xlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2xlib('ii','F1',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'F1';
        iv = ix_('X2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        pbm.gscale(ig,1) = 0.1;
        [ig,ig_] = s2xlib('ii','F2',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'F2';
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
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
        pbm.congrps = find(ismember(gtype,{'<=','==','>='}));
        [pb.cnames{1:pb.m}] = deal(cnames{pbm.congrps});
        pb.nob = ngrp-pb.m;
        pbm.objgrps = find(strcmp(gtype,'<>'));
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        pbm.gconst(ig_('F2')) = 0.6;
        pb.xlower = zeros(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,'X1'))
            pb.x0(ix_('X1'),1) = 1.0;
        else
            pb.y0(find(pbm.congrps==ig_('X1')),1) = 1.0;
        end
        if(isKey(ix_,'X2'))
            pb.x0(ix_('X2'),1) = 5.0;
        else
            pb.y0(find(pbm.congrps==ig_('X2')),1) = 5.0;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2xlib( 'ii', 'eSQR',iet_);
        elftv{it}{1} = 'X';
        [it,iet_] = s2xlib( 'ii', 'eSSQR',iet_);
        elftv{it}{1} = 'X';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = 'E1';
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        vname = 'X1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E2';
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSSQR';
        ielftype(ie) = iet_('eSSQR');
        vname = 'X2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        ig = ig_('F1');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('F2');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E2');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 6.4;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'NOR2-MN-2-2';
        varargout{1} = pb;
        varargout{2} = pbm;

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eSQR'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)^2;
        if(nargout>1)
            g_(1,1) = EV_(1)+EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0;
                varargout{3} = H_;
            end
        end

    case 'eSSQR'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = (EV_(1)-0.5)^2;
        if(nargout>1)
            g_(1,1) = EV_(1)+EV_(1)-1.0;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0;
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    case {'fx','fgx','fgHx','cx','cJx','cJHx','cIx','cIJx','cIJHx','cIJxv','fHxv',...
          'cJxv','Lxy','Lgxy','LgHxy','LIxy','LIgxy','LIgHxy','LHxyv','LIHxyv'}

        if(isfield(pbm,'name')&&strcmp(pbm.name,name))
            pbm.has_globs = [0,0];
            [varargout{1:max(1,nargout)}] = s2xlib(action,pbm,varargin{:});
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

