function varargout = PRICE4B(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : PRICE4B
%    *********
% 
%    SCIPY global optimization benchmark example PRICE4
% 
%    Fit: (2x_1^2 x_2 - x_2^3,6x_1-x_2^2+x_2) +  e = 0
% 
%    version with box-constrained feasible region
% 
%    Source:  Problem from the SCIPY benchmark set
%      https://github.com/scipy/scipy/tree/master/benchmarks/ ...
%              benchmarks/go_benchmark_functions
% 
%    SIF input: Nick Gould, July 2021
% 
%    classification = 'SBR2-MN-2-0'
% 
%    Number of data values
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'PRICE4B';

switch(action)

    case 'setup'

    pb.name      = 'PRICE4B';
    pb.sifpbname = 'PRICE4B';
    pbm.name     = 'PRICE4B';
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
        gtype{ig} = '<>';
        [ig,ig_] = s2xlib('ii','F2',ig_);
        gtype{ig} = '<>';
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 6.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 6.0;
        end
        iv = ix_('X2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        pb.xlower = zeros(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -50.0*ones(pb.n,1);
        pb.xupper = 50.0*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.x0(ix_('X1'),1) = 1.0;
        pb.x0(ix_('X2'),1) = 5.0;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2xlib( 'ii', 'eSQR',iet_);
        elftv{it}{1} = 'X';
        [it,iet_] = s2xlib( 'ii', 'eCUBE',iet_);
        elftv{it}{1} = 'X';
        [it,iet_] = s2xlib( 'ii', 'eCUBEL',iet_);
        elftv{it}{1} = 'X1';
        elftv{it}{2} = 'X2';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = 'E12';
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eCUBEL';
        ielftype(ie) = iet_('eCUBEL');
        vname = 'X1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,-50.0,50.0,[]);
        posev = find(strcmp('X1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,-50.0,50.0,[]);
        posev = find(strcmp('X2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E1';
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eCUBE';
        ielftype(ie) = iet_('eCUBE');
        vname = 'X2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,-50.0,50.0,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E2';
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        vname = 'X2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,-50.0,50.0,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2xlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for ig = 1:ngrp
            pbm.grftype{ig} = 'gL2';
        end
        ig = ig_('F1');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E12');
        pbm.grelw{ig}(posel) = 2.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E1');
        pbm.grelw{ig}(posel) = -1.0;
        ig = ig_('F2');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E2');
        pbm.grelw{ig}(posel) = -1.0;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'SBR2-MN-2-0';
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

    case 'eCUBE'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)^3;
        if(nargout>1)
            g_(1,1) = 3.0*EV_(1)^2;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 6.0*EV_(1);
                varargout{3} = H_;
            end
        end

    case 'eCUBEL'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(2)*EV_(1)^3;
        if(nargout>1)
            g_(1,1) = 3.0*EV_(2)*EV_(1)^2;
            g_(2,1) = EV_(1)^3;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 6.0*EV_(2)*EV_(1);
                H_(1,2) = 3.0*EV_(1)^2;
                H_(2,1) = H_(1,2);
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
