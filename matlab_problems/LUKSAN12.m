function varargout = LUKSAN12(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : LUKSAN12
%    *********
% 
%    Problem 12 (chained and modified HS47) in the paper
% 
%      L. Luksan
%      Hybrid methods in large sparse nonlinear least squares
%      J. Optimization Theory & Applications 89(3) 575-595 (1996)
% 
%    SIF input: Nick Gould, June 2017.
% 
%    classification = 'NOR2-AN-V-V'
% 
%   seed for dimensions
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'LUKSAN12';

switch(action)

    case 'setup'

    pb.name      = 'LUKSAN12';
    pb.sifpbname = 'LUKSAN12';
    pbm.name     = 'LUKSAN12';
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        v_('S') = 32;
        v_('N') = 3*v_('S');
        v_('N') = 2+v_('N');
        v_('M') = 6*v_('S');
        v_('1') = 1;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2xlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        v_('I') = 1;
        v_('K') = 1;
        for J=v_('1'):v_('S')
            v_('K+1') = 1+v_('K');
            v_('K+2') = 2+v_('K');
            v_('K+3') = 3+v_('K');
            v_('K+4') = 4+v_('K');
            v_('K+5') = 5+v_('K');
            v_('I+1') = 1+v_('I');
            v_('I+2') = 2+v_('I');
            [ig,ig_] = s2xlib('ii',['E',int2str(round(v_('K')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['E',int2str(round(v_('K')))];
            iv = ix_(['X',int2str(round(v_('I+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -10.0e0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -10.0e0;
            end
            [ig,ig_] = s2xlib('ii',['E',int2str(round(v_('K+1')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['E',int2str(round(v_('K+1')))];
            iv = ix_(['X',int2str(round(v_('I+2')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0e0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0e0;
            end
            [ig,ig_] = s2xlib('ii',['E',int2str(round(v_('K+2')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['E',int2str(round(v_('K+2')))];
            [ig,ig_] = s2xlib('ii',['E',int2str(round(v_('K+3')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['E',int2str(round(v_('K+3')))];
            [ig,ig_] = s2xlib('ii',['E',int2str(round(v_('K+4')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['E',int2str(round(v_('K+4')))];
            [ig,ig_] = s2xlib('ii',['E',int2str(round(v_('K+5')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['E',int2str(round(v_('K+5')))];
            iv = ix_(['X',int2str(round(v_('I+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0e0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0e0;
            end
            v_('I') = 3+v_('I');
            v_('K') = 6+v_('K');
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
        v_('K') = 1;
        for J=v_('1'):v_('S')
            v_('K+1') = 1+v_('K');
            v_('K+4') = 4+v_('K');
            v_('K+5') = 5+v_('K');
            pbm.gconst(ig_(['E',int2str(round(v_('K+1')))])) = 1.0;
            pbm.gconst(ig_(['E',int2str(round(v_('K+4')))])) = 10.0;
            pbm.gconst(ig_(['E',int2str(round(v_('K+5')))])) = 20.0;
            v_('K') = 6+v_('K');
        end
        pb.xlower = zeros(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        for I=v_('1'):v_('N')
            pb.x0(ix_(['X',int2str(I)]),1) = -1.0;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2xlib( 'ii', 'eE1',iet_);
        elftv{it}{1} = 'X0';
        [it,iet_] = s2xlib( 'ii', 'eE3',iet_);
        elftv{it}{1} = 'X3';
        [it,iet_] = s2xlib( 'ii', 'eE4',iet_);
        elftv{it}{1} = 'X4';
        [it,iet_] = s2xlib( 'ii', 'eE5',iet_);
        elftv{it}{1} = 'X0';
        elftv{it}{2} = 'X3';
        [it,iet_] = s2xlib( 'ii', 'eF5',iet_);
        elftv{it}{1} = 'X3';
        elftv{it}{2} = 'X4';
        [it,iet_] = s2xlib( 'ii', 'eE6',iet_);
        elftv{it}{1} = 'X2';
        elftv{it}{2} = 'X3';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        v_('I') = 1;
        v_('K') = 1;
        for J=v_('1'):v_('S')
            v_('K+2') = 2+v_('K');
            v_('K+3') = 3+v_('K');
            v_('K+4') = 4+v_('K');
            v_('K+5') = 5+v_('K');
            v_('I+2') = 2+v_('I');
            v_('I+3') = 3+v_('I');
            v_('I+4') = 4+v_('I');
            ename = ['E',int2str(round(v_('K')))];
            [ie,ie_] = s2xlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eE1';
            ielftype(ie) = iet_('eE1');
            ename = ['E',int2str(round(v_('K')))];
            [ie,ie_] = s2xlib('ii',ename,ie_);
            vname = ['X',int2str(round(v_('I')))];
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X0',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['E',int2str(round(v_('K+2')))];
            [ie,ie_] = s2xlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eE3';
            ielftype(ie) = iet_('eE3');
            ename = ['E',int2str(round(v_('K+2')))];
            [ie,ie_] = s2xlib('ii',ename,ie_);
            vname = ['X',int2str(round(v_('I+3')))];
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X3',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['E',int2str(round(v_('K+3')))];
            [ie,ie_] = s2xlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eE4';
            ielftype(ie) = iet_('eE4');
            ename = ['E',int2str(round(v_('K+3')))];
            [ie,ie_] = s2xlib('ii',ename,ie_);
            vname = ['X',int2str(round(v_('I+4')))];
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X4',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['E',int2str(round(v_('K+4')))];
            [ie,ie_] = s2xlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eE5';
            ielftype(ie) = iet_('eE5');
            ename = ['E',int2str(round(v_('K+4')))];
            [ie,ie_] = s2xlib('ii',ename,ie_);
            vname = ['X',int2str(round(v_('I')))];
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X0',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['E',int2str(round(v_('K+4')))];
            [ie,ie_] = s2xlib('ii',ename,ie_);
            vname = ['X',int2str(round(v_('I+3')))];
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X3',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['F',int2str(round(v_('K+4')))];
            [ie,ie_] = s2xlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eF5';
            ielftype(ie) = iet_('eF5');
            ename = ['F',int2str(round(v_('K+4')))];
            [ie,ie_] = s2xlib('ii',ename,ie_);
            vname = ['X',int2str(round(v_('I+3')))];
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X3',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['F',int2str(round(v_('K+4')))];
            [ie,ie_] = s2xlib('ii',ename,ie_);
            vname = ['X',int2str(round(v_('I+4')))];
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X4',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['E',int2str(round(v_('K+5')))];
            [ie,ie_] = s2xlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eE6';
            ielftype(ie) = iet_('eE6');
            ename = ['E',int2str(round(v_('K+5')))];
            [ie,ie_] = s2xlib('ii',ename,ie_);
            vname = ['X',int2str(round(v_('I+2')))];
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['E',int2str(round(v_('K+5')))];
            [ie,ie_] = s2xlib('ii',ename,ie_);
            vname = ['X',int2str(round(v_('I+3')))];
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X3',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            v_('I') = 3+v_('I');
            v_('K') = 6+v_('K');
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        v_('K') = 1;
        for J=v_('1'):v_('S')
            v_('K+2') = 2+v_('K');
            v_('K+3') = 3+v_('K');
            v_('K+4') = 4+v_('K');
            v_('K+5') = 5+v_('K');
            ig = ig_(['E',int2str(round(v_('K')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('K')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            ig = ig_(['E',int2str(round(v_('K+2')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('K+2')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            ig = ig_(['E',int2str(round(v_('K+3')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('K+3')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            ig = ig_(['E',int2str(round(v_('K+4')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('K+4')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['F',int2str(round(v_('K+4')))]);
            pbm.grelw{ig}(posel) = 1.;
            ig = ig_(['E',int2str(round(v_('K+5')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('K+5')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            v_('K') = 6+v_('K');
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'NOR2-AN-V-V';
        varargout{1} = pb;
        varargout{2} = pbm;

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eE1'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = 10.0e0*EV_(1)^2;
        if(nargout>1)
            g_(1,1) = 20.0e0*EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 20.0e0;
                varargout{3} = H_;
            end
        end

    case 'eE3'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = (EV_(1)-1.0e0)^2;
        if(nargout>1)
            g_(1,1) = 2.0e0*(EV_(1)-1.0e0);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0e0;
                varargout{3} = H_;
            end
        end

    case 'eE4'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = (EV_(1)-1.0e0)^3;
        if(nargout>1)
            g_(1,1) = 3.0e0*(EV_(1)-1.0e0)^2;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 6.0e0*(EV_(1)-1.0e0);
                varargout{3} = H_;
            end
        end

    case 'eE5'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(2)*EV_(1)*EV_(1);
        if(nargout>1)
            g_(1,1) = 2.0e0*EV_(2)*EV_(1);
            g_(2,1) = EV_(1)*EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 2.0e0*EV_(2);
                H_(1,2) = 2.0e0*EV_(1);
                H_(2,1) = H_(1,2);
                varargout{3} = H_;
            end
        end

    case 'eF5'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(1,2);
        U_(1,1) = U_(1,1)+1;
        U_(1,2) = U_(1,2)-1;
        IV_(1) = U_(1,:)*EV_;
        varargout{1} = sin(IV_(1));
        if(nargout>1)
            g_(1,1) = cos(IV_(1));
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_(1,1) = -sin(IV_(1));
                varargout{3} = U_.'*H_*U_;
            end
        end

    case 'eE6'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = (EV_(1)^4)*(EV_(2)^2);
        if(nargout>1)
            g_(1,1) = 4.0e0*(EV_(1)^3)*(EV_(2)^2);
            g_(2,1) = 2.0e0*(EV_(1)^4)*EV_(2);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 12.0e0*(EV_(1)^2)*(EV_(2)^2);
                H_(1,2) = 8.0e0*(EV_(1)^3)*EV_(2);
                H_(2,1) = H_(1,2);
                H_(2,2) = 2.0e0*(EV_(1)^4);
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
