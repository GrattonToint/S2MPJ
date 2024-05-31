function varargout = LUKVLI11(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : LUKVLI11
%    *********
% 
%    Source: Problem 5.11, the chained HS46 problem, 
%    due to L. Luksan and J. Vlcek,
%    "Sparse and partially separable test problems for 
%    unconstrained and equality constrained optimization",
%    Technical Report 767, Inst. Computer Science, Academy of Sciences
%    of the Czech Republic, 182 07 Prague, Czech Republic, 1999
% 
%    Equality constraints changed to inequalities
% 
%    SIF input: Nick Gould, April 2001
% 
%    classification = 'OOR2-AY-V-V'
% 
%    some useful parameters, including N, the number of variables.
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'LUKVLI11';

switch(action)

    case 'setup'

    pb.name      = 'LUKVLI11';
    pb.sifpbname = 'LUKVLI11';
    pbm.name     = 'LUKVLI11';
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        if(nargin<2)
            v_('N') = 8;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
%       Alternative values for the SIF file parameters:
% IE N                   98             $-PARAMETER
% IE N                   998            $-PARAMETER
% IE N                   9998           $-PARAMETER
% IE N                   99998          $-PARAMETER
        v_('1') = 1;
        v_('2') = 2;
        v_('3') = 3;
        v_('4') = 4;
        v_('5') = 5;
        v_('N-2') = -2+v_('N');
        v_('(N-2)/3') = fix(v_('N-2')/v_('3'));
        v_('NC') = v_('2')*v_('(N-2)/3');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2xlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('(N-2)/3')
            v_('I-1') = -1+I;
            v_('J') = v_('3')*v_('I-1');
            v_('J+1') = 1+v_('J');
            v_('J+2') = 2+v_('J');
            v_('J+3') = 3+v_('J');
            v_('J+4') = 4+v_('J');
            v_('J+5') = 5+v_('J');
            [ig,ig_] = s2xlib('ii',['OBJ1',int2str(I)],ig_);
            gtype{ig} = '<>';
            iv = ix_(['X',int2str(round(v_('J+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['X',int2str(round(v_('J+2')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] = s2xlib('ii',['OBJ2',int2str(I)],ig_);
            gtype{ig} = '<>';
            iv = ix_(['X',int2str(round(v_('J+3')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            [ig,ig_] = s2xlib('ii',['OBJ3',int2str(I)],ig_);
            gtype{ig} = '<>';
            iv = ix_(['X',int2str(round(v_('J+4')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            [ig,ig_] = s2xlib('ii',['OBJ4',int2str(I)],ig_);
            gtype{ig} = '<>';
            iv = ix_(['X',int2str(round(v_('J+5')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
        end
        for K=v_('1'):v_('2'):v_('NC')
            v_('K+1') = 1+K;
            [ig,ig_] = s2xlib('ii',['C',int2str(K)],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['C',int2str(K)];
            [ig,ig_] = s2xlib('ii',['C',int2str(round(v_('K+1')))],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['C',int2str(round(v_('K+1')))];
            iv = ix_(['X',int2str(round(v_('K+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
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
        for I=v_('1'):v_('(N-2)/3')
            pbm.gconst(ig_(['OBJ2',int2str(I)])) = 1.0;
            pbm.gconst(ig_(['OBJ3',int2str(I)])) = 1.0;
            pbm.gconst(ig_(['OBJ4',int2str(I)])) = 1.0;
        end
        for K=v_('1'):v_('2'):v_('NC')
            v_('K+1') = 1+K;
            pbm.gconst(ig_(['C',int2str(K)])) = 1.0;
            pbm.gconst(ig_(['C',int2str(round(v_('K+1')))])) = 2.0;
        end
        pb.xlower = zeros(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        for I=v_('1'):v_('3'):v_('N')
            pb.x0(ix_(['X',int2str(I)]),1) = 2.0;
        end
        for I=v_('2'):v_('3'):v_('N')
            pb.x0(ix_(['X',int2str(I)]),1) = 1.5;
        end
        for I=v_('3'):v_('3'):v_('N')
            pb.x0(ix_(['X',int2str(I)]),1) = 0.5;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2xlib( 'ii', 'eC21',iet_);
        elftv{it}{1} = 'V';
        elftv{it}{2} = 'W';
        [it,iet_] = s2xlib( 'ii', 'eS',iet_);
        elftv{it}{1} = 'V';
        elftv{it}{2} = 'W';
        [it,iet_] = s2xlib( 'ii', 'eC42',iet_);
        elftv{it}{1} = 'V';
        elftv{it}{2} = 'W';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for K=v_('1'):v_('2'):v_('NC')
            v_('K+1') = 1+K;
            v_('K+2') = 2+K;
            v_('K+3') = 3+K;
            v_('K+4') = 4+K;
            ename = ['EA',int2str(K)];
            [ie,ie_] = s2xlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eC21';
            ielftype(ie) = iet_('eC21');
            vname = ['X',int2str(K)];
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('K+3')))];
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('W',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['EB',int2str(K)];
            [ie,ie_] = s2xlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eS';
            ielftype(ie) = iet_('eS');
            vname = ['X',int2str(round(v_('K+3')))];
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('K+4')))];
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('W',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['E',int2str(round(v_('K+1')))];
            [ie,ie_] = s2xlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eC21';
            ielftype(ie) = iet_('eC21');
            ename = ['E',int2str(round(v_('K+1')))];
            [ie,ie_] = s2xlib('ii',ename,ie_);
            vname = ['X',int2str(round(v_('K+2')))];
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['E',int2str(round(v_('K+1')))];
            [ie,ie_] = s2xlib('ii',ename,ie_);
            vname = ['X',int2str(round(v_('K+3')))];
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('W',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2xlib('ii','gL2',igt_);
        [it,igt_] = s2xlib('ii','gL4',igt_);
        [it,igt_] = s2xlib('ii','gL6',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('(N-2)/3')
            ig = ig_(['OBJ1',int2str(I)]);
            pbm.grftype{ig} = 'gL2';
            ig = ig_(['OBJ2',int2str(I)]);
            pbm.grftype{ig} = 'gL2';
            ig = ig_(['OBJ3',int2str(I)]);
            pbm.grftype{ig} = 'gL4';
            ig = ig_(['OBJ4',int2str(I)]);
            pbm.grftype{ig} = 'gL6';
        end
        for K=v_('1'):v_('2'):v_('NC')
            v_('K+1') = 1+K;
            ig = ig_(['C',int2str(K)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['EA',int2str(K)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['EB',int2str(K)]);
            pbm.grelw{ig}(posel) = 1.;
            ig = ig_(['C',int2str(round(v_('K+1')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('K+1')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'OOR2-AY-V-V';
        varargout{1} = pb;
        varargout{2} = pbm;

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eC21'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(1)*EV_(2);
        if(nargout>1)
            g_(1,1) = 2.0*EV_(1)*EV_(2);
            g_(2,1) = EV_(1)*EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 2.0*EV_(2);
                H_(1,2) = 2.0*EV_(1);
                H_(2,1) = H_(1,2);
                varargout{3} = H_;
            end
        end

    case 'eS'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(1,2);
        U_(1,1) = U_(1,1)+1;
        U_(1,2) = U_(1,2)-1;
        IV_(1) = U_(1,:)*EV_;
        SIND = sin(IV_(1));
        varargout{1} = SIND;
        if(nargout>1)
            g_(1,1) = cos(IV_(1));
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_(1,1) = -SIND;
                varargout{3} = U_.'*H_*U_;
            end
        end

    case 'eC42'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)^4*EV_(2)^2;
        if(nargout>1)
            g_(1,1) = 4.0*EV_(1)^3*EV_(2)^2;
            g_(2,1) = 2.0*EV_(1)^4*EV_(2);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 12.0*EV_(1)^2*EV_(2)^2;
                H_(1,2) = 8.0*EV_(1)^3*EV_(2);
                H_(2,1) = H_(1,2);
                H_(2,2) = 2.0*EV_(1)^4;
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

    case 'gL4'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = GVAR_^4;
        if(nargout>1)
            g_ = 4.0*GVAR_^3;
            varargout{2} = g_;
            if(nargout>2)
                H_ = 12.0*GVAR_^2;
                varargout{3} = H_;
            end
        end

    case 'gL6'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = GVAR_^6;
        if(nargout>1)
            g_ = 6.0*GVAR_^5;
            varargout{2} = g_;
            if(nargout>2)
                H_ = 30.0*GVAR_^4;
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

