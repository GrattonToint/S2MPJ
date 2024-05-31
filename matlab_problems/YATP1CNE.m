function varargout = YATP1CNE(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : YATP1CNE
%    *********
% 
%    Yet another test problem involving double pseudo-stochastic constraints
%    on a square matrix. If the matrix dimension is N, the number of
%    variables is equal to  N**2 + 2* N. The equations are
%    x_{ij}^3 - A x_{ij}^2 - ( y_i + z_j ) ( x_{ij}cos(x_{ij} - sin(x_{ij}) ) = 0
%                                                     (i,j = 1, ..., N )
%    \sum_i^N sin(x_{ij}) / x_{ij} = 1                  (j = 1,..., N)
%    \sum_j^N sin(x_{ij}) / x_{ij} = 1                  (i = 1,..., N)
%    The problem is non convex.
% 
%    Source:
%    a late evening idea by Ph. Toint
% 
%    SIF input: Ph. Toint, June 2003.
%               corrected Nick Gould, March 2019
% 
%    classification = 'NOR2-AN-V-V'
% 
%    The dimension of the matrix
% 
%    nonlinear equation version
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'YATP1CNE';

switch(action)

    case 'setup'

    pb.name      = 'YATP1CNE';
    pb.sifpbname = 'YATP1CNE';
    pbm.name     = 'YATP1CNE';
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        if(nargin<2)
            v_('N') = 5;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
%       Alternative values for the SIF file parameters:
% IE N                   10             $-PARAMETER n = 120
% IE N                   50             $-PARAMETER n = 2600
% IE N                   100            $-PARAMETER n = 10200
% IE N                   200            $-PARAMETER n = 40400
% IE N                   350            $-PARAMETER n = 123200
        v_('A') = 10.0;
        v_('1') = 1;
        v_('-A') = -1.0*v_('A');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            for J=v_('1'):v_('N')
                [iv,ix_] = s2xlib('ii',['X',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['X',int2str(I),',',int2str(J)];
            end
        end
        for I=v_('1'):v_('N')
            [iv,ix_] = s2xlib('ii',['Y',int2str(I)],ix_);
            pb.xnames{iv} = ['Y',int2str(I)];
            [iv,ix_] = s2xlib('ii',['Z',int2str(I)],ix_);
            pb.xnames{iv} = ['Z',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('N')
            for J=v_('1'):v_('N')
                [ig,ig_] = s2xlib('ii',['E',int2str(I),',',int2str(J)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['E',int2str(I),',',int2str(J)];
            end
        end
        for I=v_('1'):v_('N')
            [ig,ig_] = s2xlib('ii',['ER',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['ER',int2str(I)];
            [ig,ig_] = s2xlib('ii',['EC',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['EC',int2str(I)];
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
        for I=v_('1'):v_('N')
            pbm.gconst(ig_(['ER',int2str(I)])) = 1.0;
            pbm.gconst(ig_(['EC',int2str(I)])) = 1.0;
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
            for J=v_('1'):v_('N')
                if(isKey(ix_,['X',int2str(I),',',int2str(J)]))
                    pb.x0(ix_(['X',int2str(I),',',int2str(J)]),1) = 6.0;
                else
                    pb.y0(find(pbm.congrps==ig_(['X',int2str(I),',',int2str(J)])),1) = 6.0;
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2xlib( 'ii', 'eSQ',iet_);
        elftv{it}{1} = 'X';
        [it,iet_] = s2xlib( 'ii', 'eCB',iet_);
        elftv{it}{1} = 'X';
        [it,iet_] = s2xlib( 'ii', 'eLXC',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        elftv{it}{3} = 'Z';
        [it,iet_] = s2xlib( 'ii', 'eLS',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        elftv{it}{3} = 'Z';
        [it,iet_] = s2xlib( 'ii', 'eRAT',iet_);
        elftv{it}{1} = 'X';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('N')
            for J=v_('1'):v_('N')
                ename = ['CB',int2str(I),',',int2str(J)];
                [ie,ie_] = s2xlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eCB';
                ielftype(ie) = iet_('eCB');
                vname = ['X',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('X',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['SQ',int2str(I),',',int2str(J)];
                [ie,ie_] = s2xlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eSQ';
                ielftype(ie) = iet_('eSQ');
                vname = ['X',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('X',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['DC',int2str(I),',',int2str(J)];
                [ie,ie_] = s2xlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eLXC';
                ielftype(ie) = iet_('eLXC');
                vname = ['X',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('X',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['Y',int2str(I)];
                [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('Y',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['Z',int2str(J)];
                [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('Z',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['DS',int2str(I),',',int2str(J)];
                [ie,ie_] = s2xlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eLS';
                ielftype(ie) = iet_('eLS');
                vname = ['X',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('X',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['Y',int2str(I)];
                [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('Y',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['Z',int2str(J)];
                [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('Z',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['SX',int2str(I),',',int2str(J)];
                [ie,ie_] = s2xlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eRAT';
                ielftype(ie) = iet_('eRAT');
                vname = ['X',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('X',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('N')
            for J=v_('1'):v_('N')
                ig = ig_(['E',int2str(I),',',int2str(J)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['CB',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.0;
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['SQ',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-A');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['DC',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = -1.0;
                posel = posel+1;
                pbm.grelt{ig}(posel) = ie_(['DS',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = 1.0;
            end
        end
        for I=v_('1'):v_('N')
            for J=v_('1'):v_('N')
                ig = ig_(['ER',int2str(I)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['SX',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
                ig = ig_(['EC',int2str(J)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['SX',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
            end
        end
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

    case 'eSQ'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(1);
        if(nargout>1)
            g_(1,1) = EV_(1)+EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0e0;
                varargout{3} = H_;
            end
        end

    case 'eCB'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(1)*EV_(1);
        if(nargout>1)
            g_(1,1) = 3.0e0*EV_(1)*EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 6.0e0*EV_(1);
                varargout{3} = H_;
            end
        end

    case 'eLXC'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(2,3);
        U_(1,1) = U_(1,1)+1;
        U_(2,2) = U_(2,2)+1;
        U_(2,3) = U_(2,3)+1;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        C = cos(IV_(1));
        S = sin(IV_(1));
        varargout{1} = IV_(2)*IV_(1)*C;
        if(nargout>1)
            g_(2,1) = IV_(1)*C;
            g_(1,1) = IV_(2)*(C-IV_(1)*S);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(2,1) = C-IV_(1)*S;
                H_(1,2) = H_(2,1);
                H_(1,1) = -IV_(2)*(S+S+IV_(1)*C);
                varargout{3} = U_.'*H_*U_;
            end
        end

    case 'eLS'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(2,3);
        U_(1,1) = U_(1,1)+1;
        U_(2,2) = U_(2,2)+1;
        U_(2,3) = U_(2,3)+1;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        C = cos(IV_(1));
        S = sin(IV_(1));
        varargout{1} = IV_(2)*S;
        if(nargout>1)
            g_(2,1) = S;
            g_(1,1) = IV_(2)*C;
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(2,1) = C;
                H_(1,2) = H_(2,1);
                H_(1,1) = -IV_(2)*S;
                varargout{3} = U_.'*H_*U_;
            end
        end

    case 'eRAT'

        EV_  = varargin{1};
        iel_ = varargin{2};
        C = cos(EV_(1));
        S = sin(EV_(1));
        varargout{1} = S/EV_(1);
        if(nargout>1)
            g_(1,1) = (C-S/EV_(1))/EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = -S/EV_(1)-(C+C)/EV_(1)^2+(S+S)/EV_(1)^3;
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

