function varargout = EIGENBLS(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : EIGENBLS
%    --------
% 
%    Solving symmetric eigenvalue problems as systems of
%    nonlinear equations.
% 
%    The problem is, given a symmetric matrix A, to find an orthogonal
%    matrix Q and diagonal matrix D such that A = Q(T) D Q.
% 
%    Example B: a tridiagonal matrix with diagonals 2 and off diagonals -1
% 
%    Source:  An idea by Nick Gould
% 
%             Least-squares version
% 
%    SIF input: Nick Gould, Nov 1992.
% 
%    classification = 'SUR2-AN-V-0'
% 
%    The dimension of the matrix.
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'EIGENBLS';

switch(action)

    case 'setup'

    pb.name      = 'EIGENBLS';
    pb.sifpbname = 'EIGENBLS';
    pbm.name     = 'EIGENBLS';
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        if(nargin<2)
            v_('N') = 2;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
%       Alternative values for the SIF file parameters:
% IE N                   10             $-PARAMETER     original value
% IE N                   50             $-PARAMETER
        v_('1') = 1;
        v_('2') = 2;
        v_('N-1') = -1+v_('N');
        v_(['A',int2str(round(v_('1'))),',',int2str(round(v_('1')))]) = 2.0;
        for J=v_('2'):v_('N')
            v_('J-1') = -1+J;
            v_('J-2') = -2+J;
            for I=v_('1'):v_('J-2')
                v_(['A',int2str(I),',',int2str(J)]) = 0.0;
            end
            v_(['A',int2str(round(v_('J-1'))),',',int2str(J)]) = -1.0;
            v_(['A',int2str(J),',',int2str(J)]) = 2.0;
        end
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for J=v_('1'):v_('N')
            [iv,ix_] = s2xlib('ii',['D',int2str(J)],ix_);
            pb.xnames{iv} = ['D',int2str(J)];
            for I=v_('1'):v_('N')
                [iv,ix_] = s2xlib('ii',['Q',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['Q',int2str(I),',',int2str(J)];
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for J=v_('1'):v_('N')
            for I=v_('1'):J
                [ig,ig_] = s2xlib('ii',['E',int2str(I),',',int2str(J)],ig_);
                gtype{ig} = '<>';
                [ig,ig_] = s2xlib('ii',['O',int2str(I),',',int2str(J)],ig_);
                gtype{ig} = '<>';
            end
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        for J=v_('1'):v_('N')
            pbm.gconst(ig_(['O',int2str(J),',',int2str(J)])) = 1.0;
            for I=v_('1'):J
                pbm.gconst(ig_(['E',int2str(I),',',int2str(J)])) =...
                      v_(['A',int2str(I),',',int2str(J)]);
            end
        end
        pb.xlower = zeros(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 0.0*ones(pb.n,1);
        for J=v_('1'):v_('N')
            pb.x0(ix_(['D',int2str(J)]),1) = 1.0;
            pb.x0(ix_(['Q',int2str(J),',',int2str(J)]),1) = 1.0;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2xlib( 'ii', 'en2PROD',iet_);
        elftv{it}{1} = 'Q1';
        elftv{it}{2} = 'Q2';
        [it,iet_] = s2xlib( 'ii', 'en3PROD',iet_);
        elftv{it}{1} = 'Q1';
        elftv{it}{2} = 'Q2';
        elftv{it}{3} = 'D';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for J=v_('1'):v_('N')
            for I=v_('1'):J
                for K=v_('1'):v_('N')
                    ename = ['E',int2str(I),',',int2str(J),',',int2str(K)];
                    [ie,ie_] = s2xlib('ii',ename,ie_);
                    pbm.elftype{ie} = 'en3PROD';
                    ielftype(ie) = iet_('en3PROD');
                    vname = ['Q',int2str(K),',',int2str(I)];
                    [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],0.0);
                    posev = find(strcmp('Q1',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                    vname = ['Q',int2str(K),',',int2str(J)];
                    [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],0.0);
                    posev = find(strcmp('Q2',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                    vname = ['D',int2str(K)];
                    [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],0.0);
                    posev = find(strcmp('D',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                    ename = ['O',int2str(I),',',int2str(J),',',int2str(K)];
                    [ie,ie_] = s2xlib('ii',ename,ie_);
                    pbm.elftype{ie} = 'en2PROD';
                    ielftype(ie) = iet_('en2PROD');
                    vname = ['Q',int2str(K),',',int2str(I)];
                    [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],0.0);
                    posev = find(strcmp('Q1',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                    vname = ['Q',int2str(K),',',int2str(J)];
                    [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],0.0);
                    posev = find(strcmp('Q2',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2xlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for ig = 1:ngrp
            pbm.grftype{ig} = 'gL2';
        end
        for J=v_('1'):v_('N')
            for I=v_('1'):J
                for K=v_('1'):v_('N')
                    ig = ig_(['E',int2str(I),',',int2str(J)]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) = ie_(['E',int2str(I),',',int2str(J),',',int2str(K)]);
                    pbm.grelw{ig}(posel) = 1.;
                    ig = ig_(['O',int2str(I),',',int2str(J)]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) = ie_(['O',int2str(I),',',int2str(J),',',int2str(K)]);
                    pbm.grelw{ig}(posel) = 1.;
                end
            end
        end
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'SUR2-AN-V-0';
        varargout{1} = pb;
        varargout{2} = pbm;

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'en2PROD'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(2);
        if(nargout>1)
            g_(1,1) = EV_(2);
            g_(2,1) = EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = 1.0e+0;
                H_(2,1) = H_(1,2);
                varargout{3} = H_;
            end
        end

    case 'en3PROD'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(2)*EV_(3);
        if(nargout>1)
            g_(1,1) = EV_(2)*EV_(3);
            g_(2,1) = EV_(1)*EV_(3);
            g_(3,1) = EV_(1)*EV_(2);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,2) = EV_(3);
                H_(2,1) = H_(1,2);
                H_(1,3) = EV_(2);
                H_(3,1) = H_(1,3);
                H_(2,3) = EV_(1);
                H_(3,2) = H_(2,3);
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
                H_ = 2.0e+0;
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
