function varargout = S368(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : S368
%    *********
% 
%    Wolfe's problem.
% 
%    Source:
%    P. Wolfe,
%    "Explicit solution of an optimization problem",
%    Mathematical Programming 2, 258-260, 1972.
% 
%    SIF input: Nick Gould, Oct 1992.
% 
%    See also Schittkowski #368 (for N = 8)
% 
%    classification = 'OBR2-MN-V-0'
% 
%    The number of variables is N.
% 
%       Alternative values for the SIF file parameters:
% IE N                   8              $-PARAMETER Schittkowski #368
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'S368';

switch(action)

    case 'setup'

    pb.name      = 'S368';
    pb.sifpbname = 'S368';
    pbm.name     = 'S368';
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        if(nargin<2)
            v_('N') = 10;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
% IE N                   50             $-PARAMETER
% IE N                   100            $-PARAMETER
        v_('1') = 1;
        v_('N+1') = 1+v_('N');
        v_('RN+1') = v_('N+1');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2xlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for J=v_('1'):v_('N')
            for I=v_('1'):v_('N')
                [ig,ig_] = s2xlib('ii',['M',int2str(I),',',int2str(J)],ig_);
                gtype{ig} = '<>';
                [ig,ig_] = s2xlib('ii',['P',int2str(I),',',int2str(J)],ig_);
                gtype{ig} = '<>';
            end
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        pb.xlower = zeros(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xupper = 1.0*ones(pb.n,1);
        pb.xlower = -Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        for I=v_('1'):v_('N')
            v_('RI') = I;
            v_('T') = v_('RI')/v_('RN+1');
            pb.x0(ix_(['X',int2str(I)]),1) = v_('T');
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2xlib( 'ii', 'ePRODM',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        [it,iet_] = s2xlib( 'ii', 'ePRODP',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for J=v_('1'):v_('N')
            for I=v_('1'):v_('N')
                ename = ['M',int2str(I),',',int2str(J)];
                [ie,ie_] = s2xlib('ii',ename,ie_);
                pbm.elftype{ie} = 'ePRODM';
                ielftype(ie) = iet_('ePRODM');
                vname = ['X',int2str(I)];
                [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],1.0,[]);
                posev = find(strcmp('X',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(J)];
                [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],1.0,[]);
                posev = find(strcmp('Y',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['P',int2str(I),',',int2str(J)];
                [ie,ie_] = s2xlib('ii',ename,ie_);
                pbm.elftype{ie} = 'ePRODP';
                ielftype(ie) = iet_('ePRODP');
                vname = ['X',int2str(I)];
                [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],1.0,[]);
                posev = find(strcmp('X',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(J)];
                [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],1.0,[]);
                posev = find(strcmp('Y',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for J=v_('1'):v_('N')
            for I=v_('1'):v_('N')
                ig = ig_(['M',int2str(I),',',int2str(J)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['M',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = 1.;
                ig = ig_(['P',int2str(I),',',int2str(J)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['P',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = 1.;
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'OBR2-MN-V-0';
        varargout{1} = pb;
        varargout{2} = pbm;

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'ePRODM'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = -EV_(1)^2*EV_(2)^4;
        if(nargout>1)
            g_(1,1) = -2.0e+0*EV_(1)*EV_(2)^4;
            g_(2,1) = -4.0e+0*EV_(1)^2*EV_(2)^3;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = -2.0e+0*EV_(2)^4;
                H_(2,1) = -8.0e+0*EV_(1)*EV_(2)^3;
                H_(1,2) = H_(2,1);
                H_(2,2) = -1.2e+1*EV_(1)^2*EV_(2)^2;
                varargout{3} = H_;
            end
        end

    case 'ePRODP'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)^3*EV_(2)^3;
        if(nargout>1)
            g_(1,1) = 3.0e+0*EV_(1)^2*EV_(2)^3;
            g_(2,1) = 3.0e+0*EV_(1)^3*EV_(2)^2;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 6.0e+0*EV_(1)*EV_(2)^3;
                H_(2,1) = 9.0e+0*EV_(1)^2*EV_(2)^2;
                H_(1,2) = H_(2,1);
                H_(2,2) = 6.0e+0*EV_(1)^3*EV_(2);
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
