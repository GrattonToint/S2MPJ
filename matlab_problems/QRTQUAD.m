function varargout = QRTQUAD(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : QRTQUAD
%    *********
% 
%    An unconstrained problem involving elements with different degree of
%    nonlinearity (quartics and quadratics).  It contains N elements the first
%    M of which are quartic.  The Hessian is tridiagonal for the first M
%    elements and arrow-head for the remaining ones.  The first M variables
%    have lower and upper bounds.
% 
%    SIF input: Ph. Toint, June 1992.
%               minor correction by Ph. Shott, Jan 1995.
%               correction by S. Gratton & Ph. Toint, May 2024
% 
%    classification = 'C-COBR2-MN-V-0'
% 
%       Alternative values for the SIF file parameters:
% IE N                   12             $-PARAMETER
% IE N                   120            $-PARAMETER     original value
% IE N                   1200           $-PARAMETER
% IE N                   5000           $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'QRTQUAD';

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
        v_  = containers.Map('KeyType','char', 'ValueType', 'double');
        ix_ = containers.Map('KeyType','char', 'ValueType', 'double');
        ig_ = containers.Map('KeyType','char', 'ValueType', 'double');
        if(nargs<1)
            v_('N') = 12;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
% IE M                   6              $-PARAMETER
% IE M                   10             $-PARAMETER     original value
% IE M                   100            $-PARAMETER
% IE M                   1100           $-PARAMETER
        if(nargs<2)
            v_('M') = 100;  %  SIF file default value
        else
            v_('M') = varargin{2};
        end
        v_('0') = 0;
        v_('1') = 1;
        v_('2') = 2;
        v_('M+1') = 1+v_('M');
        v_('N-1') = -1+v_('N');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('N')
            v_('RI') = I;
            v_('C') = - 10.0*v_('RI');
            [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
            gtype{ig} = '<>';
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('C')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('C');
            end
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = ix_.Count;
        ngrp   = ig_.Count;
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        for I=v_('1'):v_('N')
            pb.xlower(ix_(['X',int2str(I)]),1) = 0.0;
            pb.xupper(ix_(['X',int2str(I)])) = 10.0;
        end
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eQUART',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        elftp{it}{1} = 'P';
        [it,iet_] = s2mpjlib( 'ii', 'eQUAD',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        v_('RM') = v_('M');
        for I=v_('1'):v_('M')
            v_('RI') = I;
            v_('C') = v_('RI')/v_('RM');
            v_('I+1') = I+v_('1');
            ename = ['E',int2str(I)];
            [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eQUART';
            ielftype(ie) = iet_('eQUART');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('I+1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('P',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('C');
        end
        for I=v_('M+1'):v_('N-1')
            ename = ['E',int2str(I)];
            [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
            if(newelt)
                pbm.elftype{ie} = 'eQUAD';
                ielftype(ie) = iet_('eQUAD');
            end
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('N')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('N-1')
            ig = ig_('OBJ');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
%    Solution
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'C-COBR2-MN-V-0';
        pb.x0          = zeros(pb.n,1);
        pbm.objderlvl = 2;
        pb.objderlvl = pbm.objderlvl;
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

    case 'eQUART'

        EV_  = varargin{1};
        iel_ = varargin{2};
        XY = EV_(1)*EV_(2);
        varargout{1} = pbm.elpar{iel_}(1)*XY^4;
        if(nargout>1)
            g_(1,1) = pbm.elpar{iel_}(1)*EV_(2)*4.0*XY^3;
            g_(2,1) = pbm.elpar{iel_}(1)*EV_(1)*4.0*XY^3;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 12.0*pbm.elpar{iel_}(1)*(EV_(2)^2)*(XY^2);
                H_(2,2) = 12.0*pbm.elpar{iel_}(1)*(EV_(1)^2)*(XY^2);
                H_(1,2) = 4.0*pbm.elpar{iel_}(1)*(XY^3)+12.0*pbm.elpar{iel_}(1)*EV_(2)*...
                     EV_(1)*(XY^2);
                H_(2,1) = H_(1,2);
                varargout{3} = H_;
            end
        end

    case 'eQUAD'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = 4.0*EV_(1)*EV_(1)+2.0*EV_(2)*EV_(2)+EV_(1)*EV_(2);
        if(nargout>1)
            g_(1,1) = 8.0*EV_(1)+EV_(2);
            g_(2,1) = 4.0*EV_(2)+EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 8.0;
                H_(1,2) = 1.0;
                H_(2,1) = H_(1,2);
                H_(2,2) = 4.0;
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

