function varargout = HETmZ(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HETmZ
%    *********
% 
%    A nonlinear programming formulation of a discretization of
%    a nonlinear Chebychev problem.
% 
%    The problem is
% 
%        min  max || phi(x,w) ||, for all w in the interval I.
%         x    w
% 
%    I is discretized, and the problem solved over the
%    discrete points.
% 
%    Nonlinear programming formulation
%        min   u     s.t.  u - phi >= 0, u + phi >= 0
%        x,u
% 
%    Specific problem: I = [-1,1]
%    phi(x,w) = (1-w^2) - (0.5 x^2 - 2 x w )
% 
%    Source: R. Hettich and P. Znecke,
%    "Numerische Methoden der Approximation und semi-infiniter optimierung",
%    Teubener Studienbucher, 1982
% 
%    SIF input: Nick Gould, February, 1994.
% 
%    classification = 'LQR2-AN-2-V'
% 
%    Discretization
% 
% IE M                   2
% IE M                   100
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HETmZ';

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
        v_('M') = 500;
        v_('LOWER') = -1.0;
        v_('UPPER') = 1.0;
        v_('ONE') = 1.0;
        v_('0') = 0;
        v_('1') = 1;
        v_('DIFF') = v_('UPPER')-v_('LOWER');
        v_('RM') = v_('M');
        v_('H') = v_('DIFF')/v_('RM');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2mpjlib('ii','U',ix_);
        pb.xnames{iv} = 'U';
        [iv,ix_] = s2mpjlib('ii','X',ix_);
        pb.xnames{iv} = 'X';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        iv = ix_('U');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        for I=v_('0'):v_('M')
            v_('RI') = I;
            v_('W') = v_('RI')*v_('H');
            v_('W') = v_('W')+v_('LOWER');
            v_('2W') = 2.0*v_('W');
            v_('-2W') = -1.0*v_('2W');
            [ig,ig_] = s2mpjlib('ii',['LO',int2str(I)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['LO',int2str(I)];
            iv = ix_('U');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_('X');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('2W')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('2W');
            end
            [ig,ig_] = s2mpjlib('ii',['UP',int2str(I)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['UP',int2str(I)];
            iv = ix_('U');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_('X');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('-2W')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('-2W');
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
        pbm.congrps = [ legrps, eqgrps, gegrps ];
        [pb.cnames{1:pb.m}] = deal(cnames{pbm.congrps});
        pb.nob = ngrp-pb.m;
        pbm.objgrps = find(strcmp(gtype,'<>'));
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        for I=v_('0'):v_('M')
            v_('RI') = I;
            v_('W') = v_('RI')*v_('H');
            v_('W') = v_('W')+v_('LOWER');
            v_('W**2') = v_('W')*v_('W');
            v_('1-W**2') = 1.0-v_('W**2');
            v_('-1+W**2') = -1.0*v_('1-W**2');
            pbm.gconst(ig_(['LO',int2str(I)])) = v_('-1+W**2');
            pbm.gconst(ig_(['UP',int2str(I)])) = v_('1-W**2');
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eQUAD',iet_);
        elftv{it}{1} = 'X';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = 'QUAD';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eQUAD';
        ielftype(ie) = iet_('eQUAD');
        vname = 'X';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('0'):v_('M')
            ig = ig_(['LO',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_('QUAD');
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = -1.0;
            ig = ig_(['UP',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_('QUAD');
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'LQR2-AN-2-V';
        pb.x0          = zeros(pb.n,1);
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb,  pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm, pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eQUAD'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = 0.5*EV_(1)*EV_(1);
        if(nargout>1)
            g_(1,1) = EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 1.0;
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

