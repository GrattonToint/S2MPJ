function varargout = OET3(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : OET3
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
%    Specific problem: I = [0,1]
%    phi(x,w) = sin(w) - (x1 + x2 w + x3 w^2 )
% 
%    Source: K. Oettershagen
%    "Ein superlinear konvergenter algorithmus zur losung 
%     semi-infiniter optimierungsproblem",
%     Ph.D thesis, Bonn University, 1982
% 
%    SIF input: Nick Gould, February, 1994.
% 
%    classification = 'LLR2-AN-4-V'
% 
%    Discretization
% 
% IE M                   2
% IE M                   100
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'OET3';

switch(action)

    case 'setup'

    pb.name      = 'OET3';
    pb.sifpbname = 'OET3';
    pbm.name     = 'OET3';
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        v_('M') = 500;
        v_('LOWER') = 0.0;
        v_('UPPER') = 1.0;
        v_('0') = 0;
        v_('DIFF') = v_('UPPER')-v_('LOWER');
        v_('RM') = v_('M');
        v_('H') = v_('DIFF')/v_('RM');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2xlib('ii','U',ix_);
        pb.xnames{iv} = 'U';
        [iv,ix_] = s2xlib('ii','X1',ix_);
        pb.xnames{iv} = 'X1';
        [iv,ix_] = s2xlib('ii','X2',ix_);
        pb.xnames{iv} = 'X2';
        [iv,ix_] = s2xlib('ii','X3',ix_);
        pb.xnames{iv} = 'X3';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2xlib('ii','OBJ',ig_);
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
            v_('-W') = -1.0*v_('W');
            v_('W**2') = v_('W')*v_('W');
            v_('-W**2') = -1.0*v_('W**2');
            [ig,ig_] = s2xlib('ii',['LO',int2str(I)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['LO',int2str(I)];
            iv = ix_('U');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_('X1');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            iv = ix_('X2');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('-W')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('-W');
            end
            iv = ix_('X3');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('-W**2')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('-W**2');
            end
            [ig,ig_] = s2xlib('ii',['UP',int2str(I)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['UP',int2str(I)];
            iv = ix_('U');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_('X1');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_('X2');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('W')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('W');
            end
            iv = ix_('X3');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('W**2')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('W**2');
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
        for I=v_('0'):v_('M')
            v_('RI') = I;
            v_('W') = v_('RI')*v_('H');
            v_('W') = v_('W')+v_('LOWER');
            v_('SINW') = sin(v_('W'));
            v_('-SINW') = -1.0*v_('SINW');
            pbm.gconst(ig_(['LO',int2str(I)])) = v_('-SINW');
            pbm.gconst(ig_(['UP',int2str(I)])) = v_('SINW');
        end
        pb.xlower = zeros(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons   = [1:length(pbm.congrps)];
        pb.pbclass = 'LLR2-AN-4-V';
        pb.x0          = zeros(pb.n,1);
        varargout{1} = pb;
        varargout{2} = pbm;

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

