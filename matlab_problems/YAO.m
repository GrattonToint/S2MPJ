function varargout = YAO(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem:
%    ********
%    A linear least-sqaure problem with k-convex constraints
%       min (1/2) || f(t) - x ||^2
%    subject to the constraints
%       _ 
%       \/_k  x  >=  0,
%    where  f(t) and  x  are vectors in (n+k)-dimensional space.
% 
%    We choose f(t) = sin(t), x(1) >= 0.08 and fix x(n+i) = 0
% 
%    SIF input: Aixiang Yao, Virginia Tech., May 1995
%               modifications by Nick Gould
% 
%    classification = 'QLR2-AN-V-V'
% 
%   Number of discretization points
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'YAO';

switch(action)

    case 'setup'

    pb.name      = 'YAO';
    pb.sifpbname = 'YAO';
    pbm.name     = 'YAO';
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        if(nargin<2)
            v_('P') = 20;  %  SIF file default value
        else
            v_('P') = varargin{1};
        end
%       Alternative values for the SIF file parameters:
% IE P                   200            $-PARAMETER
% IE P                   2000           $-PARAMETER
        if(nargin<3)
            v_('k') = 2;  %  SIF file default value
        else
            v_('k') = varargin{2};
        end
% IE k                   3              $-PARAMETER
% IE k                   4              $-PARAMETER
        v_('1') = 1;
        v_('2') = 2;
        v_('P+1') = v_('P')+v_('1');
        v_('P+k') = v_('P')+v_('k');
        v_('RP') = v_('P+k');
        v_('OVP') = 1.0/v_('RP');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for i=v_('1'):v_('P+k')
            [iv,ix_] = s2xlib('ii',['X',int2str(i)],ix_);
            pb.xnames{iv} = ['X',int2str(i)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for i=v_('1'):v_('P+k')
            [ig,ig_] = s2xlib('ii',['S',int2str(i)],ig_);
            gtype{ig} = '<>';
            iv = ix_(['X',int2str(i)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            pbm.gscale(ig,1) = 2.0;
        end
        for i=v_('1'):v_('P')
            v_('i+1') = 1+i;
            [ig,ig_] = s2xlib('ii',['B',int2str(i)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['B',int2str(i)];
            iv = ix_(['X',int2str(i)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['X',int2str(round(v_('i+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -2.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -2.0;
            end
            v_('i+2') = 2+i;
            iv = ix_(['X',int2str(round(v_('i+2')))]);
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
        for i=v_('1'):v_('P+k')
            v_('Ri') = i;
            v_('iOVP') = v_('Ri')*v_('OVP');
            v_('SINI') = sin(v_('iOVP'));
            pbm.gconst(ig_(['S',int2str(i)])) = v_('SINI');
        end
        pb.xlower = zeros(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower(ix_(['X',int2str(round(v_('1')))]),1) = 0.08;
        for i=v_('P+1'):v_('P+k')
            pb.xlower(ix_(['X',int2str(i)]),1) = 0.0;
            pb.xupper(ix_(['X',int2str(i)]),1) = 0.0;
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2xlib('ii','gSQ',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for i=v_('1'):v_('P+k')
            ig = ig_(['S',int2str(i)]);
            pbm.grftype{ig} = 'gSQ';
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons   = [1:length(pbm.congrps)];
        pb.pbclass = 'QLR2-AN-V-V';
        pb.x0          = zeros(pb.n,1);
        varargout{1} = pb;
        varargout{2} = pbm;

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gSQ'

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

