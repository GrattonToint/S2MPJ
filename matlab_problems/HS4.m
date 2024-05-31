function varargout = HS4(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HS4
%    *********
% 
%    Source: problem 4 in
%    W. Hock and K. Schittkowski,
%    "Test examples for nonlinear programming codes",
%    Lectures Notes in Economics and Mathematical Systems 187, Springer
%    Verlag, Heidelberg, 1981.
% 
%    SIF input: A.R. Conn March 1990
% 
%    classification = 'OBR2-AN-2-0'
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HS4';

switch(action)

    case 'setup'

    pb.name      = 'HS4';
    pb.sifpbname = 'HS4';
    pbm.name     = 'HS4';
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2xlib('ii','X1',ix_);
        pb.xnames{iv} = 'X1';
        [iv,ix_] = s2xlib('ii','X2',ix_);
        pb.xnames{iv} = 'X2';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2xlib('ii','G1',ig_);
        gtype{ig} = '<>';
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        pbm.gscale(ig,1) = 3.0;
        [ig,ig_] = s2xlib('ii','G2',ig_);
        gtype{ig} = '<>';
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
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        pbm.gconst(ig_('G1')) = -1.0;
        pb.xlower = zeros(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xlower(ix_('X1'),1) = 1.0;
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.x0(ix_('X1'),1) = 1.125;
        pb.x0(ix_('X2'),1) = 0.125;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2xlib('ii','gCUBE',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        ig = ig_('G1');
        pbm.grftype{ig} = 'gCUBE';
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 2.66;
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'OBR2-AN-2-0';
        varargout{1} = pb;
        varargout{2} = pbm;

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gCUBE'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = GVAR_*GVAR_*GVAR_;
        if(nargout>1)
            g_ = 3.0*GVAR_*GVAR_;
            varargout{2} = g_;
            if(nargout>2)
                H_ = 6.0*GVAR_;
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
