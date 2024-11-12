function varargout = ROSSIMP1(action,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Problem : ROSSIMP1
%    The ever famous 2 variables Rosenbrock "banana valley" problem.
%    This version is constructed from ROSENBR.m, removing all (here)
%    unnecessary statements automatically generated by S2MPJ.
%    classification = 'C-CSUR2-AN-2-0'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
persistent pbm;
name = 'ROSSIMP1';
switch(action)
    case 'setup'
        pb.name      = name;
        pbm.name     = name;
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(2,2);
        pbm.A(1,2) = 1.0;
        pbm.gscale(1) = 0.01;
        pbm.A(2,1) = 1.0;
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = 2;
        ngrp   = 2;
        pbm.objgrps = [ 1, 2 ];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst =[ 0.0; 1 ];
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = [ -Inf; -Inf ];
        pb.xupper = [ +Inf; +Inf ];
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = [ -1.2; 1.0 ];
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        pbm.elftype{1} = 'eSQ';
        pbm.elvar{1}   = [ 1 ];
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grftype = { 'gL2', 'gL2' };
        pbm.grelt{1} = [ 1 ];
        pbm.grelw{1} = [ -1.0 ];
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'SUR2-AN-2-0';
        varargout{1} = pb;
        varargout{2} = pbm;
    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%
    case 'eSQ'
        EV_  = varargin{1};
        varargout{1} = EV_(1)*EV_(1);
        if(nargout>1)
            g_(1,1) = EV_(1)+EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0;
                varargout{3} = H_;
            end
        end
    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%
    case 'gL2'
        GVAR_ = varargin{1};
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


