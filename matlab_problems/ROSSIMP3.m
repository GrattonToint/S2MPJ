function varargout = ROSSIMP3(action,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Problem : ROSSIMP3
%    The ever famous 2 variables Rosenbrock "banana valley" problem
%    This version uses  1 trivial group and 2 nonlinear elements.
%    classification = 'C-CSUR2-AN-2-0'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
persistent pbm;
name = 'ROSSIMP3';
switch(action)
    case 'setup'
        pb.name      = name;
        pbm.name     = name;
        pb.n = 2;
        pb.m = 0;
        pbm.objgrps  = [ 1 ];
        pb.xlower    = [ -Inf; -Inf ];
        pb.xupper    = [ +Inf; +Inf ];
        pb.x0        = [ -1.2,; 1.0 ];
        pbm.elftype  = { 'EROSNB'};
        pbm.elvar    = { [ 1, 2 ] };
        pbm.grelt{1} = [ 1 ];
        pb.objlower  = 0.0;
        pb.pbclass   = 'SUR2-AN-2-0';
        varargout{1} = pb;
        varargout{2} = pbm;
    case 'EROSNB'    %  The ever famous "banana" function!
        x  = varargin{1};
        varargout{1} = 100.0*(x(2)-x(1)*x(1))^2+(x(1)-1.0)^2;
        if(nargout>1)
            g_(1,1) = -400.0*(x(2)-x(1)*x(1))*x(1)+2.0*(x(1)-1.0);
            g_(2,1) = 200.0*(x(2)-x(1)*x(1));
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 1200.0*x(1)*x(1)+2.0;
                H_(1,2) = -400.0;
                H_(2,2) =  200.0;
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


