function varargout = SCURLY30(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : SCURLY30
%    --------
% 
%    A banded function with semi-bandwidth 10 and
%    negative curvature near the starting point.
%    NB: scaled version of CURLY30
% 
%    Source: Nick Gould
% 
%    SIF input: Nick Gould, September 1997.
% 
%    classification = 'OUR2-AN-V-0'
% 
%    Number of variables
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'SCURLY30';

switch(action)

    case 'setup'

    pb.name      = 'SCURLY30';
    pb.sifpbname = 'SCURLY30';
    pbm.name     = 'SCURLY30';
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        if(nargin<2)
            v_('N') = 35;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
%       Alternative values for the SIF file parameters:
% IE N                   100            $-PARAMETER
% IE N                   1000           $-PARAMETER     original value
% IE N                   10000          $-PARAMETER
% IE N                   100000         $-PARAMETER
% IE N                   1000000        $-PARAMETER
        v_('K') = 30;
        v_('SCAL') = 12.0;
        v_('1') = 1;
        v_('N-K') = v_('N')-v_('K');
        v_('N-K+1') = 1+v_('N-K');
        v_('RN') = v_('N');
        v_('RN-1') = -1+v_('RN');
        v_('RN+1') = 1+v_('RN');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            v_('I-1') = -1+I;
            v_('RI-1') = v_('I-1');
            v_('RAT') = v_('RI-1')/v_('RN-1');
            v_('ARG') = v_('RAT')*v_('SCAL');
            v_(['S',int2str(I)]) = exp(v_('ARG'));
            [iv,ix_] = s2xlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('N-K')
            v_('I+K') = I+v_('K');
            for J=I:v_('I+K')
                [ig,ig_] = s2xlib('ii',['Q',int2str(I)],ig_);
                gtype{ig} = '<>';
                iv = ix_(['X',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_(['S',int2str(J)])+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_(['S',int2str(J)]);
                end
            end
        end
        for I=v_('N-K+1'):v_('N')
            for J=I:v_('N')
                [ig,ig_] = s2xlib('ii',['Q',int2str(I)],ig_);
                gtype{ig} = '<>';
                iv = ix_(['X',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_(['S',int2str(J)])+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_(['S',int2str(J)]);
                end
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
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        for I=v_('1'):v_('N')
            v_('RI') = I;
            v_('T') = v_('RI')/v_('RN+1');
            v_('T') = 0.0001*v_('T');
            v_('T') = v_('T')*v_(['S',int2str(I)]);
            pb.x0(ix_(['X',int2str(I)]),1) = v_('T');
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2xlib('ii','gP4',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for ig = 1:ngrp
            pbm.grftype{ig} = 'gP4';
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'OUR2-AN-V-0';
        varargout{1} = pb;
        varargout{2} = pbm;

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gP4'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        APB = 2.0e+1;
        varargout{1} = GVAR_*(GVAR_*(GVAR_^2-APB)-1.0e-1);
        if(nargout>1)
            g_ = 2.0e+0*GVAR_*(2.0e+0*GVAR_^2-APB)-1.0e-1;
            varargout{2} = g_;
            if(nargout>2)
                H_ = 1.2e+1*GVAR_^2-2.0e+0*APB;
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
