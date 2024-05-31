function varargout = CLPLATEB(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : CLPLATEB
%    *********
% 
%    The clamped plate problem (Strang, Nocedal, Dax)
%    The problem comes from the discretization the following problem
%    in mechanics:  a plate is clamped on one edge and loaded on the
%    opposite side.  The plate is the unit square.
% 
%    The plate is clamped on its lower edge, by fixing the
%    corresponding variables to zero.
% 
%    In this version of the problem, the weight WGHT is distributed
%    equally along the upper edge, introducing a symmetry with respect
%    to the vertical axis.
% 
%    Source:
%    J. Nocedal,
%    "Solving large nonlinear systems of equations arising in mechanics",
%    Proceedings of the Cocoyoc Numerical Analysis Conference, Mexico,
%    pp. 132-141, 1981.
% 
%    SIF input: Ph. Toint, Dec 1989.
% 
%    classification = 'OXR2-MN-V-0'
% 
%    P is the number of points in one side of the unit square
%    The number of variables is P*P, of which (P-1)*(P-1) are free.
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'CLPLATEB';

switch(action)

    case 'setup'

    pb.name      = 'CLPLATEB';
    pb.sifpbname = 'CLPLATEB';
    pbm.name     = 'CLPLATEB';
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        if(nargin<2)
            v_('P') = 4;  %  SIF file default value
        else
            v_('P') = varargin{1};
        end
%       Alternative values for the SIF file parameters:
% IE P                   7              $-PARAMETER n = 49    original value
% IE P                   10             $-PARAMETER n = 100
% IE P                   23             $-PARAMETER n = 529
% IE P                   32             $-PARAMETER n = 1024
% IE P                   71             $-PARAMETER n = 5041
        v_('WGHT') = -0.1;
        v_('1') = 1;
        v_('2') = 2;
        v_('P-1') = -1+v_('P');
        v_('RP-1') = v_('P-1');
        v_('1/P-1') = 1.0/v_('RP-1');
        v_('DISW') = v_('WGHT')*v_('1/P-1');
        v_('P2') = v_('P')*v_('P');
        v_('RP2') = v_('P2');
        v_('HP2') = 0.5*v_('RP2');
        v_('1/HP2') = 1.0/v_('HP2');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for J=v_('1'):v_('P')
            for I=v_('1'):v_('P')
                [iv,ix_] = s2xlib('ii',['X',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['X',int2str(I),',',int2str(J)];
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('2'):v_('P')
            v_('I-1') = -1+I;
            for J=v_('2'):v_('P')
                v_('J-1') = -1+J;
                [ig,ig_] = s2xlib('ii',['A',int2str(I),',',int2str(J)],ig_);
                gtype{ig} = '<>';
                pbm.gscale(ig,1) = 2.0;
                iv = ix_(['X',int2str(I),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
                iv = ix_(['X',int2str(I),',',int2str(round(v_('J-1')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
                [ig,ig_] = s2xlib('ii',['B',int2str(I),',',int2str(J)],ig_);
                gtype{ig} = '<>';
                pbm.gscale(ig,1) = 2.0;
                iv = ix_(['X',int2str(I),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
                iv = ix_(['X',int2str(round(v_('I-1'))),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
                [ig,ig_] = s2xlib('ii',['C',int2str(I),',',int2str(J)],ig_);
                gtype{ig} = '<>';
                pbm.gscale(ig,1) = v_('1/HP2');
                iv = ix_(['X',int2str(I),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
                iv = ix_(['X',int2str(I),',',int2str(round(v_('J-1')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
                [ig,ig_] = s2xlib('ii',['D',int2str(I),',',int2str(J)],ig_);
                gtype{ig} = '<>';
                pbm.gscale(ig,1) = v_('1/HP2');
                iv = ix_(['X',int2str(I),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
                iv = ix_(['X',int2str(round(v_('I-1'))),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
            end
        end
        for J=v_('1'):v_('P')
            [ig,ig_] = s2xlib('ii','W',ig_);
            gtype{ig} = '<>';
            iv = ix_(['X',int2str(round(v_('P'))),',',int2str(J)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('DISW')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('DISW');
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
        for J=v_('1'):v_('P')
            pb.xlower(ix_(['X',int2str(round(v_('1'))),',',int2str(J)]),1) = 0.0;
            pb.xupper(ix_(['X',int2str(round(v_('1'))),',',int2str(J)]),1) = 0.0;
        end
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 0.0*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2xlib('ii','gL2',igt_);
        [it,igt_] = s2xlib('ii','gL4',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('2'):v_('P')
            for J=v_('2'):v_('P')
                ig = ig_(['A',int2str(I),',',int2str(J)]);
                pbm.grftype{ig} = 'gL2';
                ig = ig_(['B',int2str(I),',',int2str(J)]);
                pbm.grftype{ig} = 'gL2';
                ig = ig_(['C',int2str(I),',',int2str(J)]);
                pbm.grftype{ig} = 'gL4';
                ig = ig_(['D',int2str(I),',',int2str(J)]);
                pbm.grftype{ig} = 'gL4';
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'OXR2-MN-V-0';
        varargout{1} = pb;
        varargout{2} = pbm;

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gL2'

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

    case 'gL4'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = GVAR_^4;
        if(nargout>1)
            g_ = 4.0*GVAR_^3;
            varargout{2} = g_;
            if(nargout>2)
                H_ = 12.0*GVAR_^2;
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

