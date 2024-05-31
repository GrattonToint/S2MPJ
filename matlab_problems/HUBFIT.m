function varargout = HUBFIT(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HUBFIT
%    *********
%    Variable dimension full rank linear problem
%    An elementary fit using the Huber loss function
% 
%    Source:
%    A.R. Conn, N. Gould and Ph.L. Toint,
%    "The LANCELOT User's Manual",
%    Dept of Maths, FUNDP, 1991.
% 
%    SIF input: Ph. Toint, Jan 1991.
% 
%    classification = 'OLR2-AN-2-1'
% 
%    Data points
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HUBFIT';

switch(action)

    case 'setup'

    pb.name      = 'HUBFIT';
    pb.sifpbname = 'HUBFIT';
    pbm.name     = 'HUBFIT';
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        v_('X1') = 0.1;
        v_('X2') = 0.3;
        v_('X3') = 0.5;
        v_('X4') = 0.7;
        v_('X5') = 0.9;
        v_('Y1') = 0.25;
        v_('Y2') = 0.3;
        v_('Y3') = 0.625;
        v_('Y4') = 0.701;
        v_('Y5') = 1.0;
        v_('C') = 0.85;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2xlib('ii','a',ix_);
        pb.xnames{iv} = 'a';
        [iv,ix_] = s2xlib('ii','b',ix_);
        pb.xnames{iv} = 'b';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2xlib('ii','Obj1',ig_);
        gtype{ig} = '<>';
        iv = ix_('a');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('X1')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('X1');
        end
        iv = ix_('b');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        pbm.gscale(ig,1) = 2.0;
        [ig,ig_] = s2xlib('ii','Obj2',ig_);
        gtype{ig} = '<>';
        iv = ix_('a');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('X2')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('X2');
        end
        iv = ix_('b');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        pbm.gscale(ig,1) = 2.0;
        [ig,ig_] = s2xlib('ii','Obj3',ig_);
        gtype{ig} = '<>';
        iv = ix_('a');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('X3')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('X3');
        end
        iv = ix_('b');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        pbm.gscale(ig,1) = 2.0;
        [ig,ig_] = s2xlib('ii','Obj4',ig_);
        gtype{ig} = '<>';
        iv = ix_('a');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('X4')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('X4');
        end
        iv = ix_('b');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        pbm.gscale(ig,1) = 2.0;
        [ig,ig_] = s2xlib('ii','Obj5',ig_);
        gtype{ig} = '<>';
        iv = ix_('a');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('X5')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('X5');
        end
        iv = ix_('b');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        pbm.gscale(ig,1) = 2.0;
        [ig,ig_] = s2xlib('ii','Cons',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'Cons';
        iv = ix_('a');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('b');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
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
        pbm.gconst(ig_('Obj1')) = v_('Y1');
        pbm.gconst(ig_('Obj2')) = v_('Y2');
        pbm.gconst(ig_('Obj3')) = v_('Y3');
        pbm.gconst(ig_('Obj4')) = v_('Y4');
        pbm.gconst(ig_('Obj5')) = v_('Y5');
        pbm.gconst(ig_('Cons')) = v_('C');
        pb.xlower = zeros(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xlower(ix_('b')) = -Inf;
        pb.xupper(ix_('b'),1) = +Inf;
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2xlib('ii','gHUBER',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        ig = ig_('Obj1');
        pbm.grftype{ig} = 'gHUBER';
        ig = ig_('Obj2');
        pbm.grftype{ig} = 'gHUBER';
        ig = ig_('Obj3');
        pbm.grftype{ig} = 'gHUBER';
        ig = ig_('Obj4');
        pbm.grftype{ig} = 'gHUBER';
        ig = ig_('Obj5');
        pbm.grftype{ig} = 'gHUBER';
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons   = [1:length(pbm.congrps)];
        pb.pbclass = 'OLR2-AN-2-1';
        pb.x0          = zeros(pb.n,1);
        varargout{1} = pb;
        varargout{2} = pbm;

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gHUBER'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        HUBERK = 1.5;
        ABSA = abs(GVAR_);
        OUT = ABSA>HUBERK;
        if(OUT)
            FF = HUBERK*ABSA-0.5*HUBERK*HUBERK;
        end
        NEGOUT = OUT&&(GVAR_<0.0);
        POSOUT = OUT&&(GVAR_>=0.0);
        if(POSOUT)
            GG = HUBERK;
        end
        if(NEGOUT)
            GG = -HUBERK;
        end
        if(OUT)
            HH = 0.0;
        end
        if(~OUT)
            FF = 0.5*ABSA*ABSA;
        end
        if(~OUT)
            GG = GVAR_;
        end
        if(~OUT)
            HH = 1.0;
        end
        varargout{1} = FF;
        if(nargout>1)
            g_ = GG;
            varargout{2} = g_;
            if(nargout>2)
                H_ = HH;
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
