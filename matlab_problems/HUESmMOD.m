function varargout = HUESmMOD(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HUESmMOD
%    *********
% 
%    Source: An inverse problem from astronomy,
%    reformulated as a convex quadratic program by
%    S. P. Hestis, SIAM Review 34 (1992) pp. 642-647.
% 
%    SIF input: Nick Gould, January 1993.
%    improvements by: Ruediger Franke (Ruediger.Franke@RZ.TU-Ilmenau.DE)
% 
%    classification = 'QLR2-MN-V-V'
% 
%    Number of variables
% 
%       Alternative values for the SIF file parameters:
% IE K                   10             $-PARAMETER
% IE K                   100            $-PARAMETER
% IE K                   1000           $-PARAMETER    original value
% IE K                   5000           $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HUESmMOD';

switch(action)

    case 'setup'

        pb.name      = name;
        pbm.name     = name;
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        if(nargin<2)
            v_('K') = 10;  %  SIF file default value
        else
            v_('K') = varargin{1};
        end
% IE K                   10000          $-PARAMETER
        v_('1') = 1;
        v_('RANGE') = 1.0;
        v_('3.0') = 3.0;
        v_('5.0') = 5.0;
        v_('RK') = v_('K');
        v_('DELTAX') = v_('RANGE')/v_('RK');
        v_('DELTAX2') = v_('DELTAX')*v_('DELTAX');
        v_('DELTAX3') = v_('DELTAX2')*v_('DELTAX');
        v_('DELTAX5') = v_('DELTAX3')*v_('DELTAX2');
        v_('DELTAX3/3') = v_('DELTAX3')/v_('3.0');
        v_('DELTAX5/5') = v_('DELTAX5')/v_('5.0');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('K')
            [iv,ix_] = s2mpjlib('ii',['M',int2str(I)],ix_);
            pb.xnames{iv} = ['M',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('K')
            [ig,ig_] = s2mpjlib('ii',['OBJ',int2str(I)],ig_);
            gtype{ig} = '<>';
        end
        for I=v_('1'):v_('K')
            v_('I-1') = -1+I;
            v_('RI') = I;
            v_('RI2') = v_('RI')*v_('RI');
            v_('RI3') = v_('RI2')*v_('RI');
            v_('RI5') = v_('RI3')*v_('RI2');
            v_('RI-1') = v_('I-1');
            v_('RI-12') = v_('RI-1')*v_('RI-1');
            v_('RI-13') = v_('RI-12')*v_('RI-1');
            v_('RI-15') = v_('RI-13')*v_('RI-12');
            v_('DIFF3') = v_('RI3')-v_('RI-13');
            v_('DIFF5') = v_('RI5')-v_('RI-15');
            v_('COEFF1') = v_('DIFF3')*v_('DELTAX3/3');
            v_('COEFF2') = v_('DIFF5')*v_('DELTAX5/5');
            [ig,ig_] = s2mpjlib('ii','E1',ig_);
            gtype{ig}  = '==';
            cnames{ig} = 'E1';
            iv = ix_(['M',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('COEFF1')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('COEFF1');
            end
            [ig,ig_] = s2mpjlib('ii','E2',ig_);
            gtype{ig}  = '==';
            cnames{ig} = 'E2';
            iv = ix_(['M',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('COEFF2')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('COEFF2');
            end
        end
        v_('RK') = v_('K');
        v_('1/RK') = 1.0/v_('RK');
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
        pbm.gconst(ig_('E1')) = 1835.2;
        pbm.gconst(ig_('E2')) = 909.8;
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 1.0*ones(pb.n,1);
        pb.y0 = 1.0*ones(pb.m,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSQ',iet_);
        elftv{it}{1} = 'U1';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('K')
            ename = ['E',int2str(I)];
            [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
            if(newelt)
                pbm.elftype{ie} = 'eSQ';
                ielftype(ie) = iet_('eSQ');
            end
            vname = ['M',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('U1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('K')
            ig = ig_(['OBJ',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('1/RK');
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
%    Solution
% LO SOLTN               0.0
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'QLR2-MN-V-V';
        varargout{1} = pb;
        varargout{2} = pbm;
% **********************
%  SET UP THE FUNCTION *
%  AND RANGE ROUTINES  *
% **********************

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eSQ'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(1);
        if(nargout>1)
            g_(1,1) = EV_(1)+EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0;
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

