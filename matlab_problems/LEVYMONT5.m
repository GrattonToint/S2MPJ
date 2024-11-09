function varargout = LEVYMONT5(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : LEVYMONT5
%    *********
%    A global optimization example due to Levy & Montalvo 
%    This problem is one of the parameterised set LEVYMONT5-LEVYMONT10
% 
%    Source:  problem 5 in
% 
%    A. V. Levy and A. Montalvo
%    "The Tunneling Algorithm for the Global Minimization of Functions"
%    SIAM J. Sci. Stat. Comp. 6(1) 1985 15:29 
%    https://doi.org/10.1137/0906002
% 
%    SIF input: Nick Gould, August 2021
% 
%    classification = 'C-CSBR2-AY-2-0'
% 
%    N is the number of variables
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'LEVYMONT5';

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
        v_  = containers.Map('KeyType','char', 'ValueType', 'double');
        ix_ = containers.Map('KeyType','char', 'ValueType', 'double');
        ig_ = containers.Map('KeyType','char', 'ValueType', 'double');
        v_('N') = 2;
        v_('A') = 1.0;
        v_('K') = 10.0;
        v_('L') = 0.25;
        v_('C') = 0.75;
        v_('1') = 1;
        v_('2') = 2;
        v_('PI/4') = atan(1.0);
        v_('PI') = 4.0*v_('PI/4');
        v_('RN') = v_('N');
        v_('A-C') = v_('A')-v_('C');
        v_('PI/N') = v_('PI')/v_('RN');
        v_('KPI/N') = v_('K')*v_('PI/N');
        v_('ROOTKPI/N') = sqrt(v_('KPI/N'));
        v_('N/PI') = v_('RN')/v_('PI');
        v_('N/KPI') = v_('N/PI')/v_('K');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('N')
            [ig,ig_] = s2mpjlib('ii',['Q',int2str(I)],ig_);
            gtype{ig} = '<>';
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('L')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('L');
            end
            pbm.gscale(ig,1) = v_('N/PI');
            [ig,ig_] = s2mpjlib('ii',['N',int2str(I)],ig_);
            gtype{ig} = '<>';
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = ix_.Count;
        ngrp   = ig_.Count;
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        for I=v_('1'):v_('N')
            pbm.gconst(ig_(['Q',int2str(I)])) = v_('A-C');
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -10.0*ones(pb.n,1);
        pb.xupper = 10.0*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 8.0*ones(pb.n,1);
        pb.x0(ix_('X1'),1) = -8.0;
        pb.x0(ix_('X2'),1) = 8.0;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eS2',iet_);
        elftv{it}{1} = 'X';
        elftp{it}{1} = 'L';
        elftp{it}{2} = 'C';
        [it,iet_] = s2mpjlib( 'ii', 'ePS2',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Z';
        elftp{it}{1} = 'L';
        elftp{it}{2} = 'C';
        elftp{it}{3} = 'A';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        ename = ['E',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eS2';
        ielftype(ie) = iet_('eS2');
        ename = ['E',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-10.0,10.0,8.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        [~,posep] = ismember('L',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('L');
        ename = ['E',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        [~,posep] = ismember('C',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('C');
        for I=v_('2'):v_('N')
            v_('I-1') = I-v_('1');
            ename = ['E',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'ePS2';
            ielftype(ie) = iet_('ePS2');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-10.0,10.0,8.0);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('I-1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-10.0,10.0,8.0);
            posev = find(strcmp('Z',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('L',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('L');
            [~,posep] = ismember('C',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('C');
            [~,posep] = ismember('A',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('A');
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('N')
            ig = ig_(['Q',int2str(I)]);
            pbm.grftype{ig} = 'gL2';
            ig = ig_(['N',int2str(I)]);
            pbm.grftype{ig} = 'gL2';
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(I)]);
            pbm.grelw{ig}(posel) = v_('ROOTKPI/N');
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
%    Solution
% LO SOLTN               0.0
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'C-CSBR2-AY-2-0';
        pbm.objderlvl = 2;
        pb.objderlvl = pbm.objderlvl;
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb, pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm,pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end

% **********************
%  SET UP THE FUNCTION *
%  AND RANGE ROUTINES  *
% **********************

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'e_globs'

        pbm = varargin{1};
        pbm.efpar = [];
        pbm.efpar(1) = 4.0*atan(1.0e0);
        varargout{1} = pbm;

    case 'eS2'

        EV_  = varargin{1};
        iel_ = varargin{2};
        PIL = pbm.efpar(1)*pbm.elpar{iel_}(1);
        V = PIL*EV_(1)+pbm.efpar(1)*pbm.elpar{iel_}(2);
        SINV = sin(V);
        COSV = cos(V);
        varargout{1} = SINV;
        if(nargout>1)
            g_(1,1) = PIL*COSV;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = -PIL*PIL*SINV;
                varargout{3} = H_;
            end
        end

    case 'ePS2'

        EV_  = varargin{1};
        iel_ = varargin{2};
        PIL = pbm.efpar(1)*pbm.elpar{iel_}(1);
        U = pbm.elpar{iel_}(1)*EV_(2)+pbm.elpar{iel_}(2)-pbm.elpar{iel_}(3);
        V = PIL*EV_(1)+pbm.efpar(1)*pbm.elpar{iel_}(2);
        SINV = sin(V);
        COSV = cos(V);
        varargout{1} = U*SINV;
        if(nargout>1)
            g_(1,1) = PIL*U*COSV;
            g_(2,1) = pbm.elpar{iel_}(1)*SINV;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = -PIL*PIL*U*SINV;
                H_(1,2) = pbm.elpar{iel_}(1)*PIL*COSV;
                H_(2,1) = H_(1,2);
                H_(2,2) = 0.0;
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gL2'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = GVAR_*GVAR_;
        if(nargout>1)
            g_ = 2.0*GVAR_;
            varargout{2} = g_;
            if(nargout>2)
                H_ = 2.0;
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    case {'fx','fgx','fgHx','cx','cJx','cJHx','cIx','cIJx','cIJHx','cIJxv','fHxv',...
          'cJxv','cJtxv','cIJtxv','Lxy','Lgxy','LgHxy','LIxy','LIgxy','LIgHxy',...
          'LHxyv','LIHxyv'}

        if(isfield(pbm,'name')&&strcmp(pbm.name,name))
            pbm.has_globs = [1,0];
            [varargout{1:max(1,nargout)}] = s2mpjlib(action,pbm,varargin{:});
        else
            disp(['ERROR: please run ',name,' with action = setup'])
            [varargout{1:nargout}] = deal(NaN);
        end

    otherwise
        disp([' ERROR: action ',action,' unavailable for problem ',name,'.m'])
    end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

