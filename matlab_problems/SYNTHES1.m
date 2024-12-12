function varargout = SYNTHES1(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : SYNTHES1
%    *********
% 
%    Source: Test problem 1 (Synthesis of processing system) in 
%    M. Duran & I.E. Grossmann,
%    "An outer approximation algorithm for a class of mixed integer nonlinear
%     programs", Mathematical Programming 36, pp. 307-339, 1986.
% 
%    SIF input: S. Leyffer, October 1997
% 
%    classification = 'C-COOR2-AN-6-6'
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'SYNTHES1';

switch(action)

    case {'setup','setup_redprec'}

        if(strcmp(action,'setup_redprec'))
            if(isfield(pbm,'ndigs'))
                rmfield(pbm,'ndigs');
            end
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
        v_('1') = 1;
        v_('3') = 3;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for I=v_('1'):v_('3')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        for I=v_('1'):v_('3')
            [iv,ix_] = s2mpjlib('ii',['Y',int2str(I)],ix_);
            pb.xnames{iv} = ['Y',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Y1');
        valA(end+1) = 5.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Y2');
        valA(end+1) = 6.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Y3');
        valA(end+1) = 8.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1');
        valA(end+1) = 10.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X3');
        valA(end+1) = -7.0;
        [ig,ig_] = s2mpjlib('ii','N1',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'N1';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X3');
        valA(end+1) = -0.8;
        [ig,ig_] = s2mpjlib('ii','N2',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'N2';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X3');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Y3');
        valA(end+1) = -2.0;
        [ig,ig_] = s2mpjlib('ii','L3',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'L3';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','L4',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'L4';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Y1');
        valA(end+1) = -2.0;
        [ig,ig_] = s2mpjlib('ii','L5',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'L5';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Y2');
        valA(end+1) = -2.0;
        [ig,ig_] = s2mpjlib('ii','L6',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'L6';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Y1');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Y2');
        valA(end+1) = 1.0;
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = ix_.Count;
        ngrp   = ig_.Count;
        legrps = find(strcmp(gtype,'<='));
        eqgrps = find(strcmp(gtype,'=='));
        gegrps = find(strcmp(gtype,'>='));
        pb.nle = length(legrps);
        pb.neq = length(eqgrps);
        pb.nge = length(gegrps);
        pb.m   = pb.nle+pb.neq+pb.nge;
        pbm.congrps = [ legrps, eqgrps, gegrps ];
        [pb.cnames{1:pb.m}] = deal(cnames{pbm.congrps});
        pb.nob = ngrp-pb.m;
        pbm.objgrps = find(strcmp(gtype,'<>'));
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        pbm.gconst(ig_('OBJ')) = -10.0;
        pbm.gconst(ig_('N2')) = -2.0;
        pbm.gconst(ig_('L6')) = 1.0;
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xupper(ix_('X1')) = 2.0;
        pb.xupper(ix_('X2')) = 2.0;
        pb.xupper(ix_('X3')) = 1.0;
        pb.xupper(ix_('Y1')) = 1.0;
        pb.xupper(ix_('Y2')) = 1.0;
        pb.xupper(ix_('Y3')) = 1.0;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eLOGXP1',iet_);
        elftv{it}{1} = 'X';
        [it,iet_] = s2mpjlib( 'ii', 'eLOGDIFF',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = 'LOGX2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eLOGXP1';
        ielftype(ie) = iet_('eLOGXP1');
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'LOGX1X2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eLOGDIFF';
        ielftype(ie) = iet_('eLOGDIFF');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        ig = ig_('OBJ');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('LOGX2');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -18.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('LOGX1X2');
        pbm.grelw{ig}(posel) = -19.2;
        ig = ig_('N1');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('LOGX2');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 0.8;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('LOGX1X2');
        pbm.grelw{ig}(posel) = 0.96;
        ig = ig_('N2');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('LOGX2');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('LOGX1X2');
        pbm.grelw{ig}(posel) = 1.2;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-COOR2-AN-6-6';
        pb.x0          = zeros(pb.n,1);
        pbm.objderlvl = 2;
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2];
        pb.conderlvl  = pbm.conderlvl;
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

    case 'eLOGXP1'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = log(EV_(1)+1.0);
        if(nargout>1)
            g_(1,1) = 1.0/(EV_(1)+1.0);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = -1.0/(EV_(1)+1.0)^2;
                varargout{3} = H_;
            end
        end

    case 'eLOGDIFF'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = log(EV_(1)-EV_(2)+1.0);
        DX = 1.0/(EV_(1)-EV_(2)+1.0);
        DXDX = -1.0/(EV_(1)-EV_(2)+1.0)^2;
        if(nargout>1)
            g_(1,1) = DX;
            g_(2,1) = -DX;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = DXDX;
                H_(1,2) = -DXDX;
                H_(2,1) = H_(1,2);
                H_(2,2) = DXDX;
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    case {'fx','fgx','fgHx','cx','cJx','cJHx','cIx','cIJx','cIJHx','cIJxv','fHxv',...
          'cJxv','cJtxv','cIJtxv','Lxy','Lgxy','LgHxy','LIxy','LIgxy','LIgHxy',...
          'LHxyv','LIHxyv'}

        if(isfield(pbm,'name')&&strcmp(pbm.name,name))
            pbm.has_globs = [0,0];
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

