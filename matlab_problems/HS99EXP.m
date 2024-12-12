function varargout = HS99EXP(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HS99EXP
%    *********
% 
%    Source: an expanded form of problem 99 in
%    W. Hock and K. Schittkowski,
%    "Test examples for nonlinear programming codes",
%    Lectures Notes in Economics and Mathematical Systems 187, Springer
%    Verlag, Heidelberg, 1981.
% 
%    SIF input: Ph. Toint, April 1991.
% 
%    classification = 'C-COOR2-AN-31-21'
% 
%    Problem data
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HS99EXP';

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
        v_('T1') = 0.0;
        v_('T2') = 25.0;
        v_('T3') = 50.0;
        v_('T4') = 100.0;
        v_('T5') = 150.0;
        v_('T6') = 200.0;
        v_('T7') = 290.0;
        v_('T8') = 380.0;
        v_('A1') = 0.0;
        v_('A2') = 50.0;
        v_('A3') = 50.0;
        v_('A4') = 75.0;
        v_('A5') = 75.0;
        v_('A6') = 75.0;
        v_('A7') = 100.0;
        v_('A8') = 100.0;
        v_('B') = 32.0;
        v_('1') = 1;
        v_('2') = 2;
        v_('7') = 7;
        v_('8') = 8;
        for I=v_('2'):v_('8')
            v_('I-1') = -1+I;
            v_(['DT',int2str(I)]) =...
                  v_(['T',int2str(I)])-v_(['T',int2str(round(v_('I-1')))]);
            v_('DTISQ') = v_(['DT',int2str(I)])*v_(['DT',int2str(I)]);
            v_(['DT',int2str(I)]) = 0.5*v_('DTISQ');
        end
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for I=v_('1'):v_('7')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
            [iv,ix_] = s2mpjlib('ii',['R',int2str(I)],ix_);
            pb.xnames{iv} = ['R',int2str(I)];
            [iv,ix_] = s2mpjlib('ii',['Q',int2str(I)],ix_);
            pb.xnames{iv} = ['Q',int2str(I)];
            [iv,ix_] = s2mpjlib('ii',['S',int2str(I)],ix_);
            pb.xnames{iv} = ['S',int2str(I)];
        end
        [iv,ix_] = s2mpjlib('ii',['R',int2str(round(v_('8')))],ix_);
        pb.xnames{iv} = ['R',int2str(round(v_('8')))];
        [iv,ix_] = s2mpjlib('ii',['Q',int2str(round(v_('8')))],ix_);
        pb.xnames{iv} = ['Q',int2str(round(v_('8')))];
        [iv,ix_] = s2mpjlib('ii',['S',int2str(round(v_('8')))],ix_);
        pb.xnames{iv} = ['S',int2str(round(v_('8')))];
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        for I=v_('2'):v_('8')
            v_('I-1') = -1+I;
            [ig,ig_] = s2mpjlib('ii',['R',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['R',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['R',int2str(I)]);
            valA(end+1) = -1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['R',int2str(round(v_('I-1')))]);
            valA(end+1) = 1.0;
            [ig,ig_] = s2mpjlib('ii',['Q',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['Q',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['Q',int2str(I)]);
            valA(end+1) = -1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['Q',int2str(round(v_('I-1')))]);
            valA(end+1) = 1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['S',int2str(round(v_('I-1')))]);
            valA(end+1) = v_(['DT',int2str(I)]);
            [ig,ig_] = s2mpjlib('ii',['S',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['S',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['S',int2str(I)]);
            valA(end+1) = -1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['S',int2str(round(v_('I-1')))]);
            valA(end+1) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['R',int2str(round(v_('8')))]);
        valA(end+1) = 1.0;
        pbm.gscale(ig,1) = -1.0;
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
        for I=v_('2'):v_('7')
            v_('RHS') = v_(['DT',int2str(I)])*v_('B');
            pbm.gconst(ig_(['Q',int2str(I)])) = v_('RHS');
            v_('RHS') = v_(['DT',int2str(I)])*v_('B');
            pbm.gconst(ig_(['S',int2str(I)])) = v_('RHS');
        end
        pbm.gconst(ig_(['Q',int2str(round(v_('8')))])) = 100000.0;
        pbm.gconst(ig_(['S',int2str(round(v_('8')))])) = 1000.0;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower(ix_(['R',int2str(round(v_('1')))]),1) = 0.0;
        pb.xupper(ix_(['R',int2str(round(v_('1')))]),1) = 0.0;
        pb.xlower(ix_(['Q',int2str(round(v_('1')))]),1) = 0.0;
        pb.xupper(ix_(['Q',int2str(round(v_('1')))]),1) = 0.0;
        pb.xlower(ix_(['S',int2str(round(v_('1')))]),1) = 0.0;
        pb.xupper(ix_(['S',int2str(round(v_('1')))]),1) = 0.0;
        for I=v_('1'):v_('7')
            pb.xlower(ix_(['X',int2str(I)]),1) = 0.0;
            pb.xupper(ix_(['X',int2str(I)])) = 1.58;
        end
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        for I=v_('1'):v_('7')
            pb.x0(ix_(['X',int2str(I)]),1) = 0.5;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSN',iet_);
        elftv{it}{1} = 'X';
        [it,iet_] = s2mpjlib( 'ii', 'eCS',iet_);
        elftv{it}{1} = 'X';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('7')
            ename = ['SNX',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSN';
            ielftype(ie) = iet_('eSN');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['CSX',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eCS';
            ielftype(ie) = iet_('eCS');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        ig = ig_('OBJ');
        pbm.grftype{ig} = 'gL2';
        for I=v_('2'):v_('8')
            v_('I-1') = -1+I;
            v_('W') = v_(['A',int2str(I)])*v_(['DT',int2str(I)]);
            ig = ig_(['R',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['CSX',int2str(round(v_('I-1')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('W');
            ig = ig_(['S',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['SNX',int2str(round(v_('I-1')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('W');
            v_('W') = v_(['A',int2str(I)])*v_(['DT',int2str(I)]);
            ig = ig_(['Q',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['SNX',int2str(round(v_('I-1')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('W');
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               -831079892.0
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-COOR2-AN-31-21';
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

    case 'eSN'

        EV_  = varargin{1};
        iel_ = varargin{2};
        SNX = sin(EV_(1));
        varargout{1} = SNX;
        if(nargout>1)
            g_(1,1) = cos(EV_(1));
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = -SNX;
                varargout{3} = H_;
            end
        end

    case 'eCS'

        EV_  = varargin{1};
        iel_ = varargin{2};
        CSX = cos(EV_(1));
        varargout{1} = CSX;
        if(nargout>1)
            g_(1,1) = -sin(EV_(1));
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = -CSX;
                varargout{3} = H_;
            end
        end

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

