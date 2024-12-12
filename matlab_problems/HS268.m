function varargout = HS268(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HS268
%    *********
% 
%    A quadratic programming problem.
% 
%    Source:
%    K. Schittkowski
%    "More Test Examples for Nonlinear Programming Codes"
%    Springer Verlag, Berlin, Lecture notes in economics and 
%    mathematical systems, volume 282, 1987
% 
%    SIF input: Michel Bierlaire and Annick Sartenaer, October 1992.
%               minor correction by Ph. Shott, Jan 1995.
% 
%    classification = 'C-CQLR2-AN-5-5'
% 
%   the number of functions
% 
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HS268';

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
        v_('D1,1') = 10197.0;
        v_('D1,2') = -12454.0;
        v_('D1,3') = -1013.0;
        v_('D1,4') = 1948.0;
        v_('D1,5') = 329.0;
        v_('D2,1') = -12454.0;
        v_('D2,2') = 20909.0;
        v_('D2,3') = -1733.0;
        v_('D2,4') = -4914.0;
        v_('D2,5') = -186.0;
        v_('D3,1') = -1013.0;
        v_('D3,2') = -1733.0;
        v_('D3,3') = 1755.0;
        v_('D3,4') = 1089.0;
        v_('D3,5') = -174.0;
        v_('D4,1') = 1948.0;
        v_('D4,2') = -4914.0;
        v_('D4,3') = 1089.0;
        v_('D4,4') = 1515.0;
        v_('D4,5') = -22.0;
        v_('D5,1') = 329.0;
        v_('D5,2') = -186.0;
        v_('D5,3') = -174.0;
        v_('D5,4') = -22.0;
        v_('D5,5') = 27.0;
        v_('B1') = -9170.0;
        v_('B2') = 17099.0;
        v_('B3') = -2271.0;
        v_('B4') = -4336.0;
        v_('B5') = -43.0;
        v_('1') = 1;
        v_('5') = 5;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for I=v_('1'):v_('5')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        [ig,ig_] = s2mpjlib('ii','NONL',ig_);
        gtype{ig} = '<>';
        for I=v_('1'):v_('5')
            [ig,ig_] = s2mpjlib('ii','LINEAR',ig_);
            gtype{ig} = '<>';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(I)]);
            valA(end+1) = v_(['B',int2str(I)]);
        end
        [ig,ig_] = s2mpjlib('ii','LINEAR',ig_);
        gtype{ig} = '<>';
        pbm.gscale(ig,1) = -0.5;
        [ig,ig_] = s2mpjlib('ii','C1',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'C1';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X3');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X4');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X5');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','C2',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'C2';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1');
        valA(end+1) = 10.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2');
        valA(end+1) = 10.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X3');
        valA(end+1) = -3.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X4');
        valA(end+1) = 5.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X5');
        valA(end+1) = 4.0;
        [ig,ig_] = s2mpjlib('ii','C3',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'C3';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1');
        valA(end+1) = -8.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X3');
        valA(end+1) = -2.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X4');
        valA(end+1) = -5.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X5');
        valA(end+1) = 3.0;
        [ig,ig_] = s2mpjlib('ii','C4',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'C4';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1');
        valA(end+1) = 8.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X3');
        valA(end+1) = 2.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X4');
        valA(end+1) = 5.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X5');
        valA(end+1) = -3.0;
        [ig,ig_] = s2mpjlib('ii','C5',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'C5';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1');
        valA(end+1) = -4.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2');
        valA(end+1) = -2.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X3');
        valA(end+1) = 3.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X4');
        valA(end+1) = -5.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X5');
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
        pbm.gconst(ig_('C1')) = -5.0;
        pbm.gconst(ig_('C2')) = 20.0;
        pbm.gconst(ig_('C3')) = -40.0;
        pbm.gconst(ig_('C4')) = 11.0;
        pbm.gconst(ig_('C5')) = -30.0;
        pbm.gconst(ig_('NONL')) = -14463.0;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 1.0*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'ePROD',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        elftp{it}{1} = 'D';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('1'):v_('5')
            for J=v_('1'):v_('5')
                ename = ['E',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'ePROD';
                ielftype(ie) = iet_('ePROD');
                vname = ['X',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
                posev = find(strcmp('X',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
                posev = find(strcmp('Y',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                [~,posep] = ismember('D',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['D',int2str(I),',',int2str(J)]);
            end
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('5')
            for J=v_('1'):v_('5')
                ig = ig_('NONL');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['E',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
            end
        end
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CQLR2-AN-5-5';
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


    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'ePROD'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = pbm.elpar{iel_}(1)*EV_(1)*EV_(2);
        if(nargout>1)
            g_(1,1) = pbm.elpar{iel_}(1)*EV_(2);
            g_(2,1) = pbm.elpar{iel_}(1)*EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 0.0;
                H_(2,2) = 0.0;
                H_(1,2) = pbm.elpar{iel_}(1);
                H_(2,1) = H_(1,2);
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

