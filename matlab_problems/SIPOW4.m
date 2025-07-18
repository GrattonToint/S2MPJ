function varargout = SIPOW4(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : SIPOW4
%    *********
% 
%    This is a discretization of a one sided approximation problem of
%    approximating the function xi * xi * eta by a linear polynomial
%    on the boundary of a circle (xi - 0.5)**2 + (eta - 0.5)**2 = 0.5
% 
%    Source: problem 4 in
%    M. J. D. Powell,
%    "Log barrier methods for semi-infinite programming calculations"
%    Numerical Analysis Report DAMTP 1992/NA11, U. of Cambridge, UK.
% 
%    SIF input: A. R. Conn and Nick Gould, August 1993
% 
%    classification = 'C-CLLR2-AN-4-V'
% 
%    Problem variants: they are identified by the values of M (even)
% 
% IE M                   20 
% IE M                   100 
% IE M                   500 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 21 VI 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'SIPOW4';

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
        v_('M') = 2000;
        v_('1') = 1;
        v_('2') = 2;
        v_('M/2') = fix(v_('M')/v_('2'));
        v_('M/2+1') = 1+v_('M/2');
        v_('RM') = v_('M');
        v_('1/RM') = 1.0/v_('RM');
        v_('ONE') = 1.0;
        v_('HALF') = 0.5;
        v_('ROOTHALF') = sqrt(v_('HALF'));
        v_('PI/4') = atan(v_('ONE'));
        v_('2PI') = 8.0*v_('PI/4');
        v_('2PI/M') = v_('2PI')*v_('1/RM');
        for J=v_('1'):v_('M/2')
            v_('RJ') = J;
            v_('THETA') = v_('RJ')*v_('2PI/M');
            v_('PI/4-T') = v_('PI/4')-v_('THETA');
            v_('COS') = cos(v_('PI/4-T'));
            v_('SIN') = sin(v_('PI/4-T'));
            v_('RTC') = v_('COS')*v_('ROOTHALF');
            v_('RTS') = v_('SIN')*v_('ROOTHALF');
            v_('-RTC') = -1.0*v_('RTC');
            v_('-RTS') = -1.0*v_('RTS');
            v_(['XI',int2str(J)]) = v_('HALF')+v_('-RTC');
            v_(['ETA',int2str(J)]) = v_('HALF')+v_('-RTS');
        end
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        [iv,ix_] = s2mpjlib('ii','X1',ix_);
        pb.xnames{iv} = 'X1';
        [iv,ix_] = s2mpjlib('ii','X2',ix_);
        pb.xnames{iv} = 'X2';
        [iv,ix_] = s2mpjlib('ii','X3',ix_);
        pb.xnames{iv} = 'X3';
        [iv,ix_] = s2mpjlib('ii','X4',ix_);
        pb.xnames{iv} = 'X4';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X4');
        valA(end+1) = 1.0;
        for J=v_('1'):v_('M/2')
            [ig,ig_] = s2mpjlib('ii',['C',int2str(J)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['C',int2str(J)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_('X1');
            valA(end+1) = 1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_('X4');
            valA(end+1) = 1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_('X2');
            valA(end+1) = v_(['XI',int2str(J)]);
            irA(end+1)  = ig;
            icA(end+1)  = ix_('X3');
            valA(end+1) = v_(['ETA',int2str(J)]);
        end
        for J=v_('1'):v_('M/2')
            v_('J+') = v_('M/2')+J;
            [ig,ig_] = s2mpjlib('ii',['C',int2str(round(v_('J+')))],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['C',int2str(round(v_('J+')))];
            irA(end+1)  = ig;
            icA(end+1)  = ix_('X1');
            valA(end+1) = 1.0;
            [ig,ig_] = s2mpjlib('ii',['C',int2str(round(v_('J+')))],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['C',int2str(round(v_('J+')))];
            irA(end+1)  = ig;
            icA(end+1)  = ix_('X2');
            valA(end+1) = v_(['XI',int2str(J)]);
            [ig,ig_] = s2mpjlib('ii',['C',int2str(round(v_('J+')))],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['C',int2str(round(v_('J+')))];
            irA(end+1)  = ig;
            icA(end+1)  = ix_('X3');
            valA(end+1) = v_(['ETA',int2str(J)]);
        end
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
        for J=v_('1'):v_('M/2')
            v_('J+') = v_('M/2')+J;
            v_('XIXI') = v_(['XI',int2str(J)])*v_(['XI',int2str(J)]);
            v_('XIXIETA') = v_('XIXI')*v_(['ETA',int2str(J)]);
            pbm.gconst(ig_(['C',int2str(J)])) = v_('XIXIETA');
            pbm.gconst(ig_(['C',int2str(round(v_('J+')))])) = v_('XIXIETA');
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,'X1'))
            pb.x0(ix_('X1'),1) = -0.1;
        else
            pb.y0(find(pbm.congrps==ig_('X1')),1) = -0.1;
        end
        if(isKey(ix_,'X2'))
            pb.x0(ix_('X2'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('X2')),1) = 0.0;
        end
        if(isKey(ix_,'X3'))
            pb.x0(ix_('X3'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('X3')),1) = 0.0;
        end
        if(isKey(ix_,'X4'))
            pb.x0(ix_('X4'),1) = 1.2;
        else
            pb.y0(find(pbm.congrps==ig_('X4')),1) = 1.2;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
% LO SOLUTION            2.0704432D-1 ! m = 20
% LO SOLUTION            2.6110334D-1 ! m = 100
% LO SOLUTION            2.7060094D-1 ! m = 500
% LO SOLUTION            2.7236200D-1 ! m = 2000
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons   = [1:length(pbm.congrps)];
        pb.pbclass = 'C-CLLR2-AN-4-V';
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

