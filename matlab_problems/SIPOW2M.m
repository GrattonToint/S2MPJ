function varargout = SIPOW2M(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : SIPOW2M
%    *********
% 
%    This is a discretization of a semi-infinite programming problem, of
%    minimizing the variable x_2 within a circle of radius 1. The circle
%    is replaced by a discrete set of equally-spaced supporting tangents.   
%    The symmetry in SIPOW1.SIF is imposed by replacing those constraints
%    by an alternative set.
% 
%    A modification proposed by Powell, section 6.
% 
%    Source: problem 2 - modified - in
%    M. J. D. Powell,
%    "Log barrier methods for semi-infinite programming calculations"
%    Numerical Analysis Report DAMTP 1992/NA11, U. of Cambridge, UK.
% 
%    SIF input: A. R. Conn and Nick Gould, August 1993
% 
%    classification = 'C-CLLR2-AN-2-V'
% 
%    Problem variants: they are identified by the values of M (even)
% 
% IE M                   20 
% IE M                   100 
% IE M                   500 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'SIPOW2M';

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
        v_('PI/4') = atan(v_('ONE'));
        v_('4PI') = 16.0*v_('PI/4');
        v_('4PI/M') = v_('4PI')*v_('1/RM');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        [iv,ix_] = s2mpjlib('ii','X1',ix_);
        pb.xnames{iv} = 'X1';
        [iv,ix_] = s2mpjlib('ii','X2',ix_);
        pb.xnames{iv} = 'X2';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2');
        valA(end+1) = 1.0;
        for J=v_('1'):v_('M/2')
            v_('RJ') = J;
            v_('RJ+1/2') = 0.5+v_('RJ');
            v_('4PIJ+/M') = v_('4PI/M')*v_('RJ+1/2');
            v_('COS') = cos(v_('4PIJ+/M'));
            v_('SIN') = sin(v_('4PIJ+/M'));
            [ig,ig_] = s2mpjlib('ii',['C',int2str(J)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['C',int2str(J)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_('X1');
            valA(end+1) = v_('COS');
            irA(end+1)  = ig;
            icA(end+1)  = ix_('X2');
            valA(end+1) = v_('SIN');
        end
        for J=v_('M/2+1'):v_('M')
            [ig,ig_] = s2mpjlib('ii',['C',int2str(J)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['C',int2str(J)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_('X1');
            valA(end+1) = 1.0;
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
        for J=v_('1'):v_('M')
            pbm.gconst(ig_(['C',int2str(J)])) = -1.0;
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,'X1'))
            pb.x0(ix_('X1'),1) = 0.8;
        else
            pb.y0(find(pbm.congrps==ig_('X1')),1) = 0.8;
        end
        if(isKey(ix_,'X2'))
            pb.x0(ix_('X2'),1) = 0.5;
        else
            pb.y0(find(pbm.congrps==ig_('X2')),1) = 0.5;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
% LO SOLUTION            -1.0
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons   = [1:length(pbm.congrps)];
        pb.pbclass = 'C-CLLR2-AN-2-V';
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

