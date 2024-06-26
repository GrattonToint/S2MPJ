function varargout = READING2(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : READING2
%    *********
% 
%    A linear optimal control problem from Nancy Nichols
%    with a given initial condition.
%    This problem arises in tide modelling.
% 
%    Source:
%    S. Lyle and N.K. Nichols,
%    "Numerical Methods for Optimal Control Problems with State Constraints",
%    Numerical Analysis Report 8/91, Dept of Mathematics, 
%    University of Reading, UK.
% 
%    SIF input: Nick Gould, July 1991.
% 
%    classification = 'LLR2-MN-V-V'
% 
%    Number of discretized points in [0,1]
% 
%       Alternative values for the SIF file parameters:
% IE N                   2              $-PARAMETER
% IE N                   50             $-PARAMETER
% IE N                   100            $-PARAMETER     original value
% IE N                   500            $-PARAMETER
% IE N                   1000           $-PARAMETER
% IE N                   2000           $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'READING2';

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
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        if(nargs<1)
            v_('N') = 10;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
% IE N                   5000           $-PARAMETER
        v_('PI') = 3.1415926535;
        v_('2PI') = 2.0*v_('PI');
        v_('PI**2') = v_('PI')*v_('PI');
        v_('8PI**2') = 8.0*v_('PI**2');
        v_('1/8PI**2') = 1.0/v_('8PI**2');
        v_('A') = 0.07716;
        v_('1/A') = 1.0/v_('A');
        v_('1/2A') = 2.0*v_('1/A');
        v_('2A') = 2.0*v_('A');
        v_('-2A') = -1.0*v_('2A');
        v_('-1/2A') = 1.0/v_('-2A');
        v_('N-1') = -1+v_('N');
        v_('RN') = v_('N');
        v_('H') = 1.0/v_('RN');
        v_('2/H') = 2.0*v_('RN');
        v_('H/2') = 0.5*v_('H');
        v_('-H/2') = -1.0*v_('H/2');
        v_('1/H') = 1.0*v_('RN');
        v_('-1/H') = -1.0*v_('RN');
        v_('H/8PI**2') = v_('1/8PI**2')*v_('H');
        v_('0') = 0;
        v_('1') = 1;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('0'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X1u',int2str(I)],ix_);
            pb.xnames{iv} = ['X1u',int2str(I)];
            [iv,ix_] = s2mpjlib('ii',['X2u',int2str(I)],ix_);
            pb.xnames{iv} = ['X2u',int2str(I)];
            [iv,ix_] = s2mpjlib('ii',['U',int2str(I)],ix_);
            pb.xnames{iv} = ['U',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('N')
            v_('RI') = I;
            v_('TI') = v_('RI')*v_('H');
            v_('2PITI') = v_('2PI')*v_('TI');
            v_('CTI') = cos(v_('2PITI'));
            v_('-CCTI') = v_('CTI')*v_('-H/2');
            v_('I-1') = -1+I;
            v_('RI-1') = v_('I-1');
            v_('TI-1') = v_('RI-1')*v_('H');
            v_('2PITI-1') = v_('2PI')*v_('TI-1');
            v_('CTI-1') = cos(v_('2PITI-1'));
            v_('-CCTI-1') = v_('CTI-1')*v_('-H/2');
            [ig,ig_] = s2mpjlib('ii','COST',ig_);
            gtype{ig} = '<>';
            iv = ix_(['X1u',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('-CCTI')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('-CCTI');
            end
            iv = ix_(['X1u',int2str(round(v_('I-1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('-CCTI-1')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('-CCTI-1');
            end
            iv = ix_(['U',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('H/8PI**2')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('H/8PI**2');
            end
            iv = ix_(['U',int2str(round(v_('I-1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('H/8PI**2')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('H/8PI**2');
            end
        end
        for I=v_('1'):v_('N')
            v_('I-1') = -1+I;
            [ig,ig_] = s2mpjlib('ii',['C1u',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['C1u',int2str(I)];
            iv = ix_(['X1u',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('1/H')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('1/H');
            end
            iv = ix_(['X1u',int2str(round(v_('I-1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('-1/H')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('-1/H');
            end
            iv = ix_(['X2u',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -0.5+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -0.5;
            end
            iv = ix_(['X2u',int2str(round(v_('I-1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -0.5+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -0.5;
            end
            [ig,ig_] = s2mpjlib('ii',['C2u',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['C2u',int2str(I)];
            iv = ix_(['X2u',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('1/H')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('1/H');
            end
            iv = ix_(['X2u',int2str(round(v_('I-1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('-1/H')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('-1/H');
            end
            iv = ix_(['U',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -0.5+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -0.5;
            end
            iv = ix_(['U',int2str(round(v_('I-1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -0.5+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -0.5;
            end
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
        pbm.congrps = [ legrps, eqgrps, gegrps ];
        [pb.cnames{1:pb.m}] = deal(cnames{pbm.congrps});
        pb.nob = ngrp-pb.m;
        pbm.objgrps = find(strcmp(gtype,'<>'));
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xlower(ix_(['X1u',int2str(round(v_('0')))]),1) = 0.0;
        pb.xupper(ix_(['X1u',int2str(round(v_('0')))]),1) = 0.0;
        pb.xlower(ix_(['X2u',int2str(round(v_('0')))]),1) = 0.0;
        pb.xupper(ix_(['X2u',int2str(round(v_('0')))]),1) = 0.0;
        for I=v_('1'):v_('N')
            pb.xlower(ix_(['X1u',int2str(I)])) = -Inf;
            pb.xupper(ix_(['X1u',int2str(I)]),1) = +Inf;
            pb.xlower(ix_(['X2u',int2str(I)]),1) = -0.125;
            pb.xupper(ix_(['X2u',int2str(I)])) = 0.125;
        end
        for I=v_('0'):v_('N')
            pb.xlower(ix_(['U',int2str(I)]),1) = -1.0;
            pb.xupper(ix_(['U',int2str(I)])) = 1.0;
        end
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons   = [1:length(pbm.congrps)];
        pb.pbclass = 'LLR2-MN-V-V';
        pb.x0          = zeros(pb.n,1);
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb,  pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm, pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
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

