function varargout = DIAGIQB(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : DIAGIQB
%    *********
%    A variable dimension indefinite quadratic problem
%    with eigenvalues clusetered towards the bottom of the spectrum
% 
%    lambda_i = i^2/n - n/2 + 1/n, i = 1, ... , n
% 
%    Source: simple test for GALAHAD gltr/glrt
% 
%    SIF input: Nick Gould, Feb 2019, corrected May 2024
% 
%    classification = 'C-CQBR2-AN-V-0'
% 
%    Number of variables (variable)
% 
%       Alternative values for the SIF file parameters:
% IE N                   10             $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'DIAGIQB';

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
        if(nargs<1)
            v_('N') = 10;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
% IE N                   50             $-PARAMETER
% IE N                   100            $-PARAMETER
% IE N                   500            $-PARAMETER
% IE N                   1000           $-PARAMETER     original value
% IE N                   5000           $-PARAMETER
% IE N                   10000          $-PARAMETER
% IE N                   100000         $-PARAMETER
% IE N                   1000000        $-PARAMETER
% IE N                   10             $-PARAMETER
        v_('1') = 1;
        v_('ONE') = 1.0;
        v_('TWO') = 2.0;
        v_('RN') = v_('N');
        v_('RN/2') = v_('RN')/v_('TWO');
        v_('1/RN') = v_('ONE')/v_('RN');
        v_('SHIFT') = v_('1/RN')-v_('RN/2');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('N')
            [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
            gtype{ig} = '<>';
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = ix_.Count;
        ngrp   = ig_.Count;
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -100000.0*ones(pb.n,1);
        pb.xupper = 1000000.0*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 1.0*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% QUADRATIC %%%%%%%%%%%%%%%%%%%
        pbm.H = sparse( pb.n, pb.n );
        for I=v_('1'):v_('N')
            v_('RI') = I;
            v_('RI2') = v_('RI')*v_('RI');
            v_('RI2/RN') = v_('RI2')/v_('RN');
            v_('H') = v_('RI2/RN')+v_('SHIFT');
            ix1 = ix_(['X',int2str(I)]);
            ix2 = ix_(['X',int2str(I)]);
            pbm.H(ix1,ix2) = v_('H')+pbm.H(ix1,ix2);
            pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
%    Solution
% LO SOLTN               0.0
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'C-CQBR2-AN-V-0';
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

