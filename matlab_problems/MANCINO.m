function varargout = MANCINO(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : MANCINO
%    *********
% 
%    Mancino's function with variable dimension.
% 
%    Source:
%    E. Spedicato,
%    "Computational experience with quasi-Newton algorithms for
%    minimization problems of moderate size",
%    Report N-175, CISE, Milano, 1975.
% 
%    See also Buckley #51 (p. 72), Schittkowski #391 (for N = 30)
% 
%    SIF input: Ph. Toint, Dec 1989.
%               correction by Ph. Shott, January, 1995.
%               correction by S. Gratton & Ph. Toint, May 2024
% 
%    classification = 'C-CSUR2-AN-V-0'
% 
%    The definitions
%      s_{i,j} = \sin \log v_{i,j}   and s_{i,j} = \cos \log v_{i,j}
%    have been used.  It seems that the additional exponent ALPHA
%    in Buckley is a typo.
% 
%    Number of variables
% 
%       Alternative values for the SIF file parameters:
% IE N                   10             $-PARAMETER
% IE N                   20             $-PARAMETER
% IE N                   30             $-PARAMETER Schittkowski #391
% IE N                   50             $-PARAMETER
% IE N                   100            $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'MANCINO';

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
        if(nargs<1)
            v_('N') = 10;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
        if(nargs<2)
            v_('ALPHA') = 5;  %  SIF file default value
        else
            v_('ALPHA') = varargin{2};
        end
        if(nargs<3)
            v_('BETA') = 14.0;  %  SIF file default value
        else
            v_('BETA') = varargin{3};
        end
        if(nargs<4)
            v_('GAMMA') = 3;  %  SIF file default value
        else
            v_('GAMMA') = varargin{4};
        end
        v_('RALPHA') = v_('ALPHA');
        v_('RN') = v_('N');
        v_('N-1') = -1+v_('N');
        v_('RN-1') = v_('N-1');
        v_('N-1SQ') = v_('RN-1')*v_('RN-1');
        v_('BETAN') = v_('BETA')*v_('RN');
        v_('BETAN2') = v_('BETAN')*v_('BETAN');
        v_('AL+1') = 1.0+v_('RALPHA');
        v_('A1SQ') = v_('AL+1')*v_('AL+1');
        v_('F0') = v_('A1SQ')*v_('N-1SQ');
        v_('F1') = -1.0*v_('F0');
        v_('F2') = v_('BETAN2')+v_('F1');
        v_('F3') = 1.0/v_('F2');
        v_('F4') = v_('BETAN')*v_('F3');
        v_('A') = -1.0*v_('F4');
        v_('-N/2') = -0.5*v_('RN');
        v_('1') = 1;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        for I=v_('1'):v_('N')
            [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
            gtype{ig} = '<>';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(I)]);
            valA(end+1) = v_('BETAN');
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = ix_.Count;
        ngrp   = ig_.Count;
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        for I=v_('1'):v_('N')
            v_('RI') = I;
            v_('I-N/2') = v_('RI')+v_('-N/2');
            v_('CI') = 1.0;
            for J=v_('1'):v_('GAMMA')
                v_('CI') = v_('CI')*v_('I-N/2');
            end
            pbm.gconst(ig_(['G',int2str(I)])) = v_('CI');
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        for I=v_('1'):v_('N')
            v_('I-1') = -1+I;
            v_('RI') = I;
            v_('H') = 0.0;
            for J=v_('1'):v_('I-1')
                v_('RJ') = J;
                v_('1/J') = 1.0/v_('RJ');
                v_('I/J') = v_('RI')*v_('1/J');
                v_('SQI/J') = sqrt(v_('I/J'));
                v_('LIJ') = log(v_('SQI/J'));
                v_('SIJ') = sin(v_('LIJ'));
                v_('CIJ') = cos(v_('LIJ'));
                v_('SA') = 1.0;
                v_('CA') = 1.0;
                for K=v_('1'):v_('ALPHA')
                    v_('SA') = v_('SA')*v_('SIJ');
                    v_('CA') = v_('CA')*v_('CIJ');
                end
                v_('SCA') = v_('SA')+v_('CA');
                v_('HIJ') = v_('SQI/J')*v_('SCA');
                v_('H') = v_('H')+v_('HIJ');
            end
            v_('I+1') = 1+I;
            for J=v_('I+1'):v_('N')
                v_('RJ') = J;
                v_('1/J') = 1.0/v_('RJ');
                v_('I/J') = v_('RI')*v_('1/J');
                v_('SQI/J') = sqrt(v_('I/J'));
                v_('LIJ') = log(v_('SQI/J'));
                v_('SIJ') = sin(v_('LIJ'));
                v_('CIJ') = cos(v_('LIJ'));
                v_('SA') = 1.0;
                v_('CA') = 1.0;
                for K=v_('1'):v_('ALPHA')
                    v_('SA') = v_('SA')*v_('SIJ');
                    v_('CA') = v_('CA')*v_('CIJ');
                end
                v_('SCA') = v_('SA')+v_('CA');
                v_('HIJ') = v_('SQI/J')*v_('SCA');
                v_('H') = v_('H')+v_('HIJ');
            end
            v_('I-N/2') = v_('RI')+v_('-N/2');
            v_('CI') = 1.0;
            for J=v_('1'):v_('GAMMA')
                v_('CI') = v_('CI')*v_('I-N/2');
            end
            v_('TMP') = v_('H')+v_('CI');
            v_('XI0') = v_('TMP')*v_('A');
            pb.x0(ix_(['X',int2str(I)]),1) = v_('XI0');
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eMANC',iet_);
        elftv{it}{1} = 'X';
        elftp{it}{1} = 'II';
        elftp{it}{2} = 'JJ';
        elftp{it}{3} = 'AL';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('1'):v_('N')
            v_('RI') = I;
            v_('I-1') = -1+I;
            for J=v_('1'):v_('I-1')
                v_('RJ') = J;
                ename = ['E',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eMANC';
                ielftype(ie) = iet_('eMANC');
                vname = ['X',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('X',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                [~,posep] = ismember('II',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_('RI');
                [~,posep] = ismember('JJ',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_('RJ');
                [~,posep] = ismember('AL',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_('RALPHA');
            end
            v_('I+1') = 1+I;
            for J=v_('I+1'):v_('N')
                v_('RJ') = J;
                ename = ['E',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eMANC';
                ielftype(ie) = iet_('eMANC');
                vname = ['X',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('X',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                [~,posep] = ismember('II',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_('RI');
                [~,posep] = ismember('JJ',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_('RJ');
                [~,posep] = ismember('AL',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_('RALPHA');
            end
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('N')
            v_('I-1') = -1+I;
            ig = ig_(['G',int2str(I)]);
            pbm.grftype{ig} = 'gL2';
            for J=v_('1'):v_('I-1')
                ig = ig_(['G',int2str(I)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['E',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = 1.;
            end
            v_('I+1') = 1+I;
            for J=v_('I+1'):v_('N')
                ig = ig_(['G',int2str(I)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['E',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = 1.;
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Least square problems are bounded below by zero
        pb.objlower = 0.0;
%    Solution
% LO SOLTN               0.0
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'C-CSUR2-AN-V-0';
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

    case 'eMANC'

        EV_  = varargin{1};
        iel_ = varargin{2};
        IAL = pbm.elpar{iel_}(3);
        IA1 = IAL-1;
        A2 = pbm.elpar{iel_}(3)-2.0;
        IA2 = IAL-2;
        IA3 = IAL-3;
        INVIJ = EV_(1)*EV_(1)+pbm.elpar{iel_}(1)/pbm.elpar{iel_}(2);
        VIJ = sqrt(INVIJ);
        V2 = VIJ*VIJ;
        DVIJ = EV_(1)/VIJ;
        LIJ = log(VIJ);
        SIJ = sin(LIJ);
        CIJ = cos(LIJ);
        DSDX = CIJ*DVIJ/VIJ;
        DCDX = -SIJ*DVIJ/VIJ;
        SUMAL = SIJ^IAL+CIJ^IAL;
        DSUMAL = pbm.elpar{iel_}(3)*(DSDX*SIJ^IA1+DCDX*CIJ^IA1);
        SCIJ = SIJ*CIJ;
        DSCIJ = SIJ*DCDX+DSDX*CIJ;
        SAL = SIJ^IA2-CIJ^IA2;
        DSAL = A2*(DSDX*SIJ^IA3-DCDX*CIJ^IA3);
        B = SUMAL+pbm.elpar{iel_}(3)*SCIJ*SAL;
        DBDX = DSUMAL+pbm.elpar{iel_}(3)*(DSCIJ*SAL+SCIJ*DSAL);
        varargout{1} = VIJ*SUMAL;
        if(nargout>1)
            g_(1,1) = EV_(1)*B/VIJ;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = (B+EV_(1)*DBDX)/VIJ-B*EV_(1)*DVIJ/V2;
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

