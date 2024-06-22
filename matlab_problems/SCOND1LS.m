function varargout = SCOND1LS(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : SCOND1LS
%    *********
% 
%    The semiconductor problem by Rheinboldt, using a finite difference
%    approximation.
%    This is the least-squares version of problem SEMICON1.
% 
%    Source: problem 10 in
%    J.J. More',
%    "A collection of nonlinear model problems"
%    Proceedings of the AMS-SIAM Summer seminar on the Computational
%    Solution of Nonlinear Systems of Equations, Colorado, 1988.
%    Argonne National Laboratory MCS-P60-0289, 1989.
% 
%    SIF input: Ph. Toint, Dec 1989.
% 
%    classification = 'SBR2-AN-V-V'
% 
%    N  = Number of discretized point inside the interval [a, b]
%    LN = Index of the last negative discretization point
%         (the interest is in the negative part)
% 
%       Alternative values for the SIF file parameters:
% IE N                   10             $-PARAMETER     original value
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'SCOND1LS';

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
        if(nargs<2)
            v_('LN') = 9;  %  SIF file default value
        else
            v_('LN') = varargin{2};
        end
% IE N                   50             $-PARAMETER
% IE LN                  45             $-PARAMETER
% IE N                   100            $-PARAMETER
% IE LN                  90             $-PARAMETER
% IE N                   500            $-PARAMETER
% IE LN                  450            $-PARAMETER
% IE N                   1000           $-PARAMETER
% IE LN                  900            $-PARAMETER
% IE N                   5000           $-PARAMETER
% IE LN                  4500           $-PARAMETER
        if(nargs<3)
            v_('LAMBDA') = 1.0;  %  SIF file default value
        else
            v_('LAMBDA') = varargin{3};
        end
        v_('A') = -0.00009;
        v_('B') = 0.00001;
        v_('UA') = 0.0;
        v_('UB') = 700.0;
        v_('CA') = 1.0e12;
        v_('CB') = 1.0e13;
        v_('BETA') = 40.0;
        v_('LN+1') = 1+v_('LN');
        v_('N+1') = 1+v_('N');
        v_('-A') = -1.0*v_('A');
        v_('B-A') = v_('B')+v_('-A');
        v_('RN+1') = v_('N+1');
        v_('TMP') = 1.0/v_('RN+1');
        v_('H') = v_('B-A')*v_('TMP');
        v_('H2') = v_('H')*v_('H');
        v_('LB') = v_('LAMBDA')*v_('BETA');
        v_('H2CA') = v_('H2')*v_('CA');
        v_('H2CB') = v_('H2')*v_('CB');
        v_('LH2CA') = v_('LAMBDA')*v_('H2CA');
        v_('LH2CB') = v_('LAMBDA')*v_('H2CB');
        v_('LUA') = v_('LAMBDA')*v_('UA');
        v_('LUB') = v_('LAMBDA')*v_('UB');
        v_('ULW') = -5.0+v_('LUA');
        v_('UUP') = 5.0+v_('LUB');
        v_('-LB') = -1.0*v_('LB');
        v_('-LUB') = -1.0*v_('LUB');
        v_('-LH2CB') = -1.0*v_('LH2CB');
        v_('0') = 0;
        v_('1') = 1;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('0'):v_('N+1')
            [iv,ix_] = s2mpjlib('ii',['U',int2str(I)],ix_);
            pb.xnames{iv} = ['U',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('N')
            v_('I+1') = 1+I;
            v_('I-1') = -1+I;
            [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
            gtype{ig} = '<>';
            iv = ix_(['U',int2str(round(v_('I-1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['U',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -2.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -2.0;
            end
            iv = ix_(['U',int2str(round(v_('I+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        for I=v_('1'):v_('LN')
            pbm.gconst(ig_(['G',int2str(I)])) = v_('LH2CA');
        end
        for I=v_('LN+1'):v_('N')
            pbm.gconst(ig_(['G',int2str(I)])) = v_('-LH2CB');
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xupper = v_('UUP')*ones(pb.n,1);
        pb.xlower = v_('ULW')*ones(pb.n,1);
        pb.xlower(ix_(['U',int2str(round(v_('0')))]),1) = v_('LUA');
        pb.xupper(ix_(['U',int2str(round(v_('0')))]),1) = v_('LUA');
        pb.xlower(ix_(['U',int2str(round(v_('N+1')))]),1) = v_('LUB');
        pb.xupper(ix_(['U',int2str(round(v_('N+1')))]),1) = v_('LUB');
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 0.0*ones(pb.n,1);
        pb.x0(ix_(['U',int2str(round(v_('0')))]),1) = v_('LUA');
        pb.x0(ix_(['U',int2str(round(v_('N+1')))]),1) = v_('LUB');
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eWE1',iet_);
        elftv{it}{1} = 'X';
        elftp{it}{1} = 'LAC';
        elftp{it}{2} = 'LAB';
        elftp{it}{3} = 'LU';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('1'):v_('N')
            ename = ['EA',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eWE1';
            ielftype(ie) = iet_('eWE1');
            vname = ['U',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,v_('ULW'),v_('UUP'),0.0);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('LAC',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('LH2CA');
            [~,posep] = ismember('LAB',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('-LB');
            [~,posep] = ismember('LU',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('LUA');
            ename = ['EB',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eWE1';
            ielftype(ie) = iet_('eWE1');
            vname = ['U',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,v_('ULW'),v_('UUP'),0.0);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('LAC',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('-LH2CB');
            [~,posep] = ismember('LAB',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('LB');
            [~,posep] = ismember('LU',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('LUB');
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2mpjlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for ig = 1:ngrp
            pbm.grftype{ig} = 'gL2';
        end
        for I=v_('1'):v_('N')
            ig = ig_(['G',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['EA',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['EB',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'SBR2-AN-V-V';
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb,  pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm, pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eWE1'

        EV_  = varargin{1};
        iel_ = varargin{2};
        FVAL =...
              pbm.elpar{iel_}(1)*exp(pbm.elpar{iel_}(2)*(EV_(1)-pbm.elpar{iel_}(3)));
        varargout{1} = FVAL;
        if(nargout>1)
            g_(1,1) = pbm.elpar{iel_}(2)*FVAL;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = pbm.elpar{iel_}(2)*pbm.elpar{iel_}(2)*FVAL;
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

