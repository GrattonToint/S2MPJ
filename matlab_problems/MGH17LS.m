function varargout = MGH17LS(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : MGH17LS
%    *********
% 
%    NIST Data fitting problem MGH17.
% 
%    Fit: y = b1 + b2*exp[-x*b4] + b3*exp[-x*b5] + e
% 
%    Source:  Problem from the NIST nonlinear regression test set
%      http://www.itl.nist.gov/div898/strd/nls/nls_main.shtml
% 
%    Reference: Osborne, M. R. (1972).  
%     Some aspects of nonlinear least squares calculations.  
%     In Numerical Methods for Nonlinear Optimization, Lootsma (Ed).  
%     New York, NY:  Academic Press, pp. 171-189.
% 
%    SIF input: Nick Gould and Tyrone Rees, Oct 2015
% 
%    classification = 'SUR2-MN-5-0'
% 
%    Number of data values
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'MGH17LS';

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
        v_('M') = 33;
        v_('N') = 5;
        v_('1') = 1;
        v_('X1') = 0.0E+0;
        v_('X2') = 1.0E+1;
        v_('X3') = 2.0E+1;
        v_('X4') = 3.0E+1;
        v_('X5') = 4.0E+1;
        v_('X6') = 5.0E+1;
        v_('X7') = 6.0E+1;
        v_('X8') = 7.0E+1;
        v_('X9') = 8.0E+1;
        v_('X10') = 9.0E+1;
        v_('X11') = 1.0E+2;
        v_('X12') = 1.1E+2;
        v_('X13') = 1.2E+2;
        v_('X14') = 1.3E+2;
        v_('X15') = 1.4E+2;
        v_('X16') = 1.5E+2;
        v_('X17') = 1.6E+2;
        v_('X18') = 1.7E+2;
        v_('X19') = 1.8E+2;
        v_('X20') = 1.9E+2;
        v_('X21') = 2.0E+2;
        v_('X22') = 2.1E+2;
        v_('X23') = 2.2E+2;
        v_('X24') = 2.3E+2;
        v_('X25') = 2.4E+2;
        v_('X26') = 2.5E+2;
        v_('X27') = 2.6E+2;
        v_('X28') = 2.7E+2;
        v_('X29') = 2.8E+2;
        v_('X30') = 2.9E+2;
        v_('X31') = 3.0E+2;
        v_('X32') = 3.1E+2;
        v_('X33') = 3.2E+2;
        v_('Y1') = 8.44E-1;
        v_('Y2') = 9.08E-1;
        v_('Y3') = 9.32E-1;
        v_('Y4') = 9.36E-1;
        v_('Y5') = 9.25E-1;
        v_('Y6') = 9.08E-1;
        v_('Y7') = 8.81E-1;
        v_('Y8') = 8.50E-1;
        v_('Y9') = 8.18E-1;
        v_('Y10') = 7.84E-1;
        v_('Y11') = 7.51E-1;
        v_('Y12') = 7.18E-1;
        v_('Y13') = 6.85E-1;
        v_('Y14') = 6.58E-1;
        v_('Y15') = 6.28E-1;
        v_('Y16') = 6.03E-1;
        v_('Y17') = 5.80E-1;
        v_('Y18') = 5.58E-1;
        v_('Y19') = 5.38E-1;
        v_('Y20') = 5.22E-1;
        v_('Y21') = 5.06E-1;
        v_('Y22') = 4.90E-1;
        v_('Y23') = 4.78E-1;
        v_('Y24') = 4.67E-1;
        v_('Y25') = 4.57E-1;
        v_('Y26') = 4.48E-1;
        v_('Y27') = 4.38E-1;
        v_('Y28') = 4.31E-1;
        v_('Y29') = 4.24E-1;
        v_('Y30') = 4.20E-1;
        v_('Y31') = 4.14E-1;
        v_('Y32') = 4.11E-1;
        v_('Y33') = 4.06E-1;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['B',int2str(I)],ix_);
            pb.xnames{iv} = ['B',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('M')
            [ig,ig_] = s2mpjlib('ii',['F',int2str(I)],ig_);
            gtype{ig} = '<>';
            iv = ix_('B1');
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
        for I=v_('1'):v_('M')
            pbm.gconst(ig_(['F',int2str(I)])) = v_(['Y',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.x0(ix_('B1'),1) = 50.0;
        pb.x0(ix_('B2'),1) = 150.0;
        pb.x0(ix_('B3'),1) = -100.0;
        pb.x0(ix_('B4'),1) = 1.0;
        pb.x0(ix_('B5'),1) = 2.0;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eE2',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftp{it}{1} = 'X';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('1'):v_('M')
            ename = ['EA',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eE2';
            ielftype(ie) = iet_('eE2');
            vname = 'B2';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'B4';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('X',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['X',int2str(I)]);
            ename = ['EB',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eE2';
            ielftype(ie) = iet_('eE2');
            vname = 'B3';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'B5';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('X',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['X',int2str(I)]);
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
        for I=v_('1'):v_('M')
            ig = ig_(['F',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['EA',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['EB',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Least square problems are bounded below by zero
        pb.objlower = 0.0;
%    Solution
% LO SOLTN               
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'SUR2-MN-5-0';
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb,  pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm, pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end
% **********************
%  SET UP THE FUNCTION *
%  AND RANGE ROUTINES  *
% **********************

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eE2'

        EV_  = varargin{1};
        iel_ = varargin{2};
        E = exp(-EV_(2)*pbm.elpar{iel_}(1));
        varargout{1} = EV_(1)*E;
        if(nargout>1)
            g_(1,1) = E;
            g_(2,1) = -EV_(1)*pbm.elpar{iel_}(1)*E;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = -pbm.elpar{iel_}(1)*E;
                H_(2,1) = H_(1,2);
                H_(2,2) = EV_(1)*pbm.elpar{iel_}(1)*pbm.elpar{iel_}(1)*E;
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

