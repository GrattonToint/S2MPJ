function varargout = RAT43LS(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : RAT43LS
%    *********
% 
%    NIST Data fitting problem RAT43.
% 
%    Fit: y = b1 / ((1+exp[b2-b3*x])**(1/b4)) + e
% 
%    Source:  Problem from the NIST nonlinear regression test set
%      http://www.itl.nist.gov/div898/strd/nls/nls_main.shtml
% 
%    Reference: Ratkowsky, D.A. (1983).  
%      Nonlinear Regression Modeling.
%      New York, NY:  Marcel Dekker, pp. 62 and 88.
% 
%    SIF input: Nick Gould and Tyrone Rees, Oct 2015
% 
%    classification = 'SUR2-MN-4-0'
% 
%    Number of data values
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'RAT43LS';

switch(action)

    case 'setup'

    pb.name      = 'RAT43LS';
    pb.sifpbname = 'RAT43LS';
    pbm.name     = 'RAT43LS';
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        v_('M') = 15;
        v_('N') = 4;
        v_('1') = 1;
        v_('X1') = 1.0;
        v_('X2') = 2.0;
        v_('X3') = 3.0;
        v_('X4') = 4.0;
        v_('X5') = 5.0;
        v_('X6') = 6.0;
        v_('X7') = 7.0;
        v_('X8') = 8.0;
        v_('X9') = 9.0;
        v_('X10') = 10.0;
        v_('X11') = 11.0;
        v_('X12') = 12.0;
        v_('X13') = 13.0;
        v_('X14') = 14.0;
        v_('X15') = 15.0;
        v_('Y1') = 16.08;
        v_('Y2') = 33.83;
        v_('Y3') = 65.80;
        v_('Y4') = 97.20;
        v_('Y5') = 191.55;
        v_('Y6') = 326.20;
        v_('Y7') = 386.87;
        v_('Y8') = 520.53;
        v_('Y9') = 590.03;
        v_('Y10') = 651.92;
        v_('Y11') = 724.93;
        v_('Y12') = 699.56;
        v_('Y13') = 689.96;
        v_('Y14') = 637.56;
        v_('Y15') = 717.41;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2xlib('ii',['B',int2str(I)],ix_);
            pb.xnames{iv} = ['B',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('M')
            [ig,ig_] = s2xlib('ii',['F',int2str(I)],ig_);
            gtype{ig} = '<>';
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
        pb.xlower = zeros(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.x0(ix_('B1'),1) = 100.0;
        pb.x0(ix_('B2'),1) = 10.0;
        pb.x0(ix_('B3'),1) = 1.0;
        pb.x0(ix_('B4'),1) = 1.0;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2xlib( 'ii', 'eE',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftv{it}{3} = 'V3';
        elftv{it}{4} = 'V4';
        elftp{it}{1} = 'X';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('1'):v_('M')
            ename = ['E',int2str(I)];
            [ie,ie_] = s2xlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eE';
            ielftype(ie) = iet_('eE');
            vname = 'B1';
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'B2';
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'B3';
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V3',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'B4';
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V4',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('X',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['X',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2xlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for ig = 1:ngrp
            pbm.grftype{ig} = 'gL2';
        end
        for I=v_('1'):v_('M')
            ig = ig_(['F',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'SUR2-MN-4-0';
        varargout{1} = pb;
        varargout{2} = pbm;

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eE'

        EV_  = varargin{1};
        iel_ = varargin{2};
        V2MV3X = EV_(2)-EV_(3)*pbm.elpar{iel_}(1);
        V4INV = 1.0/EV_(4);
        V4INVP = V4INV+1.0;
        E = exp(V2MV3X);
        E2 = E*E;
        EP1 = E+1.0;
        EP1L = log(EP1);
        EP14 = EP1^V4INV;
        EP14P1 = EP1^V4INVP;
        EP14P2 = EP1^(V4INV+2.0);
        VE = EV_(4)*EP14P1;
        VE2 = EV_(4)*EP14P2;
        V42EPP = EP14*EV_(4)^2;
        V42EP2 = EP14P1*EV_(4)^2;
        V42EP3 = EP14P1*EV_(4)^3;
        varargout{1} = EV_(1)/EP14;
        if(nargout>1)
            g_(1,1) = 1.0/EP14;
            g_(2,1) = -EV_(1)*E/VE;
            g_(3,1) = EV_(1)*pbm.elpar{iel_}(1)*E/VE;
            g_(4,1) = EV_(1)*EP1L/V42EPP;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(4,4);
                H_(1,2) = -E/VE;
                H_(2,1) = H_(1,2);
                H_(1,3) = pbm.elpar{iel_}(1)*E/VE;
                H_(3,1) = H_(1,3);
                H_(1,4) = EP1L/V42EPP;
                H_(4,1) = H_(1,4);
                H_(2,2) = EV_(1)*(E2*V4INVP/VE2-E/VE);
                H_(2,3) = EV_(1)*pbm.elpar{iel_}(1)*(E/VE-E2*V4INVP/VE2);
                H_(3,2) = H_(2,3);
                H_(2,4) = EV_(1)*E*(1.0/V42EP2-EP1L/V42EP3);
                H_(4,2) = H_(2,4);
                H_(3,3) = EV_(1)*pbm.elpar{iel_}(1)^2*(E2*V4INVP/VE2-E/VE);
                H_(3,4) = EV_(1)*pbm.elpar{iel_}(1)*E*(EP1L/V42EP3-1.0/V42EP2);
                H_(4,3) = H_(3,4);
                H_(4,4) = (EV_(1)/EP14)*(EP1L^2/EV_(4)^4-2.0*EP1L/EV_(4)^3);
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
            [varargout{1:max(1,nargout)}] = s2xlib(action,pbm,varargin{:});
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

