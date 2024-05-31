function varargout = OSBORNEA(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : OSBORNEA
%    *********
% 
%    Osborne first problem in 5 variables.
% 
%    This function  is a nonlinear least squares with 33 groups.  Each
%    group has 2 nonlinear elements and one linear element.
% 
%    Source:  Problem 17 in
%    J.J. More', B.S. Garbow and K.E. Hillstrom,
%    "Testing Unconstrained Optimization Software",
%    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
% 
%    See alos Buckley#32 (p. 77).
% 
%    SIF input: Ph. Toint, Dec 1989.
% 
%    classification = 'SUR2-MN-5-0'
% 
%    Number of groups
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'OSBORNEA';

switch(action)

    case 'setup'

    pb.name      = 'OSBORNEA';
    pb.sifpbname = 'OSBORNEA';
    pbm.name     = 'OSBORNEA';
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        v_('M') = 33;
        v_('1') = 1;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2xlib('ii','X1',ix_);
        pb.xnames{iv} = 'X1';
        [iv,ix_] = s2xlib('ii','X2',ix_);
        pb.xnames{iv} = 'X2';
        [iv,ix_] = s2xlib('ii','X3',ix_);
        pb.xnames{iv} = 'X3';
        [iv,ix_] = s2xlib('ii','X4',ix_);
        pb.xnames{iv} = 'X4';
        [iv,ix_] = s2xlib('ii','X5',ix_);
        pb.xnames{iv} = 'X5';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('M')
            [ig,ig_] = s2xlib('ii',['G',int2str(I)],ig_);
            gtype{ig} = '<>';
            iv = ix_('X1');
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
        pbm.gconst(ig_('G1')) = 0.844;
        pbm.gconst(ig_('G2')) = 0.908;
        pbm.gconst(ig_('G3')) = 0.932;
        pbm.gconst(ig_('G4')) = 0.936;
        pbm.gconst(ig_('G5')) = 0.925;
        pbm.gconst(ig_('G6')) = 0.908;
        pbm.gconst(ig_('G7')) = 0.881;
        pbm.gconst(ig_('G8')) = 0.850;
        pbm.gconst(ig_('G9')) = 0.818;
        pbm.gconst(ig_('G10')) = 0.784;
        pbm.gconst(ig_('G11')) = 0.751;
        pbm.gconst(ig_('G12')) = 0.718;
        pbm.gconst(ig_('G13')) = 0.685;
        pbm.gconst(ig_('G14')) = 0.658;
        pbm.gconst(ig_('G15')) = 0.628;
        pbm.gconst(ig_('G16')) = 0.603;
        pbm.gconst(ig_('G17')) = 0.580;
        pbm.gconst(ig_('G18')) = 0.558;
        pbm.gconst(ig_('G19')) = 0.538;
        pbm.gconst(ig_('G20')) = 0.522;
        pbm.gconst(ig_('G21')) = 0.506;
        pbm.gconst(ig_('G22')) = 0.490;
        pbm.gconst(ig_('G23')) = 0.478;
        pbm.gconst(ig_('G24')) = 0.467;
        pbm.gconst(ig_('G25')) = 0.457;
        pbm.gconst(ig_('G26')) = 0.448;
        pbm.gconst(ig_('G27')) = 0.438;
        pbm.gconst(ig_('G28')) = 0.431;
        pbm.gconst(ig_('G29')) = 0.424;
        pbm.gconst(ig_('G30')) = 0.420;
        pbm.gconst(ig_('G31')) = 0.414;
        pbm.gconst(ig_('G32')) = 0.411;
        pbm.gconst(ig_('G33')) = 0.406;
        pb.xlower = zeros(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.x0(ix_('X1'),1) = 0.5;
        pb.x0(ix_('X2'),1) = 1.5;
        pb.x0(ix_('X3'),1) = -1.0;
        pb.x0(ix_('X4'),1) = 0.01;
        pb.x0(ix_('X5'),1) = 0.02;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2xlib( 'ii', 'ePEXP',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftp{it}{1} = 'T';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('1'):v_('M')
            v_('I-1') = -1+I;
            v_('ITI') = 10*v_('I-1');
            v_('MTI') = v_('ITI');
            v_('TI') = -1.0*v_('MTI');
            ename = ['A',int2str(I)];
            [ie,ie_] = s2xlib('ii',ename,ie_);
            pbm.elftype{ie} = 'ePEXP';
            ielftype(ie) = iet_('ePEXP');
            vname = 'X2';
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X4';
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('T',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('TI');
            ename = ['B',int2str(I)];
            [ie,ie_] = s2xlib('ii',ename,ie_);
            pbm.elftype{ie} = 'ePEXP';
            ielftype(ie) = iet_('ePEXP');
            vname = 'X3';
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X5';
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('T',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('TI');
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2xlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('M')
            ig = ig_(['G',int2str(I)]);
            pbm.grftype{ig} = 'gL2';
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['A',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['B',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'SUR2-MN-5-0';
        varargout{1} = pb;
        varargout{2} = pbm;

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'ePEXP'

        EV_  = varargin{1};
        iel_ = varargin{2};
        EXPA = exp(pbm.elpar{iel_}(1)*EV_(2));
        V1EXPA = EV_(1)*EXPA;
        varargout{1} = V1EXPA;
        if(nargout>1)
            g_(1,1) = EXPA;
            g_(2,1) = pbm.elpar{iel_}(1)*V1EXPA;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = pbm.elpar{iel_}(1)*EXPA;
                H_(2,1) = H_(1,2);
                H_(2,2) = pbm.elpar{iel_}(1)*pbm.elpar{iel_}(1)*V1EXPA;
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
