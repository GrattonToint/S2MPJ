function varargout = PALMER8E(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : PALMER8E
%    *********
% 
%    A nonlinear least squares problem
%    arising from chemical kinetics.
% 
%    model: H-N=C=Se TZVP + MP2
%    fitting Y to A0 + A2 X**2 + A4 X**4 + A6 X**6 + A8 X**8 +
%                 A10 X**10 + L * EXP( -K X**2 )
% 
%    Source:
%    M.  Palmer, Edinburgh, private communication.
% 
%    SIF input: Nick Gould, 1992.
% 
%    classification = 'SBR2-RN-8-0'
% 
%    Number of data points
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'PALMER8E';

switch(action)

    case 'setup'

    pb.name      = 'PALMER8E';
    pb.sifpbname = 'PALMER8E';
    pbm.name     = 'PALMER8E';
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        v_('M') = 23;
        v_('1') = 1;
        v_('12') = 12;
        v_('X12') = 0.000000;
        v_('X13') = 0.174533;
        v_('X14') = 0.314159;
        v_('X15') = 0.436332;
        v_('X16') = 0.514504;
        v_('X17') = 0.610865;
        v_('X18') = 0.785398;
        v_('X19') = 0.959931;
        v_('X20') = 1.134464;
        v_('X21') = 1.308997;
        v_('X22') = 1.483530;
        v_('X23') = 1.570796;
        v_('Y12') = 4.757534;
        v_('Y13') = 3.121416;
        v_('Y14') = 1.207606;
        v_('Y15') = 0.131916;
        v_('Y16') = 0.000000;
        v_('Y17') = 0.258514;
        v_('Y18') = 3.380161;
        v_('Y19') = 10.762813;
        v_('Y20') = 23.745996;
        v_('Y21') = 44.471864;
        v_('Y22') = 76.541947;
        v_('Y23') = 97.874528;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2xlib('ii','A0',ix_);
        pb.xnames{iv} = 'A0';
        [iv,ix_] = s2xlib('ii','A2',ix_);
        pb.xnames{iv} = 'A2';
        [iv,ix_] = s2xlib('ii','A4',ix_);
        pb.xnames{iv} = 'A4';
        [iv,ix_] = s2xlib('ii','A6',ix_);
        pb.xnames{iv} = 'A6';
        [iv,ix_] = s2xlib('ii','A8',ix_);
        pb.xnames{iv} = 'A8';
        [iv,ix_] = s2xlib('ii','A10',ix_);
        pb.xnames{iv} = 'A10';
        [iv,ix_] = s2xlib('ii','K',ix_);
        pb.xnames{iv} = 'K';
        [iv,ix_] = s2xlib('ii','L',ix_);
        pb.xnames{iv} = 'L';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('12'):v_('M')
            v_('XSQR') = v_(['X',int2str(I)])*v_(['X',int2str(I)]);
            v_('XQUART') = v_('XSQR')*v_('XSQR');
            v_('X**6') = v_('XSQR')*v_('XQUART');
            v_('X**8') = v_('XSQR')*v_('X**6');
            v_('X**10') = v_('XSQR')*v_('X**8');
            v_('X**12') = v_('XSQR')*v_('X**10');
            v_('X**14') = v_('XSQR')*v_('X**12');
            [ig,ig_] = s2xlib('ii',['O',int2str(I)],ig_);
            gtype{ig} = '<>';
            iv = ix_('A0');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_('A2');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('XSQR')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('XSQR');
            end
            iv = ix_('A4');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('XQUART')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('XQUART');
            end
            iv = ix_('A6');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('X**6')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('X**6');
            end
            iv = ix_('A8');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('X**8')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('X**8');
            end
            iv = ix_('A10');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('X**10')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('X**10');
            end
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        for I=v_('12'):v_('M')
            pbm.gconst(ig_(['O',int2str(I)])) = v_(['Y',int2str(I)]);
        end
        pb.xlower = zeros(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xlower(ix_('A0')) = -Inf;
        pb.xupper(ix_('A0'),1) = +Inf;
        pb.xlower(ix_('A2')) = -Inf;
        pb.xupper(ix_('A2'),1) = +Inf;
        pb.xlower(ix_('A4')) = -Inf;
        pb.xupper(ix_('A4'),1) = +Inf;
        pb.xlower(ix_('A6')) = -Inf;
        pb.xupper(ix_('A6'),1) = +Inf;
        pb.xlower(ix_('A8')) = -Inf;
        pb.xupper(ix_('A8'),1) = +Inf;
        pb.xlower(ix_('A10')) = -Inf;
        pb.xupper(ix_('A10'),1) = +Inf;
        pb.xlower(ix_('L')) = -Inf;
        pb.xupper(ix_('L'),1) = +Inf;
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 1.0*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2xlib( 'ii', 'ePROD',iet_);
        elftv{it}{1} = 'K';
        elftv{it}{2} = 'L';
        elftp{it}{1} = 'XSQR';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('12'):v_('M')
            v_('XSQR') = v_(['X',int2str(I)])*v_(['X',int2str(I)]);
            ename = ['E',int2str(I)];
            [ie,ie_] = s2xlib('ii',ename,ie_);
            pbm.elftype{ie} = 'ePROD';
            ielftype(ie) = iet_('ePROD');
            vname = 'K';
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('K',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'L';
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('L',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('XSQR',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('XSQR');
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2xlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('12'):v_('M')
            ig = ig_(['O',int2str(I)]);
            pbm.grftype{ig} = 'gL2';
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'SBR2-RN-8-0';
        varargout{1} = pb;
        varargout{2} = pbm;

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'ePROD'

        EV_  = varargin{1};
        iel_ = varargin{2};
        EXPON = exp(-EV_(1)*pbm.elpar{iel_}(1));
        varargout{1} = EV_(2)*EXPON;
        if(nargout>1)
            g_(1,1) = -pbm.elpar{iel_}(1)*EV_(2)*EXPON;
            g_(2,1) = EXPON;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = pbm.elpar{iel_}(1)*pbm.elpar{iel_}(1)*EV_(2)*EXPON;
                H_(1,2) = -pbm.elpar{iel_}(1)*EXPON;
                H_(2,1) = H_(1,2);
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
