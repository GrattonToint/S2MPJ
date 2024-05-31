function varargout = ALLINITU(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : ALLINITU
%    *********
% 
%    A problem with "all in it". Intended to verify that changes
%    to LANCELOT are safe.
% 
%    Source:
%    N. Gould, private communication.
% 
%    SIF input: Nick Gould, June 1990.
% 
%    classification = 'OUR2-AY-4-0'
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'ALLINITU';

switch(action)

    case 'setup'

    pb.name      = 'ALLINITU';
    pb.sifpbname = 'ALLINITU';
    pbm.name     = 'ALLINITU';
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
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
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2xlib('ii','FT1',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2xlib('ii','FT2',ig_);
        gtype{ig} = '<>';
        iv = ix_('X3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2xlib('ii','FT3',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2xlib('ii','FT4',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2xlib('ii','FT5',ig_);
        gtype{ig} = '<>';
        iv = ix_('X4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2xlib('ii','FT6',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2xlib('ii','FNT1',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2xlib('ii','FNT2',ig_);
        gtype{ig} = '<>';
        iv = ix_('X4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2xlib('ii','FNT3',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2xlib('ii','FNT4',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2xlib('ii','FNT5',ig_);
        gtype{ig} = '<>';
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2xlib('ii','FNT6',ig_);
        gtype{ig} = '<>';
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        pbm.gconst(ig_('FT2')) = 1.0;
        pbm.gconst(ig_('FT5')) = 3.0;
        pbm.gconst(ig_('FNT2')) = 1.0;
        pbm.gconst(ig_('FNT5')) = 4.0;
        pb.xlower = zeros(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2xlib( 'ii', 'eSQR',iet_);
        elftv{it}{1} = 'X';
        [it,iet_] = s2xlib( 'ii', 'eSQR2',iet_);
        elftv{it}{1} = 'Y';
        elftv{it}{2} = 'Z';
        [it,iet_] = s2xlib( 'ii', 'eSINSQR',iet_);
        elftv{it}{1} = 'X';
        [it,iet_] = s2xlib( 'ii', 'ePRODSQR',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = 'FT3E1';
        [ie,ie_,newelt] = s2xlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = 'X1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FT4E1';
        [ie,ie_,newelt] = s2xlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = 'X2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FT4E2';
        [ie,ie_,newelt] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR2';
        ielftype(ie) = iet_('eSQR2');
        vname = 'X3';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FT56E1';
        [ie,ie_,newelt] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSINSQR';
        ielftype(ie) = iet_('eSINSQR');
        vname = 'X3';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FT5E2';
        [ie,ie_,newelt] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePRODSQR';
        ielftype(ie) = iet_('ePRODSQR');
        vname = 'X1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FNT3E1';
        [ie,ie_,newelt] = s2xlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = 'X2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FNT4E1';
        [ie,ie_,newelt] = s2xlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = 'X3';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FNT4E2';
        [ie,ie_,newelt] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR2';
        ielftype(ie) = iet_('eSQR2');
        vname = 'X4';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FNT56E1';
        [ie,ie_,newelt] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSINSQR';
        ielftype(ie) = iet_('eSINSQR');
        vname = 'X4';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FNT5E2';
        [ie,ie_,newelt] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePRODSQR';
        ielftype(ie) = iet_('ePRODSQR');
        vname = 'X2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2xlib('ii','gTRIVIAL',igt_);
        [it,igt_] = s2xlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for ig = 1:ngrp
            pbm.grftype{ig} = 'gL2';
        end
        ig = ig_('FT1');
        pbm.grftype{ig} = 'gTRIVIAL';
        ig = ig_('FT2');
        pbm.grftype{ig} = 'gTRIVIAL';
        ig = ig_('FT3');
        pbm.grftype{ig} = 'gTRIVIAL';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FT3E1');
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('FT4');
        pbm.grftype{ig} = 'gTRIVIAL';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FT4E1');
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('FT4E2');
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('FT5');
        pbm.grftype{ig} = 'gTRIVIAL';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FT56E1');
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('FT5E2');
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('FT6');
        pbm.grftype{ig} = 'gTRIVIAL';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FT56E1');
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('FNT3');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FNT3E1');
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('FNT4');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FNT4E1');
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('FNT4E2');
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('FNT5');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FNT56E1');
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('FNT5E2');
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('FNT6');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FNT56E1');
        pbm.grelw{ig}(posel) = 1.;
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'OUR2-AY-4-0';
        pb.x0          = zeros(pb.n,1);
        varargout{1} = pb;
        varargout{2} = pbm;

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eSQR'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(1);
        if(nargout>1)
            g_(1,1) = EV_(1)+EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0;
                varargout{3} = H_;
            end
        end

    case 'eSQR2'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(1,2);
        U_(1,1) = U_(1,1)+1;
        U_(1,2) = U_(1,2)+1;
        IV_(1) = U_(1,:)*EV_;
        varargout{1} = IV_(1)*IV_(1);
        if(nargout>1)
            g_(1,1) = IV_(1)+IV_(1);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_(1,1) = 2.0;
                varargout{3} = U_.'*H_*U_;
            end
        end

    case 'eSINSQR'

        EV_  = varargin{1};
        iel_ = varargin{2};
        SINX = sin(EV_(1));
        COSX = cos(EV_(1));
        varargout{1} = SINX*SINX;
        if(nargout>1)
            g_(1,1) = 2.0*SINX*COSX;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0*(COSX*COSX-SINX*SINX);
                varargout{3} = H_;
            end
        end

    case 'ePRODSQR'

        EV_  = varargin{1};
        iel_ = varargin{2};
        XX = EV_(1)*EV_(1);
        YY = EV_(2)*EV_(2);
        varargout{1} = XX*YY;
        if(nargout>1)
            g_(1,1) = 2.0*EV_(1)*YY;
            g_(2,1) = 2.0*XX*EV_(2);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 2.0*YY;
                H_(1,2) = 4.0*EV_(1)*EV_(2);
                H_(2,1) = H_(1,2);
                H_(2,2) = 2.0*XX;
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gTRIVIAL'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = GVAR_;
        if(nargout>1)
            g_ = 1.0;
            varargout{2} = g_;
            if(nargout>2)
                H_ = 0.0;
                varargout{3} = H_;
            end
        end

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
