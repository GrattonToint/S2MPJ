function varargout = CANTILVR(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem:
%    ********
% 
%    Computation of a minimum weight cantilever consisting of 5 sections of
%    square shaped tube of given thickness.  
% 
%    Source:
%    an example in a talk by W.K. Zhang and C. Fleury, LLN, 1994.
% 
%    SIF input: Ph. Toint, November 1994
% 
%    classification = 'LOR2-MN-5-1'
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'CANTILVR';

switch(action)

    case 'setup'

    pb.name      = 'CANTILVR';
    pb.sifpbname = 'CANTILVR';
    pbm.name     = 'CANTILVR';
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
        [iv,ix_] = s2xlib('ii','X5',ix_);
        pb.xnames{iv} = 'X5';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2xlib('ii','WEIGHT',ig_);
        gtype{ig} = '<>';
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.0624+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.0624;
        end
        iv = ix_('X2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.0624+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.0624;
        end
        iv = ix_('X3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.0624+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.0624;
        end
        iv = ix_('X4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.0624+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.0624;
        end
        iv = ix_('X5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.0624+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.0624;
        end
        [ig,ig_] = s2xlib('ii','CONST',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CONST';
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
        pbm.congrps = find(ismember(gtype,{'<=','==','>='}));
        [pb.cnames{1:pb.m}] = deal(cnames{pbm.congrps});
        pb.nob = ngrp-pb.m;
        pbm.objgrps = find(strcmp(gtype,'<>'));
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        pbm.gconst(ig_('CONST')) = 1.0;
        pb.xlower = zeros(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = 0.000001*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 1.0*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2xlib( 'ii', 'eINVCUBE',iet_);
        elftv{it}{1} = 'XX';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = 'E1';
        [ie,ie_,newelt] = s2xlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eINVCUBE';
            ielftype(ie) = iet_('eINVCUBE');
        end
        vname = 'X1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,0.000001,[],1.0);
        posev = find(strcmp('XX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E2';
        [ie,ie_,newelt] = s2xlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eINVCUBE';
            ielftype(ie) = iet_('eINVCUBE');
        end
        vname = 'X2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,0.000001,[],1.0);
        posev = find(strcmp('XX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E3';
        [ie,ie_,newelt] = s2xlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eINVCUBE';
            ielftype(ie) = iet_('eINVCUBE');
        end
        vname = 'X3';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,0.000001,[],1.0);
        posev = find(strcmp('XX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E4';
        [ie,ie_,newelt] = s2xlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eINVCUBE';
            ielftype(ie) = iet_('eINVCUBE');
        end
        vname = 'X4';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,0.000001,[],1.0);
        posev = find(strcmp('XX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E5';
        [ie,ie_,newelt] = s2xlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eINVCUBE';
            ielftype(ie) = iet_('eINVCUBE');
        end
        vname = 'X5';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,0.000001,[],1.0);
        posev = find(strcmp('XX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        ig = ig_('CONST');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 61.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E2');
        pbm.grelw{ig}(posel) = 37.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E3');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 19.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E4');
        pbm.grelw{ig}(posel) = 7.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E5');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'LOR2-MN-5-1';
        varargout{1} = pb;
        varargout{2} = pbm;

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eINVCUBE'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = 1.0/EV_(1)^3;
        if(nargout>1)
            g_(1,1) = -3.0/EV_(1)^4;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 12.0/EV_(1)^5;
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
