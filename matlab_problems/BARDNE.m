function varargout = BARDNE(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : BARDNE
%    *********
%    Bard problem in 3 variables.
%    This function  is a nonlinear least squares with 15 groups.  Each
%    group has a linear and a nonlinear element. This is a nonlinear equation
%    version of problem BARD
% 
%    Source: Problem 3 in
%    J.J. More', B.S. Garbow and K.E. Hillstrom,
%    "Testing Unconstrained Optimization Software",
%    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
% 
%    See also Buckley#16.
%    SIF input: Ph. Toint, Dec 1989.
%    Modification as a set of nonlinear equations: Nick Gould, Oct 2015.
% 
%    classification = 'NOR2-AN-3-15'
% 
%    Number of groups
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'BARDNE';

switch(action)

    case 'setup'

    pb.name      = 'BARDNE';
    pb.sifpbname = 'BARDNE';
    pbm.name     = 'BARDNE';
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        v_('15') = 15;
        v_('1') = 1;
        v_('8') = 8;
        v_('9') = 9;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2xlib('ii','X1',ix_);
        pb.xnames{iv} = 'X1';
        [iv,ix_] = s2xlib('ii','X2',ix_);
        pb.xnames{iv} = 'X2';
        [iv,ix_] = s2xlib('ii','X3',ix_);
        pb.xnames{iv} = 'X3';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('15')
            [ig,ig_] = s2xlib('ii',['G',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['G',int2str(I)];
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
        pbm.gconst(ig_('G1')) = 0.14;
        pbm.gconst(ig_('G2')) = 0.18;
        pbm.gconst(ig_('G3')) = 0.22;
        pbm.gconst(ig_('G4')) = 0.25;
        pbm.gconst(ig_('G5')) = 0.29;
        pbm.gconst(ig_('G6')) = 0.32;
        pbm.gconst(ig_('G7')) = 0.35;
        pbm.gconst(ig_('G8')) = 0.39;
        pbm.gconst(ig_('G9')) = 0.37;
        pbm.gconst(ig_('G10')) = 0.58;
        pbm.gconst(ig_('G11')) = 0.73;
        pbm.gconst(ig_('G12')) = 0.96;
        pbm.gconst(ig_('G13')) = 1.34;
        pbm.gconst(ig_('G14')) = 2.10;
        pbm.gconst(ig_('G15')) = 4.39;
        pb.xlower = zeros(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 1.0*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2xlib( 'ii', 'eBD',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftp{it}{1} = 'U';
        elftp{it}{2} = 'V';
        elftp{it}{3} = 'W';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('1'):v_('8')
            v_('REALI') = I;
            v_('16-I') = 16.0-v_('REALI');
            ename = ['E',int2str(I)];
            [ie,ie_,newelt] = s2xlib('ii',ename,ie_);
            if(newelt)
                pbm.elftype{ie} = 'eBD';
                ielftype(ie) = iet_('eBD');
            end
            vname = 'X2';
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X3';
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('U',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('REALI');
            [~,posep] = ismember('V',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('16-I');
            [~,posep] = ismember('W',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('REALI');
        end
        for I=v_('9'):v_('15')
            v_('REALI') = I;
            v_('16-I') = 16.0-v_('REALI');
            ename = ['E',int2str(I)];
            [ie,ie_,newelt] = s2xlib('ii',ename,ie_);
            if(newelt)
                pbm.elftype{ie} = 'eBD';
                ielftype(ie) = iet_('eBD');
            end
            vname = 'X2';
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X3';
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('U',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('REALI');
            [~,posep] = ismember('V',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('16-I');
            [~,posep] = ismember('W',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('16-I');
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('15')
            ig = ig_(['G',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'NOR2-AN-3-15';
        varargout{1} = pb;
        varargout{2} = pbm;

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eBD'

        EV_  = varargin{1};
        iel_ = varargin{2};
        Z = pbm.elpar{iel_}(2)*EV_(1)+pbm.elpar{iel_}(3)*EV_(2);
        Z2 = Z*Z;
        Z3 = Z*Z2;
        VU = pbm.elpar{iel_}(2)*pbm.elpar{iel_}(1);
        WU = pbm.elpar{iel_}(3)*pbm.elpar{iel_}(1);
        varargout{1} = pbm.elpar{iel_}(1)/Z;
        if(nargout>1)
            g_(1,1) = -VU/Z2;
            g_(2,1) = -WU/Z2;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 2.0*pbm.elpar{iel_}(2)*VU/Z3;
                H_(1,2) = 2.0*pbm.elpar{iel_}(2)*WU/Z3;
                H_(2,1) = H_(1,2);
                H_(2,2) = 2.0*pbm.elpar{iel_}(3)*WU/Z3;
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

