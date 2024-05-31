function varargout = GROWTH(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : GROWTH
%    *********
%    GROWTH problem in 3 variables
% 
%    Fit the observed growth g(n) from Gaussian Elimination
%    with complete pivoting to a function of the form
%         U1 * n ** ( U2 + LOG(n) * U3 )
% 
%    SIF input: Nick Gould, Nov, 1991.
% 
%    classification = 'NOR2-AN-3-12'
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'GROWTH';

switch(action)

    case 'setup'

    pb.name      = 'GROWTH';
    pb.sifpbname = 'GROWTH';
    pbm.name     = 'GROWTH';
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        v_('N') = 3;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2xlib('ii','U1',ix_);
        pb.xnames{iv} = 'U1';
        [iv,ix_] = s2xlib('ii','U2',ix_);
        pb.xnames{iv} = 'U2';
        [iv,ix_] = s2xlib('ii','U3',ix_);
        pb.xnames{iv} = 'U3';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2xlib('ii','G8',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'G8';
        [ig,ig_] = s2xlib('ii','G9',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'G9';
        [ig,ig_] = s2xlib('ii','G10',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'G10';
        [ig,ig_] = s2xlib('ii','G11',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'G11';
        [ig,ig_] = s2xlib('ii','G12',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'G12';
        [ig,ig_] = s2xlib('ii','G13',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'G13';
        [ig,ig_] = s2xlib('ii','G14',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'G14';
        [ig,ig_] = s2xlib('ii','G15',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'G15';
        [ig,ig_] = s2xlib('ii','G16',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'G16';
        [ig,ig_] = s2xlib('ii','G18',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'G18';
        [ig,ig_] = s2xlib('ii','G20',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'G20';
        [ig,ig_] = s2xlib('ii','G25',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'G25';
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
        pbm.gconst(ig_('G8')) = 8.0;
        pbm.gconst(ig_('G9')) = 8.4305;
        pbm.gconst(ig_('G10')) = 9.5294;
        pbm.gconst(ig_('G11')) = 10.4627;
        pbm.gconst(ig_('G12')) = 12.0;
        pbm.gconst(ig_('G13')) = 13.0205;
        pbm.gconst(ig_('G14')) = 14.5949;
        pbm.gconst(ig_('G15')) = 16.1078;
        pbm.gconst(ig_('G16')) = 18.0596;
        pbm.gconst(ig_('G18')) = 20.4569;
        pbm.gconst(ig_('G20')) = 24.25;
        pbm.gconst(ig_('G25')) = 32.9863;
        pb.xlower = zeros(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        pb.x0(ix_('U1'),1) = 100.0;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2xlib( 'ii', 'eFIT',iet_);
        elftv{it}{1} = 'U1';
        elftv{it}{2} = 'U2';
        elftv{it}{3} = 'U3';
        elftp{it}{1} = 'RN';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        ename = 'G8';
        [ie,ie_,newelt] = s2xlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eFIT';
            ielftype(ie) = iet_('eFIT');
        end
        vname = 'U1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U3';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('RN',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 8.0;
        ename = 'G9';
        [ie,ie_,newelt] = s2xlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eFIT';
            ielftype(ie) = iet_('eFIT');
        end
        vname = 'U1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U3';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('RN',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 9.0;
        ename = 'G10';
        [ie,ie_,newelt] = s2xlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eFIT';
            ielftype(ie) = iet_('eFIT');
        end
        vname = 'U1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U3';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('RN',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 10.0;
        ename = 'G11';
        [ie,ie_,newelt] = s2xlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eFIT';
            ielftype(ie) = iet_('eFIT');
        end
        vname = 'U1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U3';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('RN',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 11.0;
        ename = 'G12';
        [ie,ie_,newelt] = s2xlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eFIT';
            ielftype(ie) = iet_('eFIT');
        end
        vname = 'U1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U3';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('RN',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 12.0;
        ename = 'G13';
        [ie,ie_,newelt] = s2xlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eFIT';
            ielftype(ie) = iet_('eFIT');
        end
        vname = 'U1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U3';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('RN',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 13.0;
        ename = 'G14';
        [ie,ie_,newelt] = s2xlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eFIT';
            ielftype(ie) = iet_('eFIT');
        end
        vname = 'U1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U3';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('RN',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 14.0;
        ename = 'G15';
        [ie,ie_,newelt] = s2xlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eFIT';
            ielftype(ie) = iet_('eFIT');
        end
        vname = 'U1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U3';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('RN',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 15.0;
        ename = 'G16';
        [ie,ie_,newelt] = s2xlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eFIT';
            ielftype(ie) = iet_('eFIT');
        end
        vname = 'U1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U3';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('RN',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 16.0;
        ename = 'G18';
        [ie,ie_,newelt] = s2xlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eFIT';
            ielftype(ie) = iet_('eFIT');
        end
        vname = 'U1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U3';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('RN',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 18.0;
        ename = 'G20';
        [ie,ie_,newelt] = s2xlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eFIT';
            ielftype(ie) = iet_('eFIT');
        end
        vname = 'U1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U3';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('RN',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 20.0;
        ename = 'G25';
        [ie,ie_,newelt] = s2xlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eFIT';
            ielftype(ie) = iet_('eFIT');
        end
        vname = 'U1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U3';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('RN',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 25.0;
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2xlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        ig = ig_('G8');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('G8');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('G9');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('G9');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('G10');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('G10');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('G11');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('G11');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('G12');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('G12');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('G13');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('G13');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('G14');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('G14');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('G15');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('G15');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('G16');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('G16');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('G18');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('G18');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('G20');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('G20');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('G25');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('G25');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'NOR2-AN-3-12';
        varargout{1} = pb;
        varargout{2} = pbm;

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eFIT'

        EV_  = varargin{1};
        iel_ = varargin{2};
        LOGRN = log(pbm.elpar{iel_}(1));
        POWER = pbm.elpar{iel_}(1)^(EV_(2)+LOGRN*EV_(3));
        varargout{1} = EV_(1)*POWER;
        if(nargout>1)
            g_(1,1) = POWER;
            g_(2,1) = EV_(1)*POWER*LOGRN;
            g_(3,1) = EV_(1)*POWER*LOGRN^2;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,1) = 0.0;
                H_(1,2) = POWER*LOGRN;
                H_(2,1) = H_(1,2);
                H_(1,3) = POWER*LOGRN^2;
                H_(3,1) = H_(1,3);
                H_(2,2) = EV_(1)*POWER*LOGRN^2;
                H_(2,3) = EV_(1)*POWER*LOGRN^3;
                H_(3,2) = H_(2,3);
                H_(3,3) = EV_(1)*POWER*LOGRN^4;
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
