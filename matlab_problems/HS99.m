function varargout = HS99(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HS99
%    *********
% 
%    Source: problem 99 in
%    W. Hock and K. Schittkowski,
%    "Test examples for nonlinear programming codes",
%    Lectures Notes in Economics and Mathematical Systems 187, Springer
%    Verlag, Heidelberg, 1981.
% 
%    SIF input: Ph. Toint, April 1991.
% 
%    classification = 'OOR2-AN-7-2'
% 
%    Constants
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HS99';

switch(action)

    case 'setup'

    pb.name      = 'HS99';
    pb.sifpbname = 'HS99';
    pbm.name     = 'HS99';
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        v_('1') = 1;
        v_('2') = 2;
        v_('7') = 7;
        v_('8') = 8;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('7')
            [iv,ix_] = s2xlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2xlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        pbm.gscale(ig,1) = -1.0;
        [ig,ig_] = s2xlib('ii','Q8E',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'Q8E';
        [ig,ig_] = s2xlib('ii','S8E',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'S8E';
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
        pbm.gconst(ig_('Q8E')) = 100000.0;
        pbm.gconst(ig_('S8E')) = 1000.0;
        pb.xlower = zeros(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xupper = 1.58*ones(pb.n,1);
        pb.xlower = -Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 0.5*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2xlib( 'ii', 'eR8T',iet_);
        elftv{it}{1} = 'X1';
        elftv{it}{2} = 'X2';
        elftv{it}{3} = 'X3';
        elftv{it}{4} = 'X4';
        elftv{it}{5} = 'X5';
        elftv{it}{6} = 'X6';
        elftv{it}{7} = 'X7';
        [it,iet_] = s2xlib( 'ii', 'eQ8T',iet_);
        elftv{it}{1} = 'X1';
        elftv{it}{2} = 'X2';
        elftv{it}{3} = 'X3';
        elftv{it}{4} = 'X4';
        elftv{it}{5} = 'X5';
        elftv{it}{6} = 'X6';
        elftv{it}{7} = 'X7';
        [it,iet_] = s2xlib( 'ii', 'eS8T',iet_);
        elftv{it}{1} = 'X1';
        elftv{it}{2} = 'X2';
        elftv{it}{3} = 'X3';
        elftv{it}{4} = 'X4';
        elftv{it}{5} = 'X5';
        elftv{it}{6} = 'X6';
        elftv{it}{7} = 'X7';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = 'R8';
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eR8T';
        ielftype(ie) = iet_('eR8T');
        vname = 'X1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],1.58,0.5);
        posev = find(strcmp('X1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],1.58,0.5);
        posev = find(strcmp('X2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],1.58,0.5);
        posev = find(strcmp('X3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],1.58,0.5);
        posev = find(strcmp('X4',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X5';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],1.58,0.5);
        posev = find(strcmp('X5',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X6';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],1.58,0.5);
        posev = find(strcmp('X6',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X7';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],1.58,0.5);
        posev = find(strcmp('X7',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'Q8';
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eQ8T';
        ielftype(ie) = iet_('eQ8T');
        vname = 'X1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],1.58,0.5);
        posev = find(strcmp('X1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],1.58,0.5);
        posev = find(strcmp('X2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],1.58,0.5);
        posev = find(strcmp('X3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],1.58,0.5);
        posev = find(strcmp('X4',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X5';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],1.58,0.5);
        posev = find(strcmp('X5',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X6';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],1.58,0.5);
        posev = find(strcmp('X6',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X7';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],1.58,0.5);
        posev = find(strcmp('X7',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'S8';
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eS8T';
        ielftype(ie) = iet_('eS8T');
        vname = 'X1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],1.58,0.5);
        posev = find(strcmp('X1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],1.58,0.5);
        posev = find(strcmp('X2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],1.58,0.5);
        posev = find(strcmp('X3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],1.58,0.5);
        posev = find(strcmp('X4',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X5';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],1.58,0.5);
        posev = find(strcmp('X5',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X6';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],1.58,0.5);
        posev = find(strcmp('X6',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X7';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],1.58,0.5);
        posev = find(strcmp('X7',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2xlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        ig = ig_('OBJ');
        pbm.grftype{ig} = 'gL2';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('R8');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('Q8E');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('Q8');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('S8E');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('S8');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'OOR2-AN-7-2';
        varargout{1} = pb;
        varargout{2} = pbm;

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'e_globs'

        pbm = varargin{1};
        pbm.efpar = [];
        pbm.efpar(1) = 50.0;
        pbm.efpar(2) = 50.0;
        pbm.efpar(3) = 75.0;
        pbm.efpar(4) = 75.0;
        pbm.efpar(5) = 75.0;
        pbm.efpar(6) = 100.0;
        pbm.efpar(7) = 100.0;
        pbm.efpar(8) = 25.0;
        pbm.efpar(9) = 25.0;
        pbm.efpar(10) = 50.0;
        pbm.efpar(11) = 50.0;
        pbm.efpar(12) = 50.0;
        pbm.efpar(13) = 90.0;
        pbm.efpar(14) = 90.0;
        pbm.efpar(15) = 32.0;
        varargout{1} = pbm;

    case 'eR8T'

        EV_  = varargin{1};
        iel_ = varargin{2};
        R2 = pbm.efpar(1)*pbm.efpar(8)*cos(EV_(1));
        R3 = pbm.efpar(2)*pbm.efpar(9)*cos(EV_(2))+R2;
        R4 = pbm.efpar(3)*pbm.efpar(10)*cos(EV_(3))+R3;
        R5 = pbm.efpar(4)*pbm.efpar(11)*cos(EV_(4))+R4;
        R6 = pbm.efpar(5)*pbm.efpar(12)*cos(EV_(5))+R5;
        R7 = pbm.efpar(6)*pbm.efpar(13)*cos(EV_(6))+R6;
        varargout{1} = pbm.efpar(7)*pbm.efpar(14)*cos(EV_(7))+R7;
        if(nargout>1)
            g_(1,1) = -pbm.efpar(1)*pbm.efpar(8)*sin(EV_(1));
            g_(2,1) = -pbm.efpar(2)*pbm.efpar(9)*sin(EV_(2));
            g_(3,1) = -pbm.efpar(3)*pbm.efpar(10)*sin(EV_(3));
            g_(4,1) = -pbm.efpar(4)*pbm.efpar(11)*sin(EV_(4));
            g_(5,1) = -pbm.efpar(5)*pbm.efpar(12)*sin(EV_(5));
            g_(6,1) = -pbm.efpar(6)*pbm.efpar(13)*sin(EV_(6));
            g_(7,1) = -pbm.efpar(7)*pbm.efpar(14)*sin(EV_(7));
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(7,7);
                H_(1,1) = -pbm.efpar(1)*pbm.efpar(8)*cos(EV_(1));
                H_(2,2) = -pbm.efpar(2)*pbm.efpar(9)*cos(EV_(2));
                H_(3,3) = -pbm.efpar(3)*pbm.efpar(10)*cos(EV_(3));
                H_(4,4) = -pbm.efpar(4)*pbm.efpar(11)*cos(EV_(4));
                H_(5,5) = -pbm.efpar(5)*pbm.efpar(12)*cos(EV_(5));
                H_(6,6) = -pbm.efpar(6)*pbm.efpar(13)*cos(EV_(6));
                H_(7,7) = -pbm.efpar(7)*pbm.efpar(14)*cos(EV_(7));
                varargout{3} = H_;
            end
        end

    case 'eS8T'

        EV_  = varargin{1};
        iel_ = varargin{2};
        S2 = pbm.efpar(8)*(pbm.efpar(1)*sin(EV_(1))-pbm.efpar(15));
        S3 = pbm.efpar(9)*(pbm.efpar(2)*sin(EV_(2))-pbm.efpar(15))+S2;
        S4 = pbm.efpar(10)*(pbm.efpar(3)*sin(EV_(3))-pbm.efpar(15))+S3;
        S5 = pbm.efpar(11)*(pbm.efpar(4)*sin(EV_(4))-pbm.efpar(15))+S4;
        S6 = pbm.efpar(12)*(pbm.efpar(5)*sin(EV_(5))-pbm.efpar(15))+S5;
        S7 = pbm.efpar(13)*(pbm.efpar(6)*sin(EV_(6))-pbm.efpar(15))+S6;
        varargout{1} = pbm.efpar(14)*(pbm.efpar(7)*sin(EV_(7))-pbm.efpar(15))+S7;
        if(nargout>1)
            g_(1,1) = pbm.efpar(1)*pbm.efpar(8)*cos(EV_(1));
            g_(2,1) = pbm.efpar(2)*pbm.efpar(9)*cos(EV_(2));
            g_(3,1) = pbm.efpar(3)*pbm.efpar(10)*cos(EV_(3));
            g_(4,1) = pbm.efpar(4)*pbm.efpar(11)*cos(EV_(4));
            g_(5,1) = pbm.efpar(5)*pbm.efpar(12)*cos(EV_(5));
            g_(6,1) = pbm.efpar(6)*pbm.efpar(13)*cos(EV_(6));
            g_(7,1) = pbm.efpar(7)*pbm.efpar(14)*cos(EV_(7));
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(7,7);
                H_(1,1) = -pbm.efpar(1)*pbm.efpar(8)*sin(EV_(1));
                H_(2,2) = -pbm.efpar(2)*pbm.efpar(9)*sin(EV_(2));
                H_(3,3) = -pbm.efpar(3)*pbm.efpar(10)*sin(EV_(3));
                H_(4,4) = -pbm.efpar(4)*pbm.efpar(11)*sin(EV_(4));
                H_(5,5) = -pbm.efpar(5)*pbm.efpar(12)*sin(EV_(5));
                H_(6,6) = -pbm.efpar(6)*pbm.efpar(13)*sin(EV_(6));
                H_(7,7) = -pbm.efpar(7)*pbm.efpar(14)*sin(EV_(7));
                varargout{3} = H_;
            end
        end

    case 'eQ8T'

        EV_  = varargin{1};
        iel_ = varargin{2};
        S2 = pbm.efpar(8)*(pbm.efpar(1)*sin(EV_(1))-pbm.efpar(15));
        S3 = pbm.efpar(9)*(pbm.efpar(2)*sin(EV_(2))-pbm.efpar(15))+S2;
        S4 = pbm.efpar(10)*(pbm.efpar(3)*sin(EV_(3))-pbm.efpar(15))+S3;
        S5 = pbm.efpar(11)*(pbm.efpar(4)*sin(EV_(4))-pbm.efpar(15))+S4;
        S6 = pbm.efpar(12)*(pbm.efpar(5)*sin(EV_(5))-pbm.efpar(15))+S5;
        S7 = pbm.efpar(13)*(pbm.efpar(6)*sin(EV_(6))-pbm.efpar(15))+S6;
        DSD1 = pbm.efpar(1)*pbm.efpar(8)*cos(EV_(1));
        DSD2 = pbm.efpar(2)*pbm.efpar(9)*cos(EV_(2));
        DSD3 = pbm.efpar(3)*pbm.efpar(10)*cos(EV_(3));
        DSD4 = pbm.efpar(4)*pbm.efpar(11)*cos(EV_(4));
        DSD5 = pbm.efpar(5)*pbm.efpar(12)*cos(EV_(5));
        DSD6 = pbm.efpar(6)*pbm.efpar(13)*cos(EV_(6));
        DSD7 = pbm.efpar(7)*pbm.efpar(14)*cos(EV_(7));
        D2SD1 = -pbm.efpar(1)*pbm.efpar(8)*sin(EV_(1));
        D2SD2 = -pbm.efpar(2)*pbm.efpar(9)*sin(EV_(2));
        D2SD3 = -pbm.efpar(3)*pbm.efpar(10)*sin(EV_(3));
        D2SD4 = -pbm.efpar(4)*pbm.efpar(11)*sin(EV_(4));
        D2SD5 = -pbm.efpar(5)*pbm.efpar(12)*sin(EV_(5));
        D2SD6 = -pbm.efpar(6)*pbm.efpar(13)*sin(EV_(6));
        D2SD7 = -pbm.efpar(7)*pbm.efpar(14)*sin(EV_(7));
        Q2 = 0.5*pbm.efpar(8)*pbm.efpar(8)*(pbm.efpar(1)*sin(EV_(1))-pbm.efpar(15));
        Q3 = 0.5*pbm.efpar(9)*pbm.efpar(9)*(pbm.efpar(2)*sin(EV_(2))-pbm.efpar(15))+...
             pbm.efpar(9)*S2+Q2;
        Q4 = 0.5*pbm.efpar(10)*pbm.efpar(10)*(pbm.efpar(3)*sin(EV_(3))-pbm.efpar(15))+...
             pbm.efpar(10)*S3+Q3;
        Q5 = 0.5*pbm.efpar(11)*pbm.efpar(11)*(pbm.efpar(4)*sin(EV_(4))-pbm.efpar(15))+...
             pbm.efpar(11)*S4+Q4;
        Q6 = 0.5*pbm.efpar(12)*pbm.efpar(12)*(pbm.efpar(5)*sin(EV_(5))-pbm.efpar(15))+...
             pbm.efpar(12)*S5+Q5;
        Q7 = 0.5*pbm.efpar(13)*pbm.efpar(13)*(pbm.efpar(6)*sin(EV_(6))-pbm.efpar(15))+...
             pbm.efpar(13)*S6+Q6;
        varargout{1} =...
              0.5*pbm.efpar(14)*pbm.efpar(14)*(pbm.efpar(7)*sin(EV_(7))-pbm.efpar(15))+pbm.efpar(14)*S7+Q7;
        if(nargout>1)
            g_(1,1) = 0.5*pbm.efpar(8)*pbm.efpar(8)*pbm.efpar(1)*cos(EV_(1))+...
                 (pbm.efpar(14)+pbm.efpar(13)+pbm.efpar(12)+pbm.efpar(11)+pbm.efpar(10)+pbm.efpar(9))*DSD1;
            g_(2,1) = 0.5*pbm.efpar(9)*pbm.efpar(9)*pbm.efpar(2)*cos(EV_(2))+...
                 (pbm.efpar(14)+pbm.efpar(13)+pbm.efpar(12)+pbm.efpar(11)+pbm.efpar(10))*DSD2;
            g_(3,1) = 0.5*pbm.efpar(10)*pbm.efpar(10)*pbm.efpar(3)*cos(EV_(3))+...
                 (pbm.efpar(14)+pbm.efpar(13)+pbm.efpar(12)+pbm.efpar(11))*DSD3;
            g_(4,1) = 0.5*pbm.efpar(11)*pbm.efpar(11)*pbm.efpar(4)*cos(EV_(4))+...
                 (pbm.efpar(14)+pbm.efpar(13)+pbm.efpar(12))*DSD4;
            g_(5,1) = 0.5*pbm.efpar(12)*pbm.efpar(12)*pbm.efpar(5)*cos(EV_(5))+...
                 (pbm.efpar(14)+pbm.efpar(13))*DSD5;
            g_(6,1) = 0.5*pbm.efpar(13)*pbm.efpar(13)*pbm.efpar(6)*cos(EV_(6))+...
                 pbm.efpar(14)*DSD6;
            g_(7,1) = 0.5*pbm.efpar(14)*pbm.efpar(14)*pbm.efpar(7)*cos(EV_(7));
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(7,7);
                H_(1,1) = -0.5*pbm.efpar(8)*pbm.efpar(8)*pbm.efpar(1)*sin(EV_(1))+...
                     (pbm.efpar(14)+pbm.efpar(13)+pbm.efpar(12)+pbm.efpar(11)+pbm.efpar(10)+pbm.efpar(9))*D2SD1;
                H_(2,2) = -0.5*pbm.efpar(9)*pbm.efpar(9)*pbm.efpar(2)*sin(EV_(2))+...
                     (pbm.efpar(14)+pbm.efpar(13)+pbm.efpar(12)+pbm.efpar(11)+pbm.efpar(10))*D2SD2;
                H_(3,3) = -0.5*pbm.efpar(10)*pbm.efpar(10)*pbm.efpar(3)*sin(EV_(3))+...
                     (pbm.efpar(14)+pbm.efpar(13)+pbm.efpar(12)+pbm.efpar(11))*D2SD3;
                H_(4,4) = -0.5*pbm.efpar(11)*pbm.efpar(11)*pbm.efpar(4)*sin(EV_(4))+...
                     (pbm.efpar(14)+pbm.efpar(13)+pbm.efpar(12))*D2SD4;
                H_(5,5) = -0.5*pbm.efpar(12)*pbm.efpar(12)*pbm.efpar(5)*sin(EV_(5))+...
                     (pbm.efpar(14)+pbm.efpar(13))*D2SD5;
                H_(6,6) = -0.5*pbm.efpar(13)*pbm.efpar(13)*pbm.efpar(6)*sin(EV_(6))+...
                     pbm.efpar(14)*D2SD6;
                H_(7,7) = -0.5*pbm.efpar(14)*pbm.efpar(14)*pbm.efpar(7)*sin(EV_(7));
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
            pbm.has_globs = [15,0];
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

