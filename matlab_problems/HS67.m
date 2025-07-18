function varargout = HS67(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HS67
%    *********
% 
%    Source: problem 67 in
%    W. Hock and K. Schittkowski,
%    "Test examples for nonlinear programming codes",
%    Lectures Notes in Economics and Mathematical Systems 187, Springer
%    Verlag, Heidelberg, 1981.
% 
%    Original Source: problem 8 in
%    A.R. Colville
%    "A comparative study on nonlinear programming"
%    IBM Scientific Center Report 320-2949, New York, 1968.
% 
%    SIF input: A.R. Conn & Nick Gould, April 1991.
% 
%    classification = 'C-COOI2-AN-3-14'
% 
%    Set useful parameters
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 21 VI 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HS67';

switch(action)

    case {'setup','setup_redprec'}

        if(strcmp(action,'setup_redprec'))
            if(isfield(pbm,'ndigs'))
                rmfield(pbm,'ndigs');
            end
            pbm.ndigs = max(1,min(15,varargin{end}));
            nargs     = nargin-2;
        else
            nargs = nargin-1;
        end
        pb.name   = name;
        pbm.name  = name;
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = containers.Map('KeyType','char', 'ValueType', 'double');
        ix_ = containers.Map('KeyType','char', 'ValueType', 'double');
        ig_ = containers.Map('KeyType','char', 'ValueType', 'double');
        v_('1') = 1;
        v_('2') = 2;
        v_('3') = 3;
        v_('4') = 4;
        v_('5') = 5;
        v_('6') = 6;
        v_('7') = 7;
        v_('8') = 8;
        v_('9') = 9;
        v_('10') = 10;
        v_('11') = 11;
        v_('12') = 12;
        v_('13') = 13;
        v_('14') = 14;
        v_(['A',int2str(round(v_('1')))]) = 0.0;
        v_(['A',int2str(round(v_('2')))]) = 0.0;
        v_(['A',int2str(round(v_('3')))]) = 85.0;
        v_(['A',int2str(round(v_('4')))]) = 90.0;
        v_(['A',int2str(round(v_('5')))]) = 3.0;
        v_(['A',int2str(round(v_('6')))]) = 0.01;
        v_(['A',int2str(round(v_('7')))]) = 145.0;
        v_(['A',int2str(round(v_('8')))]) = 5000.0;
        v_(['A',int2str(round(v_('9')))]) = 2000.0;
        v_(['A',int2str(round(v_('10')))]) = 93.0;
        v_(['A',int2str(round(v_('11')))]) = 95.0;
        v_(['A',int2str(round(v_('12')))]) = 12.0;
        v_(['A',int2str(round(v_('13')))]) = 4.0;
        v_(['A',int2str(round(v_('14')))]) = 162.0;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for I=v_('1'):v_('3')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['X',int2str(round(v_('1')))]);
        valA(end+1) = 5.04;
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['X',int2str(round(v_('2')))]);
        valA(end+1) = 0.035;
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['X',int2str(round(v_('3')))]);
        valA(end+1) = 10.0;
        for I=v_('1'):v_('7')
            [ig,ig_] = s2mpjlib('ii',['AG',int2str(I)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['AG',int2str(I)];
            [ig,ig_] = s2mpjlib('ii',['AL',int2str(I)],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['AL',int2str(I)];
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = ix_.Count;
        ngrp   = ig_.Count;
        legrps = find(strcmp(gtype,'<='));
        eqgrps = find(strcmp(gtype,'=='));
        gegrps = find(strcmp(gtype,'>='));
        pb.nle = length(legrps);
        pb.neq = length(eqgrps);
        pb.nge = length(gegrps);
        pb.m   = pb.nle+pb.neq+pb.nge;
        pbm.congrps = [ legrps, eqgrps, gegrps ];
        [pb.cnames{1:pb.m}] = deal(cnames{pbm.congrps});
        pb.nob = ngrp-pb.m;
        pbm.objgrps = find(strcmp(gtype,'<>'));
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        for I=v_('1'):v_('7')
            v_('I+7') = 7+I;
            pbm.gconst(ig_(['AG',int2str(I)])) = v_(['A',int2str(I)]);
            pbm.gconst(ig_(['AL',int2str(I)])) = v_(['A',int2str(round(v_('I+7')))]);
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = 0.00001*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xupper(ix_(['X',int2str(round(v_('1')))])) = 2000.0;
        pb.xupper(ix_(['X',int2str(round(v_('2')))])) = 16000.0;
        pb.xupper(ix_(['X',int2str(round(v_('3')))])) = 120.0;
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        pb.x0(ix_(['X',int2str(round(v_('1')))]),1) = 1745.0;
        pb.x0(ix_(['X',int2str(round(v_('2')))]),1) = 12000.0;
        pb.x0(ix_(['X',int2str(round(v_('3')))]),1) = 110.0;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eY2Y5',iet_);
        elftv{it}{1} = 'U1';
        elftv{it}{2} = 'U2';
        elftv{it}{3} = 'U3';
        [it,iet_] = s2mpjlib( 'ii', 'eY2',iet_);
        elftv{it}{1} = 'U1';
        elftv{it}{2} = 'U2';
        elftv{it}{3} = 'U3';
        [it,iet_] = s2mpjlib( 'ii', 'eY3',iet_);
        elftv{it}{1} = 'U1';
        elftv{it}{2} = 'U2';
        elftv{it}{3} = 'U3';
        [it,iet_] = s2mpjlib( 'ii', 'eY4',iet_);
        elftv{it}{1} = 'U1';
        elftv{it}{2} = 'U2';
        elftv{it}{3} = 'U3';
        [it,iet_] = s2mpjlib( 'ii', 'eY5',iet_);
        elftv{it}{1} = 'U1';
        elftv{it}{2} = 'U2';
        elftv{it}{3} = 'U3';
        [it,iet_] = s2mpjlib( 'ii', 'eY6',iet_);
        elftv{it}{1} = 'U1';
        elftv{it}{2} = 'U2';
        elftv{it}{3} = 'U3';
        [it,iet_] = s2mpjlib( 'ii', 'eY7',iet_);
        elftv{it}{1} = 'U1';
        elftv{it}{2} = 'U2';
        elftv{it}{3} = 'U3';
        [it,iet_] = s2mpjlib( 'ii', 'eY8',iet_);
        elftv{it}{1} = 'U1';
        elftv{it}{2} = 'U2';
        elftv{it}{3} = 'U3';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = 'E25';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eY2Y5';
        ielftype(ie) = iet_('eY2Y5');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eY2';
        ielftype(ie) = iet_('eY2');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E3';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eY3';
        ielftype(ie) = iet_('eY3');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E4';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eY4';
        ielftype(ie) = iet_('eY4');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E5';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eY5';
        ielftype(ie) = iet_('eY5');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E6';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eY6';
        ielftype(ie) = iet_('eY6');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E7';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eY7';
        ielftype(ie) = iet_('eY7');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E8';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eY8';
        ielftype(ie) = iet_('eY8');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        ig = ig_('OBJ');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E25');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -0.063;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E3');
        pbm.grelw{ig}(posel) = 3.36;
        ig = ig_(['AG',int2str(round(v_('1')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E2');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_(['AL',int2str(round(v_('1')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E2');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_(['AG',int2str(round(v_('2')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E3');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_(['AL',int2str(round(v_('2')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E3');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_(['AG',int2str(round(v_('3')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E4');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_(['AL',int2str(round(v_('3')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E4');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_(['AG',int2str(round(v_('4')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E5');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_(['AL',int2str(round(v_('4')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E5');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_(['AG',int2str(round(v_('5')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E6');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_(['AL',int2str(round(v_('5')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E6');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_(['AG',int2str(round(v_('6')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E7');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_(['AL',int2str(round(v_('6')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E7');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_(['AG',int2str(round(v_('7')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E8');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_(['AL',int2str(round(v_('7')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E8');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-COOI2-AN-3-14';
        pbm.objderlvl = 2;
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2];
        pb.conderlvl  = pbm.conderlvl;
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb, pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm,pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end


    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'e_globs'

        pbm = varargin{1};
        pbm.efpar = [];
        pbm.efpar(1) = .FALSE.0;
        varargout{1} = pbm;

    case 'eY2Y5'

        EV_  = varargin{1};
        iel_ = varargin{2};
        if(~pbm.efpar(1))
            DUMMY = HS67(EV_(1),EV_(2),EV_(3),Y,G,H);
        end
        EVAL = .TRUE.0;
        varargout{1} = Y(2)*Y(5);
        if(nargout>1)
            g_(1,1) = Y(2)*G(5,1)+Y(5)*G(2,1);
            g_(2,1) = Y(2)*G(5,2)+Y(5)*G(2,2);
            g_(3,1) = Y(2)*G(5,3)+Y(5)*G(2,3);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,1) = Y(2)*H(5,1)+2.0*G(5,1)*G(2,1)+Y(5)*H(2,1);
                H_(1,2) = Y(2)*H(5,2)+G(5,1)*G(2,2)+G(5,2)*G(2,1)+Y(5)*H(2,2);
                H_(2,1) = H_(1,2);
                H_(1,3) = Y(2)*H(5,3)+G(5,1)*G(2,3)+G(5,3)*G(2,1)+Y(5)*H(2,3);
                H_(3,1) = H_(1,3);
                H_(2,2) = Y(2)*H(5,4)+2.0*G(5,2)*G(2,2)+Y(5)*H(2,4);
                H_(2,3) = Y(2)*H(5,5)+G(5,2)*G(2,3)+G(5,3)*G(2,2)+Y(5)*H(2,5);
                H_(3,2) = H_(2,3);
                H_(3,3) = Y(2)*H(5,6)+2.0*G(5,3)*G(2,3)+Y(5)*H(2,6);
                varargout{3} = H_;
            end
        end

    case 'eY2'

        EV_  = varargin{1};
        iel_ = varargin{2};
        if(~EVAL)
            DUMMY = HS67(EV_(1),EV_(2),EV_(3),Y,G,H);
        end
        EVAL = .TRUE.0;
        varargout{1} = Y(2);
        if(nargout>1)
            g_(1,1) = G(2,1);
            g_(2,1) = G(2,2);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = H(2,1);
                H_(2,1) = H(2,2);
                H_(1,2) = H_(2,1);
                H_(2,2) = H(2,4);
                varargout{3} = H_;
            end
        end

    case 'eY3'

        EV_  = varargin{1};
        iel_ = varargin{2};
        if(~EVAL)
            DUMMY = HS67(EV_(1),EV_(2),EV_(3),Y,G,H);
        end
        EVAL = .TRUE.0;
        varargout{1} = Y(3);
        if(nargout>1)
            g_(1,1) = G(3,1);
            g_(2,1) = G(3,2);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = H(3,1);
                H_(2,1) = H(3,2);
                H_(1,2) = H_(2,1);
                H_(2,2) = H(3,4);
                varargout{3} = H_;
            end
        end

    case 'eY4'

        EV_  = varargin{1};
        iel_ = varargin{2};
        if(~EVAL)
            DUMMY = HS67(EV_(1),EV_(2),EV_(3),Y,G,H);
        end
        EVAL = .TRUE.0;
        varargout{1} = Y(4);
        if(nargout>1)
            g_(1,1) = G(4,1);
            g_(2,1) = G(4,2);
            g_(3,1) = G(4,3);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,1) = H(4,1);
                H_(2,1) = H(4,2);
                H_(1,2) = H_(2,1);
                H_(3,1) = H(4,3);
                H_(1,3) = H_(3,1);
                H_(2,2) = H(4,4);
                H_(3,2) = H(4,5);
                H_(2,3) = H_(3,2);
                H_(3,3) = H(4,6);
                varargout{3} = H_;
            end
        end

    case 'eY5'

        EV_  = varargin{1};
        iel_ = varargin{2};
        if(~EVAL)
            DUMMY = HS67(EV_(1),EV_(2),EV_(3),Y,G,H);
        end
        EVAL = .TRUE.0;
        varargout{1} = Y(5);
        if(nargout>1)
            g_(1,1) = G(5,1);
            g_(2,1) = G(5,2);
            g_(3,1) = G(5,3);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,1) = H(5,1);
                H_(2,1) = H(5,2);
                H_(1,2) = H_(2,1);
                H_(3,1) = H(5,3);
                H_(1,3) = H_(3,1);
                H_(2,2) = H(5,4);
                H_(3,2) = H(5,5);
                H_(2,3) = H_(3,2);
                H_(3,3) = H(5,6);
                varargout{3} = H_;
            end
        end

    case 'eY6'

        EV_  = varargin{1};
        iel_ = varargin{2};
        if(~EVAL)
            DUMMY = HS67(EV_(1),EV_(2),EV_(3),Y,G,H);
        end
        EVAL = .TRUE.0;
        varargout{1} = Y(6);
        if(nargout>1)
            g_(1,1) = G(6,1);
            g_(2,1) = G(6,2);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = H(6,1);
                H_(2,1) = H(6,2);
                H_(1,2) = H_(2,1);
                H_(2,2) = H(6,4);
                varargout{3} = H_;
            end
        end

    case 'eY7'

        EV_  = varargin{1};
        iel_ = varargin{2};
        if(~EVAL)
            DUMMY = HS67(EV_(1),EV_(2),EV_(3),Y,G,H);
        end
        EVAL = .TRUE.0;
        varargout{1} = Y(7);
        if(nargout>1)
            g_(1,1) = G(7,1);
            g_(2,1) = G(7,2);
            g_(3,1) = G(7,3);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,1) = H(7,1);
                H_(2,1) = H(7,2);
                H_(1,2) = H_(2,1);
                H_(3,1) = H(7,3);
                H_(1,3) = H_(3,1);
                H_(2,2) = H(7,4);
                H_(3,2) = H(7,5);
                H_(2,3) = H_(3,2);
                H_(3,3) = H(7,6);
                varargout{3} = H_;
            end
        end

    case 'eY8'

        EV_  = varargin{1};
        iel_ = varargin{2};
        if(~EVAL)
            DUMMY = HS67(EV_(1),EV_(2),EV_(3),Y,G,H);
        end
        varargout{1} = Y(8);
        if(nargout>1)
            g_(1,1) = G(8,1);
            g_(2,1) = G(8,2);
            g_(3,1) = G(8,3);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,1) = H(8,1);
                H_(2,1) = H(8,2);
                H_(1,2) = H_(2,1);
                H_(3,1) = H(8,3);
                H_(1,3) = H_(3,1);
                H_(2,2) = H(8,4);
                H_(3,2) = H(8,5);
                H_(2,3) = H_(3,2);
                H_(3,3) = H(8,6);
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    case {'fx','fgx','fgHx','cx','cJx','cJHx','cIx','cIJx','cIJHx','cIJxv','fHxv',...
          'cJxv','cJtxv','cIJtxv','Lxy','Lgxy','LgHxy','LIxy','LIgxy','LIgHxy',...
          'LHxyv','LIHxyv'}

        if(isfield(pbm,'name')&&strcmp(pbm.name,name))
            pbm.has_globs = [1,0];
            [varargout{1:max(1,nargout)}] = s2mpjlib(action,pbm,varargin{:});
        else
            disp(['ERROR: please run ',name,' with action = setup'])
            [varargout{1:nargout}] = deal(NaN);
        end

    otherwise
        disp([' ERROR: action ',action,' unavailable for problem ',name,'.m'])
    end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

