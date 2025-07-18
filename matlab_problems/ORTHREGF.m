function varargout = ORTHREGF(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : ORTHREGF
%    *********
% 
%    An orthogonal regression problem
% 
%    The problem is to fit (orthogonally) an torus to a
%    set of points in 3D space. This set of points is generated by
%    perturbing a first set lying exactly on a predefined torus
%    centered at the origin.
% 
%    Source:
%    M. Gulliksson,
%    "Algorithms for nonlinear Least-squares with Applications to
%    Orthogonal Regression",
%    UMINF-178.90, University of Umea, Sweden, 1990.
% 
%    SIF input: Ph. Toint, June 1990.
%               minor correction by Ph. Shott, Jan 1995.
% 
%    classification = 'C-CQOR2-AY-V-V'
% 
%    square root of the number of data points
%    (number of variables = 3 * NPTS**2 + 5 )
% 
%       Alternative values for the SIF file parameters:
% IE NPTS                5              $-PARAMETER n = 80    original value
% IE NPTS                7              $-PARAMETER n = 152
% IE NPTS                10             $-PARAMETER n = 305
% IE NPTS                15             $-PARAMETER n = 680
% IE NPTS                20             $-PARAMETER n = 1205
% IE NPTS                40             $-PARAMETER n = 4805
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 21 VI 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'ORTHREGF';

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
        if(nargs<1)
            v_('NPTS') = 5;  %  SIF file default value
        else
            v_('NPTS') = varargin{1};
        end
        v_('TP4') = 1.7;
        v_('TP5') = 0.8;
        v_('PSEED') = 237.1531;
        v_('PSIZE') = 0.2;
        v_('1') = 1;
        v_('5') = 5;
        v_('PI') = 3.1415926535;
        v_('2PI') = 2.0*v_('PI');
        v_('RNPTS') = v_('NPTS');
        v_('ICR0') = 1.0/v_('RNPTS');
        v_('INCR') = v_('ICR0')*v_('2PI');
        for I=v_('1'):v_('NPTS')
            v_('I-1') = -1+I;
            v_('RI-1') = v_('I-1');
            v_('THETA1') = v_('RI-1')*v_('INCR');
            v_('ST1') = sin(v_('THETA1'));
            v_('CT1') = cos(v_('THETA1'));
            v_('P5CT1') = v_('TP5')*v_('CT1');
            v_('P4P5CT1') = v_('TP4')+v_('P5CT1');
            v_('R3') = v_('TP5')*v_('ST1');
            for J=v_('1'):v_('NPTS')
                v_('J-1') = -1+J;
                v_('RJ-1') = v_('J-1');
                v_('THETA2') = v_('RJ-1')*v_('INCR');
                v_('ST2') = sin(v_('THETA2'));
                v_('CT2') = cos(v_('THETA2'));
                v_('R1') = v_('P4P5CT1')*v_('CT2');
                v_('R2') = v_('P4P5CT1')*v_('ST2');
                v_('XSEED') = v_('THETA2')*v_('PSEED');
                v_('SSEED') = cos(v_('XSEED'));
                v_('PER-1') = v_('PSIZE')*v_('SSEED');
                v_('PERT') = 1.0+v_('PER-1');
                v_(['XD',int2str(I),',',int2str(J)]) = v_('R1')*v_('PERT');
                v_(['YD',int2str(I),',',int2str(J)]) = v_('R2')*v_('PERT');
                v_(['ZD',int2str(I),',',int2str(J)]) = v_('R3')*v_('PERT');
            end
        end
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for I=v_('1'):v_('5')
            [iv,ix_] = s2mpjlib('ii',['P',int2str(I)],ix_);
            pb.xnames{iv} = ['P',int2str(I)];
        end
        for I=v_('1'):v_('NPTS')
            for J=v_('1'):v_('NPTS')
                [iv,ix_] = s2mpjlib('ii',['X',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['X',int2str(I),',',int2str(J)];
                [iv,ix_] = s2mpjlib('ii',['Y',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['Y',int2str(I),',',int2str(J)];
                [iv,ix_] = s2mpjlib('ii',['Z',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['Z',int2str(I),',',int2str(J)];
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        for I=v_('1'):v_('NPTS')
            for J=v_('1'):v_('NPTS')
                [ig,ig_] = s2mpjlib('ii',['OX',int2str(I),',',int2str(J)],ig_);
                gtype{ig} = '<>';
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['X',int2str(I),',',int2str(J)]);
                valA(end+1) = 1.0;
                [ig,ig_] = s2mpjlib('ii',['OY',int2str(I),',',int2str(J)],ig_);
                gtype{ig} = '<>';
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['Y',int2str(I),',',int2str(J)]);
                valA(end+1) = 1.0;
                [ig,ig_] = s2mpjlib('ii',['OZ',int2str(I),',',int2str(J)],ig_);
                gtype{ig} = '<>';
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['Z',int2str(I),',',int2str(J)]);
                valA(end+1) = 1.0;
                [ig,ig_] = s2mpjlib('ii',['A',int2str(I),',',int2str(J)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['A',int2str(I),',',int2str(J)];
            end
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
        for I=v_('1'):v_('NPTS')
            for J=v_('1'):v_('NPTS')
                pbm.gconst(ig_(['OX',int2str(I),',',int2str(J)])) =...
                      v_(['XD',int2str(I),',',int2str(J)]);
                pbm.gconst(ig_(['OY',int2str(I),',',int2str(J)])) =...
                      v_(['YD',int2str(I),',',int2str(J)]);
                pbm.gconst(ig_(['OZ',int2str(I),',',int2str(J)])) =...
                      v_(['ZD',int2str(I),',',int2str(J)]);
            end
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower(ix_('P4'),1) = 0.001;
        pb.xlower(ix_('P5'),1) = 0.001;
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,'P1'))
            pb.x0(ix_('P1'),1) = 1.0;
        else
            pb.y0(find(pbm.congrps==ig_('P1')),1) = 1.0;
        end
        if(isKey(ix_,'P2'))
            pb.x0(ix_('P2'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('P2')),1) = 0.0;
        end
        if(isKey(ix_,'P3'))
            pb.x0(ix_('P3'),1) = 1.0;
        else
            pb.y0(find(pbm.congrps==ig_('P3')),1) = 1.0;
        end
        if(isKey(ix_,'P4'))
            pb.x0(ix_('P4'),1) = 1.0;
        else
            pb.y0(find(pbm.congrps==ig_('P4')),1) = 1.0;
        end
        if(isKey(ix_,'P5'))
            pb.x0(ix_('P5'),1) = 0.5;
        else
            pb.y0(find(pbm.congrps==ig_('P5')),1) = 0.5;
        end
        for I=v_('1'):v_('NPTS')
            for J=v_('1'):v_('NPTS')
                if(isKey(ix_,['X',int2str(I),',',int2str(J)]))
                    pb.x0(ix_(['X',int2str(I),',',int2str(J)]),1) =...
                          v_(['XD',int2str(I),',',int2str(J)]);
                else
                    pb.y0(find(pbm.congrps==ig_(['X',int2str(I),',',int2str(J)])),1) =...
                          v_(['XD',int2str(I),',',int2str(J)]);
                end
                if(isKey(ix_,['Y',int2str(I),',',int2str(J)]))
                    pb.x0(ix_(['Y',int2str(I),',',int2str(J)]),1) =...
                          v_(['YD',int2str(I),',',int2str(J)]);
                else
                    pb.y0(find(pbm.congrps==ig_(['Y',int2str(I),',',int2str(J)])),1) =...
                          v_(['YD',int2str(I),',',int2str(J)]);
                end
                if(isKey(ix_,['Z',int2str(I),',',int2str(J)]))
                    pb.x0(ix_(['Z',int2str(I),',',int2str(J)]),1) =...
                          v_(['ZD',int2str(I),',',int2str(J)]);
                else
                    pb.y0(find(pbm.congrps==ig_(['Z',int2str(I),',',int2str(J)])),1) =...
                          v_(['ZD',int2str(I),',',int2str(J)]);
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eTA',iet_);
        elftv{it}{1} = 'XX';
        elftv{it}{2} = 'YY';
        elftv{it}{3} = 'A';
        elftv{it}{4} = 'B';
        elftv{it}{5} = 'C';
        [it,iet_] = s2mpjlib( 'ii', 'eISQ',iet_);
        elftv{it}{1} = 'Z';
        elftv{it}{2} = 'P';
        [it,iet_] = s2mpjlib( 'ii', 'eSQ',iet_);
        elftv{it}{1} = 'XX';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('NPTS')
            for J=v_('1'):v_('NPTS')
                ename = ['EA',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eTA';
                ielftype(ie) = iet_('eTA');
                vname = ['X',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('XX',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['Y',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('YY',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = 'P1';
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('A',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = 'P2';
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('B',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = 'P4';
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('C',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['EB',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eISQ';
                ielftype(ie) = iet_('eISQ');
                vname = ['Z',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('Z',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = 'P3';
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('P',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['EC',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eSQ';
                ielftype(ie) = iet_('eSQ');
                vname = 'P5';
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('XX',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('NPTS')
            for J=v_('1'):v_('NPTS')
                ig = ig_(['OX',int2str(I),',',int2str(J)]);
                pbm.grftype{ig} = 'gL2';
                ig = ig_(['OY',int2str(I),',',int2str(J)]);
                pbm.grftype{ig} = 'gL2';
                ig = ig_(['OZ',int2str(I),',',int2str(J)]);
                pbm.grftype{ig} = 'gL2';
                ig = ig_(['A',int2str(I),',',int2str(J)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['EA',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
                posel = posel+1;
                pbm.grelt{ig}(posel) = ie_(['EB',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = 1.;
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['EC',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = -1.0;
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
%    Solution
% LO SOLTN(5)            0.990089426
% LO SOLTN(7)            1.315031322
% LO SOLTN(10)           4.515848902
% LO SOLTN(15)           9.185538338
% LO SOLTN(20)           16.20054380
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CQOR2-AY-V-V';
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

% **********************
%  SET UP THE FUNCTION *
%  AND RANGE ROUTINES  *
% **********************

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eTA'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(3,5);
        U_(1,1) = U_(1,1)+1;
        U_(1,3) = U_(1,3)-1;
        U_(2,2) = U_(2,2)+1;
        U_(2,4) = U_(2,4)-1;
        U_(3,5) = U_(3,5)+1;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        IV_(3) = U_(3,:)*EV_;
        CCSQ = IV_(3)*IV_(3);
        CCCB = CCSQ*IV_(3);
        XXYY = IV_(1)*IV_(1)+IV_(2)*IV_(2);
        T = XXYY/CCSQ;
        DTDX = 2.0*IV_(1)/CCSQ;
        DTDY = 2.0*IV_(2)/CCSQ;
        DTDC = -2.0*XXYY/CCCB;
        D2TDX2 = 2.0/CCSQ;
        D2TDY2 = 2.0/CCSQ;
        D2TDC2 = 6.0*XXYY/(CCSQ*CCSQ);
        D2TDXC = -4.0*IV_(1)/CCCB;
        D2TDYC = -4.0*IV_(2)/CCCB;
        S = sqrt(T);
        R = 0.5/S;
        DSDX = R*DTDX;
        DSDY = R*DTDY;
        DSDC = R*DTDC;
        D2SDX2 = R*(D2TDX2-0.5*DTDX*DTDX/T);
        D2SDY2 = R*(D2TDY2-0.5*DTDY*DTDY/T);
        D2SDC2 = R*(D2TDC2-0.5*DTDC*DTDC/T);
        D2SDXY = -0.5*DTDX*DSDY/T;
        D2SDXC = R*(D2TDXC-0.5*DTDX*DTDC/T);
        D2SDYC = R*(D2TDYC-0.5*DTDY*DTDC/T);
        SS = S-1.0;
        SPS = SS+SS;
        varargout{1} = SS*SS;
        if(nargout>1)
            g_(1,1) = SPS*DSDX;
            g_(2,1) = SPS*DSDY;
            g_(3,1) = SPS*DSDC;
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,1) = SPS*D2SDX2+2.0*DSDX*DSDX;
                H_(1,2) = SPS*D2SDXY+2.0*DSDX*DSDY;
                H_(2,1) = H_(1,2);
                H_(1,3) = SPS*D2SDXC+2.0*DSDX*DSDC;
                H_(3,1) = H_(1,3);
                H_(2,2) = SPS*D2SDY2+2.0*DSDY*DSDY;
                H_(2,3) = SPS*D2SDYC+2.0*DSDY*DSDC;
                H_(3,2) = H_(2,3);
                H_(3,3) = SPS*D2SDC2+2.0*DSDC*DSDC;
                varargout{3} = U_.'*H_*U_;
            end
        end

    case 'eISQ'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(1,2);
        U_(1,1) = U_(1,1)+1;
        U_(1,2) = U_(1,2)-1;
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

    case 'eSQ'

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
          'cJxv','cJtxv','cIJtxv','Lxy','Lgxy','LgHxy','LIxy','LIgxy','LIgHxy',...
          'LHxyv','LIHxyv'}

        if(isfield(pbm,'name')&&strcmp(pbm.name,name))
            pbm.has_globs = [0,0];
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

