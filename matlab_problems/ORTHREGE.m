function varargout = ORTHREGE(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : ORTHREGE
%    *********
% 
%    An orthogonal regression problem,
% 
%    The problem is to fit (orthogonally) an elliptic helix to a
%    set of points in 3D space. This set of points is generated by
%    perturbing a first set lying exactly on a predifined helix
%    centered at the origin.
% 
%    Source:
%    M. Gulliksson,
%    "Algorithms for nonlinear Least-squares with Applications to
%    Orthogonal Regression",
%    UMINF-178.90, University of Umea, Sweden, 1990.
% 
%    SIF input: Ph. Toint, June 1990.
% 
%    classification = 'QOR2-AY-V-V'
% 
%    Number of data points
%    (number of variables = 3 NPTS + 6 )
% 
%       Alternative values for the SIF file parameters:
% IE NPTS                10             $-PARAMETER n= 36    original valu
% IE NPTS                2500           $-PARAMETER n= 7506
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'ORTHREGE';

switch(action)

    case {'setup','setup_redprec'}

        if(isfield(pbm,'ndigs'))
            rmfield(pbm,'ndigs');
        end
        if(strcmp(action,'setup_redprec'))
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
            v_('NPTS') = 10;  %  SIF file default value
        else
            v_('NPTS') = varargin{1};
        end
        v_('TP4') = 1.7;
        v_('TP5') = 0.8;
        v_('TP6') = 2.0;
        v_('PSEED') = 237.1531;
        v_('PSIZE') = 0.2;
        v_('1') = 1;
        v_('6') = 6;
        v_('PI') = 3.1415926535;
        v_('2PI') = 2.0*v_('PI');
        v_('RNPTS') = v_('NPTS');
        v_('ICR0') = 1.0/v_('RNPTS');
        v_('INCR') = v_('ICR0')*v_('2PI');
        for I=v_('1'):v_('NPTS')
            v_('I-1') = -1+I;
            v_('RI-1') = v_('I-1');
            v_('THETA') = v_('RI-1')*v_('INCR');
            v_('ST') = sin(v_('THETA'));
            v_('CT') = cos(v_('THETA'));
            v_('R1') = v_('TP4')*v_('CT');
            v_('R2') = v_('TP5')*v_('ST');
            v_('R3') = v_('TP6')*v_('THETA');
            v_('XSEED') = v_('THETA')*v_('PSEED');
            v_('SSEED') = cos(v_('XSEED'));
            v_('PER-1') = v_('PSIZE')*v_('SSEED');
            v_('PERT') = 1.0+v_('PER-1');
            v_(['XD',int2str(I)]) = v_('R1')*v_('PERT');
            v_(['YD',int2str(I)]) = v_('R2')*v_('PERT');
            v_(['ZD',int2str(I)]) = v_('R3')*v_('PERT');
        end
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('6')
            [iv,ix_] = s2mpjlib('ii',['P',int2str(I)],ix_);
            pb.xnames{iv} = ['P',int2str(I)];
        end
        for I=v_('1'):v_('NPTS')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
            [iv,ix_] = s2mpjlib('ii',['Y',int2str(I)],ix_);
            pb.xnames{iv} = ['Y',int2str(I)];
            [iv,ix_] = s2mpjlib('ii',['Z',int2str(I)],ix_);
            pb.xnames{iv} = ['Z',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('NPTS')
            [ig,ig_] = s2mpjlib('ii',['OX',int2str(I)],ig_);
            gtype{ig} = '<>';
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['OY',int2str(I)],ig_);
            gtype{ig} = '<>';
            iv = ix_(['Y',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['OZ',int2str(I)],ig_);
            gtype{ig} = '<>';
            iv = ix_(['Z',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['A',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['A',int2str(I)];
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_('P1');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['B',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['B',int2str(I)];
            iv = ix_(['Y',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_('P2');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
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
            pbm.gconst(ig_(['OX',int2str(I)])) = v_(['XD',int2str(I)]);
            pbm.gconst(ig_(['OY',int2str(I)])) = v_(['YD',int2str(I)]);
            pbm.gconst(ig_(['OZ',int2str(I)])) = v_(['ZD',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        pb.xlower(ix_('P6'),1) = 0.001;
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
            pb.x0(ix_('P5'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('P5')),1) = 0.0;
        end
        if(isKey(ix_,'P6'))
            pb.x0(ix_('P6'),1) = 0.25;
        else
            pb.y0(find(pbm.congrps==ig_('P6')),1) = 0.25;
        end
        for I=v_('1'):v_('NPTS')
            if(isKey(ix_,['X',int2str(I)]))
                pb.x0(ix_(['X',int2str(I)]),1) = v_(['XD',int2str(I)]);
            else
                pb.y0(find(pbm.congrps==ig_(['X',int2str(I)])),1) = v_(['XD',int2str(I)]);
            end
            if(isKey(ix_,['Y',int2str(I)]))
                pb.x0(ix_(['Y',int2str(I)]),1) = v_(['YD',int2str(I)]);
            else
                pb.y0(find(pbm.congrps==ig_(['Y',int2str(I)])),1) = v_(['YD',int2str(I)]);
            end
            if(isKey(ix_,['Z',int2str(I)]))
                pb.x0(ix_(['Z',int2str(I)]),1) = v_(['ZD',int2str(I)]);
            else
                pb.y0(find(pbm.congrps==ig_(['Z',int2str(I)])),1) = v_(['ZD',int2str(I)]);
            end
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eEHX',iet_);
        elftv{it}{1} = 'Z';
        elftv{it}{2} = 'A';
        elftv{it}{3} = 'B';
        elftv{it}{4} = 'C';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('NPTS')
            ename = ['EA',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eEHX';
            ielftype(ie) = iet_('eEHX');
            vname = ['Z',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('Z',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'P4';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('A',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'P3';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('B',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'P6';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('C',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['EB',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eEHX';
            ielftype(ie) = iet_('eEHX');
            vname = ['Z',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('Z',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'P5';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('A',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'P3';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('B',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'P6';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('C',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('NPTS')
            ig = ig_(['OX',int2str(I)]);
            pbm.grftype{ig} = 'gL2';
            ig = ig_(['OY',int2str(I)]);
            pbm.grftype{ig} = 'gL2';
            ig = ig_(['OZ',int2str(I)]);
            pbm.grftype{ig} = 'gL2';
            ig = ig_(['A',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['EA',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = -1.0;
            ig = ig_(['B',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['EB',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = -1.0;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
%    Solution
% LO SOLTN(10)           ???
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'QOR2-AY-V-V';
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

    case 'eEHX'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(3,4);
        U_(2,1) = U_(2,1)+1;
        U_(2,3) = U_(2,3)-1;
        U_(1,2) = U_(1,2)+1;
        U_(3,4) = U_(3,4)+1;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        IV_(3) = U_(3,:)*EV_;
        T = IV_(2)/IV_(3);
        CT = cos(T);
        ST = sin(T);
        CCSQ = IV_(3)*IV_(3);
        CCCB = CCSQ*IV_(3);
        ACT = IV_(1)*CT;
        AST = IV_(1)*ST;
        varargout{1} = ACT;
        if(nargout>1)
            g_(1,1) = CT;
            g_(2,1) = -AST/IV_(3);
            g_(3,1) = AST*IV_(2)/CCSQ;
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,2) = -ST/IV_(3);
                H_(2,1) = H_(1,2);
                H_(1,3) = ST*IV_(2)/CCSQ;
                H_(3,1) = H_(1,3);
                H_(2,2) = -ACT/CCSQ;
                H_(2,3) = (AST+IV_(2)*ACT/IV_(3))/CCSQ;
                H_(3,2) = H_(2,3);
                H_(3,3) = -IV_(2)*(IV_(2)*ACT/IV_(3)+AST+AST)/CCCB;
                varargout{3} = U_.'*H_*U_;
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
            [varargout{1:max(1,nargout)}] = s2mpjlib(action,pbm,varargin{:});
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

