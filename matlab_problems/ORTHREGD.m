function varargout = ORTHREGD(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : ORTHREGD
%    *********
% 
%    An orthogonal regression problem.
% 
%    The problem is to fit (orthogonally) a planar curve to a set of points
%    in the plane. This set of points is generated by perturbing a
%    first set lying exactly on the predefined curbe.  
%    The curve is referred to as a cardioid in the original paper, 
%    but is in fact a circle.
% 
%    Source: adapted from:
%    M. Gulliksson,
%    "Algorithms for nonlinear Least-squares with Applications to
%    Orthogonal Regression",
%    UMINF-178.90, University of Umea, Sweden, 1990.
% 
%    SIF input: Ph. Toint, Mar 1991.
%               minor correction by Ph. Shott, Jan 1995.
% 
%    classification = 'QOR2-AY-V-V'
% 
%    Number of data points
%    (number of variables = 2 NPTS + 3 )
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'ORTHREGD';

switch(action)

    case 'setup'

    pb.name      = 'ORTHREGD';
    pb.sifpbname = 'ORTHREGD';
    pbm.name     = 'ORTHREGD';
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        if(nargin<2)
            v_('NPTS') = 10;  %  SIF file default value
        else
            v_('NPTS') = varargin{1};
        end
%       Alternative values for the SIF file parameters:
% IE NPTS                50             $-PARAMETER n= 105
% IE NPTS                250            $-PARAMETER n= 505
% IE NPTS                500            $-PARAMETER n= 1005
% IE NPTS                2500           $-PARAMETER n= 5005
% IE NPTS                5000           $-PARAMETER n= 10005
        v_('TZ3') = 1.7;
        v_('PSEED') = 237.1531;
        v_('PSIZE') = 0.2;
        v_('1') = 1;
        v_('0') = 0;
        v_('PI') = 3.1415926535;
        v_('2PI') = 2.0*v_('PI');
        v_('RNPTS') = v_('NPTS');
        v_('ICR0') = 1.0/v_('RNPTS');
        v_('INCR') = v_('ICR0')*v_('2PI');
        v_('Z3SQ') = v_('TZ3')*v_('TZ3');
        v_('1+TZ3SQ') = 1.0+v_('Z3SQ');
        for I=v_('1'):v_('NPTS')
            v_('I-1') = -1+I;
            v_('RI-1') = v_('I-1');
            v_('THETA') = v_('RI-1')*v_('INCR');
            v_('ST') = sin(v_('THETA'));
            v_('CT') = cos(v_('THETA'));
            v_('FACT') = v_('1+TZ3SQ')+v_('CT');
            v_('R1') = v_('FACT')*v_('CT');
            v_('R2') = v_('FACT')*v_('ST');
            v_('XSEED') = v_('THETA')*v_('PSEED');
            v_('SSEED') = cos(v_('XSEED'));
            v_('PER-1') = v_('PSIZE')*v_('SSEED');
            v_('PERT') = 1.0+v_('PER-1');
            v_(['XD',int2str(I)]) = v_('R1')*v_('PERT');
            v_(['YD',int2str(I)]) = v_('R2')*v_('PERT');
        end
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2xlib('ii','Z1',ix_);
        pb.xnames{iv} = 'Z1';
        [iv,ix_] = s2xlib('ii','Z2',ix_);
        pb.xnames{iv} = 'Z2';
        [iv,ix_] = s2xlib('ii','Z3',ix_);
        pb.xnames{iv} = 'Z3';
        for I=v_('1'):v_('NPTS')
            [iv,ix_] = s2xlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
            [iv,ix_] = s2xlib('ii',['Y',int2str(I)],ix_);
            pb.xnames{iv} = ['Y',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('NPTS')
            [ig,ig_] = s2xlib('ii',['OX',int2str(I)],ig_);
            gtype{ig} = '<>';
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            [ig,ig_] = s2xlib('ii',['OY',int2str(I)],ig_);
            gtype{ig} = '<>';
            iv = ix_(['Y',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            [ig,ig_] = s2xlib('ii',['E',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['E',int2str(I)];
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
        for I=v_('1'):v_('NPTS')
            pbm.gconst(ig_(['OX',int2str(I)])) = v_(['XD',int2str(I)]);
            pbm.gconst(ig_(['OY',int2str(I)])) = v_(['YD',int2str(I)]);
        end
        pb.xlower = zeros(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,'Z1'))
            pb.x0(ix_('Z1'),1) = 1.0;
        else
            pb.y0(find(pbm.congrps==ig_('Z1')),1) = 1.0;
        end
        if(isKey(ix_,'Z2'))
            pb.x0(ix_('Z2'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('Z2')),1) = 0.0;
        end
        if(isKey(ix_,'Z3'))
            pb.x0(ix_('Z3'),1) = 1.0;
        else
            pb.y0(find(pbm.congrps==ig_('Z3')),1) = 1.0;
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
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2xlib( 'ii', 'eTA',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        elftv{it}{3} = 'ZA';
        elftv{it}{4} = 'ZB';
        [it,iet_] = s2xlib( 'ii', 'eTB',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        elftv{it}{3} = 'ZA';
        elftv{it}{4} = 'ZB';
        elftv{it}{5} = 'ZC';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('NPTS')
            ename = ['EA',int2str(I)];
            [ie,ie_] = s2xlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eTA';
            ielftype(ie) = iet_('eTA');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['Y',int2str(I)];
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'Z1';
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('ZA',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'Z2';
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('ZB',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['EB',int2str(I)];
            [ie,ie_] = s2xlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eTB';
            ielftype(ie) = iet_('eTB');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['Y',int2str(I)];
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'Z1';
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('ZA',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'Z2';
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('ZB',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'Z3';
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('ZC',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2xlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('NPTS')
            ig = ig_(['OX',int2str(I)]);
            pbm.grftype{ig} = 'gL2';
            ig = ig_(['OY',int2str(I)]);
            pbm.grftype{ig} = 'gL2';
            ig = ig_(['E',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['EA',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['EB',int2str(I)]);
            pbm.grelw{ig}(posel) = -1.0;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'QOR2-AY-V-V';
        varargout{1} = pb;
        varargout{2} = pbm;

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eTA'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(2,4);
        U_(1,1) = U_(1,1)+1;
        U_(1,3) = U_(1,3)-1;
        U_(2,2) = U_(2,2)+1;
        U_(2,4) = U_(2,4)-1;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        T = IV_(1)*IV_(1)+IV_(2)*IV_(2);
        varargout{1} = T*T;
        if(nargout>1)
            g_(1,1) = 4.0*T*IV_(1);
            g_(2,1) = 4.0*T*IV_(2);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 4.0*(T+2.0*IV_(1)*IV_(1));
                H_(1,2) = 8.0*IV_(1)*IV_(2);
                H_(2,1) = H_(1,2);
                H_(2,2) = 4.0*(T+2.0*IV_(2)*IV_(2));
                varargout{3} = U_.'*H_*U_;
            end
        end

    case 'eTB'

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
        T = IV_(1)*IV_(1)+IV_(2)*IV_(2);
        ZZSQ = IV_(3)*IV_(3);
        T1 = 1.0+ZZSQ;
        T1SQ = T1*T1;
        varargout{1} = T*T1SQ;
        if(nargout>1)
            g_(1,1) = 2.0*IV_(1)*T1SQ;
            g_(2,1) = 2.0*IV_(2)*T1SQ;
            g_(3,1) = 4.0*T*T1*IV_(3);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,1) = 2.0*T1SQ;
                H_(1,3) = 8.0*IV_(1)*T1*IV_(3);
                H_(3,1) = H_(1,3);
                H_(2,2) = 2.0*T1SQ;
                H_(2,3) = 8.0*IV_(2)*T1*IV_(3);
                H_(3,2) = H_(2,3);
                H_(3,3) = 4.0*T*(2.0*ZZSQ+T1);
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

