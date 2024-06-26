function varargout = INTEGREQ(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : INTEGREQ
%    *********
%    The discrete integral problem
% 
%    Source:  Problem 29 in
%    J.J. More', B.S. Garbow and K.E. Hillstrom,
%    "Testing Unconstrained Optimization Software",
%    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
% 
%    SIF input: Ph. Toint, Feb 1990.
% 
%    classification = 'NOR2-AN-V-V'
% 
%    N+2 is the number of discretization points .
%    The number of free variables is N.
% 
%       Alternative values for the SIF file parameters:
% IE N                   10             $-PARAMETER
% IE N                   50             $-PARAMETER     original value
% IE N                   100            $-PARAMETER
% IE N                   500            $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'INTEGREQ';

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
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        if(nargs<1)
            v_('N') = 50;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
        v_('0') = 0;
        v_('1') = 1;
        v_('N+1') = 1+v_('N');
        v_('RN+1') = v_('N+1');
        v_('H') = 1.0/v_('RN+1');
        v_('HALFH') = 0.5e0*v_('H');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('0'):v_('N+1')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('N')
            [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['G',int2str(I)];
            iv = ix_(['X',int2str(I)]);
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
        pbm.congrps = [ legrps, eqgrps, gegrps ];
        [pb.cnames{1:pb.m}] = deal(cnames{pbm.congrps});
        pb.nob = ngrp-pb.m;
        pbm.objgrps = find(strcmp(gtype,'<>'));
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        pb.xlower(ix_(['X',int2str(round(v_('0')))]),1) = 0.0;
        pb.xupper(ix_(['X',int2str(round(v_('0')))]),1) = 0.0;
        pb.xlower(ix_(['X',int2str(round(v_('N+1')))]),1) = 0.0;
        pb.xupper(ix_(['X',int2str(round(v_('N+1')))]),1) = 0.0;
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,['X',int2str(round(v_('0')))]))
            pb.x0(ix_(['X',int2str(round(v_('0')))]),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_(['X',int2str(round(v_('0')))])),1) = 0.0;
        end
        for I=v_('1'):v_('N')
            v_('REALI') = I;
            v_('IH') = v_('REALI')*v_('H');
            v_('IH-1') = -1.0+v_('IH');
            v_('TI') = v_('IH')*v_('IH-1');
            if(isKey(ix_,['X',int2str(I)]))
                pb.x0(ix_(['X',int2str(I)]),1) = v_('TI');
            else
                pb.y0(find(pbm.congrps==ig_(['X',int2str(I)])),1) = v_('TI');
            end
        end
        if(isKey(ix_,['X',int2str(round(v_('N+1')))]))
            pb.x0(ix_(['X',int2str(round(v_('N+1')))]),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_(['X',int2str(round(v_('N+1')))])),1) = 0.0;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eVBCUBE',iet_);
        elftv{it}{1} = 'V';
        elftp{it}{1} = 'B';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for J=v_('1'):v_('N')
            v_('REALJ') = J;
            v_('TJ') = v_('REALJ')*v_('H');
            v_('1+TJ') = 1.0+v_('TJ');
            ename = ['A',int2str(J)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eVBCUBE';
            ielftype(ie) = iet_('eVBCUBE');
            vname = ['X',int2str(J)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('B',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('1+TJ');
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('N')
            v_('REALI') = I;
            v_('TI') = v_('REALI')*v_('H');
            v_('-TI') = -1.0*v_('TI');
            v_('1-TI') = 1.0+v_('-TI');
            v_('P1') = v_('1-TI')*v_('HALFH');
            for J=v_('1'):I
                v_('REALJ') = J;
                v_('TJ') = v_('REALJ')*v_('H');
                v_('WIL') = v_('P1')*v_('TJ');
                ig = ig_(['G',int2str(I)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['A',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('WIL');
            end
            v_('I+1') = 1+I;
            v_('P2') = v_('TI')*v_('HALFH');
            for J=v_('I+1'):v_('N')
                v_('REALJ') = J;
                v_('TJ') = v_('REALJ')*v_('H');
                v_('-TJ') = -1.0*v_('TJ');
                v_('1-TJ') = 1.0+v_('-TJ');
                v_('WIU') = v_('P2')*v_('1-TJ');
                ig = ig_(['G',int2str(I)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['A',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('WIU');
            end
        end
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'NOR2-AN-V-V';
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb,  pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm, pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eVBCUBE'

        EV_  = varargin{1};
        iel_ = varargin{2};
        VPLUSB = EV_(1)+pbm.elpar{iel_}(1);
        varargout{1} = VPLUSB^3;
        if(nargout>1)
            g_(1,1) = 3.0*VPLUSB^2;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 6.0*VPLUSB;
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

