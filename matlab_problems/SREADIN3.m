function varargout = SREADIN3(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : SREADIN3
%    *********
% 
%    A nonlinear optimal control problem from Nancy Nichols
%    with a periodic boundary condition and special scaling.
%    This problem arises in tide modelling.
% 
%    SIF input: Nick Gould, July 1991.
% 
%    classification = 'OOR2-MN-V-V'
% 
%    Number of discretized points in [0,1]
% 
%       Alternative values for the SIF file parameters:
% IE N                   2              $-PARAMETER n=6, m=3
% IE N                   5              $-PARAMETER n=12, m=6
% IE N                   50             $-PARAMETER n=102, m=51
% IE N                   100            $-PARAMETER n=202, m=101   original value
% IE N                   500            $-PARAMETER n=1002, m=501
% IE N                   1000           $-PARAMETER n=2002, m=1001
% IE N                   2000           $-PARAMETER n=4002, m=2001
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'SREADIN3';

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
            v_('N') = 10;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
% IE N                   5000           $-PARAMETER n=10002, m=5001
        v_('PI') = 3.1415926535;
        v_('2PI') = 2.0*v_('PI');
        v_('A') = 0.07716;
        v_('1/A') = 1.0/v_('A');
        v_('1/2A') = 0.5*v_('1/A');
        v_('2A') = 2.0*v_('A');
        v_('-2A') = -1.0*v_('2A');
        v_('-1/2A') = 1.0/v_('-2A');
        v_('N-1') = -1+v_('N');
        v_('RN') = v_('N');
        v_('H') = 1.0/v_('RN');
        v_('2/H') = 2.0*v_('RN');
        v_('H/2') = 0.5*v_('H');
        v_('1/H') = 1.0*v_('RN');
        v_('-1/H') = -1.0*v_('RN');
        v_('0') = 0;
        v_('1') = 1;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('0'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
            [iv,ix_] = s2mpjlib('ii',['U',int2str(I)],ix_);
            pb.xnames{iv} = ['U',int2str(I)];
            pb.xscale(iv,1) = v_('RN');
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('N')
            [ig,ig_] = s2mpjlib('ii',['I',int2str(I)],ig_);
            gtype{ig} = '<>';
        end
        for I=v_('1'):v_('N')
            v_('RI') = I;
            v_('TI') = v_('RI')*v_('H');
            v_('2PITI') = v_('2PI')*v_('TI');
            v_('CTI') = cos(v_('2PITI'));
            v_('CCTI') = v_('CTI')*v_('-1/2A');
            v_('I-1') = -1+I;
            v_('RI-1') = v_('I-1');
            v_('TI-1') = v_('RI-1')*v_('H');
            v_('2PITI-1') = v_('2PI')*v_('TI-1');
            v_('CTI-1') = cos(v_('2PITI-1'));
            v_('CCTI-1') = v_('CTI-1')*v_('-1/2A');
            [ig,ig_] = s2mpjlib('ii',['C',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['C',int2str(I)];
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('1/H')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('1/H');
            end
            iv = ix_(['X',int2str(round(v_('I-1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('-1/H')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('-1/H');
            end
            iv = ix_(['U',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('CCTI')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('CCTI');
            end
            iv = ix_(['U',int2str(round(v_('I-1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('CCTI-1')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('CCTI-1');
            end
        end
        [ig,ig_] = s2mpjlib('ii','PERIOD',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'PERIOD';
        iv = ix_(['X',int2str(round(v_('0')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_(['X',int2str(round(v_('N')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        pbm.gscale(ig,1) = v_('H');
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
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        for I=v_('0'):v_('N')
            pb.xlower(ix_(['X',int2str(I)]),1) = -0.5;
            pb.xupper(ix_(['X',int2str(I)])) = 0.5;
        end
        for I=v_('0'):v_('N')
            pb.xlower(ix_(['U',int2str(I)]),1) = 0.0;
            pb.xupper(ix_(['U',int2str(I)])) = 1.0;
        end
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 0.25*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'ePROD',iet_);
        elftv{it}{1} = 'U';
        elftv{it}{2} = 'X';
        elftp{it}{1} = 'P';
        [it,iet_] = s2mpjlib( 'ii', 'eENERGY',iet_);
        elftv{it}{1} = 'U';
        elftv{it}{2} = 'X';
        elftp{it}{1} = 'T';
        elftp{it}{2} = 'HOVER2';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('0'):v_('N')
            v_('RI') = I;
            v_('TI') = v_('RI')*v_('H');
            ename = ['I',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eENERGY';
            ielftype(ie) = iet_('eENERGY');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.25);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['U',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.25);
            posev = find(strcmp('U',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('T',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('TI');
            [~,posep] = ismember('HOVER2',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('H/2');
        end
        for I=v_('0'):v_('N')
            ename = ['NC',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'ePROD';
            ielftype(ie) = iet_('ePROD');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.25);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['U',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.25);
            posev = find(strcmp('U',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('P',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('1/2A');
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('N')
            v_('I-1') = -1+I;
            ig = ig_(['I',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['I',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = -1.0;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['I',int2str(round(v_('I-1')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = -1.0;
        end
        for I=v_('1'):v_('N')
            v_('I-1') = -1+I;
            ig = ig_(['C',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['NC',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['NC',int2str(round(v_('I-1')))]);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'OOR2-MN-V-V';
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb,  pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm, pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end
% ********************
%  SET UP THE GROUPS *
%  ROUTINE           *
% ********************

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'ePROD'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = pbm.elpar{iel_}(1)*EV_(1)*EV_(2);
        if(nargout>1)
            g_(1,1) = pbm.elpar{iel_}(1)*EV_(2);
            g_(2,1) = pbm.elpar{iel_}(1)*EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(2,1) = pbm.elpar{iel_}(1);
                H_(1,2) = H_(2,1);
                varargout{3} = H_;
            end
        end

    case 'eENERGY'

        EV_  = varargin{1};
        iel_ = varargin{2};
        C2PIT = cos(2.0*3.141592653589*pbm.elpar{iel_}(1));
        varargout{1} = pbm.elpar{iel_}(2)*EV_(1)*(EV_(2)-C2PIT)^2;
        if(nargout>1)
            g_(1,1) = pbm.elpar{iel_}(2)*(EV_(2)-C2PIT)^2;
            g_(2,1) = pbm.elpar{iel_}(2)*2.0*EV_(1)*(EV_(2)-C2PIT);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(2,1) = pbm.elpar{iel_}(2)*2.0*(EV_(2)-C2PIT);
                H_(1,2) = H_(2,1);
                H_(2,2) = pbm.elpar{iel_}(2)*2.0*EV_(1);
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

