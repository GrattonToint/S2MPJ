function varargout = MODBEALENE(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : MODBEALENE
%    *********
%    A variation on Beale's problem in 2 variables
%    This is a nonlinear equation variant of MODBEALE
% 
%    Source: An adaptation by Ph. Toint of Problem 5 in
%    J.J. More', B.S. Garbow and K.E. Hillstrom,
%    "Testing Unconstrained Optimization Software",
%    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
% 
%    See also Buckley#89.
%    SIF input: Ph. Toint, Mar 2003.
%               Nick Gould (nonlinear equation version), Jan 2019
% 
%    classification = 'NOR2-AN-V-V'
% 
%    The number of variables is  2 * N/2
% 
%       Alternative values for the SIF file parameters:
% IE N/2                 1              $-PARAMETER     original value
% IE N/2                 2              $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'MODBEALENE';

switch(action)

    case 'setup'

    pb.name      = 'MODBEALENE';
    pb.sifpbname = 'MODBEALENE';
    pbm.name     = 'MODBEALENE';
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        if(nargin<2)
            v_('N/2') = 5;  %  SIF file default value
        else
            v_('N/2') = varargin{1};
        end
% IE N/2                 100            $-PARAMETER
% IE N/2                 1000           $-PARAMETER
% IE N/2                 10000          $-PARAMETER
        if(nargin<3)
            v_('ALPHA') = 50.0;  %  SIF file default value
        else
            v_('ALPHA') = varargin{2};
        end
        v_('1') = 1;
        v_('N') = v_('N/2')+v_('N/2');
        v_('N/2-1') = -1+v_('N/2');
        v_('ALPHINV') = 1.0/v_('ALPHA');
        v_('RALPHINV') = sqrt(v_('ALPHINV'));
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for J=v_('1'):v_('N')
            [iv,ix_] = s2xlib('ii',['X',int2str(J)],ix_);
            pb.xnames{iv} = ['X',int2str(J)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('N/2-1')
            v_('I-1') = -1+I;
            v_('2I-1') = v_('I-1')+v_('I-1');
            v_('J') = 1+v_('2I-1');
            v_('J+1') = 1+v_('J');
            v_('J+2') = 2+v_('J');
            [ig,ig_] = s2xlib('ii',['BA',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['BA',int2str(I)];
            [ig,ig_] = s2xlib('ii',['BB',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['BB',int2str(I)];
            [ig,ig_] = s2xlib('ii',['BC',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['BC',int2str(I)];
            [ig,ig_] = s2xlib('ii',['L',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['L',int2str(I)];
            iv = ix_(['X',int2str(round(v_('J+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 6.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 6.0;
            end
            iv = ix_(['X',int2str(round(v_('J+2')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            pbm.gscale(ig,1) = v_('RALPHINV');
        end
        [ig,ig_] = s2xlib('ii',['BA',int2str(round(v_('N/2')))],ig_);
        gtype{ig}  = '==';
        cnames{ig} = ['BA',int2str(round(v_('N/2')))];
        [ig,ig_] = s2xlib('ii',['BB',int2str(round(v_('N/2')))],ig_);
        gtype{ig}  = '==';
        cnames{ig} = ['BB',int2str(round(v_('N/2')))];
        [ig,ig_] = s2xlib('ii',['BC',int2str(round(v_('N/2')))],ig_);
        gtype{ig}  = '==';
        cnames{ig} = ['BC',int2str(round(v_('N/2')))];
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
        for I=v_('1'):v_('N/2')
            pbm.gconst(ig_(['BA',int2str(I)])) = 1.5;
            pbm.gconst(ig_(['BB',int2str(I)])) = 2.25;
            pbm.gconst(ig_(['BC',int2str(I)])) = 2.625;
        end
        pb.xlower = zeros(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 1.0*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2xlib( 'ii', 'ePRODB',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftp{it}{1} = 'POW';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('1'):v_('N/2')
            v_('I-1') = -1+I;
            v_('2I-1') = v_('I-1')+v_('I-1');
            v_('J') = 1+v_('2I-1');
            v_('J+1') = 1+v_('J');
            ename = ['AE',int2str(I)];
            [ie,ie_,newelt] = s2xlib('ii',ename,ie_);
            if(newelt)
                pbm.elftype{ie} = 'ePRODB';
                ielftype(ie) = iet_('ePRODB');
            end
            vname = ['X',int2str(round(v_('J')))];
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('J+1')))];
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('POW',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = 1.0;
            ename = ['BE',int2str(I)];
            [ie,ie_,newelt] = s2xlib('ii',ename,ie_);
            if(newelt)
                pbm.elftype{ie} = 'ePRODB';
                ielftype(ie) = iet_('ePRODB');
            end
            vname = ['X',int2str(round(v_('J')))];
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('J+1')))];
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('POW',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = 2.0;
            ename = ['CE',int2str(I)];
            [ie,ie_,newelt] = s2xlib('ii',ename,ie_);
            if(newelt)
                pbm.elftype{ie} = 'ePRODB';
                ielftype(ie) = iet_('ePRODB');
            end
            vname = ['X',int2str(round(v_('J')))];
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('J+1')))];
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('POW',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = 3.0;
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('N/2')
            ig = ig_(['BA',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['AE',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            ig = ig_(['BB',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['BE',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            ig = ig_(['BC',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['CE',int2str(I)]);
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
        pb.pbclass = 'NOR2-AN-V-V';
        varargout{1} = pb;
        varargout{2} = pbm;

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'ePRODB'

        EV_  = varargin{1};
        iel_ = varargin{2};
        T = 1.0-EV_(2)^pbm.elpar{iel_}(1);
        POWM1 = pbm.elpar{iel_}(1)-1.0;
        W = -pbm.elpar{iel_}(1)*EV_(2)^POWM1;
        varargout{1} = EV_(1)*T;
        if(nargout>1)
            g_(1,1) = T;
            g_(2,1) = EV_(1)*W;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 0.0;
                H_(1,2) = W;
                H_(2,1) = H_(1,2);
                H_(2,2) = -EV_(1)*pbm.elpar{iel_}(1)*POWM1*EV_(2)^(pbm.elpar{iel_}(1)-2.0);
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
