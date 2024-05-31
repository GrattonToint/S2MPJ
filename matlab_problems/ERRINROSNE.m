function varargout = ERRINROSNE(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : ERRINROSNE
%    --------
% 
%    A nonlinear function similar to the chained Rosenbrock
%    problem CHNROSNB.
%    This is a nonlinear equation variant of ERRINROS
% 
%    Source:
%    An error in specifying problem CHNROSNB.
%    SIF input: Ph. Toint, Sept 1990.
%               Nick Gould (nonlinear equation version), Jan 2019
% 
%    classification = 'NOR2-AN-V-V'
% 
%    Number of variables (at most 50)
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'ERRINROSNE';

switch(action)

    case 'setup'

    pb.name      = 'ERRINROSNE';
    pb.sifpbname = 'ERRINROSNE';
    pbm.name     = 'ERRINROSNE';
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        if(nargin<2)
            v_('N') = 10;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
        v_('ALPH1') = 1.25;
        v_('ALPH2') = 1.40;
        v_('ALPH3') = 2.40;
        v_('ALPH4') = 1.40;
        v_('ALPH5') = 1.75;
        v_('ALPH6') = 1.20;
        v_('ALPH7') = 2.25;
        v_('ALPH8') = 1.20;
        v_('ALPH9') = 1.00;
        v_('ALPH10') = 1.10;
        v_('ALPH11') = 1.50;
        v_('ALPH12') = 1.60;
        v_('ALPH13') = 1.25;
        v_('ALPH14') = 1.25;
        v_('ALPH15') = 1.20;
        v_('ALPH16') = 1.20;
        v_('ALPH17') = 1.40;
        v_('ALPH18') = 0.50;
        v_('ALPH19') = 0.50;
        v_('ALPH20') = 1.25;
        v_('ALPH21') = 1.80;
        v_('ALPH22') = 0.75;
        v_('ALPH23') = 1.25;
        v_('ALPH24') = 1.40;
        v_('ALPH25') = 1.60;
        v_('ALPH26') = 2.00;
        v_('ALPH27') = 1.00;
        v_('ALPH28') = 1.60;
        v_('ALPH29') = 1.25;
        v_('ALPH30') = 2.75;
        v_('ALPH31') = 1.25;
        v_('ALPH32') = 1.25;
        v_('ALPH33') = 1.25;
        v_('ALPH34') = 3.00;
        v_('ALPH35') = 1.50;
        v_('ALPH36') = 2.00;
        v_('ALPH37') = 1.25;
        v_('ALPH38') = 1.40;
        v_('ALPH39') = 1.80;
        v_('ALPH40') = 1.50;
        v_('ALPH41') = 2.20;
        v_('ALPH42') = 1.40;
        v_('ALPH43') = 1.50;
        v_('ALPH44') = 1.25;
        v_('ALPH45') = 2.00;
        v_('ALPH46') = 1.50;
        v_('ALPH47') = 1.25;
        v_('ALPH48') = 1.40;
        v_('ALPH49') = 0.60;
        v_('ALPH50') = 1.50;
        v_('1') = 1;
        v_('2') = 2;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2xlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('2'):v_('N')
            v_('I-1') = -1+I;
            [ig,ig_] = s2xlib('ii',['SQ',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['SQ',int2str(I)];
            iv = ix_(['X',int2str(round(v_('I-1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            [ig,ig_] = s2xlib('ii',['B',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['B',int2str(I)];
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
        pbm.congrps = find(ismember(gtype,{'<=','==','>='}));
        [pb.cnames{1:pb.m}] = deal(cnames{pbm.congrps});
        pb.nob = ngrp-pb.m;
        pbm.objgrps = find(strcmp(gtype,'<>'));
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        for I=v_('2'):v_('N')
            pbm.gconst(ig_(['B',int2str(I)])) = 1.0;
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        for I=v_('1'):v_('N')
            if(isKey(ix_,['X',int2str(I)]))
                pb.x0(ix_(['X',int2str(I)]),1) = -1.0;
            else
                pb.y0(find(pbm.congrps==ig_(['X',int2str(I)])),1) = -1.0;
            end
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2xlib( 'ii', 'eETYPE',iet_);
        elftv{it}{1} = 'V1';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('2'):v_('N')
            ename = ['ELA',int2str(I)];
            [ie,ie_,newelt] = s2xlib('ii',ename,ie_);
            if(newelt)
                pbm.elftype{ie} = 'eETYPE';
                ielftype(ie) = iet_('eETYPE');
            end
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('2'):v_('N')
            v_('AI2') = v_(['ALPH',int2str(I)])*v_(['ALPH',int2str(I)]);
            v_('AI') = 16.0*v_('AI2');
            ig = ig_(['SQ',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['ELA',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('AI');
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

    case 'eETYPE'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = -EV_(1)^2;
        if(nargout>1)
            g_(1,1) = -2.0*EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = -2.0;
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

