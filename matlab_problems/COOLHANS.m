function varargout = COOLHANS(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : COOLHANS
%    *********
% 
%    A problem arising from the analysis of a Cooley-Hansen economy with
%    loglinear approximation.  The problem is to solve the matrix equation
%                  A * X * X + B * X + C = 0
%    where A, B and C are known N times N matrices and X an unknown matrix
%    of matching dimension.  The instance considered here has N = 3.
% 
%    Source:
%    S. Ceria, private communication, 1995.
% 
%    SIF input: Ph. Toint, Feb 1995.
% 
%    classification = 'NQR2-RN-9-9'
% 
%    order of the matrix equation
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'COOLHANS';

switch(action)

    case 'setup'

    pb.name      = 'COOLHANS';
    pb.sifpbname = 'COOLHANS';
    pbm.name     = 'COOLHANS';
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        v_('N') = 3;
        v_('A1,1') = 0.0;
        v_('A2,1') = 0.13725e-6;
        v_('A3,1') = 0.0;
        v_('A1,2') = 0.0;
        v_('A2,2') = 937.62;
        v_('A3,2') = 0.0;
        v_('A1,3') = 0.0;
        v_('A2,3') = -42.207;
        v_('A3,3') = 0.0;
        v_('B1,1') = 0.0060893;
        v_('B2,1') = 0.13880e-6;
        v_('B3,1') = -0.13877e-6;
        v_('B1,2') = -44.292;
        v_('B2,2') = -1886.0;
        v_('B3,2') = 42.362;
        v_('B1,3') = 2.0011;
        v_('B2,3') = 42.362;
        v_('B3,3') = -2.0705;
        v_('C1,1') = 0.0;
        v_('C2,1') = 0.0;
        v_('C3,1') = 0.0;
        v_('C1,2') = 44.792;
        v_('C2,2') = 948.21;
        v_('C3,2') = -42.684;
        v_('C1,3') = 0.0;
        v_('C2,3') = 0.0;
        v_('C3,3') = 0.0;
        v_('1') = 1;
        v_('2') = 2;
        v_('3') = 3;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            for J=v_('1'):v_('N')
                [iv,ix_] = s2xlib('ii',['X',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['X',int2str(I),',',int2str(J)];
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for K=v_('1'):v_('N')
            for L=v_('1'):v_('N')
                for M=v_('1'):v_('N')
                    [ig,ig_] = s2xlib('ii',['G',int2str(K),',',int2str(L)],ig_);
                    gtype{ig}  = '==';
                    cnames{ig} = ['G',int2str(K),',',int2str(L)];
                    iv = ix_(['X',int2str(M),',',int2str(L)]);
                    if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                        pbm.A(ig,iv) = v_(['B',int2str(K),',',int2str(M)])+pbm.A(ig,iv);
                    else
                        pbm.A(ig,iv) = v_(['B',int2str(K),',',int2str(M)]);
                    end
                end
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
        for K=v_('1'):v_('N')
            for L=v_('1'):v_('N')
                v_('-C') = -1.0*v_(['C',int2str(K),',',int2str(L)]);
                pbm.gconst(ig_(['G',int2str(K),',',int2str(L)])) = v_('-C');
            end
        end
        pb.xlower = zeros(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2xlib( 'ii', 'en2PR',iet_);
        elftv{it}{1} = 'XX';
        elftv{it}{2} = 'YY';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for K=v_('1'):v_('N')
            for L=v_('1'):v_('N')
                for M=v_('1'):v_('N')
                    ename = ['E',int2str(K),',',int2str(M),',',int2str(L)];
                    [ie,ie_] = s2xlib('ii',ename,ie_);
                    pbm.elftype{ie} = 'en2PR';
                    ielftype(ie) = iet_('en2PR');
                    vname = ['X',int2str(K),',',int2str(M)];
                    [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('XX',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                    vname = ['X',int2str(M),',',int2str(L)];
                    [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('YY',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for L=v_('1'):v_('N')
            for P=v_('1'):v_('N')
                for M=v_('1'):v_('N')
                    ig = ig_(['G',int2str(round(v_('1'))),',',int2str(L)]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) = ie_(['E',int2str(P),',',int2str(M),',',int2str(L)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_(['A',int2str(round(v_('1'))),',',int2str(P)]);
                end
            end
        end
        for L=v_('1'):v_('N')
            for P=v_('1'):v_('N')
                for M=v_('1'):v_('N')
                    ig = ig_(['G',int2str(round(v_('2'))),',',int2str(L)]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) = ie_(['E',int2str(P),',',int2str(M),',',int2str(L)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_(['A',int2str(round(v_('2'))),',',int2str(P)]);
                end
            end
        end
        for L=v_('1'):v_('N')
            for P=v_('1'):v_('N')
                for M=v_('1'):v_('N')
                    ig = ig_(['G',int2str(round(v_('3'))),',',int2str(L)]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) = ie_(['E',int2str(P),',',int2str(M),',',int2str(L)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_(['A',int2str(round(v_('3'))),',',int2str(P)]);
                end
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'NQR2-RN-9-9';
        pb.x0          = zeros(pb.n,1);
        varargout{1} = pb;
        varargout{2} = pbm;

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'en2PR'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(2);
        if(nargout>1)
            g_(1,1) = EV_(2);
            g_(2,1) = EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = 1.0;
                H_(2,1) = H_(1,2);
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
