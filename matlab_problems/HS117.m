function varargout = HS117(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HS117
%    *********
% 
%    Source: problem 117 in
%    W. Hock and K. Schittkowski,
%    "Test examples for nonlinear programming codes",
%    Lectures Notes in Economics and Mathematical Systems 187, Springer
%    Verlag, Heidelberg, 1981.
% 
%    SIF input: Nick Gould, August 1991.
% 
%    classification = 'OQR2-AN-15-5'
% 
%    Number of constraints
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HS117';

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
        v_('N') = 5;
        v_('M') = 10;
        v_('M+N') = v_('M')+v_('N');
        v_('1') = 1;
        v_('E1') = -15.0;
        v_('E2') = -27.0;
        v_('E3') = -36.0;
        v_('E4') = -18.0;
        v_('E5') = -12.0;
        v_('C1,1') = 30.0;
        v_('C2,1') = -20.0;
        v_('C3,1') = -10.0;
        v_('C4,1') = 32.0;
        v_('C5,1') = -10.0;
        v_('C1,2') = -20.0;
        v_('C2,2') = 39.0;
        v_('C3,2') = -6.0;
        v_('C4,2') = -31.0;
        v_('C5,2') = 32.0;
        v_('C1,3') = -10.0;
        v_('C2,3') = -6.0;
        v_('C3,3') = 10.0;
        v_('C4,3') = -6.0;
        v_('C5,3') = -10.0;
        v_('C1,4') = 32.0;
        v_('C2,4') = -31.0;
        v_('C3,4') = -6.0;
        v_('C4,4') = 39.0;
        v_('C5,4') = -20.0;
        v_('C1,5') = -10.0;
        v_('C2,5') = 32.0;
        v_('C3,5') = -10.0;
        v_('C4,5') = -20.0;
        v_('C5,5') = 30.0;
        v_('D1') = 4.0;
        v_('D2') = 8.0;
        v_('D3') = 10.0;
        v_('D4') = 6.0;
        v_('D5') = 2.0;
        v_('A1,1') = -16.0;
        v_('A2,1') = 0.0;
        v_('A3,1') = -3.5;
        v_('A4,1') = 0.0;
        v_('A5,1') = 0.0;
        v_('A6,1') = 2.0;
        v_('A7,1') = -1.0;
        v_('A8,1') = -1.0;
        v_('A9,1') = 1.0;
        v_('A10,1') = 1.0;
        v_('A1,2') = 2.0;
        v_('A2,2') = -2.0;
        v_('A3,2') = 0.0;
        v_('A4,2') = -2.0;
        v_('A5,2') = -9.0;
        v_('A6,2') = 0.0;
        v_('A7,2') = -1.0;
        v_('A8,2') = -2.0;
        v_('A9,2') = 2.0;
        v_('A10,2') = 1.0;
        v_('A1,3') = 0.0;
        v_('A2,3') = 0.0;
        v_('A3,3') = 2.0;
        v_('A4,3') = 0.0;
        v_('A5,3') = -2.0;
        v_('A6,3') = -4.0;
        v_('A7,3') = -1.0;
        v_('A8,3') = -3.0;
        v_('A9,3') = 3.0;
        v_('A10,3') = 1.0;
        v_('A1,4') = 1.0;
        v_('A2,4') = 4.0;
        v_('A3,4') = 0.0;
        v_('A4,4') = -4.0;
        v_('A5,4') = 1.0;
        v_('A6,4') = 0.0;
        v_('A7,4') = -1.0;
        v_('A8,4') = -2.0;
        v_('A9,4') = 4.0;
        v_('A10,4') = 1.0;
        v_('A1,5') = 0.0;
        v_('A2,5') = 2.0;
        v_('A3,5') = 0.0;
        v_('A4,5') = -1.0;
        v_('A5,5') = -2.8;
        v_('A6,5') = 0.0;
        v_('A7,5') = -1.0;
        v_('A8,5') = -1.0;
        v_('A9,5') = 5.0;
        v_('A10,5') = 1.0;
        v_('B1') = -40.0;
        v_('B2') = -2.0;
        v_('B3') = -0.25;
        v_('B4') = -4.0;
        v_('B5') = -4.0;
        v_('B6') = -1.0;
        v_('B7') = -40.0;
        v_('B8') = -60.0;
        v_('B9') = 5.0;
        v_('B10') = 1.0;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('M+N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for J=v_('1'):v_('M')
            v_('-BJ') = -1.0*v_(['B',int2str(J)]);
            [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
            gtype{ig} = '<>';
            iv = ix_(['X',int2str(J)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('-BJ')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('-BJ');
            end
        end
        for J=v_('1'):v_('N')
            for K=v_('1'):v_('N')
                v_('M+K') = v_('M')+K;
                v_('2CKJ') = 2.0*v_(['C',int2str(K),',',int2str(J)]);
                [ig,ig_] = s2mpjlib('ii',['C',int2str(J)],ig_);
                gtype{ig}  = '>=';
                cnames{ig} = ['C',int2str(J)];
                iv = ix_(['X',int2str(round(v_('M+K')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('2CKJ')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('2CKJ');
                end
            end
            for K=v_('1'):v_('M')
                v_('-AKJ') = -1.0*v_(['A',int2str(K),',',int2str(J)]);
                [ig,ig_] = s2mpjlib('ii',['C',int2str(J)],ig_);
                gtype{ig}  = '>=';
                cnames{ig} = ['C',int2str(J)];
                iv = ix_(['X',int2str(K)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('-AKJ')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('-AKJ');
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
        pbm.congrps = [ legrps, eqgrps, gegrps ];
        [pb.cnames{1:pb.m}] = deal(cnames{pbm.congrps});
        pb.nob = ngrp-pb.m;
        pbm.objgrps = find(strcmp(gtype,'<>'));
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        for J=v_('1'):v_('N')
            v_('-EJ') = -1.0*v_(['E',int2str(J)]);
            pbm.gconst(ig_(['C',int2str(J)])) = v_('-EJ');
        end
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 0.001*ones(pb.n,1);
        pb.y0 = 0.001*ones(pb.m,1);
        if(isKey(ix_,'X7'))
            pb.x0(ix_('X7'),1) = 60.0;
        else
            pb.y0(find(pbm.congrps==ig_('X7')),1) = 60.0;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSQUARE',iet_);
        elftv{it}{1} = 'XJ';
        [it,iet_] = s2mpjlib( 'ii', 'eCUBE',iet_);
        elftv{it}{1} = 'XJ';
        [it,iet_] = s2mpjlib( 'ii', 'ePROD',iet_);
        elftv{it}{1} = 'XI';
        elftv{it}{2} = 'XJ';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for J=v_('1'):v_('N')
            v_('M+J') = v_('M')+J;
            ename = ['2D',int2str(J)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eCUBE';
            ielftype(ie) = iet_('eCUBE');
            vname = ['X',int2str(round(v_('M+J')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.001);
            posev = find(strcmp('XJ',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['3D',int2str(J)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQUARE';
            ielftype(ie) = iet_('eSQUARE');
            vname = ['X',int2str(round(v_('M+J')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.001);
            posev = find(strcmp('XJ',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            for K=v_('1'):v_('N')
                v_('M+K') = v_('M')+K;
                ename = ['C',int2str(K),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'ePROD';
                ielftype(ie) = iet_('ePROD');
                vname = ['X',int2str(round(v_('M+K')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.001);
                posev = find(strcmp('XI',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(round(v_('M+J')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.001);
                posev = find(strcmp('XJ',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for J=v_('1'):v_('N')
            v_('2DJ') = 2.0*v_(['D',int2str(J)]);
            ig = ig_('OBJ');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['2D',int2str(J)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('2DJ');
            v_('3DJ') = 3.0*v_(['D',int2str(J)]);
            ig = ig_(['C',int2str(J)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['3D',int2str(J)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('3DJ');
            for K=v_('1'):v_('N')
                ig = ig_('OBJ');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['C',int2str(K),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_(['C',int2str(K),',',int2str(J)]);
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               32.34867897
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'OQR2-AN-15-5';
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb,  pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm, pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end
% **********************
%  SET UP THE FUNCTION *
%  AND RANGE ROUTINES  *
% **********************

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eSQUARE'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(1);
        if(nargout>1)
            g_(1,1) = 2.0e+0*EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0e+0;
                varargout{3} = H_;
            end
        end

    case 'eCUBE'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(1)*EV_(1);
        if(nargout>1)
            g_(1,1) = 3.0e+0*EV_(1)*EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 6.0e+0*EV_(1);
                varargout{3} = H_;
            end
        end

    case 'ePROD'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(2);
        if(nargout>1)
            g_(1,1) = EV_(2);
            g_(2,1) = EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = 1.0e+0;
                H_(2,1) = H_(1,2);
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

