function varargout = OBSTCLBL(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : OBSTCLBL
%    *********
% 
%    A quadratic obstacle problem by Dembo and Tulowitzki
% 
%    The problem comes from the obstacle problem on a rectangle.
%    The rectangle is discretized into (px-1)(py-1) little rectangles. The
%    heights of the considered surface above the corners of these little
%    rectangles are the problem variables,  There are px*py of them.
% 
%    Source:
%    R. Dembo and U. Tulowitzki,
%    "On the minimization of quadratic functions subject to box
%    constraints",
%    WP 71, Yale University (new Haven, USA), 1983.
% 
%    See also More 1989 (Problem B, Starting point L)
% 
%    SIF input: Ph. Toint, Dec 1989.
% 
%    classification = 'QBR2-AY-V-0'
% 
%    PX is the number of points along the X side of the rectangle
%    PY is the number of points along the Y side of the rectangle
% 
%       Alternative values for the SIF file parameters:
% IE PX                  4              $-PARAMETER n = 16
% IE PY                  4              $-PARAMETER
% 
% IE PX                  10             $-PARAMETER n = 100     original value
% IE PY                  10             $-PARAMETER             original value
% 
% IE PX                  23             $-PARAMETER n = 529
% IE PY                  23             $-PARAMETER
% 
% IE PX                  32             $-PARAMETER n = 1024
% IE PY                  32             $-PARAMETER
% 
% IE PX                  75             $-PARAMETER n = 5625
% IE PY                  75             $-PARAMETER
% 
% IE PX                  100            $-PARAMETER n = 10000
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'OBSTCLBL';

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
            v_('PX') = 5;  %  SIF file default value
        else
            v_('PX') = varargin{1};
        end
% IE PY                  100            $-PARAMETER
        if(nargs<2)
            v_('PY') = 20;  %  SIF file default value
        else
            v_('PY') = varargin{2};
        end
% IE PX                  125            $-PARAMETER n = 15625
% IE PY                  125            $-PARAMETER
        if(nargs<3)
            v_('C') = 1.0;  %  SIF file default value
        else
            v_('C') = varargin{3};
        end
        v_('PX-1') = -1+v_('PX');
        v_('RPX-1') = v_('PX-1');
        v_('HX') = 1.0/v_('RPX-1');
        v_('1/HX') = 1.0/v_('HX');
        v_('PY-1') = -1+v_('PY');
        v_('RPY-1') = v_('PY-1');
        v_('HY') = 1.0/v_('RPY-1');
        v_('1/HY') = 1.0/v_('HY');
        v_('HXHY') = v_('HX')*v_('HY');
        v_('HX/HY') = v_('HX')*v_('1/HY');
        v_('HY/HX') = v_('HY')*v_('1/HX');
        v_('HY/4HX') = 0.25*v_('HY/HX');
        v_('HX/4HY') = 0.25*v_('HX/HY');
        v_('C0') = v_('HXHY')*v_('C');
        v_('LC') = -1.0*v_('C0');
        v_('1') = 1;
        v_('2') = 2;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for J=v_('1'):v_('PX')
            for I=v_('1'):v_('PY')
                [iv,ix_] = s2mpjlib('ii',['X',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['X',int2str(I),',',int2str(J)];
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('2'):v_('PY-1')
            for J=v_('2'):v_('PX-1')
                [ig,ig_] = s2mpjlib('ii',['G',int2str(I),',',int2str(J)],ig_);
                gtype{ig} = '<>';
                iv = ix_(['X',int2str(I),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('LC')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('LC');
                end
            end
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        for J=v_('1'):v_('PX')
            pb.xlower(ix_(['X',int2str(round(v_('1'))),',',int2str(J)]),1) = 0.0;
            pb.xupper(ix_(['X',int2str(round(v_('1'))),',',int2str(J)]),1) = 0.0;
            pb.xlower(ix_(['X',int2str(round(v_('PY'))),',',int2str(J)]),1) = 0.0;
            pb.xupper(ix_(['X',int2str(round(v_('PY'))),',',int2str(J)]),1) = 0.0;
        end
        for I=v_('2'):v_('PY-1')
            pb.xlower(ix_(['X',int2str(I),',',int2str(round(v_('PX')))]),1) = 0.0;
            pb.xupper(ix_(['X',int2str(I),',',int2str(round(v_('PX')))]),1) = 0.0;
            pb.xlower(ix_(['X',int2str(I),',',int2str(round(v_('1')))]),1) = 0.0;
            pb.xupper(ix_(['X',int2str(I),',',int2str(round(v_('1')))]),1) = 0.0;
        end
        for I=v_('2'):v_('PY-1')
            v_('I-1') = -1+I;
            v_('RI-1') = v_('I-1');
            v_('XI1') = v_('RI-1')*v_('HY');
            v_('3XI1') = 9.2*v_('XI1');
            v_('SXI1') = sin(v_('3XI1'));
            for J=v_('2'):v_('PX-1')
                v_('J-1') = -1+J;
                v_('RJ-1') = v_('J-1');
                v_('XI2') = v_('RJ-1')*v_('HX');
                v_('3XI2') = 9.3*v_('XI2');
                v_('SXI2') = sin(v_('3XI2'));
                v_('L1') = v_('SXI1')*v_('SXI2');
                v_('L2') = v_('L1')*v_('L1');
                v_('LOW') = v_('L2')*v_('L1');
                pb.xlower(ix_(['X',int2str(I),',',int2str(J)]),1) = v_('LOW');
            end
        end
        for I=v_('2'):v_('PY-1')
            v_('I-1') = -1+I;
            v_('RI-1') = v_('I-1');
            v_('XI1') = v_('RI-1')*v_('HY');
            v_('3XI1') = 9.2*v_('XI1');
            v_('SXI1') = sin(v_('3XI1'));
            for J=v_('2'):v_('PX-1')
                v_('J-1') = -1+J;
                v_('RJ-1') = v_('J-1');
                v_('XI2') = v_('RJ-1')*v_('HX');
                v_('3XI2') = 9.3*v_('XI2');
                v_('SXI2') = sin(v_('3XI2'));
                v_('L1') = v_('SXI1')*v_('SXI2');
                v_('L2') = v_('L1')*v_('L1');
                v_('UPP') = 0.02+v_('L2');
                pb.xupper(ix_(['X',int2str(I),',',int2str(J)])) = v_('UPP');
            end
        end
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        for J=v_('1'):v_('PX')
            pb.x0(ix_(['X',int2str(round(v_('1'))),',',int2str(J)]),1) = 0.0;
            pb.x0(ix_(['X',int2str(round(v_('PY'))),',',int2str(J)]),1) = 0.0;
        end
        for I=v_('2'):v_('PY-1')
            pb.x0(ix_(['X',int2str(I),',',int2str(round(v_('PX')))]),1) = 0.0;
            pb.x0(ix_(['X',int2str(I),',',int2str(round(v_('1')))]),1) = 0.0;
        end
        for I=v_('2'):v_('PY-1')
            v_('I-1') = -1+I;
            v_('RI-1') = v_('I-1');
            v_('XI1') = v_('RI-1')*v_('HY');
            v_('3XI1') = 9.2*v_('XI1');
            v_('SXI1') = sin(v_('3XI1'));
            for J=v_('2'):v_('PX-1')
                v_('J-1') = -1+J;
                v_('RJ-1') = v_('J-1');
                v_('XI2') = v_('RJ-1')*v_('HX');
                v_('3XI2') = 9.3*v_('XI2');
                v_('SXI2') = sin(v_('3XI2'));
                v_('L1') = v_('SXI1')*v_('SXI2');
                v_('L2') = v_('L1')*v_('L1');
                v_('LOW') = v_('L2')*v_('L1');
                pb.x0(ix_(['X',int2str(I),',',int2str(J)]),1) = v_('LOW');
            end
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eISQ',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('2'):v_('PY-1')
            v_('I-1') = -1+I;
            v_('I+1') = 1+I;
            for J=v_('2'):v_('PX-1')
                v_('J-1') = -1+J;
                v_('J+1') = 1+J;
                ename = ['A',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eISQ';
                ielftype(ie) = iet_('eISQ');
                vname = ['X',int2str(round(v_('I+1'))),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['B',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eISQ';
                ielftype(ie) = iet_('eISQ');
                vname = ['X',int2str(I),',',int2str(round(v_('J+1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['C',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eISQ';
                ielftype(ie) = iet_('eISQ');
                vname = ['X',int2str(round(v_('I-1'))),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['D',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eISQ';
                ielftype(ie) = iet_('eISQ');
                vname = ['X',int2str(I),',',int2str(round(v_('J-1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('2'):v_('PY-1')
            for J=v_('2'):v_('PX-1')
                ig = ig_(['G',int2str(I),',',int2str(J)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['A',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = v_('HY/4HX');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['B',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = v_('HX/4HY');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['C',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = v_('HY/4HX');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['D',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = v_('HX/4HY');
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN(4)            -0.0081108
% LO SOLTN(10)           2.87503823
% LO SOLTN(23)           6.51932527
% LO SOLTN(32)           6.88708670
% LO SOLTN(75)           ???
% LO SOLTN(100)          ???
% LO SOLTN(125)          ???
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'QBR2-AY-V-0';
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

