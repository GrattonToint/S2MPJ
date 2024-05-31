function varargout = POWELLBC(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : POWELLBC
%    --------
% 
%    A bound-constrained optimization problem to 
%    separate points within a square in the plane
% 
%    Given points p_j = ( x_2j-1 , x_2j ), i = 1, ..., p
% 
%    minimize sum_k=2^p sum_j=1^k-1 1 / || p_j - p_k ||_2
% 
%    subject to 0 <= x_i <= 1, i = 1, ..., 2p = n
% 
%    Source: 
%    M. J. D. Powell
%    Private communication (Optbridge, 2006)
% 
%    SIF input: Nick Gould, Aug 2006.
% 
%    classification = 'OBR2-AN-V-0'
% 
%    Number of points
% 
%       Alternative values for the SIF file parameters:
% IE P                   2              $-PARAMETER
% IE P                   5              $-PARAMETER
% IE P                   10             $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'POWELLBC';

switch(action)

    case 'setup'

    pb.name      = 'POWELLBC';
    pb.sifpbname = 'POWELLBC';
    pbm.name     = 'POWELLBC';
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        if(nargin<2)
            v_('P') = 12;  %  SIF file default value
        else
            v_('P') = varargin{1};
        end
% IE P                   100            $-PARAMETER
% IE P                   500            $-PARAMETER
        v_('1') = 1;
        v_('2') = 2;
        v_('N') = 2*v_('P');
        v_('RN') = v_('N');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2xlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2xlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        pb.xlower = zeros(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = 0.0*ones(pb.n,1);
        pb.xupper = 1.0*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        for I=v_('1'):v_('N')
            v_('RI') = I;
            v_('T') = v_('RI')/v_('RN');
            v_('T') = v_('T')*v_('T');
            pb.x0(ix_(['X',int2str(I)]),1) = v_('T');
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2xlib( 'ii', 'eINVNRM',iet_);
        elftv{it}{1} = 'XJ';
        elftv{it}{2} = 'YJ';
        elftv{it}{3} = 'XK';
        elftv{it}{4} = 'YK';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for K=v_('2'):v_('P')
            v_('K-1') = K-v_('1');
            v_('2K') = v_('2')*K;
            v_('2K-1') = v_('2K')-v_('1');
            for J=v_('1'):v_('K-1')
                v_('2J') = v_('2')*J;
                v_('2J-1') = v_('2J')-v_('1');
                ename = ['E',int2str(K),',',int2str(J)];
                [ie,ie_,newelt] = s2xlib('ii',ename,ie_);
                if(newelt)
                    pbm.elftype{ie} = 'eINVNRM';
                    ielftype(ie) = iet_('eINVNRM');
                end
                vname = ['X',int2str(round(v_('2J-1')))];
                [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
                posev = find(strcmp('XJ',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(round(v_('2K-1')))];
                [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
                posev = find(strcmp('XK',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(round(v_('2J')))];
                [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
                posev = find(strcmp('YJ',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(round(v_('2K')))];
                [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
                posev = find(strcmp('YK',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for K=v_('2'):v_('P')
            v_('K-1') = K-v_('1');
            for J=v_('1'):v_('K-1')
                ig = ig_('OBJ');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['E',int2str(K),',',int2str(J)]);
                pbm.grelw{ig}(posel) = 1.;
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'OBR2-AN-V-0';
        varargout{1} = pb;
        varargout{2} = pbm;

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eINVNRM'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(2,4);
        U_(1,1) = U_(1,1)+1;
        U_(1,3) = U_(1,3)-1;
        U_(2,2) = U_(2,2)+1;
        U_(2,4) = U_(2,4)-1;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        NORM = 1.0/sqrt(IV_(1)*IV_(1)+IV_(2)*IV_(2));
        varargout{1} = NORM;
        if(nargout>1)
            g_(1,1) = -IV_(1)*NORM^3;
            g_(2,1) = -IV_(2)*NORM^3;
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = (3.0*IV_(1)*IV_(1)*NORM^2-1.0)*NORM^3;
                H_(1,2) = 3.0*IV_(1)*IV_(2)*NORM^5;
                H_(2,1) = H_(1,2);
                H_(2,2) = (3.0*IV_(2)*IV_(2)*NORM^2-1.0)*NORM^3;
                varargout{3} = U_.'*H_*U_;
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

