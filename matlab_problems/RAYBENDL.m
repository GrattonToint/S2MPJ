function varargout = RAYBENDL(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem:
%    ********
% 
%    A ray bending problem.  A ray across a inhomogeneous 2D medium is
%    represented by a piecewise linear curve whose knots can be chosen.  
%    The problem is then to optimize the positions of these knots in order 
%    to obtain a ray path corresponding to the minimum travel time from 
%    source to receiver,  according to Fermat principle.
% 
%    The problem becomes harder and harder when the dimesnion increases
%    because the knots are getting closer and closer and the objective
%    has a nondifferentiable kink when two knots coincide.  The difficulty
%    is less apparent when exact second derivatives are not used.
% 
%    Source: a test example in
%    T.J. Moser, G. Nolet and R. Snieder,
%    "Ray bending revisited",
%    Bulletin of the Seism. Society of America 21(1).
% 
%    SIF input: Ph Toint, Dec 1991.
% 
%    classification = 'OXR2-MY-V-0'
% 
%    number of  knots  ( >= 4 )
%    ( n = 2( NKNOTS - 1 ) ) 
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'RAYBENDL';

switch(action)

    case 'setup'

    pb.name      = 'RAYBENDL';
    pb.sifpbname = 'RAYBENDL';
    pbm.name     = 'RAYBENDL';
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        if(nargin<2)
            v_('NKNOTS') = 4;  %  SIF file default value
        else
            v_('NKNOTS') = varargin{1};
        end
%       Alternative values for the SIF file parameters:
% IE NKNOTS              11             $-PARAMETER n = 20
% IE NKNOTS              21             $-PARAMETER n = 40     original value
% IE NKNOTS              32             $-PARAMETER n = 62
% IE NKNOTS              64             $-PARAMETER n = 126
% IE NKNOTS              512            $-PARAMETER n = 1022
% IE NKNOTS              1024           $-PARAMETER n = 2046
        v_('XSRC') = 0.0;
        v_('ZSRC') = 0.0;
        v_('XRCV') = 100.0;
        v_('ZRCV') = 100.0;
        v_('NK-1') = -1+v_('NKNOTS');
        v_('NK-2') = -2+v_('NKNOTS');
        v_('0') = 0;
        v_('1') = 1;
        v_('2') = 2;
        v_('3') = 3;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('0'):v_('NKNOTS')
            [iv,ix_] = s2xlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
            [iv,ix_] = s2xlib('ii',['Z',int2str(I)],ix_);
            pb.xnames{iv} = ['Z',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('NKNOTS')
            [ig,ig_] = s2xlib('ii',['TIME',int2str(I)],ig_);
            gtype{ig} = '<>';
            pbm.gscale(ig,1) = 2.0;
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        pb.xlower = zeros(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower(ix_(['X',int2str(round(v_('0')))]),1) = v_('XSRC');
        pb.xupper(ix_(['X',int2str(round(v_('0')))]),1) = v_('XSRC');
        pb.xlower(ix_(['Z',int2str(round(v_('0')))]),1) = v_('ZSRC');
        pb.xupper(ix_(['Z',int2str(round(v_('0')))]),1) = v_('ZSRC');
        pb.xlower(ix_(['X',int2str(round(v_('NKNOTS')))]),1) = v_('XRCV');
        pb.xupper(ix_(['X',int2str(round(v_('NKNOTS')))]),1) = v_('XRCV');
        pb.xlower(ix_(['Z',int2str(round(v_('NKNOTS')))]),1) = v_('ZRCV');
        pb.xupper(ix_(['Z',int2str(round(v_('NKNOTS')))]),1) = v_('ZRCV');
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        v_('XRANGE') = v_('XRCV')-v_('XSRC');
        v_('ZRANGE') = v_('ZRCV')-v_('ZSRC');
        v_('RKNOTS') = v_('NKNOTS');
        for I=v_('0'):v_('NKNOTS')
            v_('REALI') = I;
            v_('FRAC') = v_('REALI')/v_('RKNOTS');
            v_('XINCR') = v_('FRAC')*v_('XRANGE');
            v_('ZINCR') = v_('FRAC')*v_('ZRANGE');
            v_('XC') = v_('XSRC')+v_('XINCR');
            v_('ZC') = v_('ZSRC')+v_('ZINCR');
            pb.x0(ix_(['X',int2str(I)]),1) = v_('XC');
            pb.x0(ix_(['Z',int2str(I)]),1) = v_('ZC');
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2xlib( 'ii', 'eTT',iet_);
        elftv{it}{1} = 'X1';
        elftv{it}{2} = 'X2';
        elftv{it}{3} = 'Z1';
        elftv{it}{4} = 'Z2';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('NKNOTS')
            v_('I-1') = -1+I;
            ename = ['T',int2str(I)];
            [ie,ie_,newelt] = s2xlib('ii',ename,ie_);
            if(newelt)
                pbm.elftype{ie} = 'eTT';
                ielftype(ie) = iet_('eTT');
            end
            vname = ['X',int2str(round(v_('I-1')))];
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['Z',int2str(round(v_('I-1')))];
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('Z1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['Z',int2str(I)];
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('Z2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('NKNOTS')
            ig = ig_(['TIME',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['T',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'OXR2-MY-V-0';
        varargout{1} = pb;
        varargout{2} = pbm;

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'e_globs'

        pbm = varargin{1};
        pbm.efpar = [];
        pbm.efpar(1) = 0.01;
        varargout{1} = pbm;

    case 'eTT'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(3,4);
        U_(3,1) = U_(3,1)-1;
        U_(3,2) = U_(3,2)+1;
        U_(1,3) = U_(1,3)+1;
        U_(2,4) = U_(2,4)+1;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        IV_(3) = U_(3,:)*EV_;
        C0 = 1.0+pbm.efpar(1)*IV_(1);
        C1 = 1.0+pbm.efpar(1)*IV_(2);
        DCDZ = pbm.efpar(1);
        V = 1.0/C1+1.0/C0;
        VDZ0 = -DCDZ/(C0*C0);
        VDZ1 = -DCDZ/(C1*C1);
        VDZ0Z0 = 2.0*DCDZ*DCDZ/C0^3;
        VDZ1Z1 = 2.0*DCDZ*DCDZ/C1^3;
        DZ1 = IV_(2)-IV_(1);
        R = sqrt(IV_(3)*IV_(3)+DZ1*DZ1);
        RDX = IV_(3)/R;
        RDZ1 = DZ1/R;
        RDZ0 = -RDZ1;
        RDXDX = (1.0-IV_(3)*IV_(3)/(R*R))/R;
        RDXZ1 = -IV_(3)*DZ1/R^3;
        RDXZ0 = -RDXZ1;
        RDZ1Z1 = (1.0-DZ1*DZ1/(R*R))/R;
        RDZ0Z0 = RDZ1Z1;
        RDZ0Z1 = -RDZ1Z1;
        varargout{1} = V*R;
        if(nargout>1)
            g_(3,1) = V*RDX;
            g_(1,1) = V*RDZ0+VDZ0*R;
            g_(2,1) = V*RDZ1+VDZ1*R;
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(3,3) = V*RDXDX;
                H_(3,1) = VDZ0*RDX+V*RDXZ0;
                H_(1,3) = H_(3,1);
                H_(3,2) = VDZ1*RDX+V*RDXZ1;
                H_(2,3) = H_(3,2);
                H_(1,1) = V*RDZ0Z0+VDZ0Z0*R+2.0*VDZ0*RDZ0;
                H_(1,2) = V*RDZ0Z1+VDZ1*RDZ0+VDZ0*RDZ1;
                H_(2,1) = H_(1,2);
                H_(2,2) = V*RDZ1Z1+VDZ1Z1*R+2.0*VDZ1*RDZ1;
                varargout{3} = U_.'*H_*U_;
            end
        end

    %%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    case {'fx','fgx','fgHx','cx','cJx','cJHx','cIx','cIJx','cIJHx','cIJxv','fHxv',...
          'cJxv','Lxy','Lgxy','LgHxy','LIxy','LIgxy','LIgHxy','LHxyv','LIHxyv'}

        if(isfield(pbm,'name')&&strcmp(pbm.name,name))
            pbm.has_globs = [1,0];
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

