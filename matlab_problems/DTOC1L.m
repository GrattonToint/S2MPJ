function varargout = DTOC1L(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : DTOC1L
%    *********
% 
%    This is a discrete time optimal control (DTOC) problem.  
%    The system has N time periods, NX control variables and NY state variables.
%    The parameter mu in the original problem formulation is set to zero, 
%    yielding linear transition functions, hence the L in the problem's name.
% 
%    The problem is convex.
% 
%    Sources: problem 1 (with mu = 0) in
%    T.F. Coleman and A. Liao,
%    "An Efficient Trust Region Method for Unconstrained Discret-Time Optimal
%    Control Problems",
%    Tech. Report, ctc93tr144,  Advanced Computing Research Institute, 
%    Cornell University, 1992.
% 
%    L.Z. Liao and C.A. Shoemaker,
%    "Advantages of differential dynamic programming over Newton's method for
%    discrete-time optimal control problems",
%    Tech. Report ctc92tr97, Advanced Computing Research Institute, 
%    Cornell University, 1992.
% 
%    SIF input: Ph. Toint, August 1993
% 
%    classification = 'OLR2-AN-V-V'
% 
%    Problem variants: they are identified by the values of
%    the parameter vector ( N, NX, NY )
% 
%    The problem has (N-1)*NX+N*NY  variables (of which NY are fixed),
%    and (N-1)*NY constraints
% 
%       Alternative values for the SIF file parameters:
% IE N                   10             $-PARAMETER # periods  } original value
% IE NX                  2              $-PARAMETER # controls } n=   58, m=  36
% IE NY                  4              $-PARAMETER # states   }
% 
% IE N                   50             $-PARAMETER # periods  }
% IE NX                  2              $-PARAMETER # controls } n=  298, m= 196
% IE NY                  4              $-PARAMETER # states   }
% 
% IE N                   100            $-PARAMETER # periods  }
% IE NX                  2              $-PARAMETER # controls } n=  598, m= 396
% IE NY                  4              $-PARAMETER # states   }
% 
% IE N                   500            $-PARAMETER # periods  }
% IE NX                  2              $-PARAMETER # controls } n= 2998, m=1996
% IE NY                  4              $-PARAMETER # states   }
% 
% IE N                   1000           $-PARAMETER # periods  }
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'DTOC1L';

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
% IE NX                  2              $-PARAMETER # controls } n= 5998, m=3996
        if(nargs<2)
            v_('NX') = 2;  %  SIF file default value
        else
            v_('NX') = varargin{2};
        end
% IE NY                  4              $-PARAMETER # states   }
        if(nargs<3)
            v_('NY') = 4;  %  SIF file default value
        else
            v_('NY') = varargin{3};
        end
% IE N                   10             $-PARAMETER # periods  }
% IE NX                  5              $-PARAMETER # controls } n=  145, m=  90
% IE NY                  10             $-PARAMETER # states   }
% IE N                   50             $-PARAMETER # periods  }
% IE NX                  5              $-PARAMETER # controls } n=  745, m= 490
% IE NY                  10             $-PARAMETER # states   }
% IE N                   100            $-PARAMETER # periods  }
% IE NX                  5              $-PARAMETER # controls } n= 1495, m= 990
% IE NY                  10             $-PARAMETER # states   }
% IE N                   500            $-PARAMETER # periods  }
% IE NX                  5              $-PARAMETER # controls } n= 7495, m=4990
% IE NY                  10             $-PARAMETER # states   }
% IE N                   1000           $-PARAMETER # periods  }
% IE NX                  5              $-PARAMETER # controls } n=14995, m=9990
% IE NY                  10             $-PARAMETER # states   }
        v_('N-1') = -1+v_('N');
        v_('1') = 1;
        v_('2') = 2;
        v_('NY-1') = -1+v_('NY');
        v_('NX+NY') = v_('NX')+v_('NY');
        v_('RXY') = v_('NX+NY');
        v_('1/RXY') = 1.0/v_('RXY');
        for J=v_('1'):v_('NX')
            for I=v_('1'):v_('NY')
                v_('I-J') = I-J;
                v_('RI-J') = v_('I-J');
                v_(['B',int2str(I),',',int2str(J)]) = v_('RI-J')*v_('1/RXY');
            end
        end
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for T=v_('1'):v_('N-1')
            for I=v_('1'):v_('NX')
                [iv,ix_] = s2mpjlib('ii',['X',int2str(T),',',int2str(I)],ix_);
                pb.xnames{iv} = ['X',int2str(T),',',int2str(I)];
            end
        end
        for T=v_('1'):v_('N')
            for I=v_('1'):v_('NY')
                [iv,ix_] = s2mpjlib('ii',['Y',int2str(T),',',int2str(I)],ix_);
                pb.xnames{iv} = ['Y',int2str(T),',',int2str(I)];
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for T=v_('1'):v_('N-1')
            for I=v_('1'):v_('NX')
                [ig,ig_] = s2mpjlib('ii',['OX',int2str(T),',',int2str(I)],ig_);
                gtype{ig} = '<>';
                iv = ix_(['X',int2str(T),',',int2str(I)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
            end
        end
        for T=v_('1'):v_('N')
            for I=v_('1'):v_('NY')
                [ig,ig_] = s2mpjlib('ii',['OY',int2str(T),',',int2str(I)],ig_);
                gtype{ig} = '<>';
                iv = ix_(['Y',int2str(T),',',int2str(I)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
            end
        end
        for T=v_('1'):v_('N-1')
            v_('T+1') = 1+T;
            [ig,ig_] = s2mpjlib('ii',['TT',int2str(T),',',int2str(round(v_('1')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['TT',int2str(T),',',int2str(round(v_('1')))];
            iv = ix_(['Y',int2str(round(v_('T+1'))),',',int2str(round(v_('1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['TT',int2str(T),',',int2str(round(v_('1')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['TT',int2str(T),',',int2str(round(v_('1')))];
            iv = ix_(['Y',int2str(T),',',int2str(round(v_('1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 0.5+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 0.5;
            end
            iv = ix_(['Y',int2str(T),',',int2str(round(v_('2')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 0.25+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 0.25;
            end
            for I=v_('1'):v_('NX')
                [ig,ig_] = s2mpjlib('ii',['TT',int2str(T),',',int2str(round(v_('1')))],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['TT',int2str(T),',',int2str(round(v_('1')))];
                iv = ix_(['X',int2str(T),',',int2str(I)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) =...
                          v_(['B',int2str(round(v_('1'))),',',int2str(I)])+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_(['B',int2str(round(v_('1'))),',',int2str(I)]);
                end
            end
            for J=v_('2'):v_('NY-1')
                v_('J-1') = -1+J;
                v_('J+1') = 1+J;
                [ig,ig_] = s2mpjlib('ii',['TT',int2str(T),',',int2str(J)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['TT',int2str(T),',',int2str(J)];
                iv = ix_(['Y',int2str(round(v_('T+1'))),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
                iv = ix_(['Y',int2str(T),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 0.5+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 0.5;
                end
                iv = ix_(['Y',int2str(T),',',int2str(round(v_('J-1')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -0.25+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -0.25;
                end
                iv = ix_(['Y',int2str(T),',',int2str(round(v_('J+1')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 0.25+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 0.25;
                end
                for I=v_('1'):v_('NX')
                    [ig,ig_] = s2mpjlib('ii',['TT',int2str(T),',',int2str(J)],ig_);
                    gtype{ig}  = '==';
                    cnames{ig} = ['TT',int2str(T),',',int2str(J)];
                    iv = ix_(['X',int2str(T),',',int2str(I)]);
                    if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                        pbm.A(ig,iv) = v_(['B',int2str(J),',',int2str(I)])+pbm.A(ig,iv);
                    else
                        pbm.A(ig,iv) = v_(['B',int2str(J),',',int2str(I)]);
                    end
                end
            end
            [ig,ig_] =...
                  s2mpjlib('ii',['TT',int2str(T),',',int2str(round(v_('NY')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['TT',int2str(T),',',int2str(round(v_('NY')))];
            iv = ix_(['Y',int2str(round(v_('T+1'))),',',int2str(round(v_('NY')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] =...
                  s2mpjlib('ii',['TT',int2str(T),',',int2str(round(v_('NY')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['TT',int2str(T),',',int2str(round(v_('NY')))];
            iv = ix_(['Y',int2str(T),',',int2str(round(v_('NY')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 0.5+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 0.5;
            end
            iv = ix_(['Y',int2str(T),',',int2str(round(v_('NY-1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -0.25+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -0.25;
            end
            for I=v_('1'):v_('NX')
                [ig,ig_] =...
                      s2mpjlib('ii',['TT',int2str(T),',',int2str(round(v_('NY')))],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['TT',int2str(T),',',int2str(round(v_('NY')))];
                iv = ix_(['X',int2str(T),',',int2str(I)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_(['B',int2str(round(v_('NY'))),',',int2str(I)])+...
                         pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_(['B',int2str(round(v_('NY'))),',',int2str(I)]);
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
        for T=v_('1'):v_('N-1')
            for I=v_('1'):v_('NX')
                pbm.gconst(ig_(['OX',int2str(T),',',int2str(I)])) = -0.5;
            end
        end
        for T=v_('1'):v_('N')
            for I=v_('1'):v_('NY')
                pbm.gconst(ig_(['OY',int2str(T),',',int2str(I)])) = -0.25;
            end
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        for I=v_('1'):v_('NY')
            pb.xlower(ix_(['Y',int2str(round(v_('1'))),',',int2str(I)]),1) = 0.0;
            pb.xupper(ix_(['Y',int2str(round(v_('1'))),',',int2str(I)]),1) = 0.0;
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2mpjlib('ii','gL4',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for T=v_('1'):v_('N-1')
            for I=v_('1'):v_('NX')
                ig = ig_(['OX',int2str(T),',',int2str(I)]);
                pbm.grftype{ig} = 'gL4';
            end
        end
        for T=v_('1'):v_('N')
            for I=v_('1'):v_('NY')
                ig = ig_(['OY',int2str(T),',',int2str(I)]);
                pbm.grftype{ig} = 'gL4';
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
% LO SOLUTION(  10,2, 4) 0.0735931360
% LO SOLUTION(  50,2, 4) 0.2299411960
% LO SOLUTION( 100,2, 4) 0.4253681120
% LO SOLUTION( 500,2, 4) 1.9887794988
% LO SOLUTION(1000,2, 4) 3.9430507151
% LO SOLUTION(  10,5,10) 1.1498579294
% LO SOLUTION(  50,5,10) 6.1678479713
% LO SOLUTION( 100,5,10) 12.439954329
% LO SOLUTION( 500,5,10) 62.616843379
% LO SOLUTION(1000,5,10) 125.33793359
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons   = [1:length(pbm.congrps)];
        pb.pbclass = 'OLR2-AN-V-V';
        pb.x0          = zeros(pb.n,1);
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

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gL4'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = GVAR_^4;
        if(nargout>1)
            g_ = 4.0*GVAR_^3;
            varargout{2} = g_;
            if(nargout>2)
                H_ = 12.0*GVAR_^2;
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

