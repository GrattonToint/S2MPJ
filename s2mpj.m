%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ probname, exitc, errors ] = s2mpj( sifpbname, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   S2MPJ: a SIF decoder to Matlab, Python and Julia
%
%   The purpose of this function is to parse the content of a (correct) SIF problem file into a Matlab, Python or Julia
%   file, which can then be used for computing values of the  objective function and constraints, as well as that of
%   their derivatives, directly within Matlab/Python/Julia and without any additional installation or interfacing with
%   MEX files or Fortran programs. It also allows to compute the product of the objective function's or Lagrangian's
%   Hessian or constraints' Jacobian times a user-supplied vector. Importantly, the produced Matlab/Python file,
%   hereafter called the 'output file', is not extensive, meaning that loops are not unrolled (thereby maintaining a
%   reasonably compact description).
%
%   In what follows, the optimization problem under consideration is assumed to be of the form
%
%      min_x  f(x)
%
%   such that
%
%      xlower_i <= x_i <= xupper_i   (i=1,...,n)
%
%   and
%
%      clower_i <= c_i(x)  <= cupper_i   (i=1 ..., m).
%
%   In this description, we have that
%
%      f(x) = \sum_{i in OG} F_i[ a_i(x) ] / gscale_i
%
%   where
%
%      a_i(x) = \sum_{j in GE_i} w_{i,j} f_j(x) - A(i,:)*x - constant_i
%
%   with OG being the set of "objective groups" and GE_i the set of nonlinear elements of group i.
%   We also have that
%
%      c_i(x) = C_i[ b_i(x) ] / gscale_i
%
%   where
%
%      b_i(x) = \sum_{j in GE_i} w_{i,j} f_j(x) - A(i,:)*x - constant_i.
%
%   The functions F_i[.] and C_i[.] are the (possibly nonlinear) univariate "group functions". 
%
%   NOTE: This format for the constraints *differs* from that used in the SIF standard where they are given, depending
%         on their type, by
%
%            0 >= c_i(x) >= range_i.   ( type = <= )
%
%                 c_i(x) = 0,          ( type = == )
%
%         or
%
%            0 <= c_i(x) <= range_i,   ( type = >= )
%
%         for some types (<=, == or >=),  constants (constant_i) and ranges (range_i) (possibly) specified in the
%         problem's SIF file (see Section 2.1 of the SIF report).  Note that range_i <= 0 for constraints of type <=,
%         and range_i >=0 for constraints of type >=.
%
%   We let c(x) = [ c_1(x) ; ...; c_m(x) ]. We also consider the associated Lagrangian function
%
%      L(x,y) = f(x) + y'*c(x).
%
%   As it may convenient to restrict one's attention to a subset
%      
%      I = [ i_1, ..., i_p]  
%
%   of the constraints' indeces {1,...,m}, we also consider the 'I-restricted' variants of c(x) and L(x,y) given by
%
%      c_I(x) = [ c_{i_1}(x) ; ... ; c_{i_p}(x) ]
%
%   and
%
%      L_I(x,y) = f(x) + y_I'*c_I(x).
%
%   According to widespread convention, the gradients of f(x), L(x,y) and L_I(x,y) are column vectors while the Jacobian
%   of the constraints has its rows defined as the gradients of individual constraints. When S2MPJ produces a Matlab file,
%   gradients, Jacobians and Hessians are produced in Matlab sparse format, in Python CSR matrix format or in Julia
%   sparse matrix format.

%   METHOD OVERVIEW
%
%   The sifpbname.SIF SIF file (which must be in the Matlab path) is read one line at a time, and each line is then
%   translated, sometimes is a somewhat indirect way, into the corresponding Matlab/Python/Julia commands. The decoder
%   - first sets up the problem data structure (including the value of its constant parameters) in the 'setup' action
%     (ie. 'case' statement in the Matlab output file or __init__ method in the Python output file or the outer
%     if-then-else statement in the Julia output file);
%   - then includes other actions defining the involved nonlinear functions, as extracted from the Fortran statements in
%     the SIF nonlinear ELEMENTS and GROUPS sections;
%   - adds a final call to the S2MPJ-supplied s2mpjlib.m (for Matlab), s2mpjlib.py (for Python) or s2mpjlib.jl (for Julia)
%     that accesses the previously defined actions and data-structures to perform the required evaluation tasks.
%   The decoder then produces an output file called 'probname.m', 'probname.py' or 'probname.jl' in the current directory,
%   possibly overwriting an existing file with the same name (see the S2MPJ OPTIONS section below to modify the directory
%   where the SIF file is found and that where the output file is written).  This file is intended for direct calls
%   producing values of the objective function, constraints or Lagrangian (possibly with their derivatives), as well as
%   products of the objective function's or Lagrangian's Hessian or constraints' Jacobian times a user-supplied vector
%   (see below).
%
%   The name probname of the output file is, in most cases, identical to sifpbname, but differs from it when the sifpbname
%   string starts with a digit (it is then prefixed by 'n') or contains one or more of the characters '+', '-', '*' or '/'
%   (these characters are then replaced by 'p', 'm', 't', and 'd', respectively). This renaming is necessary to allow the
%   output file to be used as a Matlab/Julia function or as a Python class.
%
%   Within the 'setup'/__init__ action, instructions are written in the output file to define the various problem
%   parameters, as well as vectors of bounds on the variables, constraints and objective, and the variables/multipliers
%   starting values. Loops in the SIF file are directly transformed into loops in the output file. Structured sets of the
%   SIF files (variables, groups, elements) are characterized by the fact that their components may be defined using
%   different names and multi-indexes. Because most optimization codes only recognize linear structures (a vector of
%   variables, a vector of constraints), S2MPJ transforms variables, groups and elements into linear 'flat' unidimensional
%   structures. Since this transformation may depend on problem internal parameters such as loop limits, themselves
%   depending on problem input arguments only known at runtime, the output file uses the S2MPJ-supplied function
%   s2mpj_ii(...) to compute the relevant indices within 'setup'/__init__ at runtime.
%
%   The names of SIF entities have a potentially complex structure, possibly involving symbolic or numerical indices, but
%   also numerical digits within the names which may (or may not) correspond to indeces. S2MPJ attempts to sort this out in
%   the s2mpjname function, in order to assign a Matlab/Python/Julia-comptatible name to each of these entities. These
%   names are then used in a dictionaries associating names and values (v_ for parameters) or index (ix_ for variables,
%   ig_ for groups, iet_ for element types, ie_ for elements and itg_ for group types).
%
%   Note that the SIF file may also include scaling factors for the variables of the groups' linear terms.  These are read
%   by S2MPJ and are, by default, applied to the relevant coefficients during 'setup', in order to speed up later
%   evaluation calls. Again this *differs* from what is done by the historical Fortran decoder, in which variable scaling
%   is ignored entirely (this less efficient approach may nevertheless be requested as an option).
%
%   Special actions in the output file ('e@globs" and 'g@globs') are finally needed for the GLOBALS subsections in the
%   nonlinear ELEMENTS and GROUPS sections. The definitions in both these special acctions are incorporated in the problem
%   data structure just before the main call to evaluations, and then checked for within s2mpjlib (meaning that the
%   'e@glob' and 'g@globs' actions may be called there if needed).
%
%   The code uses three different structs to pass information around:
%   pbs  is the structure internal to the decoding during the 'setup'/__init__ action and is used internally by the
%        Matlab s2mpj to pass information on entities like loops, element/ group ftype(s) across the different
%        subfunctions used in this action;
%   pb   is the structure used to pass information on the problem to the user once the setup/__init__ action is
%        completed;
%   pbm  is the stucture used to pass information between the different actions of the output file (occuring on
%        successive calls) and the nonlinear functions' evaluation functions/methods in s2mpjlib.
%
%   Input:
%
%      sifpbname  : a string containing the name of the problem to be decoded (S2MPJ then reads the 'probname.SIF' file for
%                   input).
%      varargin{1}: if present, vararagin{1} allows the specification of S2MPJ options (see the section on S2MPJ OPTIONS
%                   below). The use of the suitable option is mandatory for decoding SIF files into a Python or a Julia
%                   output file.
%
%   Output:
%
%      probname: a string containing the name of the Matlab or Python output file (without the .m, .py or .jl suffix),
%                see above.
%
%      exitc   : the number of errors which occured before termination (or crash). exitc = 0 thus indicates error-free
%                execution.
%      errors  : a cell/list of length exitc, whose entries contain a brief description of the error(s) found.
%           
%      Beyond those arguments, S2MPJ typically produces an output file called 'probname.m', 'probname.py'  or 'probname.jl'
%      in the current directory, possibly overwriting an existing file with the same name. 
%
%      The Matlab output file
%      ----------------------
%
%      This Matlab output file is is a Matlab function designed to provide the following user interface. A call of the
%      form
%
%         varargout = probname( action, varargin )
%
%      produces result(s) stored in varargout according to the following choices of the 'action' argument:
%
%      action = 'setup'  ( varargin = { problem_parameters } )
%
%         This action outputs two structs called 'pb' and 'pbm'. The fields of 'pb' are defined as follows:
%
%         pb.name      is a string containing the problem's name
%         pb.sifpbname if present, is the original name of the SIF file before its modification to ensure Matlab/Python
%                      compatibility
%         pb.n         is the problem's number of variables,
%         pb.nob       is the number of objective groups
%         pb.nle       is the number of <= constraints
%         pb.neq       is the number of == constraints
%         pb.nge       is the number of >= constraints
%         pb.m         is the total number of constraints
%         pb.licons    is a integer vector containing the indeces of the linear constraints
%         pb.pbclass   is the problem's SIF classification
%         pb.x0        is the starting points for the variables
%         pb.xlower    is the real vector of lower bounds on the variables
%         pb.xupper    is the real vector of upper bounds on the variables
%         pb.xtype     is a string whose i-th position describes the type of the i-th variable
%                      'r' : the variable is real
%                      'i' : the variable is integer
%                      'b' : the variable is binary (zero-one)
%                      If 'xtype' is not a field of pb, all variables are assumed to be of the 'r' type.
%         pb.xscale    if nonempty, is a real vector containing the scaling factors to be applied by the user on the
%                      occurrences of the variables in the linear terms of the problem's groups (that is on the columns
%                      of pbm.A).  If empty, the scaling factors are applied to the relevant terms internally to the
%                      decoder without need for further user action (see the pbxscale option of s2mpj below)
%         pb.y0        is a real vector containing the starting values for the constraint's multipliers
%         pb.clower    is the real vector of lower bounds on the constraints values
%         pb.cupper    is the real vector of upper bounds on the constraints values
%         pb.objlower  is a real lower bound of the objective function's value, if known
%         pb.objupper  is a real upper bound of the objective function's value, if known
%         pb.xnames    is a cell of strings containing the name of the variables
%         pb.cnames    is a cell of strings containing the name of the constraints
%         pb.objderlvl is the level of derivatives available for the objective function
%                       (0 = function value, 1 = function value and gradient, 2 = function value, gradient and Hessian)
%                       NOTE: all CUTEst problems have objective function's value, gradient and Hessian available.
%         pb.conderlvl is the level of derivatives available for the constraints.  If of length = 1, the level is
%                       the same for all constraints.  Otherwise, it is of length pb.m and the derivative level
%                       for the i-th constraint (that is the constraint whose group number is pb.congrps(i))
%                       is given by pbm.conderlvl(i).
%                       NOTE: all CUTEst problems have value, gradient and Hessian available for all constraints.
%
%         Here, 'problem parameters' are parameters identified by a $-PARAMETER string in the SIF file. They are assigned
%         in the order where they appear in varargin, which is the same as that used in the SIF and Matlab/Python/Julia
%         files. If varargin is too short in that it does not provide a value for a parameter, the SIF default value is
%         used.
%
%         If constraints are present, they are reported in the following order:
%
%                    1:pb.nle                     : constraints of <= type,
%                pb.nle+1:pb.nle+pb.neq           : constraints of == type,
%            pb.nle+pb.neq+1:pb.nle+pb.neq+pb.nge : constraints of >= type.
%
%         This is only for convenience, but this 'ordered' feature  should not be included in a 'S2MPJ format' for problem
%         description, the type of each constraint being determined by the values the lower and upper bounds on its value.
%         This may also be superseded by setting the keepcorder flag in the S2MPJ options (see below). If the problem has
%         no constraints or bound constraints only, then lincons, y0, clower, cupper and cnames are not fields of the pb
%         struct returned on setup. pb.m is still defined in this case, but its value is 0. The fields objlower
%         and objupper may also be missing if no value is provided in the SIF file. When specific S2MPJ options are used
%         (again, see below), the fields pb.xnames, pb.cnames may be missing from the pb struct.
%
%         The fields of 'pbm' are defined as follows:
%
%         pbm.name      is a string containing the name of the problem
%         pbm.objgrps   is the list of indeces of the objective groups
%         pbm.congrps   is the list of indeces of the contraint groups
%         pbm.A         is the matrix whose lines contain the coefficients of the linear terms in the groups
%         pbm.gconst    is the list of group's constants
%         pbm.H         is the Hessian matrix specified in the QUADRATICS section
%         pbm.enames    is a cell of strings containing the names of the nonlinear elements
%         pbm.elftype   is a cell of strings containing the names of the element's ftype
%         pbm.elvar     is a cell of integer vectors containing the indeces of the elemental variables for each element
%         pbm.elpar     is a cell of real vectors containing the values of the elemental parameters for each element
%         pbm.gscale    is a real vector containing the group's scaling factors
%         pbm.grnames   is a cell of strings containing the group's names
%         pbm.grftype   is a cell of strings containing the group's ftypes
%         pbm.grelt     is a cell of integer vectors containing the indeces of the nonlinear elements occuring in each
%                       group
%         pbm.grelw     is a cell of real vectors containing the weights of the nonlinear elements occuring in each group
%         pbm.grpar     is a cell of real vectors containing the values of the parameters of the group functions
%         pbm.efpar     is a real vector containing the values of global ELEMENTS parameters, once computed in the e_globs
%                       section using their Fortran expression
%         pbm.gfpar     is a real vector containing the values of global GROUPS parameters, once computed in the g_globs
%                        section using their Fortran expression
%         pbm.objderlvl is the level of derivatives available for the objective function
%                       (0 = function value, 1 = function value and gradient, 2 = function value, gradient and Hessian)
%                       NOTE: all CUTEst problems have value, gradient and Hessian available.
%         pbm.conderlvl is the level of derivatives available for the constraints.  If of length = 1, the level is
%                       the same for all constraints.  Otherwise, it is of length pb.m and the derivative level for the i-th
%                       constraint (that is the constraint whose group number is pbm.congrps(i)) is given by pbm.conderlvl(i).
%                       NOTE: all CUTEst problems have value, gradient and Hessian available.
%         pbm.ndigs     if present and >0, indicates that reduced precision arithmetic with pbm.ndigs digits is requested
% 
%         Depending on the problem's nature, some fields may be missing from the pbm struct.  Only the fields name and
%         A are mandatory (a problem only containing those would be an homogeneous linear system). The presence of the
%         fields enames and gnames in the pbm struct depends on the specific S2MPJ options defined (see below).
%
%        A call to probname( 'setup', ...) MUST PRECEDE any call to probname with other actions (so as to setup the
%        problem's data structure).
%
%     action = 'fx'                      ( varargin = { x }
%
%         varargout{1} = the value of the objective function at x
%
%     action in { 'fgx', ['fx'] }        ( varargin = { x }
%
%         varargout{1} = the value of the objective function at x
%         varargout{2} = the vector giving the gradient of the objective function at x
%
%     action in { 'fgHx', ['fx'] }       ( varargin = { x }
%
%         varargout{1} = the value of the objective function at x
%         varargout{2} = the vector giving the gradient of the objective function at x
%         varargout{3} = the square symmetric sparse matrix giving the Hessian of the objective function at x
%
%     action = 'cx'                      ( varargin =  x )
%
%         varargout{1} = the vector of contraints' values at x
%
%     action in { 'cJx', ['cx'] }        ( varargin =  x )
%
%         varargout{1} = the vector of contraints' values at x
%         varargout{2} = the constraints' Jacobian at x
%
%     action in { 'cJHx', ['cx'] }       ( varargin =  x )
%
%         varargout{1} = the vector of contraints' values at x
%         varargout{2} = the constraints' Jacobian at x
%         varargout{3} = a cell of length m, whose j-th component is Hc_j(x), the Hessian of the i-th constraint at x
%                        (j = 1, ..., m)
%
%     action = 'cIx'                     ( varargin = ( x, I ) )
%
%         varargout{1} = the vector of values of the 'I-restricted' constraints c__I(x)
%
%     action in { 'cIJx', ['cIx'] }      ( varargin = ( x, I ) )
%
%         varargout{1} = the vector of values of the 'I-restricted' constraints c__I(x)
%         varargout{2} = the Jacobian of c_I(x) at x
%
%     action in { 'cIJHx', ['cIx'] }     ( varargin = ( x, I ) )
%
%         varargout{1} = the vector of values of the 'I-restricted' constraints c__I(x)
%         varargout{2} = the Jacobian of c_I(x) at x
%         varargout{3} = a cell of length p, whose j-th component is Hc_j(x), the Hessian of the j-th constraint at x
%                       (j= i_1, ...., i_p)
%
%     action in 'fHxv'                   ( varargin = ( x, v ) )
%
%         varargout{1} = the product of the objective function's Hessian at x times the vector v
%
%     action in 'cJxv'                   ( varargin = ( x, v ) )
%
%         varargout{1} = the product of the constraints' Jacobian at x times the vector v
%
%     action in 'cIJxv'                  ( varargin = ( x, v, I ) )
%
%         varargout{1} = the product of the 'I-restricted' Jacobian J_I at x times the vector v
%
%     action = 'Lxy'                     ( varargin = ( x, y ) )
%
%         varargout{1} = the value of the Lagrangian function L(x,y) at ( x, y )
%
%     action in { 'Lgxy', ['Lxy'] }      ( varargin = ( x, y ) )
%
%         varargout{1} = the value of the Lagrangian function L(x,y) at ( x, y )
%         varargout{2} = the vector giving the gradient wrt x of the Lagrangian at ( x, y )
%
%     action in { 'LgHxy',  ['Lxy'] }    ( varargin = ( x, y ) )
%
%         varargout{1} = the value of the Lagrangian function L(x,y) at ( x, y )
%         varargout{2} = the vector giving the gradient wrt x of the Lagrangian at ( x, y )
%         varargout{3} = the square symmetric sparse matrix giving the Hessian  wrt x of the Lagrangian function at
%                        (x, y )
%
%     action = 'LIxy'                    ( varargin = ( x, y_I, I ) )
%
%         varargout{1} = the value of the 'I-restricted' Lagrangian function L_I(x,y_I) at ( x, y_I )
%         varargout{2} = the vector giving the gradient wrt x of the 'I-restricted' Lagrangian at ( x, y_I )
%         varargout{3} = the square symmetric sparse matrix giving the Hessian  wrt x of the 'I-restricted' Lagrangian
%                        function at ( x, y_I )
%
%     action in { 'LIgxy', ['LIxy'] }    ( varargin = ( x, y_I, I ) )
%
%         varargout{1} = the value of the 'I-restricted' Lagrangian function L_I(x,y_I) at ( x, y_I )
%         varargout{2} = the vector giving the gradient wrt x of the 'I-restricted' Lagrangian at ( x, y_I )
%
%     action = { 'LIgHxy', ['LIxy'] }    ( varargin = ( x, y_I, I ) )
%
%         varargout{1} = the value of the 'I-restricted' Lagrangian function L_I(x,y_I) at ( x, y )
%         varargout{2} = the vector giving the gradient wrt x of the 'I-restricted' Lagrangian at ( x, y_I )
%         varargout{3} = the square symmetric sparse matrix giving the Hessian wrt x of the 'I-restricted' Lagrangian
%                        function at ( x, y_I )
%
%     action ='LHxyv'                   ( varargin = ( x, y, v ) )
%
%         varargout{1} = the product of the Lagragian's Hessian wrt x at ( x, y ) times the vector v
%
%     action = 'LIHxyv'                 ( varargin = ( x, y_I, v, I ) )
%
%         varargout{1} = the product of the 'I-restricted' Lagragian's Hessian wrt x at ( x, y_I ) times the vector v
%
%
%     Note that actions 'cx', 'cIx', 'Jxv', 'JIxv', 'Lxy', 'LIxy', HLxyv' or 'HLxyv' are of course meaningless when the
%     problem is unconstrained or only has bound constraints.
%
%     In the above description, the actions in square brackets [] are optional variants for action names.
%
%     The Python output file
%     ----------------------
%
%     The Python output file is a Python class whose name is that of the problem and which is derived from the parent
%     CUTEst_problem class described in the s2mpjlib.py file. Its inherits the evaluation methods of this parent class,
%     which are organized in a manner similar to actions for the Matlab output file, except that the Matlab 'setup' action
%     is replaced by a simple call to the class with the problem parameters, the setting up of the problem structure(s)
%     being then performed by the __init__ method of the class. So, if the problem is given by the probname class, setting
%     it up is achieved by the call
%
%          probname(problems parameters)
%
%     while the evaluation task corresponding to the Matlab call [outputs] = probname( action, x, [other args] ) is
%     performed by the Python call
%
%          ouputs = probname.action( x, [other args] )
%
%     for all actions (except 'setup') listed above. The real vectors resulting from the various evaluations or being
%     present in the pb and pbm structures are column-oriented numpy.ndarrays.  Real matrices are sparse.csr matrices, or
%     lists of sparse.csr matrices (for the Hessians of the constraints). Note that pb and pbm struct have the same fields
%     as those described for the Matlab case, subject to the same conditions for being present or not.
%
%     Note that all fields of the pb and pbm structs in the Matlab case are now field of the probname class.
%
%     The Julia output file
%     ---------------------
%
%     The Julia output file is a cross between the Matlab and Python ones. It is a Julia function whose interface is
%     identical to that of the Matlab output, EXCEPT that the second argument in the call must be the pbm struct. This
%     struct now has an additional field "call" so that pbm.call contains the Julia call to the output file after
%     execution with the "setup" action.
%
%   TYPICAL USE OF MATLAB OUTPUT FILES
%
%   The use of S2MPJ is typically as follows.
%   1) a SIF problem (PROBLEM.SIF, say) is decoded by the Matlab command
%
%         PROBLEM = s2mpj( 'SIFPROBLEM' );
%
%      or
%
%         PROBLEM = s2mpj( 'SIFPROBLEM', options );
%
%      where options is a struct whose content is described below in the S2MPJ OPTIONS paragraph. This assumes that the
%      PROBLEM.SIF file is in the Matlab path. A file PROBLEM.m is then produced in the current directory (potentially
%      overwriting an existing one with the same name).
%
%   2) the problem data structure is then setup, along with the starting points, bounds and other components of the pb
%      struct (see above) by issuing the command
%
%         pb = PROBLEM( 'setup', args{:} )
%
%      where args{:} is a comma separated list of problem-dependent arguments (such as problem's dimension, for instance)
%      identified by the string $-PARAMETER in the SIF file;
%
%    3) the value(s) of the problem functions (objective and constraints) at a vector x, together with values of their
%       first and second derivatives (if requested) are then computed by issuing one of the commands
%
%         fx              = PROBLEM( 'fx', x );  
%         [ fx, gx ]      = PROBLEM( 'fgx', x );             or     [ fx, gx ]      = PROBLEM( 'fx', x );
%         [ fx, gx, Hx ]  = PROBLEM( 'fgHx', x );            or     [ fx, gx, Hx ]  = PROBLEM( 'fx', x );
%
%       for the objective function and (if constraints are present)
%
%         cx              = PROBLEM( 'cx', x );
%         [ cx, Jx ]      = PROBLEM( 'cJx', x );             or     [ cx, Jx ]      = PROBLEM( 'cx', x );
%         [ cx, Jx, Hx ]  = PROBLEM( 'cJHx', x );            or     [ cx, Jx, Hx ]  = PROBLEM( 'cx', x );
%
%       for the constraints.  In this case, it is possible to restrict one's attention to a subset I of the constraints
%       (ie using the 'I-restricted' version of c(x)) by using
%
%         cix                = PROBLEM( 'cIx', x, I ); 
%         [ cIx, JIx ]       = PROBLEM( 'cIJx', x, I );      or    [ cIx, JIx ]       = PROBLEM( 'cIx', x, I );
%         [ cIx, JIx, HIx ]  = PROBLEM( 'cIJHx', x, I );     or    [ cIx, JIx, HIx ]  = PROBLEM( 'cIx', x, I );
%
%       The product of the objective function's Hessian (at x) times a user-supplied vector v can be obtained by the
%       command
%
%         Hxv = PROBLEM( 'fHxv', xx, v );
%
%       while the product of the (potentially I-restricted) constraints' Jacobian (at x) times v is computed by issuing
%       one of the commands
%
%         Jxv  = PROBLEM( 'cJxv' , x, v );
%         JIxv = PROBLEM( 'cIJxv', x, v, I );
%
%       and the product of the (potentially I restrticted) constraints' Jacobian (at x) transposed times v is 
%       computed by issuing one of the commands
%
%         Jtxv  = PROBLEM( 'cJtxv' , x, v );
%         JItxv = PROBLEM( 'cIJtxv', x, v, I );
%
%       Finally,  the value (and derivatives) of the Lagrangian function L(x,y) with respect to x times a user-supplied
%       vector v is obtained by one of the commands
%
%          Lxy                 = PROBLEM( 'Lxy', x, y )
%          [ Lxy, Lgxy ]       = PROBLEM( 'Lgxy', x, y )     or   [ Lxy, Lgxy ]       = PROBLEM( 'Lxy', x, y )
%          [ Lxy, Lgxy, LHxy ] = PROBLEM( 'LgHxy', x, y )    or   [ Lxy, Lgxy, LHxy ] = PROBLEM( 'Lxy', x, y )
%
%       while the product of the Lagrangian's Hessian times a vector v can be computed with the command
%
%          LHxyv = PROBLEM( 'HLxyv', x, y, v )
%
%       Finally, and as above for constraints, the Lagrangian may be I-restricted in the commands
%
%          [LIxy, LIgxy, LIHxy ] = PROBLEM( 'LIxy', x, y, I )
%          LIHxyv                = PROBLEM( 'LIHxyv' , x, y, v, I )
%
%   Because PROBLEM.m uses a persistent structure to pass problem structure across its various actions, such calls are
%   valid as long as the Matlab variable space is not cleared and as long as another OTHERPROBLEM.m file is not called.
%
%
%   TYPICAL USE OF PYTHON OUTPUT FILES
%
%   Again the use of the Python output file is very similar to that of the Matlab file. The various steps described above
%   for the the use of the Matlab output file are now (given the PROBNAME.SIF file) as follows.  The decoding of the SIF
%   file and production of the PROBNAME.py file containing the PROBNAME class is performed by the call
%
%       PROBNAME = s2mpj( 'SIFPBNAME', inpy )
%
%   where inpy is a struct describing the S2MPJ options (see next paragraph) with inpy.language = 'python' (other options
%   may be specified if desired.)  After importing the functions of PROBNAME module by
%
%       from PROBNAME import *
%
%   the setup of the problem is now performed by a call of the form
%
%       PB = PROBNAME( arg[:] )
%
%   where args[:] is a comma separated list of problem-dependent arguments. The subsequent evaluation tasks are then
%   obtained by calling one (or more) of
%
%       fx          = PB.fx( x )
%       fx, gx      = PB.fgx( x )
%       fx, gx, Hx  = PB.fgHx( x )
%       Hxv         = PB.fHxv( x, v )
%       cx          = PB.cx( x )
%       cx, Jx      = PB.cJx( x )
%       cx, Jx, HXs = PB.cJHx( x )
%       Jxv         = PB.cJxv( x, v )
%       cIx         = PB.cIx( x, I )
%       cIx, cIJx   = PB.cIJx( x, I )
%       cIx, cIJx, cIJHxs = PB.cIJHx( x, I )
%       cIHv        = PB.cIHxv( x, v, I )
%       Lxxy        = PB.Lxy( x, y )
%       Lxy, Lgxy   = PB.Lgxy( x, y )
%       Lxy, Lgxy, LgHxy = PB.LgHxy( x, y )
%       LHxyv       = PB.LHxyv( x, y, v )
%       LIxy        = PB.LIxy( x, y, I )
%       LIxy, LIgxy = PB.LIgxy( x, y, I )
%       LIxy, LIgxy, LIgHxy = PB.LIgHxy( x, y, I )
%       LIHxyv      = PB.LIHxyv( x, y, v, I )
%
%   
%   TYPICAL USE OF JULIA OUTPUT FILES
%
%   Finally, the typical use of the Julia output file is as follows. The PROBNAME.jl output file is produced from the
%   PROBNAME.SIF file by the call
%
%       PROBNAME = s2mpj( 'PROBNAME', injl )
%
%   where injl is a struct describing the S2MPJ options (see next paragraph) with injl.language = 'julia' (other options
%   may be specified if desired.) After including the PROBNAME functions in the Julia program using
%
%       include( "PROBNAME.jl" )
%
%   the setup of the problem is performed by a call of the form
%
%       pb, pbm = PROBNAME( "setup", arg[:] )
%
%   where args[:] is a comma separated list of problem-dependent arguments. The subsequent evaluation tasks are then
%   obtained by calling one (or more) of
%
%       fx          = PROBNAME( "fx", pbm, x )
%       fx, gx      = PROBNAME( "fgx", pbm, x )
%       fx, gx, Hx  = PROBNAME( "fgHx", pbm, x )
%       Hxv         = PROBNAME( "fHxv", pbm x, v )
%       cx          = PROBNAME( "cx", pbm, x )
%       cx, Jx      = PROBNAME( "cJx", pbm, x )
%       cx, Jx, HXs = PROBNAME( "cJHx", pbm,  x )
%       Jxv         = PROBNAME( "cJxv", pbm, x, v )
%       cIx         = PROBNAME( "cIx", pbm, x, I )
%       cIx, cIJx   = PROBNAME( "cIJx", pbm, x, I )
%       cIx, cIJx, cIJHxs = PROBNAME( "cIJHx", pbm,  x, I )
%       cIHv        = PROBNAME( "cIHxv", pbm, x, v, I )
%       Lxxy        = PROBNAME( "Lxy", pbm, x, y )
%       Lxy, Lgxy   = PROBNAME( "Lgxy", pbm, x, y )
%       Lxy, Lgxy, LgHxy = PROBNAME( "LgHxy", pbm, x, y )
%       LHxyv       = PROBNAME( "LHxyv", pbm, x, y, v )
%       LIxy        = PROBNAME( "LIxy", pbm, x, y, I )
%       LIxy, LIgxy = PROBNAME( "LIgxy(", pbm, x, y, I )
%       LIxy, LIgxy, LIgHxy = PROBNAME( "LIgHxy", pbm, x, y, I )
%       LIHxyv      = PROBNAME( "LIHxyv", pbm, x, y, v, I )
%
%   NOTE that pbm always occurs as the second input argument.
%   
%   THE S2MPJ OPTIONS
%
%   S2MPJ gives the user some access to internal options, typically to control file production during execution, save memory,
%   help debugging or remain as close as possible to the original Fortran SIF decoder. These options are accessible by
%   suitably setting the first (and only) variable input argument varargin{1}.
%
%   If varargin{1} is a struct then
%
%        * varargin{1}.language is is a string specifying the type of output produced by S2MPJ:
%               'matlab'    : a Matlab problem file named 'probname.m' (where probname is the name of the problem) is
%                             written in the current directory, possibly overwriting an existing file with the same name.
%               'python'    : a Python problem file named 'probname.py' (where probname is the name of the problem) is
%                             written in the current directory, possibly overwriting an existing file with the same name.
%               'julia'     : a Julia problem file named 'probname.jl' (where probname is the name of the problem) is
%                             written in the current directory, possibly overwriting an existing file with the same name.
%               'stdma'     : the content of a potentialMatlab problem file is printed on the standard output (no file is
%                             produced)
%               'stdpy'     : the content of a potential Python problem file is printed on the standard output (no file is
%                             produced)
%               'stdjl'     : the content of a potential Julia problem file is printed on the standard output (no file is
%                             produced)
%               (default: 'matlab')
%
%        * varargin{1}.showsiflines is a binary flag which is true iff the SIF data lines must be printed on the standard
%               output. This is intended for debugging.
%               (default: 0);
%
%        * varargin{1}.getxnames is a binary flag which is true iff the names of the variables must be returned in
%               pb.xnames on exit of the setup action
%               (default: 1)
%
%        * varargin{1}.getcnames is a binary flag which is true iff the names of the constraints must be returned in
%               pb.cnames on exit of the setup action
%               (default: 1)
%
%        * varargin{1}.getenames is a binary flag which is true iff the names of the nonlinear elementss must be returned
%               in pbm.elnames on exit of the setup action
%               (default: 0)
%
%        * varargin{1}.getgnames is a binary flag which is true iff the names of the groups elements must be returned in
%               pbm.grnames on exit of the setup action
%               (default: 0)
%
%        * varargin{1}.pbxscale is a binary flag which is true iff the variable's scaling for the columns of pbm.A
%               must be provided in pb.xscale, instead of being applied internally to the linear factors within
%               the problem represented by the output file.
%               NOTE: the historical Fortran decoder ignores the variable's scaling.  To preserve compatibility with
%                     problems decoded with it, this option should be activated (preventing internal scaling within s2mpj)
%                     ... and pb.xscale should then be ignored.
%               (default: 0)
%
%        * varargin{1}.keepcorder is a binary flag which is true iff, on exit of the setup action, the constraints must
%               not be reordered to appear in the order <=, followed by == and then by >=, but should instead appear in
%               the order in which they are defined in the SIF file.  (This may potentially be useful if the SIF file is
%               constructed to reflect the problem's structure, but otherwise often results in an constraint ordering
%               which is difficult to interpret)
%               (default: 0)
%
%        * varargin{1}.keepcformat is a binary flag which is true iff, on exit of the setup action, the constraints must
%               be specified using the SIF format (ie specifying ranges and types, see above) instead of specifying lower
%               and upper bounds (clower and cupper) on the constraints values. If set, the fields clower and cupper of
%               pb are replaced by the field 'ctypes' of (a cell whose i-th component contains a string ('<=', '==', or
%               '>=' ) defining the type of the i-th constraint) and, if defined, a field "ranges" defining the ranges of
%               the constraints. 
%               (default: 0)
%
%        * varargin{1}.writealtsets is a binary flag which is true iff alternative sets of constants, ranges, bounds,
%               starting points or objective function bounds must be written (as comments) in the output file.
%               (default: 0)
%
%        * varargin{1}.sifcomments is a binary flag which is true iff the comments appearing in the SIF file must be
%               repeated in the output file.
%               (default: 0)
%
%        * varargin{1}.disperrors is a binary flag which is true iff the error messages are to be displayed on the screen
%               as soon as the error is detected
%               (default: 1 )
%
%        * varargin{1}.sifdir is a string giving the path to the directory wherethe problem's SIF file must be read
%               (default: in the Matlab path)
%
%        * varargin{1}.outdir is a string giving the path to the directory where the output file (.m or .py) must be
%               written
%               (default: '.' )
%
%        * varargin{1}.redprec is a binary flag whic is true iff the problem's setup must allow for variable precision
%                arithmetic to be used in evaluations
%                (default: 1 )
%
%        * varargin{1}.dicttype is a string indicating the type of dictionary facility to be used by Matlab.
%                'native'   : use the dictionary facility included in the Matlab language,
%                'custom'   : use a dictionary based on containers maps.
%                (default: 'custom')
%
%        Not every of the above fields must be defined, each field being tested for presence and value individually.
%
%   If varargin{1} is not a struct or is not present, all options take their default values.
%              
%
%   LIMITATIONS:
%
%   In its current form, the following features are NOT supported:
%   1) the call to external FORTRAN subroutines;
%   2) the D data lines in the data GROUPS section;
%   3) the occurrence of blank characters within indeces (probably against the standard anyway);
%   5) the FREE FORM version of the SIF file;
%   4) ... and probably other unforeseen strange SIF and FORTRAN constructs.
%
%   In addition, S2MPJ *DOES NOT* provide a comprehensive explanation of possible errors in the SIF file (it is supposed to
%   adhere to the standard): it may thus crash on errors, or the produced output file may also crash. Should this happen,
%   running S2MPJ with varargin{1].showsiflines = 1, usually helps spotting the SIF error (or detecting the S2MPJ bug).
%
%   FURTHER IMPLEMENTATION NOTES:
%
%   1) For clear identification, all lines of the output problem file are written using the printmline and printcmline
%      functions.
%   2) Functionnally, remembering the names of elements and groups is not
%      necessary. It is however possible, for debugging purposes, to store them in the pbm struct by setting the flags
%      getenames and getgnames to 1 (see S2MPJ internal options above).
%   3) It is possible to disable the storing of variables' and constraints' names (saving some memory) by setting the
%      getxnames and/or getcnames flags to 0 in the internal options below. Should this be activated, the pb struct on
%      exit of the setup action has no xnames and/or cnames field(s).
%   4) S2MPJ makes the assumption that RANGES are only meaningful for inequality constraints. The SIF standard appears to
%      leave open the  possibility of ranges for other types of groups, but fails to mention what they can/could be used
%      for. Thus S2MPJ ignores the ranges of equality constraints or objective groups, should them be supplied in the SIF
%      file.
%   5) S2MPJ handles the constraints in the SIF format (in order to have coherent definition of constants for both the
%      objective and constraints, and only converts them to the S2MPJ format (if requested) at the end of the 'setup'  or
%      __init__ action (returning the correct information to the user in the pb struct).
%   6) The representation of nonlinear elements in the output file uses the following data structures:
%
%      ELEMENT TYPES
%
%      iet_      : elftype name -> elftype flat index
%      elftv{iet}: names of elemental variables for eltfype iet
%      elftp{iet}: names of elemental parameters for elftype iet    
%
%      pbm.elftype{iet}
%      - name
%      - invar : name of elftype internal variables
%
%      ELEMENTS
%
%      ie_             : element name -> element flat index
%
%      pbm.elftype{ie} : name of element type for element ie
%      pbm.elvar{ie}   : flat index of elemental variables for element ie                                   
%      pbm.elpar{ie}   : values of parameters for element ie
%                                             
%   7) Because the SIF standard typically imposes the DEFAULT type of entities such as constants, ranges, bounds or
%      starting point to be defined first and optionally, and because this default might not be zero (the value that
%      Matlab uses to fill-in unaffected vector entries), one needs to state the default value before any other value is
%      attached to an entity's component. This may then clash with an explicitly stated DEFAULT value, and a mechanism is
%      used in S2MPJ to avoid stating the SIF default for as long as one is certain that no other default is stated to
%      supersede it. This is achieved by storing the code line defining the SIF default in memory (in pending{}) and
%      only using it if the next line is not a DEFAULT line. A similar technique is also used for the DO loop statements
%      (they may be modified by a subsequent DI line, or not) and for the decoding of the Fortran expressions of the
%      nonlinear functions, because they may be continued using continuation lines (A+, I+, E+, F+, G+ and H+): the
%      definition code is stored (in codeg for gradients and codeH for Hessians) until a non-continuation line is found.
%
%   8) Because all CUTEst problems provide values, gradients and Hessians for all elements, the level of available
%      derivatives for the objective function and constraints is set to 2.  Should one consider decoding problems where
%      some values or derivatives are missing, this has to be implemented (by storing information when decoding ELEMENTS
%      and GROUPS).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   PROGRAMMING: S. Gratton (Python and Julia adaptations)
%                Ph. Toint  (Matlab code, Python and Julia adaptations),
%                started VI 2023,
                 this_version = '25 XI 2024';
%                Apologies in advance for the bugs!
%                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Set the default S2MPJ options.

pbs.fidma      = 0;  %  No Matlab file output for now
pbs.fidpy      = 0;  %  No Python file output for now
pbs.fidjl      = 0;  %  No Julia file output for now
getxnames      = 1;  %  Store the variables' names in pb.xnames      
getcnames      = 1;  %  Store the constraints' names in pb.cnames    
getenames      = 0;  %  Do not store the elements names in pbm.elnames      
getgnames      = 0;  %  Do not store the group names in pbm.grnames         
keepcorder     = 0;  %  Do not reorder the constraints as <=, ==, >=, but keep the order in which they appear in the
                     %  SIF file 
keepcformat    = 0;  %  Do not keep the SIF format for constraint specification (ie specifying constants and ranges)
                     %  instead of producing lower and upper bounds on the constraints values (clowerand cupper). 
writealtsets   = 0;  %  Do not write (as comments) active (i.e. not commented off)
                     %  alternate sets of constants, ranges, bounds, starting points or objective function bounds in the
                     %  output file
showsiflines   = 0;  %  Do not show the SIF data lines in the output file
sifcomments    = 0;  %  Do not include SIF comments in the output file
dispwarning    = 1;  %  Display a warning when a problem is renamed for Matlab/Python/julia compatibility, that is when
                     %  probname differs from sifpbname
pbs.disperrors = 1;  %  Display error messages asap
extxscale      = 1;  %  Apply the variable's scaling internally, instead of passing them the the user in pb.xscale
sifdir         = ''; %  Find the problem in the Matlab path
outdir         = ''; %  Write the output file in the current directory
pbs.dicttype   = 'custom'; % use a custom dictionary based on containers maps.

%  Determine if the Symbolic Math Toolbox is installed, in which case request support of reduced precision arithmetic
%  by default for Matlab output files (by setting redprec to 1).

v              = ver;
redprec        = any( strcmp( 'Symbolic Math Toolbox', {v.Name} ) );

%  Initialize the container for error messages.

pbs.errors     = {}; %  So far so good :-)

%  Possibly rename the problem name so that the output file can be used as a function name in Matlab or as a class name
%  in Python.

probname       = nlfname( sifpbname );
if ( ~strcmp( probname, sifpbname ) && dispwarning )
   disp( [ '*** WARNING: problem ', sifpbname, ' renamed to ', ...
            probname, ' for Matlab/Python/Julia compatibility' ] )
end

%  If S2MPJ options are specified by the user, use them instead of the defaults.

if ( nargin > 1 )

   if ( isstruct( varargin{1} ) )
      options = varargin{1};

      %  Simple options
      
      if ( isfield( options, 'getxnames' ) )
         getxnames = options.getxnames;
      end
      if ( isfield( options, 'getcnames' ) )
         getcnames = options.getcnames;
      end
      if ( isfield( options, 'getenames' ) )
         getenames = options.getenames;
      end
      if ( isfield( options, 'getgnames' ) )
         getgnames = options.getgnames;
      end
      if ( isfield( options, 'keepcorder' ) )
         keepcorder = options.keepcorder;
      end
      if ( isfield( options, 'keepcformat' ) )
         keepcformat = options.keepcformat;
      end
      if ( isfield( options, 'writealtsets' ) )
         writealtsets = options.writealtsets;
      end
      if ( isfield( options, 'sifcomments' ) )
         sifcomments = options.sifcomments;
      end
      if ( isfield( options, 'disperrors' ) )
         pbs.disperrors = options.disperrors;
      end
      if ( isfield( options, 'sifdir' ) )
         sifdir = options.sifdir;
      end
      if ( isfield( options, 'outdir' ) )
         outdir = options.outdir;
      end
      if ( isfield( options, 'pbxscale' ) )
         extxscale = 1;
      end
      if ( isfield( options, 'redprec' ) )
         redprec = options.redprec;
      end
      if ( isfield( options, 'dicttype' ) )
         pbs.dicttype = options.dicttype;
      end
      
      % Prepare for the desired output language.
      
      if ( isfield( options, 'language' ) )
         switch ( options.language )
         case 'matlab'
            if ( isempty( outdir ) )
                mafile = [ './', probname, '.m' ];    %  the name of the Matlab file
            else
                mafile = [ outdir, '/', probname, '.m' ];
            end
            pbs.fidma  = fopen( mafile, 'w' );
            pbs.lang   = 'matlab';                    %  Matlab code will be generated
         case 'python'
            if ( isempty( outdir ) )
               pyfile = [ './', probname, '.py' ];    %  the name of the Python file
            else
               pyfile = [ outdir, '/', probname, '.py' ];
            end
            pbs.fidpy = fopen( pyfile, 'w' );            
            pbs.lang  = 'python';                     %  Python code will be generated
         case 'julia'
            if ( isempty( outdir ) )
               jlfile = [ './', probname, '.jl' ];    %  the name of the Julia file
            else
               jlfile = [ outdir, '/', probname, '.jl' ];
            end
            pbs.fidjl = fopen( jlfile, 'w' );
            pbs.lang  = 'julia';                      %  Julia code will be generated
         case 'stdma'
            pbs.fidma = 1;
            pbs.lang  = 'matlab';                     %  Matlab code will be printed on STD
         case 'stdpy'
            pbs.fidpy = 1;
            pbs.lang  = 'python';                     %  Python code will be printed on STD
         case 'stdjl'
            pbs.fidjl = 1;
            pbs.lang  = 'julia';                      %  Julia code will be printed on STD
         end
      end

      %  See if SIF lines must be shown.

      if ( isfield( options, 'showsiflines' ) )
          showsiflines = options.showsiflines;
      end

   else
      pbs.errors{end+1} = 'ERROR: second argument of S2MPJ is not a struct (options)';
      disp( pbs.errors{end} );
      exitc = length( pbs.errors );
      return
   end
end

%  If no output stream was specified, output a Matlab file.

if ( pbs.fidma == 0 && pbs.fidpy == 0 && pbs.fidjl == 0 )
   mafile    = [ './', probname, '.m' ];      %  the name of the Matlab file
   pbs.fidma = fopen( mafile, 'w' );
   pbs.lang  = 'matlab';                      %  Matlab code will be generated
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                      The SIF file 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Open the problem SIF file.

if ( contains( sifpbname, '.SIF' ) )
   fname = sifpbname;
else
   fname = [ sifpbname, '.SIF' ];
end
if ( ~isempty( sifdir ) )
   fname = [ sifdir, '/', fname ];
end
if ( ~exist( fname, 'file' ) )
   pbs.errors{end+1} = sprintf( 'ERROR: can''t find the %s file!', fname );
   if ( pbs.disperrors )
      disp( pbs.errors{end} )
   end
   errors = pbs.errors;
   exitc  = length( errors );
   return
end
fidSIF = fopen( fname, 'r');
if ( fidSIF < 0 )
   pbs.errors{end+1} = sprintf( 'ERROR: can''t open %s', fname );
   if ( pbs.disperrors )
      disp( pbs.errors{end} )
   end
   errors = pbs.errors;
   exitc  = length( errors );
   return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                    The output file options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Define the default indentation space and the indentation level to be used in the output file.

bindent = '    ';
indlvl  = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                        Initializations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Prepare the internal loop structure.

pbs.loop    = struct([]);  % the array of struct whose elements describe the successive loops
pbs.actloop = [];          % the positions (in pbs.loop) of the currently active loops

%  Prepare a few holders.

elftype = {}; % the internal list of element types
grftype = {}; % the internal list of group types

%  Names for sets of constants, ranges, bounds on the variables, start points and objective bounds as specified in field 2
%  of the relevant SIF lines.  These listsare used to detect when a change in the current set occurs (allowing statements
%  defining sets beyond the first to be written in the output file as comments).

cname = '';
rname = '';
bname = '';
sname = '';
oname = '';

%  Positioning indicators

section         =''; % the name of the current data section
incomments      = 0; % true while in the initial comments section
indata          = 1; % true while in the data section
inobounds       = 0; % rtue while in the objective bounds section
ingroups        = 0; % true while in the nonlinear elements section
inelements      = 0; % true while in the nonlinear groups section
pbs.nline       = 0; % the current line number
varsdef         = 0; % true if the variables have been defined
grpsdef         = 0; % true if the groups (obj, lt, eq and gt) have been defined
arguments       = 0; % true if a comment indicating alternative aruments has been written

%  Logical switches

has_xscale      = 0; % true if there is a variable scaling factor different from 1  
has_bounds      = 0; % true if the problem has an explicit BOUNDS section
has_xlowdef     = 0; % true if default values for lower bounds have been explicitly defined
has_xuppdef     = 0; % true if default values for upper bounds have been explicitly defined
has_intvars     = 0; % true if the problem has integer variables
has_binvars     = 0; % true if the problem has binary variables
has_x0def       = 0; % true if the default for x0 has been explicitly specified
has_y0def       = 0; % true if the default for y0 has been explicitly specified
has_constraints = 0; % true if the problem has constraints other than bounds
has_leconstr    = 0; % true if the problem has <= constraints
has_eqconstr    = 0; % true if the problem has == constraints
has_geconstr    = 0; % true if the problem has >= constraints
has_nonlinc     = 0; % true if the problem has nonlinear constraints
has_constants   = 0; % true if group constants are explicitly specified
has_ranges      = 0; % true if group ranges are explicitly specified
has_start       = 0; % true if starting point(s) are explicitly specified
has_elpar       = 0; % true if at least one element type has a parameter
has_grpar       = 0; % true if at least one group type has a parameter
has_A           = 0; % true if the matrix A of linear terms for groups is not empty
has_H           = 0; % true if a quadratic term is defined (in the QUADRATICS section)
has_ngrp        = 0; % true if ngrp has been defined.  The idea is to define ngrp just before
                     % it is actually used (for instance in defining defaults for constants or ranges) and at the latest
                     % in the GLOBAL DIMENSIONS section of the output file
prevlineispass  = 0; % true if the previous line written in the Python output file is 'pass', which allows avoiding
                     % multiple "pass" lines

%  Prepare the dictionary of integer and real  parameters ( unquoted names --> quoted names ), including standard integer
%  parameters but also indeces of no longer active loops indeces, providing a way to remember their last value. This
%  dictionary thus provides a mechanism to extract the current (during decoding) (string) valueof an index, irrespective
%  of it being a past loop index or a standard integer parameter.

switch ( pbs.dicttype )
case 'native'
   pbs.irpdict = configureDictionary( 'string', 'string' );
case 'custom'
   pbs.irpdict = containers.Map( 'KeyType', 'char', 'ValueType', 'char' );
end

%  A translation table (dictionary) is necessary to define the correspondance between the internal names (in pbm) of
%  parameters and the names used in Fortran expressions in the nonlinear ELEMENTS and GROUPS sections. The global ELEMENTS
%  and GROUPS parameters are stored in pbm.efpar and pbm.gfpar, respectively.

switch( pbs.dicttype )
case 'native'
   globdict = configureDictionary( 'string', 'string' );  % for elements
case 'custom'
   globdict = containers.Map( 'KeyType', 'char', 'ValueType', 'char' );  % for elements
end
n_eglobs = 0;     % the number of global Fortran parameters for nonlinear elements
n_gglobs = 0;     % the number of global Fortran parameters for nonlinear elements

%  Some further initializations

pending     = {}; %  a temporary holder for a Matlab/python command defining a loop (see DO LOOPS section), or for a
                  % Fortran expression (used to cope with A+, I+ and E+ lines in the ELEMENT or GROUPS sections)
pendingkey  = ''; %  a string identifying the type of pending statements
readsifpars = {}; %  the list of already read $-PARAMETERS of the SIF file
globdims    = 1;  %  the global problem's dimension must still  be written in the output file
xlowdef     = []; %  the value of an explicit default lower bound on the variables when > -Inf
xuppdef     = []; %  the value of an explicit default upper bound on the variables when < +Inf
x0def       = []; %  the value of an explicit nonzero default starting point for the variables

%  Equivalent names for some data sections (for compatibility with MPS and other standards)

EQVARIABLES = { 'VARIABLES', 'COLUMNS' };
EQGROUPS    = { 'GROUPS', 'ROWS', 'CONSTRAINTS' };
EQCONSTANTS = { 'CONSTANTS', 'RHS', 'RHS''' };
EQQUADRATIC = { 'QUADRATIC', 'QMATRIX', 'QUADS', 'QUADOBJ', 'QSECTION', 'QMATRIX', 'HESSIAN' };

%  The list of section headers potentially following the VARIABLES and GROUPS sections. When a section in this list is met,
%  the global dimensions of the problem are computed and written before the new section is interpreted.

GLOBDIMTIME = { 'RANGES', 'START POINT', 'BOUNDS', 'OBJECT BOUND', 'ELEMENT TYPE', 'ELEMENT USES', 'GROUP TYPE', ...
                'GROUP USES', 'ENDATA' };
GLOBDIMTIME = union( GLOBDIMTIME, EQCONSTANTS );
GLOBDIMTIME = union( GLOBDIMTIME, EQQUADRATIC );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Loop over the lines of the problem file  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Open the problem's SIF file and start reading.

while ( ~feof( fidSIF ) )  %  Within the SIF file

   %  Read a data line.
   
   line      = fgetl( fidSIF );
   lline     = length( line );
   pbs.nline = pbs.nline + 1;

   if ( showsiflines )
      printmline( sprintf( ' line %d:  %s\n', pbs.nline, line ),                                   0, bindent, pbs.fidma );
      printpline( sprintf( ' line %d:  %s\n', pbs.nline, line ),                                   0, bindent, pbs.fidpy );
      printjline( sprintf( ' line %d:  %s\n', pbs.nline, line ),                                   0, bindent, pbs.fidjl );
   end
   
   %  Blank line: do nothing 

   if ( isempty( lline ) || isempty( strtrim( line ) ) )
   
      if ( incomments )
         printcmline( ' ',                                                                            bindent, pbs.fidma );
         printcpline( ' ',                                                                            bindent, pbs.fidpy );
         printcjline( ' ',                                                                            bindent, pbs.fidjl );
      end
      continue  % next SIF data line
      
   %  Initial comment lines: find the problem name and its classification
   %  (as well as potential other comments).
   
   elseif ( line(1) == '*'  )

      if ( lline > 1 )
         posc = strfind( line, 'classification' );
         if ( ~isempty( posc ) )
             switch ( pbs.lang )
             case 'matlab'
                classification = [ '''C-C',   strtrim( line( posc+15:end ) ), '''' ];
                printcmline( sprintf( '    classification = %s',classification ),                     bindent, pbs.fidma );
             case { 'python', 'julia' }
                classification = [ '"C-C',    strtrim( line( posc+15:end ) ), '"' ];
                printcpline( sprintf( '    classification = %s',classification ),                     bindent, pbs.fidpy );
                printcjline( sprintf( '    classification = %s',classification ),                     bindent, pbs.fidjl );
             end 
         elseif ( contains( line, 'Problem :' ) )
             printcmline( sprintf( '    Problem : %s', probname ),                                    bindent, pbs.fidma );
             printcpline( sprintf( '    Problem : %s', probname ),                                    bindent, pbs.fidpy );
             printcjline( sprintf( '    Problem : %s', probname ),                                    bindent, pbs.fidjl );
         elseif ( contains( line, '$-PARAMETER' ) )
             if ( arguments == 0 )
                printcmline( '       Alternative values for the SIF file parameters:' ,               bindent, pbs.fidma );
                printcpline( '           Alternative values for the SIF file parameters:' ,           bindent, pbs.fidpy );
                printcjline( '       Alternative values for the SIF file parameters:' ,               bindent, pbs.fidjl );
                arguments = 1;
             end
             printcmline( sprintf( ' %s', line( 2:end ) ),                                            bindent, pbs.fidma );
             printcpline( sprintf( ' %s', line( 2:end ) ),                                            bindent, pbs.fidpy );
             printcjline( sprintf( ' %s', line( 2:end ) ),                                            bindent, pbs.fidjl );
         elseif ( incomments || sifcomments || inobounds )
             printcmline( sprintf( ' %s', line( 2:end ) ),                                            bindent, pbs.fidma );
             printcpline( sprintf( ' %s', line( 2:end ) ),                                            bindent, pbs.fidpy );
             printcjline( sprintf( ' %s', line( 2:end ) ),                                            bindent, pbs.fidjl );
         end
      end
      continue  % next SIF data line
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %                               Write the pending statements in the output file
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %
   %  This section of the code checks whether or not some code is pending (ie waiting to be written in the output file if
   %  necessary).  This occurs when pending code is not 'final' at the moment of decoding, in that it can be modified or
   %  made irrelevant by later lines.  It is written in the output file here if still relevant or in its possibly modified
   %  form.
   
   %  1) The rest of the output file header

   if ( strcmp( pendingkey, 'header' ) )

      %  Keeping the decoded information available for later calls of the output file. This is achieved in Matlab by making
      %  the struct pbm persistent and later storing in pbm the data required  for actions different from 'setup' and
      %  passing information to s2mpjlib.

      switch ( pbs.lang )
      
      case 'matlab'
      
         printcmline( ' ',                                                                            bindent, pbs.fidma );
         printcmline( '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',   bindent, pbs.fidma );
         printcmline( ['   Translated to Matlab by S2MPJ version ', this_version ],                   bindent, pbs.fidma );
         printcmline( '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',   bindent, pbs.fidma );
         printmline( ' ',                                                                     0,      bindent, pbs.fidma );
         printmline( 'persistent pbm;',                                                       0,      bindent, pbs.fidma );
         printmline( ' ',                                                                     0,      bindent, pbs.fidma );
         printmline( sprintf( 'name = ''%s'';', probname ),                                   0,      bindent, pbs.fidma );

         %   Start the SETUP section of the output file.

         printmline( ' ',                                                                     0,      bindent, pbs.fidma );
         printmline( 'switch(action)', 0, bindent, pbs.fidma );   %  the 'action' switch
         printmline( ' ',                                                                     0,      bindent, pbs.fidma );
         if ( redprec )
            printmline( 'case {''setup'',''setup_redprec''}',                                 1,      bindent, pbs.fidma ); 
         else
            printmline( 'case ''setup''',                                                     1,      bindent, pbs.fidma );
         end         
         printmline( ' ',                                                                     0,      bindent, pbs.fidma );
         if ( redprec )
            printmline( 'if(strcmp(action,''setup_redprec''))',                               2,      bindent, pbs.fidma );
            printmline( 'if(isfield(pbm,''ndigs''))',                                         3,      bindent, pbs.fidma );
            printmline( 'rmfield(pbm,''ndigs'');',                                            4,      bindent, pbs.fidma );
            printmline( 'end',                                                                3,      bindent, pbs.fidma );
            printmline( 'pbm.ndigs = max(1,min(15,varargin{end}));',                          3     , bindent, pbs.fidma );
            printmline( 'nargs     = nargin-2;',                                              3,      bindent, pbs.fidma );
            printmline( 'else',                                                               2,      bindent, pbs.fidma );
            printmline( 'nargs = nargin-1;',                                                  3,      bindent, pbs.fidma );
            printmline( 'end',                                                                2,      bindent, pbs.fidma );
         end
         printmline( sprintf( 'pb.name   = name;' ),                                          2,      bindent, pbs.fidma );
         printmline( sprintf( 'pbm.name  = name;' ),                                          2,      bindent, pbs.fidma );
         if ( ~strcmp( sifpbname, probname ) )
            printmline( sprintf( 'pb.sifpbname  = ''%s'';', sifpbname ),                      2,      bindent, pbs.fidma );
         end
         printmline( '%%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%',                   indlvl, bindent, pbs.fidma );

         %  Initialize the dictionaries for parameters' values (v_), variables' 'flat'
         %  indeces and groups' 'flat' indeces. 

         switch( pbs.dicttype )
         case 'native'
            printmline( 'v_  = configureDictionary(''string'',''double'');',                  indlvl, bindent, pbs.fidma );
            printmline( 'ix_ = configureDictionary(''string'',''double'');',                  indlvl, bindent, pbs.fidma );
            printmline( 'ig_ = configureDictionary(''string'',''double'');',                  indlvl, bindent, pbs.fidma );
         case 'custom'
            printmline( 'v_  = containers.Map(''KeyType'',''char'', ''ValueType'', ''double'');', ...
                                                                                              indlvl, bindent, pbs.fidma );
            printmline( 'ix_ = containers.Map(''KeyType'',''char'', ''ValueType'', ''double'');', ...
                                                                                              indlvl, bindent, pbs.fidma );
            printmline( 'ig_ = containers.Map(''KeyType'',''char'', ''ValueType'', ''double'');', ...
                                                                                              indlvl, bindent, pbs.fidma );
         end
      case 'python'
      
         printcpline( ' ',                                                                            bindent, pbs.fidpy );
         printcpline( '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',   bindent, pbs.fidpy );
         printcpline( ['   Translated to Python by S2MPJ version ', this_version ],                   bindent, pbs.fidpy );
         printcpline( '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',   bindent, pbs.fidpy );
         printpline( ' ',                                                                     0,      bindent, pbs.fidpy );
         printpline( sprintf( 'name = ''%s''', probname ),                                    1,      bindent, pbs.fidpy );

         %   Start the SETUP section of the output file.

         printpline( ' ',                                                                     0,      bindent, pbs.fidpy );
         printpline( 'def __init__(self, *args): ',                                           1,      bindent, pbs.fidpy );
         printpline( 'import numpy as np',                                                    2,      bindent, pbs.fidpy );
         printpline( 'from scipy.sparse import csr_matrix',                                   2,      bindent, pbs.fidpy );
         if ( ~strcmp( sifpbname, probname ) )
            printpline( sprintf( 'self.sifpbname = ''%s''', sifpbname ),                      2,      bindent, pbs.fidpy );
         end
         printpline( 'nargin   = len(args)',                                                  2,      bindent, pbs.fidpy );
         printpline( ' ',                                                                     0,      bindent, pbs.fidpy );
         printpline( '#%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%',                   indlvl, bindent, pbs.fidpy );

         %  Initialize the dictionaries for parameters' values (v_), variables' 'flat'
         %  indeces and groups' 'flat' indeces. 

         printpline( 'v_  = {}',                                                              indlvl, bindent, pbs.fidpy );
         printpline( 'ix_ = {}',                                                              indlvl, bindent, pbs.fidpy );
         printpline( 'ig_ = {}',                                                              indlvl, bindent, pbs.fidpy );

      case 'julia'
      
         printcjline( ' ',                                                                            bindent, pbs.fidjl );
         printcjline( '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',   bindent, pbs.fidjl );
         printcjline( ['   Translated to Julia by S2MPJ version ', this_version ],                    bindent, pbs.fidjl );
         printcjline( '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',   bindent, pbs.fidjl );
         printjline( ' ',                                                                     0,      bindent, pbs.fidjl );
         printjline( sprintf( 'name = "%s"', probname ),                                      1,      bindent, pbs.fidjl );

         %   Start the SETUP section of the output file.

         printjline( ' ',                                                                     0,      bindent, pbs.fidjl );
         printjline( 'if action == "setup"',                                                  1,      bindent, pbs.fidjl );
         printjline( 'pb           = PB(name)',                                               2,      bindent, pbs.fidjl );
         printjline( 'pbm          = PBM(name)',                                              2,      bindent, pbs.fidjl );
         if ( ~strcmp( sifpbname, probname ) )
            printjline( sprintf( 'pb.sifpbname = "%s"', sifpbname ),                          2,      bindent, pbs.fidjl );
         end
         printjline( 'nargin       = length(args)',                                           2,      bindent, pbs.fidjl );
         printjline( 'pbm.call     = getfield( Main, Symbol( name ) )',                       2,      bindent, pbs.fidjl );
         printjline( ' ',                                                                     0,      bindent, pbs.fidjl );
         printjline( '#%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%',                   indlvl, bindent, pbs.fidjl );

         %  Initialize the dictionaries for parameters' values (v_), variables' 'flat'
         %  indeces and groups' 'flat' indeces. 

         printjline( 'v_  = Dict{String,Float64}();',                                         indlvl, bindent, pbs.fidjl );
         printjline( 'ix_ = Dict{String,Int}();',                                             indlvl, bindent, pbs.fidjl );
         printjline( 'ig_ = Dict{String,Int}();',                                             indlvl, bindent, pbs.fidjl );

      end

      %  No pending code yet.
      
      pendingkey = '';

   end

   % 2) The DO loop statements that are not modified by a subsequent DI specification
   
   if ( strcmp( pendingkey, 'loop') && ~strcmp( line(2:3), 'DI' ) )
      printmline( pending{1},                                                                 indlvl, bindent, pbs.fidma );
      printpline( pending{1},                                                                 indlvl, bindent, pbs.fidpy );
      printjline( pending{1},                                                                 indlvl, bindent, pbs.fidjl );
      pending        = {};
      pendingkey     = '';
      indlvl         = indlvl + 1;
      prevlineispass = 0;
   end

   % 3) The global problem's dimensions, once the VARIABLES and GROUPS sections have been decoded.
   %    Note that the number of variables may still be increased when 'nonlinear variables" are found in the ELEMENT USES
   %    section.

   if ( globdims && ismember( strtrim( line ), GLOBDIMTIME  ) )

      switch( pbs.lang )
      case 'matlab'
         printmline( '%%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%',                   indlvl, bindent, pbs.fidma );
         switch ( pbs.dicttype )
         case 'native'
            printmline( 'pb.n    = numEntries(ix_);',                                         indlvl, bindent, pbs.fidma );
            if ( ~has_ngrp )
               printmline( 'ngrp = numEntries(ig_);',                                         indlvl, bindent, pbs.fidma );
               has_ngrp = 1;
            end
         case 'custom'
            printmline( 'pb.n   = ix_.Count;',                                                indlvl, bindent, pbs.fidma );
            if ( ~has_ngrp )
               printmline( 'ngrp   = ig_.Count;',                                             indlvl, bindent, pbs.fidma );
               has_ngrp = 1;
            end
         end
         if ( has_constraints )
            printmline( 'legrps = find(strcmp(gtype,''<=''));',                               indlvl, bindent, pbs.fidma );
            printmline( 'eqgrps = find(strcmp(gtype,''==''));',                               indlvl, bindent, pbs.fidma );
            printmline( 'gegrps = find(strcmp(gtype,''>=''));',                               indlvl, bindent, pbs.fidma );
            printmline( 'pb.nle = length(legrps);',                                           indlvl, bindent, pbs.fidma );
            printmline( 'pb.neq = length(eqgrps);',                                           indlvl, bindent, pbs.fidma );
            printmline( 'pb.nge = length(gegrps);',                                           indlvl, bindent, pbs.fidma );
            printmline( 'pb.m   = pb.nle+pb.neq+pb.nge;',                                     indlvl, bindent, pbs.fidma );
            if ( keepcorder )
               printmline( 'pbm.congrps = find(ismember(gtype,{''<='',''=='',''>=''}));',     indlvl, bindent, pbs.fidma );
            else
               printmline( 'pbm.congrps = [ legrps, eqgrps, gegrps ];',                       indlvl, bindent, pbs.fidma );
            end
            if ( getcnames )
               printmline( '[pb.cnames{1:pb.m}] = deal(cnames{pbm.congrps});',                indlvl, bindent, pbs.fidma );
            end
            printmline( 'pb.nob = ngrp-pb.m;',                                                indlvl, bindent, pbs.fidma );
            printmline( 'pbm.objgrps = find(strcmp(gtype,''<>''));',                          indlvl, bindent, pbs.fidma );
         else
            printmline( 'pbm.objgrps = [1:ngrp];',                                            indlvl, bindent, pbs.fidma );
            printmline( 'pb.m        = 0;',                                                   indlvl, bindent, pbs.fidma );
          end
      case 'python'
         printpline( '#%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%',                   indlvl, bindent, pbs.fidpy );
         printpline( 'self.n   = len(ix_)',                                                   indlvl, bindent, pbs.fidpy );
         if ( ~has_ngrp )
            printpline( 'ngrp   = len(ig_)',                                                  indlvl, bindent, pbs.fidpy );
            has_ngrp = 1;
         end
         if ( has_constraints )
            printpline( 'legrps = np.where(gtype==''<='')[0]',                                indlvl, bindent, pbs.fidpy );
            printpline( 'eqgrps = np.where(gtype==''=='')[0]',                                indlvl, bindent, pbs.fidpy );
            printpline( 'gegrps = np.where(gtype==''>='')[0]',                                indlvl, bindent, pbs.fidpy );
            printpline( 'self.nle = len(legrps)',                                             indlvl, bindent, pbs.fidpy );
            printpline( 'self.neq = len(eqgrps)',                                             indlvl, bindent, pbs.fidpy );
            printpline( 'self.nge = len(gegrps)',                                             indlvl, bindent, pbs.fidpy );
            printpline( 'self.m   = self.nle+self.neq+self.nge',                              indlvl, bindent, pbs.fidpy );
            if ( keepcorder )
               printpline( 'self.congrps = np.where(gtype==''<='' or gtype==''=='' or gytype==''>='')[0])', ...
                                                                                              indlvl, bindent, pbs.fidpy );
            else
               printpline( 'self.congrps = np.concatenate((legrps,eqgrps,gegrps))',           indlvl, bindent, pbs.fidpy );
            end
            if ( getcnames )
               printpline( 'self.cnames = cnames[self.congrps]',                               indlvl, bindent, pbs.fidpy );
            end
            printpline( 'self.nob = ngrp-self.m',                                             indlvl, bindent, pbs.fidpy );
            printpline( 'self.objgrps = np.where(gtype==''<>'')[0]',                          indlvl, bindent, pbs.fidpy );
         else
            printpline( 'self.objgrps = np.arange(ngrp)',                                     indlvl, bindent, pbs.fidpy );
            printpline( 'self.m       = 0',                                                   indlvl, bindent, pbs.fidpy );
         end
      case 'julia'
         printjline( '#%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%',                   indlvl, bindent, pbs.fidjl );
         printjline( 'pb.n   = length(ix_)',                                                  indlvl, bindent, pbs.fidjl );
         if ( ~has_ngrp )
            printjline( 'ngrp   = length(ig_)',                                               indlvl, bindent, pbs.fidjl );
            has_ngrp = 1;
         end
         if ( has_constraints )
            printjline( 'legrps = findall(x->x=="<=",gtype)',                                 indlvl, bindent, pbs.fidjl );
            printjline( 'eqgrps = findall(x->x=="==",gtype)',                                 indlvl, bindent, pbs.fidjl );
            printjline( 'gegrps = findall(x->x==">=",gtype)',                                 indlvl, bindent, pbs.fidjl );
            printjline( 'pb.nle = length(legrps)',                                            indlvl, bindent, pbs.fidjl );
            printjline( 'pb.neq = length(eqgrps)',                                            indlvl, bindent, pbs.fidjl );
            printjline( 'pb.nge = length(gegrps)',                                            indlvl, bindent, pbs.fidjl );
            printjline( 'pb.m   = pb.nle+pb.neq+pb.nge',                                      indlvl, bindent, pbs.fidjl );
            if ( keepcorder )
               printjline( 'pbm.congrps = findall(x->x!="<>",gtype)',                         indlvl, bindent, pbs.fidjl );
            else
               printjline( 'pbm.congrps = [[legrps;eqgrps];gegrps]',                          indlvl, bindent, pbs.fidjl );
            end
            printjline( 'pb.nob = ngrp-pb.m',                                                 indlvl, bindent, pbs.fidjl );
            printjline( 'pbm.objgrps = findall(x->x=="<>",gtype)',                            indlvl, bindent, pbs.fidjl );
         else
            printjline( 'pbm.objgrps = collect(1:ngrp)',                                      indlvl, bindent, pbs.fidjl );
            printjline( 'pb.m        = 0',                                                    indlvl, bindent, pbs.fidjl );
         end
      
      end
      globdims = 0;
   end

   % 4) The implicit DEFAULT specifications that are not superseded by an explicit one,

   if ( strcmp( pendingkey, 'default' ) && ( length( line ) < 23 || ~strcmp( line(15:23), '''DEFAULT''' ) ) )
      if ( contains( pending{1}, 'BOUNDS' ) )
         if ( ~has_xlowdef && ~has_xuppdef )  % there is no BOUNDS section in the SIF file
            printmline( pending{1},                                                           indlvl, bindent, pbs.fidma );
            printpline( pending{1},                                                           indlvl, bindent, pbs.fidpy );
            printjline( pending{1},                                                           indlvl, bindent, pbs.fidjl );
         end
         if ( ~has_xlowdef )
            printmline( 'pb.xlower = zeros(pb.n,1);',                                         indlvl, bindent, pbs.fidma );
            printpline( 'self.xlower = np.zeros((self.n,1))',                                 indlvl, bindent, pbs.fidpy );
            printjline( 'pb.xlower = zeros(Float64,pb.n)',                                    indlvl, bindent, pbs.fidjl );
            has_xlowdef = 1;
         end
         if ( ~has_xuppdef ) 
            printmline( 'pb.xupper = Inf*ones(pb.n,1);',                                      indlvl, bindent, pbs.fidma );
            printpline( 'self.xupper = np.full((self.n,1),float(''inf''))',                   indlvl, bindent, pbs.fidpy );
            printjline( 'pb.xupper =    fill(Inf,pb.n)',                                      indlvl, bindent, pbs.fidjl );
            has_xuppdef = 1;
         end
      elseif (contains( pending{1}, 'START POINT' ) )
         if ( ~has_x0def || ( has_constraints && ~has_y0def ) )
            printmline( pending{1},                                                           indlvl, bindent, pbs.fidma );
            printpline( pending{1},                                                           indlvl, bindent, pbs.fidpy );
            printjline( pending{1},                                                           indlvl, bindent, pbs.fidjl );
         end
         if ( ~has_x0def )
            printmline( 'pb.x0(1:pb.n,1) = zeros(pb.n,1);',                                   indlvl, bindent, pbs.fidma );
            printpline( 'self.x0 = np.zeros((self.n,1))',                                     indlvl, bindent, pbs.fidpy );
            printjline( 'pb.x0 = zeros(Float64,pb.n)',                                        indlvl, bindent, pbs.fidjl );
            has_x0def = 1;
         end
         if ( has_constraints && ~has_y0def )
            printmline( 'pb.y0 = zeros(pb.m,1);',                                             indlvl, bindent, pbs.fidma );
            printpline( 'self.y0 = np.zeros((self.m,1))',                                     indlvl, bindent, pbs.fidpy );
            printjline( 'pb.y0 = zeros(Float64,pb.m)',                                        indlvl, bindent, pbs.fidjl );
            has_y0def = 1;
         end
      else
         for ii = 1:length( pending )
            printmline( pending{ii},                                                          indlvl, bindent, pbs.fidma );
            printpline( pending{ii},                                                          indlvl, bindent, pbs.fidpy );
            printjline( pending{ii},                                                          indlvl, bindent, pbs.fidjl );
         end       
      end
      pending    = {};
      pendingkey = '';
   end
   
   %  5) The Matlab/Python versions of the Fortran assignment in the ELEMENTS and GROUPS sections

   if ( strcmp( pendingkey, 'fortran' ) && ~ismember( line(2:3), {'A+', 'I+', 'E+', 'F+' } ) )
      for ii = 1:length( pending )
         printmline( pending{ii},                                                             indlvl, bindent, pbs.fidma );
         printpline( pending{ii},                                                             indlvl, bindent, pbs.fidpy );
         printjline( pending{ii},                                                             indlvl, bindent, pbs.fidjl );
      end       
      pending    = {};
      pendingkey = '';
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Headers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   %  Problem name

   if ( lline >= 4 && strcmp( line(1:4), 'NAME' ) )
   
%      probname = strtrim( line(15:end) );

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %                                  The function calling sequence
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      switch ( pbs.lang )
      case 'matlab'
         printmline( sprintf( 'function varargout = %s(action,varargin)', probname ),              0, bindent, pbs.fidma );
         printmline( ' ',                                                                          0, bindent, pbs.fidma );
         printcmline( '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',   bindent, pbs.fidma );
         printcmline( ' ',                                                                            bindent, pbs.fidma );
      case 'python'
         printpline( 'from s2mpjlib import *',                                                     0, bindent, pbs.fidpy );
         printpline( sprintf( 'class  %s(CUTEst_problem):', probname ),                            0, bindent, pbs.fidpy );
         printpline( ' ',                                                                          0, bindent, pbs.fidpy );
         printcpline( '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',   bindent, pbs.fidpy );
         printcpline( ' ', bindent, pbs.fidpy );
      case 'julia'
         printjline( sprintf( 'function %s(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)', ...
                              probname ),                                                          0, bindent, pbs.fidjl );
         printcjline( ' ',                                                                            bindent, pbs.fidjl );
         printcjline( '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',    bindent, pbs.fidjl );
         printcjline( ' ',                                                                            bindent, pbs.fidjl );
      
      end
      pendingkey = 'header';    %  remember to output the rest of the output file header

      incomments = 1;

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %                          The section headers and section specific initializations
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   %  Data groups

   elseif( ismember( strtrim( line ), EQGROUPS ) )

      %  Start the GROUPS section.

      section    = 'DGROUPS';
      incomments = 0;
      grpsdef    = 1;
      switch( pbs.lang )
      case 'matlab'
         printmline( '%%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%',                  indlvl, bindent, pbs.fidma );
      case 'python'
         printpline( '#%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%',                  indlvl, bindent, pbs.fidpy );
         printpline( 'self.gscale  = np.array([])',                                           indlvl, bindent, pbs.fidpy );
         printpline( 'self.grnames = np.array([])',                                           indlvl, bindent, pbs.fidpy );
         printpline( 'cnames       = np.array([])',                                           indlvl, bindent, pbs.fidpy );
         printpline( 'self.cnames  = np.array([])',                                           indlvl, bindent, pbs.fidpy );
         printpline( 'gtype        = np.array([])',                                           indlvl, bindent, pbs.fidpy );
      case 'julia'
         printjline( '#%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%',                  indlvl, bindent, pbs.fidjl );
         printjline( 'gtype = String[]',                                                      indlvl, bindent, pbs.fidjl );
      end

      %  Initialize the vectors irA, icA and valA, if not already done.
      
      if ( ~varsdef )
         initA( indlvl, bindent, pbs );
      end   
      prevgname = '';
      
      consineq = {}; % the names of the inequality constraints (for which a range is meaningful)
      
   %  Variables

   elseif( ismember( strtrim( line ), EQVARIABLES ) )
           
      %  Start the VARIABLES section.

      section    = 'VARIABLES';
      incomments = 0;
      varsdef    = 1;
      incomments = 0;
      printmline( '%%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%',                     indlvl, bindent, pbs.fidma );
      printpline( '#%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%',                     indlvl, bindent, pbs.fidpy );
      printjline( '#%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%',                     indlvl, bindent, pbs.fidjl );
      if ( getxnames )
      
         % a holder for the variable's names
         
         printmline( 'pb.xnames = {};',                                                       indlvl, bindent, pbs.fidma );
         printpline( 'self.xnames = np.array([])',                                            indlvl, bindent, pbs.fidpy );
      end

      % a holder for the variable's scalings
      
      printpline( 'self.xscale = np.array([])',                                               indlvl, bindent, pbs.fidpy );
      printjline( 'pb.xscale = Float64[]',                                                    indlvl, bindent, pbs.fidjl );
      
       % holders for the indeces of the integer and binary variables
       
      printpline( 'intvars   = np.array([])',                                                 indlvl, bindent, pbs.fidpy );
      printpline( 'binvars   = np.array([])',                                                 indlvl, bindent, pbs.fidpy );

      printjline( 'intvars = Int64[]',                                                        indlvl, bindent, pbs.fidjl );
      printjline( 'binvars = Int64[]',                                                        indlvl, bindent, pbs.fidjl );
      if ( grpsdef )
         switch ( pbs.dicttype )
         case 'native'
            printmline( 'ngrp   = numEntries(ig_);',                                          indlvl, bindent, pbs.fidma );
         case 'custom'
            printmline( 'ngrp   = ig_.Count;',                                                indlvl, bindent, pbs.fidma );
         end
         printpline( 'ngrp   = len(ig_)',                                                     indlvl, bindent, pbs.fidpy );
         printjline( 'ngrp   = length(ig_)',                                                  indlvl, bindent, pbs.fidjl );
         has_ngrp = 1;
      end

      %  Initialize the vectors irA, icA and valA, if not already done.
      
      if ( ~grpsdef )
         initA( indlvl, bindent, pbs );
      end   
       
   %  Group's constants
   
   elseif( ismember( strtrim( line ), EQCONSTANTS ) )
  
      section    = 'CONSTANTS';
      switch( pbs.lang )
      case 'matlab'
         pending{1}   = '%%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%';
         pending{end+1} = sprintf( 'pbm.gconst = zeros(ngrp,1);' );
      case 'python'
         pending{1}   = '#%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%';
         pending{end+1} = sprintf( 'self.gconst = np.zeros((ngrp,1))' );
      case 'julia'
         pending{1}   = '#%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%';
         pending{end+1} = sprintf( 'pbm.gconst = zeros(Float64,ngrp)' );
      end
      pendingkey    = 'default';

   %  Groups' ranges
   
   elseif ( lline >= 6 && strcmp( line(1:6), 'RANGES' ) ) 

      section    = 'RANGES';
      switch( pbs.lang )
      case 'matlab'
         pending{1}    = '%%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%';
      case 'python'
         pending{1}    = '#%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%';
         pending{end+1}    = 'grange = np.full((ngrp,1),None)';
      case 'julia'
         pending{1}    = '#%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%';
         pending{end+1}    = 'grange = Vector{Float64}(undef,ngrp)';
      end
      if( has_leconstr )
         switch( pbs.lang )
         case 'matlab'
            pending{end+1} = sprintf( 'grange(legrps,1) = Inf*ones(pb.nle,1);' );
         case 'python'
            pending{end+1} = sprintf( 'grange[legrps] = np.full((self.nle,1),float(''inf''))' );
         case 'julia'
            pending{end+1} = sprintf( 'grange[legrps,1] = fill(Inf,pb.nle)' );
         end
      end
      if( has_geconstr )
         switch ( pbs.lang )
         case 'matlab'
            pending{end+1} = sprintf( 'grange(gegrps,1) = Inf*ones(pb.nge,1);' );
         case 'python'
            pending{end+1} = sprintf( 'grange[gegrps] = np.full((self.nge,1),float(''inf''))' );
         case 'julia'
            pending{end+1} = sprintf( 'grange[gegrps,1] = fill(Inf,pb.nge)' );
         end
      end
      pendingkey     = 'default';
      prevlineispass = 0;
      
   %  Bounds on the variables
   
   elseif ( lline >= 6 && strcmp( line(1:6), 'BOUNDS' ) )

      section    = 'BOUNDS';
      switch ( pbs.lang )
      case 'matlab'
         pending{1} = '%%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%';
      case 'python'
         pending{1} = '#%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%';
      case 'julia'
         pending{1} = '#%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%';
      end
      pendingkey     = 'default';
      MPSbounds      = 0; %  the number of satisfied conditions for MPS-style bounds on the variables
      prevlineispass = 0;
      has_bounds     = 1;

   %  Starting values for the variables
      
   elseif ( lline >= 11 && strcmp( line(1:11), 'START POINT' ) )

      section    = 'STARTPOINT';
      switch ( pbs.lang )
      case 'matlab'
         pending{1} = '%%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%';
      case 'python'
         pending{1} = '#%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%';
      case 'julia'
         pending{1} = '#%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%';
      end
      pendingkey     = 'default';
      prevlineispass = 0;
      
   %  Objective quadratic term

   elseif( ismember( strtrim( line ), EQQUADRATIC ) )

      section    = 'QUADRATIC';
      switch( pbs.lang )
      case 'matlab'
         printmline( '%%%%%%%%%%%%%%%%%%%%% QUADRATIC %%%%%%%%%%%%%%%%%%%',                   indlvl, bindent, pbs.fidma );
         printmline( 'irH  = [];',                                                            indlvl, bindent, pbs.fidma );
         printmline( 'icH  = [];',                                                            indlvl, bindent, pbs.fidma );
         printmline( 'valH = [];',                                                            indlvl, bindent, pbs.fidma );
      case 'python'
         printpline( '#%%%%%%%%%%%%%%%%%%%% QUADRATIC %%%%%%%%%%%%%%%%%%%',                   indlvl, bindent, pbs.fidpy );
         printpline( 'irH  = np.array([],dtype=int)',                                         indlvl, bindent, pbs.fidpy );
         printpline( 'icH  = np.array([],dtype=int)',                                         indlvl, bindent, pbs.fidpy );
         printpline( 'valH = np.array([],dtype=float)',                                       indlvl, bindent, pbs.fidpy );
      case 'julia'
         printjline( '#%%%%%%%%%%%%%%%%%%%% QUADRATIC %%%%%%%%%%%%%%%%%%%',                   indlvl, bindent, pbs.fidjl );
         printjline( 'irH  = Int64[]',                                                        indlvl, bindent, pbs.fidjl );
         printjline( 'icH  = Int64[]',                                                        indlvl, bindent, pbs.fidjl );
         printjline( 'valH = Float64[]',                                                      indlvl, bindent, pbs.fidjl );
      end
      
   %  Nonlinear element types

   elseif ( lline >= 12 && strcmp( line(1:12), 'ELEMENT TYPE' ) )

      section    = 'ELFTYPE';
      switch( pbs.lang )
      case 'matlab'
         printmline( '%%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%',                   indlvl, bindent, pbs.fidma );
         switch( pbs.dicttype )
         case 'native'
            printmline( 'iet_ = configureDictionary(''string'',''double'');',                 indlvl, bindent, pbs.fidma );
         case 'custom'
            printmline( 'iet_ = containers.Map(''KeyType'', ''char'', ''ValueType'',''double'');', ...
                                                                                              indlvl, bindent, pbs.fidma );
         end
      case 'python'
         printpline( '#%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%',                   indlvl, bindent, pbs.fidpy );
         printpline( 'iet_  = {}',                                                            indlvl, bindent, pbs.fidpy );
         printpline( 'elftv = []',                                                            indlvl, bindent, pbs.fidpy );
      case 'julia'
         printjline( '#%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%',                   indlvl, bindent, pbs.fidjl );
         printjline( 'iet_  = Dict{String,Int}()',                                            indlvl, bindent, pbs.fidjl );
         printjline( 'elftv = Vector{Vector{String}}()',                                      indlvl, bindent, pbs.fidjl );
      end
      
      prevetype = '';
      
   %  Nonlinear elements

   elseif ( lline >= 12 && strcmp( line(1:12), 'ELEMENT USES' ) ) 

      %  ie_ is the dictionnary associating the falt index of an element to its run-time name.
      %  elftype{} is the internal structure defining an element type (name, elvar, epar)
      %  ielftype(ie) gives the index of the elftype of element whose flat index is ie, so that
      %  elftype{ielftype(ie)} is the structure of element with flatindex ie.
      
      section    = 'ELUSES';
      switch( pbs.lang )
      case 'matlab'
         printmline( '%%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%',                   indlvl, bindent, pbs.fidma );
         switch ( pbs.dicttype )
         case 'native'
            printmline( 'ie_ = configureDictionary(''string'',''double'');',                  indlvl, bindent, pbs.fidma );
         case 'custom'
            printmline( 'ie_ = containers.Map(''KeyType'',''char'',''ValueType'',''double'');',indlvl, bindent, pbs.fidma);
         end
         if ( getenames )
            printmline( 'pbm.enames  = {};',                                                  indlvl, bindent, pbs.fidma );
         end
         printmline( 'pbm.elftype = {};',                                                     indlvl, bindent, pbs.fidma );
         printmline( 'ielftype    = [];',                                                     indlvl, bindent, pbs.fidma );
         printmline( 'pbm.elvar   = {};',                                                     indlvl, bindent, pbs.fidma );
         if ( has_elpar )
            printmline( 'pbm.elpar   = {};',                                                  indlvl, bindent, pbs.fidma );
         end
      case 'python'
         printpline( '#%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%',                   indlvl, bindent, pbs.fidpy );
         printpline( 'ie_ = {}',                                                              indlvl, bindent, pbs.fidpy );
         if ( getenames )
            printpline( 'self.enames  = np.array([])',                                        indlvl, bindent, pbs.fidpy );
         end
         printpline( 'self.elftype = np.array([])',                                           indlvl, bindent, pbs.fidpy );
         printpline( 'ielftype     = np.array([])',                                           indlvl, bindent, pbs.fidpy );
         printpline( 'self.elvar   = []',                                                     indlvl, bindent, pbs.fidpy );
         if ( has_elpar )
            printpline( 'self.elpar   = []',                                                  indlvl, bindent, pbs.fidpy );
         end
      case 'julia'
         printjline( '#%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%',                   indlvl, bindent, pbs.fidjl );
         printjline( 'ie_      = Dict{String,Int}()',                                         indlvl, bindent, pbs.fidjl );
         printjline( 'ielftype = Vector{Int64}()',                                            indlvl, bindent, pbs.fidjl );
      end
      
      defelftype  = '';  % nonempty if an explicit default ftype is declared for elements
      prevename   = '';
      
   %  Nonlinear group types

   elseif ( lline >= 10 && strcmp( line(1:10), 'GROUP TYPE' ) ) 

      section    = 'GRFTYPE';
      switch( pbs.lang )
      case 'matlab'
         printmline( '%%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%',                   indlvl, bindent, pbs.fidma );
         switch (pbs.dicttype )
         case 'native'
            printmline( 'igt_ = configureDictionary(''string'',''double'');',                 indlvl, bindent, pbs.fidma );
         case 'custom'
            printmline( 'igt_ = containers.Map(''KeyType'',''char'',''ValueType'',''double'');', ...
                                                                                              indlvl, bindent, pbs.fidma );
         end
      case 'python'
         printpline( '#%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%',                   indlvl, bindent, pbs.fidpy );
         printpline( 'igt_ = {}',                                                             indlvl, bindent, pbs.fidpy );
      case 'julia'
         printjline( '#%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%',                   indlvl, bindent, pbs.fidjl );
         printjline( 'igt_ = Dict{String,Int}()',                                             indlvl, bindent, pbs.fidjl );
      end
      
   %  Assignment of elements to groups

   elseif ( lline >= 10 && strcmp( line(1:10), 'GROUP USES' ) ) 

      section    = 'GRUSES';
      switch ( pbs.lang )
      case 'matlab'
         printmline( '%%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%',                   indlvl, bindent, pbs.fidma );
         printmline( '[pbm.grelt{1:ngrp}] = deal([]);',                                       indlvl, bindent, pbs.fidma );
         printmline( 'nlc = [];',                                                             indlvl, bindent, pbs.fidma );
      case 'python'
         printpline( '#%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%',               indlvl,     bindent, pbs.fidpy );
         printpline( 'self.grelt   = []' ,                                                indlvl,     bindent, pbs.fidpy );
         printpline( 'for ig in np.arange(0,ngrp):',                                      indlvl,     bindent, pbs.fidpy )
         printpline( 'self.grelt.append(np.array([]))',                                   indlvl + 1, bindent, pbs.fidpy );
         printpline( 'self.grftype = np.array([])',                                       indlvl,     bindent, pbs.fidpy );
         printpline( 'self.grelw   = []',                                                 indlvl,     bindent, pbs.fidpy );
         printpline( 'nlc         = np.array([])',                                        indlvl,     bindent, pbs.fidpy );
         if ( has_grpar )
            printpline( 'self.grpar   = []',                                              indlvl,     bindent, pbs.fidpy );
         end
      case 'julia'
         printjline( '#%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%',               indlvl,     bindent, pbs.fidjl );
         printjline( 'for ig in 1:ngrp',                                                  indlvl,     bindent, pbs.fidjl );
         printjline( 'arrset(pbm.grelt,ig,Int64[])',                                      indlvl + 1, bindent, pbs.fidjl );
         printjline( 'end',                                                               indlvl,     bindent, pbs.fidjl );
         printjline( 'nlc = Int64[]',                                                     indlvl,     bindent, pbs.fidjl );
      end
      prevgname = '';

   %  Known bounds on the objective function's values
   
   elseif ( lline >= 12 && strcmp( line(1:12), 'OBJECT BOUND' ) ) 

      section    = 'OBOUND';
      printmline( '%%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%',                      indlvl, bindent, pbs.fidma );
      printpline( '#%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%',                      indlvl, bindent, pbs.fidpy );
      printjline( '#%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%',                      indlvl, bindent, pbs.fidjl );
      prevlineispass = 0;
      inobounds      = 1;

   %  Nonlinear elements...

   elseif ( lline >= 8 && strcmp( line(1:8), 'ELEMENTS' ) )
   
      printmline( ' ',                                                                    0,          bindent, pbs.fidma );
      printmline( '%%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%',                  indlvl - 1, bindent, pbs.fidma );
      printpline( ' ',                                                                    0,          bindent, pbs.fidpy );
      printpline( '#%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%',                  indlvl - 1, bindent, pbs.fidpy );
      printjline( ' ',                                                                    0,          bindent, pbs.fidjl );
      printjline( '#%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%',                  indlvl - 1, bindent, pbs.fidjl );
      inelements = 1;
      codeg      = {}; % the Matlab holder for code describing the element's gradient
      codeH      = {}; % the Matlab holder for code describing the element's Hessian
      inobounds  = 0;
      
   %  ... declarations of names...
      
   elseif ( lline >= 11 && strcmp( line(1:11), 'TEMPORARIES' ) && inelements ) 

      section = 'ETEMPS';

   %  ... assignemt of parameters...
   
   elseif ( lline >= 7 && strcmp( line(1:7), 'GLOBALS' ) && inelements ) 

      section = 'EGLOBS';
      
      %  There are global elements parameters: open a specific action in the output file

      switch( pbs.lang )
      case 'matlab'
         printmline( ' ',                                                                     0,      bindent, pbs.fidma );
         printmline( 'case ''e_globs''',                                                      1,      bindent, pbs.fidma );
         printmline( ' ',                                                                     0,      bindent, pbs.fidma );
         printmline( 'pbm = varargin{1};',                                                    2,      bindent, pbs.fidma );
      case 'python'
         printpline( ' ',                                                                     0,      bindent, pbs.fidpy );
         printpline( '@staticmethod',                                                         1,      bindent, pbs.fidpy );
         printpline( 'def e_globs(self):',                                                    1,      bindent, pbs.fidpy );
         printpline( ' ',                                                                     0,      bindent, pbs.fidpy );
         printpline( 'import numpy as np',                                                    2,      bindent, pbs.fidpy );
      case 'julia'
         printjline( ' ',                                                                     0,      bindent, pbs.fidjl );
         printjline( 'elseif action == "e_globs"',                                            1,      bindent, pbs.fidjl );
         printjline( ' ',                                                                     0,      bindent, pbs.fidjl );
         printjline( 'pbm = args[1]',                                                         2,      bindent, pbs.fidjl );
      end

      
   elseif ( lline >= 11 && strcmp( line(1:11), 'INDIVIDUALS' ) && inelements ) 

   %  ... and computation of the relevant values

      if ( n_eglobs )
      
         % Close the eglobs action.
          
         printmline( 'varargout{1} = pbm;',                                                   indlvl, bindent, pbs.fidma );
         printpline( 'return pbm',                                                            indlvl, bindent, pbs.fidpy );
         printjline( 'return pbm',                                                            indlvl, bindent, pbs.fidjl );
      end
      
      section  = 'EINDIVS';
      todefine = []; % 1 if the ith internal variable must still be defined, 0 if it has been defined.
                     % Empty if no internal variables.
                     % Needed here because conclude_nlf is called before any element has been seen
      
   %  Nonlinear groups

   elseif ( lline >= 6 && strcmp( line(1:6), 'GROUPS' ) )

      switch ( pbs.lang )
      case 'matlab'
         printmline( ' ',                                                                 0,          bindent, pbs.fidma );
         printmline( '%%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%',              indlvl - 1, bindent, pbs.fidma );      case 'python'
         printpline( ' ',                                                                 0,          bindent, pbs.fidpy );
         printpline( '#%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%',              indlvl - 1, bindent, pbs.fidpy );      case 'julia'
         printjline( ' ',                                                                 0,          bindent, pbs.fidjl );
         printjline( '#%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%',              indlvl - 1, bindent, pbs.fidjl );
      end
      ingroups  = 1;

      %  Reset the global dictionary for ELEMENTS/GROUPS
      
      switch (pbs.dicttype )
      case 'native'
         globdict  = configureDictionary( 'string', 'string' ); 
      case 'custom'
         globdict  = containers.Map('KeyType','char','ValueType','char' );
      end
      
      codeg     = {};                  % the Matlab holder for code describing the group function's gradient
      codeH     = {};                  % the Matlab holder for code describing the group function's Hessian
      inobounds = 0;
      
   %  ... declarations of names...
      
   elseif ( lline >= 11 && strcmp( line(1:11), 'TEMPORARIES' ) && ingroups ) 

      section = 'GTEMPS';

   %  ... assignemt of parameters...
   
   elseif ( lline >= 7 && strcmp( line(1:7), 'GLOBALS' ) && ingroups ) 

      section = 'GGLOBS';
      
      %  There are global group parameters: open a specific action in the output file

      switch ( pbs.lang )
      case 'matlab'
         printmline( ' ',                                                                     0,      bindent, pbs.fidma );
         printmline( 'case ''g_globs''',                                                      1,      bindent, pbs.fidma );
         printmline( ' ',                                                                     0,      bindent, pbs.fidma );
         printmline( 'pbm = varargin{1};',                                                    2,      bindent, pbs.fidma );
      case 'python'
         printpline( ' ',                                                                     0,      bindent, pbs.fidpy );
         printpline( '@staticmethod',                                                         1,      bindent, pbs.fidpy );
         printpline( 'def g_globs(self):',                                                    1,      bindent, pbs.fidpy );
         printpline( ' ',                                                                     0,      bindent, pbs.fidpy );
      case 'julia'
         printjline( 'elseif action == "g_globs"',                                            1,      bindent, pbs.fidjl );
         printjline( ' ',                                                                     0,      bindent, pbs.fidjl );
         printjline( 'pbm = args[1]',                                                         2,      bindent, pbs.fidjl );
      end

   %  ... and computation of the relevant values
   
   elseif ( lline >= 11 && strcmp( line(1:11), 'INDIVIDUALS' ) && ingroups ) 

      %  Output the pending instructions in the output file, if any.

      if ( n_gglobs )

         % Close the gglobs action.
           
         printmline( 'varargout{1} = pbm;',                                                   indlvl, bindent, pbs.fidma );
         printpline( 'return pbm',                                                            indlvl, bindent, pbs.fidpy );
         printjline( 'return pbm',                                                            indlvl, bindent, pbs.fidjl );
      end
      section = 'GINDIVS';

   %  End of problem data

   elseif( lline >= 6 && strcmp( line(1:6), 'ENDATA' ) )

      if ( indata )
         indata = 0;

         if ( has_A || has_H)
         
            printmline( '%%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%',                indlvl, bindent, pbs.fidma );
            printpline( '#%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%',                indlvl, bindent, pbs.fidpy );
            printjline( '#%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%',                indlvl, bindent, pbs.fidjl );

            %  Build A from the vectors irA, icA and valA.

            if ( has_A )
               printmline( 'pbm.A = sparse(irA,icA,valA,ngrp,pb.n);',                         indlvl, bindent, pbs.fidma );
               printpline( 'self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))',       indlvl, bindent, pbs.fidpy );
               printjline( 'pbm.A = sparse(irA,icA,valA,ngrp,pb.n)',                          indlvl, bindent, pbs.fidjl );
            end

            %  Build H from the vectors irH, icH and valH.

            if ( has_H )
               printmline( 'pbm.H = sparse(irH,icH,valH,pb.n,pb.n);',                         indlvl, bindent, pbs.fidma );
               printpline( 'self.H = csr_matrix((valH,(irH,icH)),shape=(self.n,self.n))',     indlvl, bindent, pbs.fidpy );
               printjline( 'pbm.H = sparse(irH,icH,valH,pb.n,pb.n)',                          indlvl, bindent, pbs.fidjl );
            end
         end

         %  The decoding of the SETUP section is (nearly) complete. Conclude it by writing the last missing bits in the
         %  output file.
         
         %  Missing defaults

         printmline( '%%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%',                   indlvl, bindent, pbs.fidma );
         printpline( '#%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%',                   indlvl, bindent, pbs.fidpy );
         printjline( '#%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%',                   indlvl, bindent, pbs.fidjl );

         %  1) missing default ftype for groups.  Note this cannot happen for elements because
         %     an empty ELEMENT USES simply indicates there are no elements (and hence no defaults ftype for them).

         %  2) Missing variable's types

         if ( has_intvars || has_binvars )
            printmline( 'pb.xtype = repmat(''r'',1,pb.n);',                                   indlvl, bindent, pbs.fidma );
            printpline( 'self.xtype = np.full(self.n,''r'')',                                 indlvl, bindent, pbs.fidpy );
            if ( has_intvars )
               printmline( 'pb.xtype(intvars) = ''i'';',                                      indlvl, bindent, pbs.fidma );
               printpline( 'self.xtype[intvars] = np.full(len(intvars),''i'')',               indlvl, bindent, pbs.fidpy );
               printjline( 'pb.xtype[intvars] = fill("i",length(intvars))',                   indlvl, bindent, pbs.fidjl );
            end
            if ( has_binvars ) 
               printmline( 'pb.xtype(binvars) = ''b'';',                                      indlvl, bindent, pbs.fidma );
               printpline( 'self.xtype[binvars] = np.full(len(binvars),''b'')',               indlvl, bindent, pbs.fidpy );
               printjline( 'pb.xtype[binvars] = fill("b",length(binvars))',                   indlvl, bindent, pbs.fidjl );
            end
         end

         %  3) Missing defaults for bounds on the variables

         if ( ~has_bounds )
            switch ( pbs.lang )
            case 'matlab'
               printmline( 'pb.xlower = zeros(pb.n,1);',                                      indlvl, bindent, pbs.fidma );
               printmline( 'pb.xupper = +Inf*ones(pb.n,1);',                                  indlvl, bindent, pbs.fidma );
            case 'python'
               printpline( 'self.xlower = np.zeros((self.n,1))',                              indlvl, bindent, pbs.fidpy );
               printpline( 'self.xupper = np.full((self.n,1),+float(''Inf''))',               indlvl, bindent, pbs.fidpy );
            case 'julia'
               printjline( 'pb.xlower = zeros(Float64,pb.n)',                                 indlvl, bindent, pbs.fidjl );
               printjline( 'pb.xupper =    fill(Inf,pb.n)',                                   indlvl, bindent, pbs.fidjl );
            end   
         end

         %  Now define the constraints right-hand sides.
         
         if ( has_constraints )

            %  1) In the old SIF format
            
            if ( keepcformat )
            
               if ( has_ranges )
                  printmline( 'pb.ranges = grange(pbm.congrps);',                             indlvl, bindent, pbs.fidma );
                  printpline( 'self.ranges = grange[self.congrps]',                           indlvl, bindent, pbs.fidpy );
                  printjline( 'pb.ranges = grange[pbm.congrps]',                              indlvl, bindent, pbs.fidjl );
               end
               printmline( 'pb.ctypes = gtype(pbm.congrps);',                                 indlvl, bindent, pbs.fidma );
               printpline( 'self.ctypes = np.array([])',                                      indlvl, bindent, pbs.fidpy );
               printpline( 'self.ctypes = gtype[self.congrps]',                               indlvl, bindent, pbs.fidpy );
               printjline( 'pb.ctypes = gtype[pbm.congrps]',                                  indlvl, bindent, pbs.fidjl );
               
            %  2) In the better(?) S2MPJ format: computate the lower and upper bounds on constraints values. 

            else
            
               printmline( '%%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%',             indlvl, bindent, pbs.fidma );
               %
               printpline( '#%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%',             indlvl, bindent, pbs.fidpy );
               printpline( 'self.clower = np.full((self.m,1),-float(''Inf''))',               indlvl, bindent, pbs.fidpy );
               printpline( 'self.cupper = np.full((self.m,1),+float(''Inf''))',               indlvl, bindent, pbs.fidpy );
               %
               printjline( '#%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%',             indlvl, bindent, pbs.fidjl );
               printjline( 'pb.clower = -1*fill(Inf,pb.m)',                                   indlvl, bindent, pbs.fidjl );
               printjline( 'pb.cupper =    fill(Inf,pb.m)',                                   indlvl, bindent, pbs.fidjl );

               %   Find the indeces of the constraints with TRIVIAL group type.

               % 2a) <= constraints
               
               if ( has_leconstr )

                  %  Set the <= contraint's lower bound.
                  
                  if ( has_ranges )
                     printmline( 'pb.clower(1:pb.nle) = grange(legrps);',                     indlvl, bindent, pbs.fidma );
                     printpline( 'self.clower[np.arange(self.nle)] = grange[legrps]',         indlvl, bindent, pbs.fidpy );
                     printjline( 'pb.clower[1:pb.nle] = grange[legrps]',                      indlvl, bindent, pbs.fidjl );
                  else
                     printmline( 'pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);',                indlvl, bindent, pbs.fidma );
                  end

                  %  Set the <= constraint's upper bound
                   
                  printmline( 'pb.cupper(1:pb.nle) = zeros(pb.nle,1);',                       indlvl, bindent, pbs.fidma );
                  printpline( 'self.cupper[np.arange(self.nle)] = np.zeros((self.nle,1))',    indlvl, bindent, pbs.fidpy );
                  printjline( 'pb.cupper[1:pb.nle] = zeros(Float64,pb.nle)',                  indlvl, bindent, pbs.fidjl );
               end

               %  2b) == constraints: set both lower and upper bounds to the constant.
               
               if ( has_eqconstr )
               
                  printmline( 'pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);',         indlvl, bindent, pbs.fidma );
                  printmline( 'pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);',         indlvl, bindent, pbs.fidma );
                  %
                  printpline( 'self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))', ...
                                                                                              indlvl, bindent, pbs.fidpy );
                  printpline( 'self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))', ...
                                                                                              indlvl, bindent, pbs.fidpy );
                  %
                  printjline( 'pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)',    indlvl, bindent, pbs.fidjl );
                  printjline( 'pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)',    indlvl, bindent, pbs.fidjl );
               end

               %  2c) >= constraints
               
               if ( has_geconstr )
               
                  %  Set the >= constraint's lower bound.
                  
                  printmline( 'pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);',           indlvl, bindent, pbs.fidma );
                  printpline( 'self.clower[np.arange(self.nle+self.neq,self.m)] = np.zeros((self.nge,1))',  ...
                                                                                              indlvl, bindent, pbs.fidpy );
                  printjline( 'pb.clower[pb.nle+pb.neq+1:pb.m] = zeros(Float64,pb.nge)',      indlvl, bindent, pbs.fidjl );

                  %  Set the >= constraint's upper bound.
                  
                  if ( has_ranges )
                     printmline( 'pb.cupper(1:pb.nge) = grange(gegrps);',                     indlvl, bindent, pbs.fidma );
                     printpline( 'self.cupper[np.arange(self.nge)] = grange[gegrps]',         indlvl, bindent, pbs.fidpy );
                     printjline( 'pb.cupper[1:pb.nge] = grange[gegrps]',                      indlvl, bindent, pbs.fidjl );
                  else
                     printmline( 'pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);',                indlvl, bindent, pbs.fidma );
                     printjline( 'pb.cupper[1:pb.nge] = fill(Inf,pb.nge)',                    indlvl, bindent, pbs.fidjl );
                  end
               end
            end
         end

         %  Perform variable scaling if requested
         
         if ( has_A )

            if ( has_xscale )
               if ( ~extxscale )
                  switch ( pbs.lang )
                  case 'matlab'
                     printmline( '%%%%%%%%%%%%%%%% VARIABLES'' SCALING %%%%%%%%%%%%%%%',      indlvl, bindent, pbs.fidma );
                     printmline( 'sA2 = size(pbm.A,2);',                                      indlvl, bindent, pbs.fidma );
                     printmline( 'for j = 1:min([sA2,pb.n,length(pb.xscale)])',               indlvl, bindent, pbs.fidma );
                     printmline( '    if ( pb.xscale(j) ~= 0.0 && pb.xscale(j) ~= 1.0 )',     indlvl, bindent, pbs.fidma );
                     printmline( '        for i = find(pbm.A(:,j))',                          indlvl, bindent, pbs.fidma );
                     printmline( '              pbm.A(i,j) = pbm.A(i,j)/pb.xscale(j);',       indlvl, bindent, pbs.fidma );
                     printmline( '        end',                                               indlvl, bindent, pbs.fidma );
                     printmline( '    end',                                                   indlvl, bindent, pbs.fidma );
                     printmline( 'end',                                                       indlvl, bindent, pbs.fidma );
                     printmline( 'pb.xscale = [];',                                           indlvl, bindent, pbs.fidma );
                  case 'python'
                     printpline( '#%%%%%%%%%%%%%%% VARIABLES'' SCALING %%%%%%%%%%%%%%%',      indlvl, bindent, pbs.fidpy );
                     printpline( 'lxs = len(self.xscale);',                                   indlvl, bindent, pbs.fidpy );
                     printpline( 'for j in np.arange(0,min(sA2,self.n,lxs)):',                indlvl, bindent, pbs.fidpy );
                     printpline( 'if not self.xscale[j] is None and self.xscale[j] != 0.0 and self.xscale[j] != 1.0:', ...
                                                                                          indlvl + 1, bindent, pbs.fidpy );
                     printpline( '        for i in np.where(self.A[:,j]!=0)[0]:',             indlvl, bindent, pbs.fidpy );
                     printpline( '              self.A[i,j] = self.A[i,j]/self.xscale[j]',    indlvl, bindent, pbs.fidpy );
                     printpline( 'self.xscale = []',                                          indlvl, bindent, pbs.fidpy );
                  case 'julia'
                     printjline( '#%%%%%%%%%%%%%%% VARIABLES'' SCALING %%%%%%%%%%%%%%%',      indlvl, bindent, pbs.fidjl );
                     printjline( 'sA2 = size(pbm.A,2);',                                      indlvl, bindent, pbs.fidjl );
                     printjline( 'lxs = length(pb.xscale);',                                  indlvl, bindent, pbs.fidjl );
                     printjline( 'for j = 1:min(sA2,pb.n,lxs)',                               indlvl, bindent, pbs.fidjl );
                     printjline( '    if pb.xscale[j] != 0.0 && pb.xscale[j] != 1.0',         indlvl, bindent, pbs.fidjl );
                     printjline( '        for i in findall(x->x!=0,pbm.A[:,j])',              indlvl, bindent, pbs.fidjl );
                     printjline( '              pbm.A[i,j] = pbm.A[i,j]/pb.xscale[j]',        indlvl, bindent, pbs.fidjl );
                     printjline( '        end',                                               indlvl, bindent, pbs.fidjl );
                     printjline( '    end',                                                   indlvl, bindent, pbs.fidjl );
                     printjline( 'end',                                                       indlvl, bindent, pbs.fidjl );
                     printjline( 'pb.xscale = []',                                            indlvl, bindent, pbs.fidjl );
                  end
               end
            end
         end

         %  Finally assign the output values which have not been assigned yet.

         printmline( '%%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%',                   indlvl, bindent, pbs.fidma );
         printpline( '#%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%',                   indlvl, bindent, pbs.fidpy );
         printjline( '#%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%',                   indlvl, bindent, pbs.fidjl );

         %  1) Indeces of linear constraints, if relevant
         
         if ( has_constraints )
            if ( has_nonlinc )
               printmline( '[~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);',                ...
                                                                                              indlvl, bindent, pbs.fidma );
               printpline( 'self.lincons = np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0]', ...
                                                                                              indlvl, bindent, pbs.fidpy );
               printjline( 'pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)',             ...
                                                                                              indlvl, bindent, pbs.fidjl );
            else
               printmline( 'pb.lincons   = [1:length(pbm.congrps)];',                         indlvl, bindent, pbs.fidma );
               printpline( 'self.lincons   = np.arange(len(self.congrps))',                   indlvl, bindent, pbs.fidpy );
               printjline( 'pb.lincons   = collect(1:length(pbm.congrps))',                   indlvl, bindent, pbs.fidjl );
            end
         end

         %  2) Problem classification

         printmline( sprintf( 'pb.pbclass = %s;', classification ),                           indlvl, bindent, pbs.fidma );        
         printpline( sprintf( 'self.pbclass   = %s',  classification ),                       indlvl, bindent, pbs.fidpy );        
         printjline( sprintf( 'pb.pbclass = %s',  classification ),                           indlvl, bindent, pbs.fidjl );

         if ( ~has_start )
             printmline( 'pb.x0          = zeros(pb.n,1);',                                   indlvl, bindent, pbs.fidma );
             printpline( 'self.x0        = np.zeros((self.n,1))',                             indlvl, bindent, pbs.fidpy );
             printjline( 'pb.x0          = zeros(Float64,pb.n)',                              indlvl, bindent, pbs.fidjl );
         end

         %  3) level of available derivatives
         %     Note that first and second derivatives are availabvle for all CUTEst problems.

         printmline( 'pbm.objderlvl = 2;',                                                    indlvl, bindent, pbs.fidma );
         printmline( 'pb.objderlvl = pbm.objderlvl;',                                         indlvl, bindent, pbs.fidma );
         printpline( 'self.objderlvl = 2',                                                    indlvl, bindent, pbs.fidpy );
         printjline( 'pbm.objderlvl = 2',                                                     indlvl, bindent, pbs.fidjl );
         printjline( 'pb.objderlvl = pbm.objderlvl;',                                         indlvl, bindent, pbs.fidjl );
         if ( has_constraints )
             printmline( 'pbm.conderlvl = [2];',                                              indlvl, bindent, pbs.fidma );
             printmline( 'pb.conderlvl  = pbm.conderlvl;',                                    indlvl, bindent, pbs.fidma );
             printpline( 'self.conderlvl = [2]',                                              indlvl, bindent, pbs.fidpy );
             printjline( 'pbm.conderlvl = [2]',                                               indlvl, bindent, pbs.fidjl );
             printjline( 'pb.conderlvl  = pbm.conderlvl;',                                    indlvl, bindent, pbs.fidjl );
         end

         %  4) Reduced precision option
         
         switch( pbs.lang )
         case 'matlab'

            %  Prepare for variable-precision evaluations by converting pb and pbm to VPA.

            if ( redprec )
               printmline( '%%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%',        indlvl,     bindent, pbs.fidma );
               printmline( 'if(strcmp(action,''setup_redprec''))',                        indlvl,     bindent, pbs.fidma );
               printmline( 'varargout{1} = s2mpjlib(''convert'',pb, pbm.ndigs);',         indlvl + 1, bindent, pbs.fidma );
               printmline( 'varargout{2} = s2mpjlib(''convert'',pbm,pbm.ndigs);',         indlvl + 1, bindent, pbs.fidma );
               printmline( 'else',                                                        indlvl,     bindent, pbs.fidma );
               printmline( 'varargout{1} = pb;',                                          indlvl + 1, bindent, pbs.fidma );
               printmline( 'varargout{2} = pbm;',                                         indlvl + 1, bindent, pbs.fidma );
               printmline( 'end',                                                         indlvl,     bindent, pbs.fidma );
            else
               printmline( 'varargout{1} = pb;',                                          indlvl,     bindent, pbs.fidma );
               printmline( 'varargout{2} = pbm;',                                         indlvl,     bindent, pbs.fidma );
            end
         case 'python'
         case 'julia'
            printjline( 'return pb, pbm',                                                     indlvl, bindent, pbs.fidjl );
         end
         printmline( ' ',                                                                     0,      bindent, pbs.fidma );
         printpline( ' ',                                                                     0,      bindent, pbs.fidpy );
         printjline( ' ',                                                                     0,      bindent, pbs.fidjl );

      elseif ( inelements )
      
         %  Conclude the code for the last element by writing pending code  for derivatives.
            
         conclude_nlf( codeg, codeH, isempty( todefine ),  indlvl, bindent, pbs );
         inelements = 0;

      else % ingroups
       
          %  Conclude the code for the last group by writing pending code for derivatives.
            
         conclude_nlf( codeg, codeH, 1, indlvl, bindent, pbs );
         ingroups = 0;
 
      end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Data lines   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   elseif ( indata )

      incomments = 0;
      
      %  Parse the data line into its up to 6 fields.

      [ f, nf ]    = parsedataline( line );
      is_param_def = 1;         %  first attempt to define a parameter or a loop

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  DO LOOPS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      %  The Matlab/Python instruction defining the loop is not immediately written in the output file as soon as defined
      %  by a DO data line, because it may have to be modified by an immediately subsequent DI data line.  It is thus kept
      %  in memory (in pending{1}) until the next dataline which is not a DI line is read.
      
      switch ( f{1} )

      %  Start of new loop
      
      case 'DO'

         %  Create a new loop with those bounds, making sure not to use pbs.loop before it is created.

         index         = s2mpjname( '', f{2}, pbs );    %  loop indeces
         iename        = s2mpjname( '', f{5}, pbs );
         iloop         = length( pbs.loop ) + 1;        %  internal loop index
         pbs.actloop   = [ pbs.actloop, iloop ];        %  update the list of active loops
         pbs.loop{iloop}.index = index(2:end-1);        %  populate the internal loop struct
         pbs.loop{iloop}.iend  = iename(2:end-1);
         switch ( pbs.lang )
         case 'matlab'
            pending{1} = sprintf( 'for %s=%s:%s', ...
                         index(2:end-1), s2mpjvalue('', f{3}, pbs ), s2mpjvalue( '', f{5}, pbs ) );
         case 'python'
            pending{1} = sprintf( 'for %s in range(int(%s),int(%s)+1):', ...
                         index(2:end-1), s2mpjvalue('', f{3}, pbs ), s2mpjvalue( '', f{5}, pbs ) );
         case 'julia'
            pending{1} = sprintf( 'for %s = Int64(%s):Int64(%s)', ...
                         index(2:end-1), s2mpjvalue('', f{3}, pbs ), s2mpjvalue( '', f{5}, pbs ) );
         end
         pendingkey    = 'loop';
         prevename     = '';
         prevgname     = '';
         continue      %  get the next SIF data line

      %  Loop increment specification: if met, modify the pending loop statement and write it in the output file.

      case 'DI' 
         switch( pbs.lang )
         case 'matlab'
            printmline( replace( pending{1}, ':', sprintf( ':%s:', s2mpjvalue( '', f{3}, pbs ) ) ), ...
                                                                                              indlvl, bindent, pbs.fidma );
         case 'python'
            er = strfind( pending{1}, '):' );
            printpline( sprintf( '%s,int(%s)):', pending{1}(1:er-1), s2mpjvalue( '', f{3}, pbs ) ), ...
                                                                                              indlvl, bindent, pbs.fidpy );
         case 'julia'
            printjline( replace( pending{1}, ':', sprintf( ':Int64(%s):', s2mpjvalue( '', f{3}, pbs ) ) ), ...
                                                                                              indlvl, bindent, pbs.fidjl );
         end
         pending        = {};
         pendingkey     = '';
         indlvl         = indlvl + 1;
         prevlineispass = 0;
         continue   %  get the next SIF data line

      %  End of a loop: once the loop is completed, write a Matlab/Python instruction defining an integer parameter
      %  holding the last value of the loop index. Also adjust the output file's indentation.

      case 'OD'

         indlvl      = indlvl - 1;
         printmline( 'end',                                                                   indlvl, bindent, pbs.fidma );
         printjline( 'end',                                                                   indlvl, bindent, pbs.fidjl );
         iloop = pbs.actloop(end);
         switch( pbs.lang )
         case { 'matlab', 'python' }
            pbs.irpdict( pbs.loop{iloop}.index ) =  [ '''', pbs.loop{iloop}.iend, '''' ];
         case 'julia'
            pbs.irpdict( pbs.loop{iloop}.index ) =  [ '"', pbs.loop{iloop}.iend, '"' ];
         end
         pbs.actloop = pbs.actloop(1:end-1);
         prevename = '';
         prevgname = '';
         continue    %  get the next SIF data line

      %  End of all loops

      case 'ND'
      
         for i = length( pbs.actloop ):-1:1
            indlvl      = indlvl - 1;
            printmline( 'end',                                                                indlvl, bindent, pbs.fidma );
            printjline( 'end',                                                                indlvl, bindent, pbs.fidjl );
            iloop       = pbs.actloop(end);
            switch( pbs.lang )
            case { 'matlab', 'python' }
               pbs.irpdict( pbs.loop{iloop}.index ) =  [ '''', pbs.loop{iloop}.iend, '''' ];
            case 'julia'
               pbs.irpdict( pbs.loop{iloop}.index ) =  [ '"', pbs.loop{iloop}.iend, '"' ];
            end
            pbs.actloop = pbs.actloop(1:end-1);
         end
         prevename = '';
         prevgname = '';
         continue       %  get the next SIF data line

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      %  A definition line is written for each of the various data line types defining parameters (f{2}), each time
      %  screening the name of the relevant parameters using s2mpjname.

      case { 'IE', 'RE', 'AE' }

         %  Detect SIF file parameter and transform them into output file arguments provided enough of them are present in
         %  the call.
         
         name2 = s2mpjname( f{1}(1), f{2}, pbs );
         if ( nf >= 4  && contains( line, '$-PARAMETER' ) )
            if ( ~ismember( name2, readsifpars ) )
               readsifpars{end+1} = name2;
               ninputs = length( readsifpars );
            end
            switch( pbs.lang )
            case 'matlab'

               %  If variable precision is requested, save the last input (the number of digits) and test for a shorter
               %  list  of input arguments.
             
               printmline( sprintf( 'if(nargs<%s)', int2str( ninputs ) ),                 indlvl,     bindent, pbs.fidma );
               if ( strcmp( f{1}, 'IE' ) )
                  printmline( sprintf( 'v_(%s) = %s;  %%  SIF file default value', ...
                                       name2, round( d2e(f{4}) ) ),                       indlvl + 1, bindent, pbs.fidma );
               else
                  printmline( sprintf( 'v_(%s) = %s;  %%  SIF file default value', ...
                                       name2, d2e(f{4}) ),                                indlvl + 1, bindent, pbs.fidma );
               end
               printmline( 'else',                                                        indlvl,     bindent, pbs.fidma );
               evalstr = [ 'varargin{', int2str( ninputs ), '}' ];
               printmline( sprintf( 'v_(%s) = %s;', name2, evalstr ),                     indlvl + 1, bindent, pbs.fidma );
               printmline( 'end',                                                         indlvl,     bindent, pbs.fidma );
            case 'python'
               printpline( sprintf( 'if nargin<%s:', int2str( ninputs ) ),                indlvl,     bindent, pbs.fidpy );
               if ( strcmp( f{1}, 'IE' ) )
                  printpline( sprintf( 'v_[%s] = int(%s);  #  SIF file default value', ...
                                       name2, round( d2e(f{4}) ) ),                       indlvl + 1, bindent, pbs.fidpy );
               else
                  printpline( sprintf( 'v_[%s] = float(%s);  #  SIF file default value', ...
                                       name2, d2e(f{4}) ),                                indlvl + 1, bindent, pbs.fidpy );
               end
               printpline( 'else:', indlvl, bindent, pbs.fidpy );
               evalstr = [ 'args[', int2str( ninputs-1 ), ']' ];
               if  ( strcmp( f{1}, 'IE' ) )
                  printpline( sprintf( 'v_[%s] = int(%s)', name2, evalstr ),              indlvl + 1, bindent, pbs.fidpy );
               else
                  printpline( sprintf( 'v_[%s] = float(%s)', name2, evalstr ),            indlvl + 1, bindent, pbs.fidpy );
               end               
            case 'julia'
               printjline( sprintf( 'if nargin<%s', int2str( ninputs ) ),                 indlvl,     bindent, pbs.fidjl );
               if ( strcmp( f{1}, 'IE' ) )
                  printjline( sprintf( 'v_[%s] = Int64(%s);  #  SIF file default value', ...
                                       name2, round( d2e(f{4}) ) ),                       indlvl + 1, bindent, pbs.fidjl );
               else
                  printjline( sprintf( 'v_[%s] = Float64(%s);  #  SIF file default value', ...
                                       name2, d2e(f{4}) ),                                indlvl + 1, bindent, pbs.fidjl );
               end
               printjline( 'else', indlvl, bindent, pbs.fidjl );
               evalstr = [ 'args[', int2str( ninputs ), ']' ];
               if  ( strcmp( f{1}, 'IE' ) )
                  printjline( sprintf( 'v_[%s] = Int64(%s);', name2, evalstr ),           indlvl + 1, bindent, pbs.fidjl );
               else
                  printjline( sprintf( 'v_[%s] = Float64(%s);', name2, evalstr ),         indlvl + 1, bindent, pbs.fidjl );
               end
               printjline( 'end',                                                         indlvl,     bindent, pbs.fidjl );
            end
            
            
         %  Standard parameter assignment
         
         else
            printmline( sprintf( 'v_(%s) = %s;', name2, d2e( f{4} ) ),                        indlvl, bindent, pbs.fidma );
            printpline( sprintf( 'v_[%s] = %s',  name2, d2e( f{4} ) ),                        indlvl, bindent, pbs.fidpy );
            printjline( sprintf( 'v_[%s] = %s',  name2, d2e( f{4} ) ),                        indlvl, bindent, pbs.fidjl );
         end

      case 'IR'

         name2  = s2mpjname ( f{1}(1), f{2}, pbs );
         value3 = s2mpjvalue( f{1}(1), f{3}, pbs );
         printmline( sprintf( 'v_(%s) = fix(%s);',        name2, value3 ),                    indlvl, bindent, pbs.fidma );
         printpline( sprintf( 'v_[%s] = int(np.fix(%s))', name2, value3 ),                    indlvl, bindent, pbs.fidpy );
         printjline( sprintf( 'v_[%s] = trunc(Int,%s)',   name2, value3 ),                    indlvl, bindent, pbs.fidjl );

      case { 'RI', 'AI' }

         [ name2, isidx ]  = s2mpjname ( f{1}(1), f{2}, pbs );
         if ( isidx )
            switch ( pbs.lang )
            case 'matlab'
               name2 = [ '''', name2, '''' ];
            case { 'python', 'julia' }
               name2 = [ '"', name2, '"'];
            end
         end
         value3 = s2mpjvalue( f{1}(1), f{3}, pbs );
         printmline( sprintf( 'v_(%s) = %s;',          name2, value3 ),                       indlvl, bindent, pbs.fidma );
         printpline( sprintf( 'v_[%s] = float(%s)'  ,  name2, value3 ),                       indlvl, bindent, pbs.fidpy );
         printjline( sprintf( 'v_[%s] = Float64(%s)',  name2, value3 ),                       indlvl, bindent, pbs.fidjl );
         
      case { 'IA', 'RA', 'AA' }

         name2  = s2mpjname ( f{1}(1), f{2}, pbs );
         value3 = s2mpjvalue( f{1}(1), f{3}, pbs );
         printmline( sprintf( 'v_(%s) = %s+%s;', name2, d2e( f{4} ), value3 ),                indlvl, bindent, pbs.fidma );
         printpline( sprintf( 'v_[%s] = %s+%s',  name2, d2e( f{4} ), value3 ),                indlvl, bindent, pbs.fidpy );
         printjline( sprintf( 'v_[%s] = %s+%s',  name2, d2e( f{4} ), value3 ),                indlvl, bindent, pbs.fidjl );

      case { 'IS', 'RS', 'AS' }

         name2  = s2mpjname ( f{1}(1), f{2}, pbs );
         value3 = s2mpjvalue( f{1}(1), f{3}, pbs );
         printmline( sprintf( 'v_(%s) = %s-%s;', name2, d2e( f{4} ), value3 ),                indlvl, bindent, pbs.fidma );
         printpline( sprintf( 'v_[%s] = %s-%s',  name2, d2e( f{4} ), value3 ),                indlvl, bindent, pbs.fidpy );
         printjline( sprintf( 'v_[%s] = %s-%s',  name2, d2e( f{4} ), value3 ),                indlvl, bindent, pbs.fidjl );

      case { 'IM', 'RM', 'AM' }

         name2  = s2mpjname ( f{1}(1), f{2}, pbs );
         value3 = s2mpjvalue( f{1}(1), f{3}, pbs );
         printmline( sprintf( 'v_(%s) = %s*%s;', name2, d2e( f{4} ), value3 ),                indlvl, bindent, pbs.fidma );
         printpline( sprintf( 'v_[%s] = %s*%s',  name2, d2e( f{4} ), value3 ),                indlvl, bindent, pbs.fidpy );
         printjline( sprintf( 'v_[%s] = %s*%s',  name2, d2e( f{4} ), value3 ),                indlvl, bindent, pbs.fidjl );

      case 'ID'

         name2  = s2mpjname ( f{1}(1), f{2}, pbs );
         value3 = s2mpjvalue( f{1}(1), f{3}, pbs );
         printmline( sprintf( 'v_(%s) = fix(%s/%s);',        name2, d2e( f{4} ), value3 ),    indlvl, bindent, pbs.fidma );
         printpline( sprintf( 'v_[%s] = int(np.fix(%s/%s))', name2, d2e( f{4} ), value3 ),    indlvl, bindent, pbs.fidpy );
         printjline( sprintf( 'v_[%s] = trunc(Int,(%s/%s))', name2, d2e( f{4} ), value3 ),    indlvl, bindent, pbs.fidjl );

      case { 'RD', 'AD' }

         name2  = s2mpjname ( f{1}(1), f{2}, pbs );
         value3 = s2mpjvalue( f{1}(1), f{3}, pbs );
         printmline( sprintf( 'v_(%s) = %s/%s;', name2, d2e( f{4} ), value3 ),                indlvl, bindent, pbs.fidma );
         printpline( sprintf( 'v_[%s] = %s/%s',  name2, d2e( f{4} ), value3 ),                indlvl, bindent, pbs.fidpy );
         printjline( sprintf( 'v_[%s] = %s/%s',  name2, d2e( f{4} ), value3 ),                indlvl, bindent, pbs.fidjl );

      case { 'I=', 'R=', 'A=' }

         name2  = s2mpjname ( f{1}(1), f{2}, pbs );
         value3 = s2mpjvalue( f{1}(1), f{3}, pbs );
         printmline( sprintf( 'v_(%s) = %s;', name2, value3 ),                                indlvl, bindent, pbs.fidma );
         printpline( sprintf( 'v_[%s] = %s',  name2, value3 ),                                indlvl, bindent, pbs.fidpy );
         printjline( sprintf( 'v_[%s] = %s',  name2, value3 ),                                indlvl, bindent, pbs.fidjl );

      case { 'I+', 'R+', 'A+' }
       
         name2  = s2mpjname ( f{1}(1), f{2}, pbs );
         value3 = s2mpjvalue( f{1}(1), f{3}, pbs );
         value5 = s2mpjvalue( f{1}(1), f{5}, pbs );
         printmline( sprintf( 'v_(%s) = %s+%s;', name2, value3, value5 ),                     indlvl, bindent, pbs.fidma );
         printpline( sprintf( 'v_[%s] = %s+%s',  name2, value3, value5 ),                     indlvl, bindent, pbs.fidpy );
         printjline( sprintf( 'v_[%s] = %s+%s',  name2, value3, value5 ),                     indlvl, bindent, pbs.fidjl );

      case { 'I-', 'R-', 'A-' }

         name2  = s2mpjname ( f{1}(1), f{2}, pbs );
         value3 = s2mpjvalue( f{1}(1), f{3}, pbs );
         value5 = s2mpjvalue( f{1}(1), f{5}, pbs );
         printmline( sprintf( 'v_(%s) = %s-%s;', name2, value3, value5 ),                     indlvl, bindent, pbs.fidma );
         printpline( sprintf( 'v_[%s] = %s-%s',  name2, value3, value5 ),                     indlvl, bindent, pbs.fidpy );
         printjline( sprintf( 'v_[%s] = %s-%s',  name2, value3, value5 ),                     indlvl, bindent, pbs.fidjl );

      case { 'I*', 'R*', 'A*' }

         name2  = s2mpjname ( f{1}(1), f{2}, pbs );
         value3 = s2mpjvalue( f{1}(1), f{3}, pbs );
         value5 = s2mpjvalue( f{1}(1), f{5}, pbs );
         printmline( sprintf( 'v_(%s) = %s*%s;', name2, value3, value5 ),                     indlvl, bindent, pbs.fidma );
         printpline( sprintf( 'v_[%s] = %s*%s',  name2, value3, value5 ),                     indlvl, bindent, pbs.fidpy );
         printjline( sprintf( 'v_[%s] = %s*%s',  name2, value3, value5 ),                     indlvl, bindent, pbs.fidjl );

      case 'I/'

         name2  = s2mpjname ( f{1}(1), f{2}, pbs );
         value3 = s2mpjvalue( f{1}(1), f{3}, pbs );
         value5 = s2mpjvalue( f{1}(1), f{5}, pbs );
         printmline( sprintf( 'v_(%s) = fix(%s/%s);',        name2, value3, value5 ),         indlvl, bindent, pbs.fidma );
         printpline( sprintf( 'v_[%s] = int(np.fix(%s/%s))', name2, value3, value5 ),         indlvl, bindent, pbs.fidpy );
         printjline( sprintf( 'v_[%s] = trunc(Int,(%s/%s))', name2, value3, value5 ),         indlvl, bindent, pbs.fidjl );

      case { 'R/', 'A/' }

         name2  = s2mpjname ( f{1}(1), f{2}, pbs );
         value3 = s2mpjvalue( f{1}(1), f{3}, pbs );
         value5 = s2mpjvalue( f{1}(1), f{5}, pbs );
         printmline( sprintf( 'v_(%s) = %s/%s;', name2, value3, value5 ),                     indlvl, bindent, pbs.fidma );
         printpline( sprintf( 'v_[%s] = %s/%s',  name2, value3, value5 ),                     indlvl, bindent, pbs.fidpy );
         printjline( sprintf( 'v_[%s] = %s/%s',  name2, value3, value5 ),                     indlvl, bindent, pbs.fidjl );

      case { 'RF', 'AF', 'R(', 'A(' }

         switch ( pbs.lang )
         case { 'matlab', 'julia' }
            switch( f{3} )
            case 'ABS'
               ff = 'abs';
            case 'SQRT'
               ff = 'sqrt';
            case 'EXP'
               ff = 'exp';
            case 'LOG'
               ff = 'log';
            case 'LOG10'
               ff = 'log10';
            case 'SIN'
               ff = 'sin';
            case 'COS'
               ff = 'cos';
            case 'TAN'
               ff = 'tan';
            case 'ARCSIN'
               ff = 'asin';
            case 'ARCCOS'
               ff = 'acos';
            case 'ARCTAN'
               ff = 'atan';
            case 'HYPSIN'
               ff = 'sinh';
            case 'HYPCOS'
               ff = 'cosh';
            case 'HYPTAN'
               ff = 'tanh';
            end
         case 'python'
            switch( f{3} )
            case 'ABS'
               ff = 'np.absolute';
            case 'SQRT'
               ff = 'np.sqrt';
            case 'EXP'
               ff = 'np.exp';
            case 'LOG'
               ff = 'np.log';
            case 'LOG10'
               ff = 'np.log10';
            case 'SIN'
               ff = 'np.sin';
            case 'COS'
               ff = 'np.cos';
            case 'TAN'
               ff = 'np.tan';
            case 'ARCSIN'
               ff = 'np.arcsin';
            case 'ARCCOS'
               ff = 'np.arccos';
            case 'ARCTAN'
               ff = 'np.arctan';
            case 'HYPSIN'
               ff = 'np.sinh';
            case 'HYPCOS'
               ff = 'np.cosh';
            case 'HYPTAN'
               ff = 'np.tanh';
            end
         end
         name2 = s2mpjname( f{1}(1), f{2}, pbs );
         switch( f{1} )
         case { 'RF', 'AF' }
            printmline( sprintf( 'v_(%s) = %s(%s);',name2, ff, d2e(f{4}) ),                   indlvl, bindent, pbs.fidma );
            printpline( sprintf( 'v_[%s] = %s(%s)', name2, ff, d2e(f{4}) ),                   indlvl, bindent, pbs.fidpy );
            printjline( sprintf( 'v_[%s] = %s(%s)', name2, ff, d2e(f{4}) ),                   indlvl, bindent, pbs.fidjl );
         case { 'R(', 'A(' }
            value5 = s2mpjvalue( f{1}(1), f{5}, pbs );
            printmline( sprintf( 'v_(%s) = %s(%s);',name2, ff, value5 ),                      indlvl, bindent, pbs.fidma );
            printpline( sprintf( 'v_[%s] = %s(%s)', name2, ff, value5 ),                      indlvl, bindent, pbs.fidpy );
            printjline( sprintf( 'v_[%s] = %s(%s)', name2, ff, value5 ),                      indlvl, bindent, pbs.fidjl );
         end
         
      otherwise

         is_param_def = 0;  %  The line does not involve parameters.
         
      end

      %  Store the names of the integer and real paramaters.
      
      if ( is_param_def && ( f{1}(1) == 'I' || f{1}(1) == 'R' ) )
         pbs.irpdict( name2(2:end-1) ) = name2;
      end

      %  If a paremeter or a loop was set, get the next SIF line.

      if ( is_param_def )
          continue
      end
   
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  STANDARD DATA LINES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      switch ( section )
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   VARIABLES   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %
      %  This section uses the following logical indicators:
      %  has_xscale       true iff the problem has scaling factors for the variables that are different from 1.0 or 0.0
      %                   (when numerical).
      %                   NOTE: no control is made on the value of the scaling factor if it is a real parameter
      %  has_intvars      true iff the problem has integer variables
      %  has_binvars      true iff the problem has binary (O/1) variables
      
      case 'VARIABLES'

         vname = s2mpjname( '', f{2}, pbs );
         printmline( sprintf( '[iv,ix_] = s2mpjlib(''ii'',%s,ix_);', vname ),                 indlvl, bindent, pbs.fidma );
         printpline( sprintf( '[iv,ix_,_] = s2mpj_ii(%s,ix_)', vname ),                       indlvl, bindent, pbs.fidpy );
         printjline( sprintf( 'iv,ix_,_ = s2mpj_ii(%s,ix_)', vname ),                         indlvl, bindent, pbs.fidjl );
         if ( getxnames )
            printmline( sprintf( 'pb.xnames{iv} = %s;', vname ),                              indlvl, bindent, pbs.fidma );
            printpline( sprintf( 'self.xnames=arrset(self.xnames,iv,%s)',  vname ),           indlvl, bindent, pbs.fidpy );
            printjline( sprintf( 'arrset(pb.xnames,iv,%s)',  vname ),                         indlvl, bindent, pbs.fidjl );
         end
         
         if ( nf > 2 )
            switch( f{3} )
            case '''SCALE'''
               if ( ~isempty( f{1} ) && f{1}(1) == 'Z' )
                  printmline( sprintf( 'pb.xscale(iv,1) = %s;',            getv1( f{1}, f, pbs ) ), ...
                                                                                              indlvl, bindent, pbs.fidma );
                  printpline( sprintf( 'self.xscale = arrset(self.xscale,iv,%s)', getv1( f{1}, f, pbs ) ), ...
                                                                                              indlvl, bindent, pbs.fidpy );
                  printjline( sprintf( 'arrset(pb.xscale,iv,%s)',     getv1( f{1}, f, pbs ) ), ...
                                                                                              indlvl, bindent, pbs.fidjl );
                  has_xscale = 1;                  
               else
                  xsc = d2e( f{4} );
                  if ( abs( str2double( xsc ) - 1 ) > 1.e-10 && abs( str2double( xsc) ) > 1.0e-15 )
                     printmline( sprintf( 'pb.xscale(iv,1) = %s;',               xsc),        indlvl, bindent, pbs.fidma );
                     printpline( sprintf( 'self.xscale = arrset(self.xscale,iv,%s)', xsc),    indlvl, bindent, pbs.fidpy );
                     printjline( sprintf( 'arrset(pb.xscale,iv,%s)',             xsc),        indlvl, bindent, pbs.fidjl );
                     has_xscale = 1;
                  end
               end

            case '''INTEGER'''
               if ( has_intvars )
                  printmline( 'intvars(end+1) = iv;',                                         indlvl, bindent, pbs.fidma );
               else
                  printmline( 'intvars(1) = iv;',                                             indlvl, bindent, pbs.fidma );
               end
               printpline( 'intvars = intvars.append(iv)',                                    indlvl, bindent, pbs.fidpy );
               printjline( 'arrset(intvars,length(intvars)+1,iv)',                            indlvl, bindent, pbs.fidjl );
               has_intvars = 1;
            case '''ZERO-ONE'''
               if ( has-binvars )
                  printmline( 'binvars(end+1) = iv;',                                         indlvl, bindent, pbs.fidma );
               else
                  printmline( 'binvars(1) = iv;',                                             indlvl, bindent, pbs.fidma );
               end
               printpline( 'binvars = binvars.append(iv)',                                    indlvl, bindent, pbs.fidpy );
               printjline( 'arrset(binvars,length(binvars)+1,iv)',                            indlvl, bindent, pbs.fidjl );
               has_binvars = 1;
            otherwise

               %  If the groups have already been defined, check for the presence of linear variables in them.
               
               if ( grpsdef && nf > 3 )     %  Variable f{2} occurs linearly in group f{3}
                  has_A = 1;
                  gname =  s2mpjname('',f{3},pbs);
                  add2A( 'row', gname, getv1( f{1}, f, pbs ), indlvl, bindent, pbs )
                  if ( nf > 5 )            %   Also in group f{5}
                     gname = s2mpjname('',f{5},pbs);
                     add2A( 'row', gname, getv2( f ), indlvl, bindent, pbs )
                  end
               end
            end
         end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   DATA GROUPS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %
      %  This section uses the following logical indicators:
      %  has_consttaints    true iff the problem has constraints other than bounds on the variables
      %  has_leconstr       true iff the problem has <= constraints
      %  has_eqconstr       true iff the problem has == constraints
      %  has_geconstr       true iff the problem has >= constraints

      case 'DGROUPS'

         if ( isempty( f{1} ) )
            pbs.errors{end+1} = sprintf( 'ERROR in line %d:  empty field 1.', pbs.nline );
            if ( pbs.disperrors )
               disp( pbs.errors{end} )
            end
         else
            [ ~, act ] = formact( f{1} );
            [ gname, ~, explindex ]  = s2mpjname( '', f{2}, pbs);

            %  Find the flat index of the group, unless theis group is the same as the previous one (for which ig is still
            %  defined). This requires the name of the group to be the same and not involve an index which is not that of
            %  an active loop (explindex = 0), because this index may have been modified since the group ig has been
            %  defined.
            
            if ( ~strcmp( f{2}, prevgname ) || explindex )
               printmline( sprintf( '[ig,ig_] = s2mpjlib(''ii'',%s,ig_);', gname ),           indlvl, bindent, pbs.fidma );
               printpline( sprintf( '[ig,ig_,_] = s2mpj_ii(%s,ig_)'      , gname ),           indlvl, bindent, pbs.fidpy );
               printjline( sprintf( 'ig,ig_,_ = s2mpj_ii(%s,ig_)'        , gname ),           indlvl, bindent, pbs.fidjl );
               if ( getgnames )
                  printmline( sprintf( 'pbm.grnames{ig} = %s;',                     gname ),  indlvl, bindent, pbs.fidma );
                  printpline( sprintf( 'self.grnames = arrset(self.grnames,ig,%s)', gname ),  indlvl, bindent, pbs.fidpy );
                  printjline( sprintf( 'arrset(pbm.grnames,ig,%s)',                 gname ),  indlvl, bindent, pbs.fidjl );
               end
               prevgname = f{2};
               ghastype  = 0;
            end
            
            %  Find and remember the group's type.

            if ( ~ghastype )
               switch ( act )
               case 'N'
                  printmline( 'gtype{ig} = ''<>'';',                                          indlvl, bindent, pbs.fidma );
                  printpline( 'gtype = arrset(gtype,ig,''<>'')',                              indlvl, bindent, pbs.fidpy );
                  printjline( 'arrset(gtype,ig,"<>")',                                        indlvl, bindent, pbs.fidjl );
               case 'L'
                  printmline( 'gtype{ig}  = ''<='';',                                         indlvl, bindent, pbs.fidma );
                  printpline( 'gtype = arrset(gtype,ig,''<='')',                              indlvl, bindent, pbs.fidpy );
                  printjline( 'arrset(gtype,ig,"<=")',                                        indlvl, bindent, pbs.fidjl );
                  if ( getcnames )
                     printmline( sprintf( 'cnames{ig} = %s;',              gname ),           indlvl, bindent, pbs.fidma );
                     printpline( sprintf( 'cnames = arrset(cnames,ig,%s)', gname ),           indlvl, bindent, pbs.fidpy );
                     printjline( sprintf( 'arrset(pb.cnames,ig,%s)',       gname ),           indlvl, bindent, pbs.fidjl );
                  end
                  has_constraints = 1;
                  has_leconstr    = 1;
                  consineq{end+1} = gname;
               case 'E'
                  printmline( 'gtype{ig}  = ''=='';',                                         indlvl, bindent, pbs.fidma );
                  printpline( 'gtype = arrset(gtype,ig,''=='')',                              indlvl, bindent, pbs.fidpy );
                  printjline( 'arrset(gtype,ig,"==")',                                        indlvl, bindent, pbs.fidjl );
                  if ( getcnames )
                     printmline( sprintf( 'cnames{ig} = %s;',              gname ),           indlvl, bindent, pbs.fidma );
                     printpline( sprintf( 'cnames = arrset(cnames,ig,%s)', gname ),           indlvl, bindent, pbs.fidpy );
                     printjline( sprintf( 'arrset(pb.cnames,ig,%s)',       gname ),           indlvl, bindent, pbs.fidjl );
                  end
                  has_constraints = 1;
                  has_eqconstr    = 1;
               case 'G'
                  printmline( 'gtype{ig}  = ''>='';',                                         indlvl, bindent, pbs.fidma );
                  printpline( 'gtype = arrset(gtype,ig,''>='')',                              indlvl, bindent, pbs.fidpy );
                  printjline( 'arrset(gtype,ig,">=")',                                        indlvl, bindent, pbs.fidjl );
                  consineq{end+1} = gname;
                  if ( getcnames )
                     printmline( sprintf( 'cnames{ig} = %s;',              gname ),           indlvl, bindent, pbs.fidma );
                     printpline( sprintf( 'cnames = arrset(cnames,ig,%s)', gname ),           indlvl, bindent, pbs.fidpy );
                     printjline( sprintf( 'arrset(pb.cnames,ig,%s)',       gname ),           indlvl, bindent, pbs.fidjl );
                  end
                  has_constraints = 1;
                  has_geconstr    = 1;
               case 'D'
                  pbs.errors{end+1} = ' *** UNSUPPORTED: the D data lines in the data GROUP section is not supported ***';
                  if ( pbs.disperrors )
                     disp( pbs.errors{end} )
                  end
                  errors = pbs.errors;
                  exitc  = length( errors );
                  return
               otherwise
                  pbs.errors{end+1} = sprintf( 'ERROR in line %d: unknown field 1.', pbs.nline );
                  if ( pbs.disperrors )
                     disp( pbs.errors{end} )
                  end
               end
               ghastype = 1;
            end
            
            %   If variables have already been defined, check for the occurence of groups' linear terms.
            
            if ( varsdef )
         
               %  Find the group's linear part, if any, or the group's scaling

               if ( nf > 2 )
                  switch ( f{3} )
                  case '''SCALE'''      %   Scaling
                      if ( ~isempty( f{1} ) && f{1}(1) == 'Z' )
                         printmline( sprintf( 'pbm.gscale(ig,1) = %s;', getv1( f{1}, f, pbs ) ), ...
                                                                                              indlvl, bindent, pbs.fidma );
                         printpline( sprintf( 'self.gscale = arrset(self.gscale,ig,float(%s))',  ...
                                              getv1( f{1}, f, pbs ) ),                        indlvl, bindent, pbs.fidpy );
                         printjline( sprintf( 'arrset(pbm.gscale,ig,Float64(%s))', getv1( f{1}, f, pbs ) ), ...
                                                                                              indlvl, bindent, pbs.fidjl );
                      else
                         gsc = d2e( f{4} );
                         if ( abs( str2double( gsc ) - 1 ) > 1.e-10 )
                            printmline( sprintf( 'pbm.gscale(ig,1) = %s;',  gsc ),            indlvl, bindent, pbs.fidma );
                            printpline( sprintf( 'self.gscale = arrset(self.gscale,ig,float(%s))', gsc ), ...
                                                                                              indlvl, bindent, pbs.fidpy );
                            printjline( sprintf( 'arrset(pbm.gscale,ig,Float64(%s))', gsc ),  indlvl, bindent, pbs.fidjl );
                         end
                      end
                  otherwise             %  Variable f{3} occurs linearly in group f{2} of index ig
                     has_A = 1;
                     vname = s2mpjname('',f{3},pbs);
                     add2A( 'col', vname, getv1( f{1}, f, pbs ), indlvl, bindent, pbs )
                  end
                  if ( nf > 5 )         %   Also variable f{5}
                     vname = s2mpjname('',f{5},pbs);
                     add2A( 'col', vname, getv2( f ), indlvl, bindent, pbs )
                  end
               end
            elseif ( nf > 2 && strcmp( f{3}, '''SCALE''' ) )
               printmline( sprintf( 'pbm.gscale(ig,1) = %s;',                getv1(f{1},f,pbs) ), ...
                                                                                              indlvl, bindent, pbs.fidma );
               printpline( sprintf( 'self.gscale = arrset(self.gscale,ig,float(%s))', getv1(f{1},f,pbs) ), ...
                                                                                              indlvl, bindent, pbs.fidpy );
               printjline( sprintf( 'arrset(pbm.gscale,ig,Float64(%s))',              getv1(f{1},f,pbs) ), ...
                                                                                              indlvl, bindent, pbs.fidjl );
            end
         end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  CONSTANTS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      case 'CONSTANTS'

         %  If requested, write the definitions for alternate sets of constants as comments in the output file.
         
         if ( ~isempty( cname ) )
            if ( ~strcmp( f{2}, cname ) )
               if ( writealtsets )
                  printcmline( sprintf( '%s', line(2:end) ), bindent, pbs.fidma );
                  if( ~prevlineispass )
                     printpline(  sprintf( 'pass # %s', line(2:end) ),                        indlvl, bindent, pbs.fidpy );
                     prevlineispass = 1;
                  end
                  printcjline( sprintf( '%s', line(2:end) ),                                          bindent, pbs.fidjl );
               else
                  if ( ~prevlineispass )
                     printpline( 'pass',                                                      indlvl, bindent, pbs.fidpy );
                     prevlineispass = 1;
                  end
               end
               continue
            end
         else
            cname = f{2};
         end

         %  Define the group's constants.
         
         form = formact( f{1} );
         if ( strcmp( f{3}, '''DEFAULT''' ) )
            switch( pbs.lang )
            case 'matlab'
               printmline( '%%%%%%%%%%%%%%%%%%%  CONSTANTS %%%%%%%%%%%%%%%%%%%',              indlvl, bindent, pbs.fidma );
               printmline( sprintf( 'pbm.gconst = %s*ones(ngrp,1);', getv1(f{1},f,pbs) ),     indlvl, bindent, pbs.fidma );
            case 'python'
               printpline( '#%%%%%%%%%%%%%%%%%%  CONSTANTS %%%%%%%%%%%%%%%%%%%',              indlvl, bindent, pbs.fidpy );
               printpline( sprintf( 'self.gconst = np.full((ngrp,1),%s)', getv1(f{1},f,pbs) ),indlvl, bindent, pbs.fidpy );
            case 'julia'
               printjline( '#%%%%%%%%%%%%%%%%%%  CONSTANTS %%%%%%%%%%%%%%%%%%%',              indlvl, bindent, pbs.fidjl );
               printjline( sprintf( 'pbm.gconst = fill(%s,ngrp)', getv1(f{1},f,pbs) ),        indlvl, bindent, pbs.fidjl );
            end
            pending    = {};
            pendingkey = '';
         else
            gname = s2mpjname( '', f{3}, pbs );
            printmline( sprintf( 'pbm.gconst(ig_(%s)) = %s;', gname, getv1(f{1},f,pbs) ),     indlvl, bindent, pbs.fidma );
            printpline( sprintf( 'self.gconst = arrset(self.gconst,ig_[%s],float(%s))', gname, getv1(f{1},f,pbs) ), ...
                                                                                              indlvl, bindent, pbs.fidpy );
            printjline( sprintf( 'pbm.gconst[ig_[%s]] = Float64(%s)', gname, getv1(f{1},f,pbs) ),      ...
                                                                                              indlvl, bindent, pbs.fidjl );
            if ( nf > 5 )
               gname = s2mpjname( '', f{5}, pbs );
               printmline( sprintf( 'pbm.gconst(ig_(%s)) = %s;', gname, getv2(f) ),           indlvl, bindent, pbs.fidma );
               printpline( sprintf( 'self.gconst = arrset(self.gconst,ig_[%s],float(%s))', gname, getv2(f) ), ...
                                                                                              indlvl, bindent, pbs.fidpy );
               printjline( sprintf( 'pbm.gconst[ig_[%s]] = Float64(%s)', gname, getv2(f) ),   indlvl, bindent, pbs.fidjl ); 
            end
         end
         has_constants = 1;
         

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  RANGES  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      %  Remember that RANGES are only meaningful for inequality constraints.  They are thus ignored in the (unlikely)
      %  case where they would be specified for groups occuring in the objective function or equality constraints.
      
      case 'RANGES'

         %  If requested, write the definitions for alternate sets of ranges as comments in the output file.
         
         if ( ~isempty( rname ) )
            if ( ~strcmp( f{2}, rname ) )
               if ( writealtsets )
                  printcmline( sprintf( '%s', line(2:end) ),                                          bindent, pbs.fidma );
                  if ( ~prevlineispass )
                     printpline(  sprintf( 'pass # %s', line(2:end) ),                        indlvl, bindent, pbs.fidpy );
                     prevlineispass = 1;
                  end
                  printcjline( sprintf( '%s', line(2:end) ),                                          bindent, pbs.fidjl );
               else
                  if ( ~prevlineispass )
                     printpline( 'pass',                                                      indlvl, bindent, pbs.fidpy );
                     prevlineispass = 1;
                  end
               end
               continue
            end
         else
            rname = f{2};
         end

         %  Define the group's ranges.
         
         form = formact( f{1} );
         if ( strcmp( f{3}, '''DEFAULT''' ) )
            switch ( pbs.lang )
            case 'matlab'
               printmline( '%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%' ,             indlvl, bindent, pbs.fidma );
               printmline( sprintf( 'grange(legrps,1) = %s;', getv1(f{1},f,pbs) ),            indlvl, bindent, pbs.fidma );
               printmline( sprintf( 'grange(gegrps,1) = %s',  getv1(f{1},f,pbs) ),            indlvl, bindent, pbs.fidma );
            case 'python'
               printpline( '#%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%' ,             indlvl, bindent, pbs.fidpy );
               printpline( sprintf( 'grange = arrset(grange,legrps,np.full(self.nle,float(%s)))',   ...
                           getv1( f{1}, f, pbs ) ),                                           indlvl, bindent, pbs.fidpy );
               printpline( sprintf( 'grange = arrset(grange,gegrps,np.full(self.nge,float(%s)))',   ...
                           getv1( f{1}, f, pbs ) ),                                           indlvl, bindent, pbs.fidpy );
            case 'julia'
               printjline( '#%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%' ,             indlvl, bindent, pbs.fidjl );
               printjline( sprintf( 'grange[legrps] = fill(%s,pb.nle)', getv1( f{1}, f, pbs ) ),  ...
                                                                                              indlvl, bindent, pbs.fidjl );
               printjline( sprintf( 'grange[gegrps] = fill(%s,pb.nge)', getv1( f{1}, f, pbs ) ),  ...
                                                                                              indlvl, bindent, pbs.fidjl );
            end
            pending    = {};
            pendingkey = '';
         else
            gname = s2mpjname( '', f{3}, pbs );
            if ( ismember( gname, consineq ) )  %  Check this group is an inequality constraint
               printmline( sprintf( 'grange(ig_(%s)) = %s;', gname, getv1(f{1},f,pbs) ),      indlvl, bindent, pbs.fidma );
               printpline( sprintf( 'grange = arrset(grange,ig_[%s],float(%s))',  gname, getv1( f{1}, f, pbs ) ), ...
                                                                                              indlvl, bindent, pbs.fidpy );
               printjline( sprintf( 'arrset(grange,ig_[%s],Float64(%s))',  gname, getv1( f{1}, f, pbs ) ),     ...
                                                                                              indlvl, bindent, pbs.fidjl );
               if ( nf > 5 )
                  gname = s2mpjname( '', f{5}, pbs );
                  if ( ismember( gname, consineq ) )  %  Check this group is an inequality constraint
                     printmline( sprintf( 'grange(ig_(%s)) = %s;',              gname, getv2(f) ), ...
                                                                                              indlvl, bindent, pbs.fidma );
                     printpline( sprintf( 'grange = arrset(grange,ig_[%s],float(%s))', gname, getv2(f) ), ...
                                                                                              indlvl, bindent, pbs.fidpy );
                     printjline( sprintf( 'arrset(grange,ig_[%s],Float64(%s))', gname, getv2(f) ),     ...
                                                                                              indlvl, bindent, pbs.fidjl );
                  end
               end
            end
         end
         has_ranges = 1;
         
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  BOUNDS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      %  The special cas of MPS bounds (where upper and lower bounds are simulaneously changed to conform the the MPS
      %  standard) is detected and handled in this section. This section uses various logical indicators:
      %  has_lower/has_upper      indicates whether or not teh problem involves explicit lower/upper bounds on the
      %                           variables;
      %  has_xlowdef/has_xuppdef  indicates whether or not an explicit default value has been defined for lower/upper
      %                           bounds.
      %  MPSbounds                indicates whether or not the default define "MPS variables", for which logical rules do
      %                           not apply in to comply with the (historical) MPS standard (see the SIF report).
      
      case 'BOUNDS'

         %  If requested, write the definitions for alternate sets of bounds as comments in the output file.
         
         if ( ~isempty( bname ) )
            if ( ~strcmp( f{2}, bname ) )
               if ( writealtsets )
                  printcmline( sprintf( '%s', line(2:end) ),                                          bindent, pbs.fidma );
                  if ( ~prevlineispass )
                     printpline(  sprintf( 'pass # %s', line(2:end) ),                        indlvl, bindent, pbs.fidpy );
                     prevlineispass = 1;
                  end
                  printcjline( sprintf( '%s', line(2:end) ),                                          bindent, pbs.fidjl );
               else
                  if ( ~prevlineispass )
                     printpline( 'pass',                                                      indlvl, bindent, pbs.fidpy );
                     prevlineispass = 1;
                  end
               end
               continue
            end
         else
            bname = f{2};
         end

        %  The defaults

         if ( strcmp( f{3}, '''DEFAULT''' ) )
            if ( ~has_xlowdef && ~has_xuppdef )
               printmline( '%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%',              indlvl, bindent, pbs.fidma );
               printpline( '#%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%',              indlvl, bindent, pbs.fidpy );
               printjline( '#%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%',              indlvl, bindent, pbs.fidjl );
            end
            switch ( f{1} )
            case { 'FR', 'XR', 'ZR' }
               switch ( pbs.lang )
               case 'matlab'
                  printmline( 'pb.xlower = -Inf*ones(pb.n,1);',                               indlvl, bindent, pbs.fidma );
                  printmline( 'pb.xupper = +Inf*ones(pb.n,1);',                               indlvl, bindent, pbs.fidma );
               case 'python'
                  printpline( 'self.xlower = np.full((self.n,1),-float(''Inf''))',            indlvl, bindent, pbs.fidpy );
                  printpline( 'self.xupper = np.full((self.n,1),+float(''Inf''))',            indlvl, bindent, pbs.fidpy );
               case 'julia'
                  printjline( 'pb.xlower = -1*fill(Inf,pb.n)',                                indlvl, bindent, pbs.fidjl );
                  printjline( 'pb.xupper =    fill(Inf,pb.n)',                                indlvl, bindent, pbs.fidjl );
               end
               has_xlowdef = 1;
               has_xuppdef = 1;
               pending     = {};
               pendingkey  = '';
            case { 'MI', 'XM' }
               printmline( 'pb.xlower = -Inf*ones(pb.n,1);',                                  indlvl, bindent, pbs.fidma );
               printpline( 'self.xlower =  np.full((self.n,1),-float(''Inf''))',              indlvl, bindent, pbs.fidpy );
               printjline( 'pb.xlower =  -1*fill(Inf,pb.n)',                                  indlvl, bindent, pbs.fidjl );
               has_xlowdef = 1;
            case { 'PL', 'XP' }
               printmline( 'pb.xupper = +Inf*ones(pb.n,1);',                                  indlvl, bindent, pbs.fidma );
               printpline( 'self.xupper = np.full((self.n,1),+float(''Inf''))',               indlvl, bindent, pbs.fidpy );
               printjline( 'pb.xupper =    fill(Inf,pb.n)',                                   indlvl, bindent, pbs.fidjl );
               MPSbounds   = MPSbounds + 1;
               has_xuppdef = 1;
            case { 'LO', 'XL', 'ZL' }
               xlowdef = getv1( f{1}, f, pbs );
               printmline( sprintf( 'pb.xlower = %s*ones(pb.n,1);', xlowdef ),                indlvl, bindent, pbs.fidma );
               printpline( sprintf( 'self.xlower = np.full((self.n,1),%s)', xlowdef ),        indlvl, bindent, pbs.fidpy );
               printjline( sprintf( 'pb.xlower = fill(%s,pb.n)', xlowdef ),                   indlvl, bindent, pbs.fidjl );
               if ( abs( str2double( f{3} ) ) < 1.e-10 )
                  MPSbounds = MPSbounds + 1;
               end
               has_xlowdef = 1;
            case { 'FX', 'XX', 'ZX' }
               xlowdef = getv1( f{1}, f, pbs );
               xuppdef = xlowdef;
               switch ( pbs.lang )
               case 'matlab'
                  printmline( sprintf( 'pb.xlower = %s*ones(pb.n,1);', xlowdef ),             indlvl, bindent, pbs.fidma );
                  printmline( sprintf( 'pb.xupper = %s*ones(pb.n,1);', xuppdef ),             indlvl, bindent, pbs.fidma );
               case 'python'
                  printpline( sprintf( 'self.xlower = np.full((self.n,1),%s)', xlowdef ),     indlvl, bindent, pbs.fidpy );
                  printpline( sprintf( 'self.xupper = np.full((self.n,1),%s)', xuppdef ),     indlvl, bindent, pbs.fidpy );
               case 'julia'
                  printjline( sprintf( 'pb.xlower = fill(%s,pb.n)', xlowdef ),                indlvl, bindent, pbs.fidjl );
                  printjline( sprintf( 'pb.xupper = fill(%s,pb.n)', xuppdef ),                indlvl, bindent, pbs.fidjl );
               end
               has_xlowdef = 1;
               has_xuppdef = 1;
            case { 'UP', 'XU', 'ZU' }
               xuppdef = getv1( f{1}, f, pbs );
               printmline( sprintf( 'pb.xupper = %s*ones(pb.n,1);', xuppdef ),                indlvl, bindent, pbs.fidma );
               printpline( sprintf( 'self.xupper = np.full((self.n,1),%s)', xuppdef ),        indlvl, bindent, pbs.fidpy );
               printjline( sprintf( 'pb.xupper = fill(%s,pb.n)', xuppdef ),                   indlvl, bindent, pbs.fidjl );
               has_xuppdef = 1;
            end
            pending = {};
            if ( ~has_xlowdef )
               switch ( pbs.lang )
               case 'matlab'
                  pending{1} = 'pb.xlower = zeros(pb.n,1);';
               case 'python'
                  pending{1} = 'self.xlower = np.zeros((self.n,1))';
               case 'julia'
                  pending{1} = 'pb.xlower = zeros(Float64,pb.n)';
               end
               pendingkey = 'default';
            elseif ( ~has_xuppdef )
               switch( pbs.lang )
               case 'matlab'
                  pending{1} = 'pb.xupper = +Inf*ones(pb.n,1);';
               case 'python'
                  pending{1} = 'self.xupper = np.full((self.n,1),+float(''inf''))';
               case 'julia'
                  pending{1} = 'pb.xupper = fill(Inf,pb.n)';
               end
               pendingkey = 'default';
            else
               pending    = {};
               pendingkey = '';
            end

         %  The specific bounds
         
         else

            vname = s2mpjname('',f{3},pbs);
            switch ( f{1} )
            case { 'FR', 'XR', 'ZR' }
               switch ( pbs.lang )
               case 'matlab'
                  printmline( sprintf( 'pb.xlower(ix_(%s)) = -Inf;', vname ),                 indlvl, bindent, pbs.fidma );
                  printmline( sprintf( 'pb.xupper(ix_(%s),1) = +Inf;', vname ),               indlvl, bindent, pbs.fidma );
               case 'python'
                  printpline( sprintf( 'self.xlower[ix_[%s]] = -float(''Inf'')', vname ),     indlvl, bindent, pbs.fidpy );
                  printpline( sprintf( 'self.xupper[ix_[%s]] = +float(''Inf'')', vname ),     indlvl, bindent, pbs.fidpy );
               case 'julia'
                  printjline( sprintf( 'pb.xlower[ix_[%s]] = -Inf', vname ),                  indlvl, bindent, pbs.fidjl );
                  printjline( sprintf( 'pb.xupper[ix_[%s]] = +Inf', vname ),                  indlvl, bindent, pbs.fidjl );
               end
            case { 'MI', 'XM' }  
               printmline( sprintf( 'pb.xlower(ix_(%s),1) = -Inf;', vname ),                  indlvl, bindent, pbs.fidma );
               printpline( sprintf( 'self.xlower = arrset(self.xlower,ix_[%s],-float(''Inf''))', vname ), ...
                                                                                              indlvl, bindent, pbs.fidpy );
               printjline( sprintf( 'pb.xlower[ix_[%s]] = -Inf', vname ),                     indlvl, bindent, pbs.fidjl );
               if ( MPSbounds >= 2 )
                  printmline( sprintf( 'pb.xlower(ix_(%s),1) = 0;', vname ),                  indlvl, bindent, pbs.fidma );
                  printpline( sprintf( 'self.xlower[ix_[%s]] = 0.', vname ),                  indlvl, bindent, pbs.fidpy );
                  printjline( sprintf( 'pb.xlower[ix_[%s]] = 0.',   vname ),                  indlvl, bindent, pbs.fidjl );
               end
            case { 'PL', 'XP' }
               printmline( sprintf( 'pb.xupper(ix_(%s)) = +Inf;', vname ),                    indlvl, bindent, pbs.fidma );
               printpline( sprintf( 'self.xupper[ix_[%s]] = +float(''Inf'')', vname ),        indlvl, bindent, pbs.fidpy );
               printjline( sprintf( 'pb.xupper[ix_[%s]] = +Inf', vname ),                     indlvl, bindent, pbs.fidjl );
            case { 'LO', 'XL', 'ZL' }
               printmline( sprintf( 'pb.xlower(ix_(%s),1) = %s;', vname, getv1(f{1},f,pbs) ), indlvl, bindent, pbs.fidma );
               printpline( sprintf( 'self.xlower[ix_[%s]] = %s', vname, getv1(f{1},f,pbs) ),  indlvl, bindent, pbs.fidpy );
               printjline( sprintf( 'pb.xlower[ix_[%s]] = %s', vname, getv1(f{1},f,pbs) ),    indlvl, bindent, pbs.fidjl );
            case { 'FX', 'XX', 'ZX' }
               evalstr = getv1(f{1},f,pbs);
               switch ( pbs.lang )
               case 'matlab'
                  printmline( sprintf( 'pb.xlower(ix_(%s),1) = %s;', vname, evalstr ),        indlvl, bindent, pbs.fidma );
                  printmline( sprintf( 'pb.xupper(ix_(%s),1) = %s;', vname, evalstr ),        indlvl, bindent, pbs.fidma );
               case 'python'
                  printpline( sprintf( 'self.xlower[ix_[%s]] = %s', vname, evalstr ),         indlvl, bindent, pbs.fidpy );
                  printpline( sprintf( 'self.xupper[ix_[%s]] = %s', vname, evalstr ),         indlvl, bindent, pbs.fidpy );
               case 'julia'
                  printjline( sprintf( 'pb.xlower[ix_[%s]] = %s', vname, evalstr ),           indlvl, bindent, pbs.fidjl );
                  printjline( sprintf( 'pb.xupper[ix_[%s]] = %s', vname, evalstr ),           indlvl, bindent, pbs.fidjl );
               end
            case { 'UP', 'XU' }
               printmline( sprintf( 'pb.xupper(ix_(%s)) = %s;', vname, getv1(f{1},f,pbs) ),   indlvl, bindent, pbs.fidma );
               printpline( sprintf( 'self.xupper[ix_[%s]] = %s', vname, getv1(f{1},f,pbs) ),  indlvl, bindent, pbs.fidpy );
               printjline( sprintf( 'pb.xupper[ix_[%s]] = %s',  vname, getv1(f{1},f,pbs) ),   indlvl, bindent, pbs.fidjl );
               if ( MPSbounds >= 2 && abs( str2double( f{3} ) ) < 1.e-10 )
                  printmline( sprintf( 'pb.xlower(ix_(%s),1) = -Inf;', vname ),               indlvl, bindent, pbs.fidma );
                  printpline( sprintf( 'self.xlower[ix_[%s]] = -float(''Inf'')', vname ),     indlvl, bindent, pbs.fidpy );
                  printjline( sprintf( 'pb.xlower[ix_[%s]] = -Inf', vname ),                  indlvl, bindent, pbs.fidjl );
               end
            case 'ZU'
               printmline( sprintf( 'pb.xupper(ix_(%s)) = %s;', vname, getv1(f{1},f,pbs) ),   indlvl, bindent, pbs.fidma );
               printpline( sprintf( 'self.xupper[ix_[%s]] = %s', vname, getv1(f{1},f,pbs) ),  indlvl, bindent, pbs.fidpy );
               printjline( sprintf( 'pb.xupper[ix_[%s]] = %s',  vname, getv1(f{1},f,pbs) ),   indlvl, bindent, pbs.fidjl );
               has_xupper = 1;
            end
         end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% START POINT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      case 'STARTPOINT'

         %  If requested, write the definitions for alternate sets of starting points as comments in the output file.
         
         if ( ~isempty( sname ) )
            if ( ~strcmp( f{2}, sname ) )
               if ( writealtsets )
                  printcmline( sprintf( '%s', line(2:end) ), bindent, pbs.fidma );
                  if ( ~prevlineispass )
                     printpline(  sprintf( 'pass # %s ', line(2:end) ),                       indlvl, bindent, pbs.fidpy );
                     prevlineispass = 1.
                  end
                  printcjline( sprintf( '%s', line(2:end) ), bindent, pbs.fidjl );
               else
                  if ( ~prevlineispass )
                     printpline( 'pass', indlvl, bindent, pbs.fidpy );
                     prevlineispass = 1;
                  end
               end
               continue
            end
         else
            sname = f{2};
         end
         
         %  The defaults

         if ( strcmp( f{3}, '''DEFAULT''' ) )
            printmline( '%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%',                 indlvl, bindent, pbs.fidma );
            printpline( '#%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%',                 indlvl, bindent, pbs.fidpy );
            printjline( '#%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%',                 indlvl, bindent, pbs.fidjl );
            switch ( f{1} )
            case  { 'V', 'XV', 'ZV' }
               x0def = getv1( f{1}, f, pbs );
               printmline( sprintf( 'pb.x0 = %s*ones(pb.n,1);', x0def ),                      indlvl, bindent, pbs.fidma );
               printpline( sprintf( 'self.x0 = np.full((self.n,1),float(%s))', x0def ),       indlvl, bindent, pbs.fidpy );
               printjline( sprintf( 'pb.x0 = fill(Float64(%s),pb.n)', x0def ),                indlvl, bindent, pbs.fidjl );
               has_x0def = 1;
            case  { 'X', '', 'Z' }
               x0def = getv1( f{1}, f, pbs );
               printmline( sprintf( 'pb.x0 = %s*ones(pb.n,1);', x0def ),                      indlvl, bindent, pbs.fidma );
               printpline( sprintf( 'self.x0 = np.full((self.n,1),float(%s))', x0def ),       indlvl, bindent, pbs.fidpy );
               printjline( sprintf( 'pb.x0 = fill(Float64(%s),pb.n)', x0def ),                indlvl, bindent, pbs.fidjl );
               has_x0def = 1;
               if ( has_constraints )
                  printmline( sprintf( 'pb.y0 = %s*ones(pb.m,1);', x0def ),                   indlvl, bindent, pbs.fidma );
                  printpline( sprintf( 'self.y0 = np.full((self.m,1),float(%s))', x0def ),    indlvl, bindent, pbs.fidpy );
                  printjline( sprintf( 'pb.y0 = fill(Float64(%s),pb.m)', x0def ),             indlvl, bindent, pbs.fidjl );
                  has_y0def = 1;
               end
            case { 'M', 'XM', 'ZM' }
               printmline( sprintf( 'pb.y0 = %s*ones(pb.m,1);', getv1(f{1},f,pbs) ),          indlvl, bindent, pbs.fidma );
               printpline( sprintf( 'self.y0 = np.full((self.m,1),float(%s))', getv1(f{1},f,pbs) ),  ...
                                                                                              indlvl, bindent, pbs.fidpy );
               printjline( sprintf( 'pb.y0 = fill(Float64(%s),pb.m)', getv1(f{1},f,pbs) ),    indlvl, bindent, pbs.fidjl );
               has_y0def = 1;
            end
            pending    = {};
            pendingkey = '';

         %  The specific values
         
         else
 
            switch ( f{1} )
            case { 'V', 'XV', 'ZV' }                              %  variables
               vname = s2mpjname('',f{3},pbs);
               printmline( sprintf( 'pb.x0(ix_(%s),1) = %s;',       vname, getv1(f{1},f,pbs ) ), ...
                                                                                              indlvl, bindent, pbs.fidma );
               printpline( sprintf( 'self.x0[ix_[%s]] = float(%s)', vname, getv1(f{1},f,pbs ) ), ...
                                                                                              indlvl, bindent, pbs.fidpy );
               printjline( sprintf( 'pb.x0[ix_[%s]] = Float64(%s)', vname, getv1(f{1},f,pbs ) ), ...
                                                                                              indlvl, bindent, pbs.fidjl );
               if ( nf > 5 )
                  vname = s2mpjname('',f{5},pbs);
                  printmline( sprintf( 'pb.x0(ix_(%s),1) = %s;',       vname, getv2(f) ),     indlvl, bindent, pbs.fidma );
                  printpline( sprintf( 'self.x0[ix_[%s]] = float(%s)', vname, getv2(f) ),     indlvl, bindent, pbs.fidpy );
                  printjline( sprintf( 'pb.x0[ix_[%s]] = Float64(%s)', vname, getv2(f) ),     indlvl, bindent, pbs.fidjl );
               end
            case { 'M', 'XM', 'ZM' }                              %  multipliers
               gname = s2mpjname('',f{3},pbs);
               printmline( sprintf( 'pb.y0(find(self.congrps==ig_(%s,1))) = %s;', gname, getv1(f{1},f,pbs ) ),        ...
                                                                                              indlvl, bindent, pbs.fidma );
               printpline( sprintf( 'self.y0 = arrset(self.y0,np.where(self.congrps.data==ig_[%s])[0],float(%s))',    ...
                           gname, getv1(f{1},f,pbs ) ),                                       indlvl, bindent, pbs.fidpy );
               printjline( sprintf( 'pb.y0[findall(x->x==ig[%s],pbm.congrps)] = Float64(%s)',                        ...
                           gname, getv1(f{1},f,pbs ) ),                                       indlvl, bindent, pbs.fidjl );
               if ( nf > 5 )
                  gname = s2mpjname('',f{5},pbs);
                  printmline( sprintf( 'pb.y0(findfirst(pbm.congrps==ig_(%s)),1) = %s;', gname, getv2(f) ),          ...    
                                                                                              indlvl, bindent, pbs.fidma );
                  printpline( sprintf( 'self.y0 = arrset(self.y0,np.where(self.congrps.data==ig_[%s])[0],float(%s))', ...
                              gname, getv2(f) ),                                              indlvl, bindent, pbs.fidpy );
                  printjline( sprintf( 'pb.y0[findfirst(x->x==ig_[%s],pbm.congrps)] = Float64(%s)',                  ...
                              gname, getv2(f) ),                                              indlvl, bindent, pbs.fidjl );
               end
            case { '', 'X', 'Z' }                                 %  mix of variables and multipliers
               evalstr = getv1(f{1},f,pbs);
               vname   = s2mpjname('',f{3},pbs);
               if ( has_constraints )
                  switch ( pbs.lang )
                  case 'matlab'
                     printmline( sprintf( 'if(isKey(ix_,%s))', vname ),                   indlvl,     bindent, pbs.fidma );
                     printmline( sprintf( 'pb.x0(ix_(%s),1) = %s;', vname, evalstr ),     indlvl + 1, bindent, pbs.fidma );
                     printmline( 'else',                                                  indlvl,     bindent, pbs.fidma );
                     gname = s2mpjname('',f{3},pbs);
                     printmline( sprintf( 'pb.y0(find(pbm.congrps==ig_(%s)),1) = %s;', gname, evalstr ), ... 
                                                                                          indlvl + 1, bindent, pbs.fidma );
                     printmline( 'end',                                                   indlvl,     bindent, pbs.fidma );
                  case 'python'
                     printpline( sprintf( 'if(%s in ix_):',      vname ),                 indlvl,     bindent, pbs.fidpy );
                     printpline( sprintf( 'self.x0[ix_[%s]] = float(%s)', vname, evalstr ),  ...
                                                                                          indlvl + 1, bindent, pbs.fidpy );
                     printpline( 'else:',                                                 indlvl,     bindent, pbs.fidpy );
                     gname = s2mpjname('',f{3},pbs);
                     printpline( sprintf( 'self.y0 = arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_[%s]),float(%s))', ...
                                          gname, evalstr ),                               indlvl + 1, bindent, pbs.fidpy );
                  case 'julia'
                     printjline( sprintf( 'if haskey(ix_,%s)',              vname ),      indlvl,     bindent, pbs.fidjl );
                     printjline( sprintf( 'pb.x0[ix_[%s]] = Float64(%s)',vname, evalstr ),indlvl + 1, bindent, pbs.fidjl );
                     printjline( 'else',                                                  indlvl,     bindent, pbs.fidjl );
                     gname = s2mpjname('',f{3},pbs);
                     printjline( sprintf( 'pb.y0[findfirst(x->x==ig_[%s],pbm.congrps)] = Float64(%s)', ...
                                          gname, evalstr ),                               indlvl + 1, bindent, pbs.fidjl );
                     printjline( 'end',                                                   indlvl,     bindent, pbs.fidjl );
                  end
               else
                  printmline( sprintf( 'pb.x0(ix_(%s),1) = %s;',       vname, evalstr ),  indlvl,     bindent, pbs.fidma );
                  printpline( sprintf( 'self.x0[ix_[%s]] = float(%s)', vname, evalstr ),  indlvl,     bindent, pbs.fidpy );
                  printjline( sprintf( 'pb.x0[ix_[%s]] = Float64(%s)', vname, evalstr ),  indlvl,     bindent, pbs.fidjl );
               end
               if ( nf > 5 )
                  evalstr  = getv2(f);
                  vname    = s2mpjname('',f{5},pbs);
                  if ( has_constraints )
                     switch ( pbs.lang )
                     case 'matlab'
                        printmline( sprintf( 'if(isKey(ix_,%s))', vname),                 indlvl,     bindent, pbs.fidma );
                        printmline( sprintf( 'pb.x0(ix_(%s),1) = %s;', vname, evalstr ),  indlvl + 1, bindent, pbs.fidma );
                        printmline( 'else',                                               indlvl,     bindent, pbs.fidma );
                        gname = s2mpjname('',f{5},pbs);
                        printmline( sprintf( 'pb.y0(find(pbm.congrps==ig(%s)),1) = %s;', gname, evalstr ),     ...
                                                                                          indlvl + 1, bindent, pbs.fidma );
                        printmline( 'end',                                                indlvl,     bindent, pbs.fidma );
                     case 'python'
                        printpline( sprintf( 'if(%s in ix_):', vname),                    indlvl,     bindent, pbs.fidpy );
                        printpline( sprintf( 'self.x0[ix_[%s]] = float(%s)', vname, evalstr ),  ...
                                                                                          indlvl + 1, bindent, pbs.fidpy );
                        printpline( 'else:',                                              indlvl,     bindent, pbs.fidpy );
                        gname = s2mpjname('',f{5},pbs);
                        printpline( sprintf( 'self.y0 = arrset(self.y0,np.where(self.congrps==ig_[%s])[0],float(%s))', ...
                                          gname, evalstr ),                               indlvl + 1, bindent, pbs.fidpy );
                     case 'julia'
                        printjline( sprintf( 'if haskey(ix_,%s)', vname),                 indlvl,     bindent, pbs.fidjl );
                        printjline( sprintf( 'pb.x0[ix_[%s]] = Float64(%s)', vname, evalstr ), ...
                                                                                          indlvl + 1, bindent, pbs.fidjl );
                        printjline( 'else',                                               indlvl,     bindent, pbs.fidjl );
                        gname = s2mpjname('',f{5},pbs);
                        printjline( sprintf( 'pb.y0[findfirst(x->x==ig_[%s],pbm.congrps)] = Float64(%s)',        ...
                                          gname, evalstr ),                               indlvl + 1, bindent, pbs.fidjl );
                        printjline( 'end',                                                indlvl,     bindent, pbs.fidjl );
                     end
                  else
                     printmline( sprintf( 'pb.x0(ix_(%s),1) = %s;', vname, evalstr ),     indlvl,     bindent, pbs.fidma );
                     printpline( sprintf( 'self.y0 = arrset(self.x0,ix_[%s],float(%s))', vname, evalstr ), ...
                                                                                          indlvl,     bindent, pbs.fidpy );
                     printjline( sprintf( 'pb.x0[ix_[%s]] = Float64(%s)', vname, evalstr ),  ...
                                                                                          indlvl,     bindent, pbs.fidjl );
                  end
               end
            end
         end
         has_start = 1;

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  QUADRATIC  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      case 'QUADRATIC'

          vname1 = s2mpjname( '', f{2}, pbs );
          vname2 = s2mpjname( '', f{3}, pbs );
          switch( f{1} )
          case { '', 'X' }
             add2H( vname1, vname2, getv1( f{1}, f, pbs ), indlvl, bindent, pbs );
             if ( nf > 5 )
                vname2 = s2mpjname( '', f{5}, pbs );
                add2H( vname1, vname2, getv2( f ), indlvl, bindent, pbs );
             end
          case 'Z'
             add2H( vname1, vname2, getv1( f{1}, f, pbs ), indlvl, bindent, pbs );
          end
          has_H = 1;
          
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  ELFTYPE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      %  Element types are defined both in the S2MPJ execution and in the output file.  The first is used in the analysis
      %  of the ELEMENTS (INDIVIDUALS) section while the second is need to assign ftype to elements. Note that the latter
      %  does not need to inform the output file of the existence of internal variables at this stage.
      
      case 'ELFTYPE'

         eftname =  [ 'e', nlfname( f{2} ) ];
         pos = ismemstr( eftname, elftype, 'name' );
         if ( ~pos )
            pos = length( elftype)+1;
            elftype{pos}.name = eftname;
         end

         switch ( f{1} )

         %  Elemental variables
         
         case 'EV'
            if ( ~strcmp( prevetype, f{2} ) )
               printmline( sprintf( '[it,iet_] = s2mpjlib( ''ii'', ''%s'',iet_);', eftname ), ...
                                                                                              indlvl, bindent, pbs.fidma );
               printpline( sprintf( '[it,iet_,_] = s2mpj_ii( ''%s'', iet_)', eftname ),       indlvl, bindent, pbs.fidpy );
               printjline( sprintf( 'it,iet_,_ = s2mpj_ii( "%s", iet_)' ,    eftname ),       indlvl, bindent, pbs.fidjl );
               prevetype = f{2};
            end                           
            if ( isfield( elftype{pos}, 'elvar' ) )
               posv = length( elftype{pos}.elvar ) + 1;
               elftype{pos}.elvar{posv} = f{3};
               printmline( sprintf( 'elftv{it}{%d} = ''%s'';',            posv,   f{3} ),     indlvl, bindent, pbs.fidma );
               printpline( sprintf( 'elftv = loaset(elftv,it,%d,''%s'')', posv-1, f{3} ),     indlvl, bindent, pbs.fidpy );
               printjline( sprintf( 'loaset(elftv,it,%d,"%s")',           posv,   f{3} ),     indlvl, bindent, pbs.fidjl );
            else
               elftype{pos}.elvar{1} = f{3};
               printmline( sprintf( 'elftv{it}{1} = ''%s'';',            f{3} ),              indlvl, bindent, pbs.fidma );
               printpline( sprintf( 'elftv = loaset(elftv,it,0,''%s'')', f{3} ),              indlvl, bindent, pbs.fidpy );
               printjline( sprintf( 'loaset(elftv,it,1,"%s")',           f{3} ),              indlvl, bindent, pbs.fidjl );
            end
            if ( nf > 4 )
               posv = length( elftype{pos}.elvar ) + 1;
               elftype{pos}.elvar{posv} = f{5};
               printmline( sprintf( 'elftv{it}{%d} = ''%s'';',            posv,   f{5} ),     indlvl, bindent, pbs.fidma );
               printpline( sprintf( 'elftv = loaset(elftv,it,%d,''%s'')', posv-1, f{5} ),     indlvl, bindent, pbs.fidpy );
               printjline( sprintf( 'loaset(elftv,it,%d,"%s")',           posv,   f{5} ),     indlvl, bindent, pbs.fidjl );
            end

         %  Internal variables
         
         case 'IV'
            if (  isfield( elftype{pos}, 'invar' ) )
               elftype{pos}.invar{end+1} = f{3};
            else
               elftype{pos}.invar{1} = f{3};
            end     
            if ( nf > 4 )
               elftype{pos}.invar{end+1} = f{5};
            end

         %  Element's parameters
         
         case 'EP'

            %  Have elemental parameters been defined for the current type? If not, initialize their list.
            
            if ( ~has_elpar )
               printpline( 'elftp = []',                                                      indlvl, bindent, pbs.fidpy );
               printjline( 'elftp = Vector{Vector{String}}()',                                indlvl, bindent, pbs.fidjl );
            end

            %   Is this the same element type as the previous one? If not, find its flat index and define the associated
            %   list of elemental variables.
            
            if ( ~strcmp( prevetype, f{2} ) )
               printmline( sprintf( '[it,iet_] = s2mpjlib(''ii'',''%s'',iet_);', eftname ),   indlvl, bindent, pbs.fidma );
               printpline( sprintf( '[it,iet_,_] = s2mpj_ii(''%s'',iet_)', eftname ),         indlvl, bindent, pbs.fidpy );
               printjline( sprintf( 'it,iet_,_ = s2mpj_ii("%s",iet_)', eftname ),             indlvl, bindent, pbs.fidjl ); 
               prevetype = f{2};
            end
            if ( isfield( elftype{pos}, 'elpar' ) )
               posp = length( elftype{pos}.elpar ) + 1;
               elftype{pos}.elpar{posp} = f{3};
               printmline( sprintf( 'elftp{it}{%d} = ''%s'';',    posp,   f{3} ),             indlvl, bindent, pbs.fidma );
               printpline( sprintf( 'elftp = loaset(elftp,it,%d,''%s'')', posp-1, f{3} ),     indlvl, bindent, pbs.fidpy );
               printjline( sprintf( 'loaset(elftp,it,%d,"%s")', posp, f{3} ),                 indlvl, bindent, pbs.fidjl );
            else
               elftype{pos}.elpar{1} = f{3};
               printmline( sprintf( 'elftp{it}{1} = ''%s'';',    f{3} ),                      indlvl, bindent, pbs.fidma );
               printpline( sprintf( 'elftp = loaset(elftp,it,0,''%s'')', f{3} ),              indlvl, bindent, pbs.fidpy );
               printjline( sprintf( 'loaset(elftp,it,1,"%s")', f{3} ),                        indlvl, bindent, pbs.fidjl );
            end
            if ( nf > 4 )
               posp = length( elftype{pos}.elpar ) + 1;
               elftype{pos}.elpar{posp} = f{5};
               printmline( sprintf( 'elftp{it}{%d} = ''%s'';',    posp,   f{5} ),             indlvl, bindent, pbs.fidma );
               printpline( sprintf( 'elftp = loaset(elftp,it,%d,''%s'')', posp-1, f{5} ),     indlvl, bindent, pbs.fidpy );
               printjline( sprintf( 'loaset(elftp,it,%d,"%s")', posp, f{5} ),                 indlvl, bindent, pbs.fidjl );
            end
            has_elpar = 1;
         end


      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  ELUSES  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      %  The complexity of this section is caused by the fact that the SIF standard does not require the ftype of an
      %  element to explicitly defined before elemental variables or parameters are assigned to it other than by a
      %  'DEFAULT' statement.  However, the total number of elements is unknown when a 'DEFAULT' line is read,  preventing
      %  a blanket assignment of ftype (as is the case for group ftypes). Thus, if a new element is met, its type must be
      %  defined, either in an explicit T or XT line or by its 'DEFAULT' value, before its is actually used. In the latter
      %  case, a (non-empty) default must have been defined (unless there is an error in the SIF file). Unfortunately,
      %  this requires checking if the considered element is new (newelt=1) in the output file whenever a 'DEFAULT' type
      %  has been defined. If the element is not new (newelt=0), then its elftype is already known and is retrieved from
      %  pbm.elftype.
      
      case 'ELUSES'

         [ ~, act ] = formact( f{1} );
         if ( strcmp( f{2}, '''DEFAULT''' ) )
            switch ( act )
            case 'T'
               defelftype = [ 'e', nlfname(f{3}) ];
            end
            continue
         end

         %  Is this element the same as the previous one? If not, find its flat index, and possibly store its name. This
         %  is also the case if the element's name depends on an index which is not that of an active loop (such an index
         %  may be explicitly modified between the previous element declaration and the previous one. This detected in
         %  s2mpjname and the flag explindex = 1 is then returned.

         [ ename,~,explindex] = s2mpjname( '', f{2}, pbs);
         if ( ~strcmp( f{2}, prevename ) || explindex ) 
            printmline( sprintf( 'ename = %s;', ename ),                                      indlvl, bindent, pbs.fidma );
            printpline( sprintf( 'ename = %s',  ename ),                                      indlvl, bindent, pbs.fidpy );
            printjline( sprintf( 'ename = %s',  ename ),                                      indlvl, bindent, pbs.fidjl );
            if ( ~isempty( defelftype ) )
               printmline( '[ie,ie_,newelt] = s2mpjlib(''ii'',ename,ie_);',                   indlvl, bindent, pbs.fidma );
               printpline( '[ie,ie_,newelt] = s2mpj_ii(ename,ie_)',                           indlvl, bindent, pbs.fidpy );
               printjline( 'ie,ie_,newelt = s2mpj_ii(ename,ie_)',                             indlvl, bindent, pbs.fidjl );
            else
               printmline( '[ie,ie_] = s2mpjlib(''ii'',ename,ie_);',                          indlvl, bindent, pbs.fidma );
               printpline( '[ie,ie_,_] = s2mpj_ii(ename,ie_)',                                indlvl, bindent, pbs.fidpy );
               printjline( 'ie,ie_,_  = s2mpj_ii(ename,ie_)',                                 indlvl, bindent, pbs.fidjl );
            end
            if ( getenames )        
               printmline( 'pbm.enames{ie} = ename;',                                         indlvl, bindent, pbs.fidma );
               printpline( 'self.enames = arrset(self.enames,ie,ename)',                      indlvl, bindent, pbs.fidpy );
               printjline( 'arrset(pbm.enames,ie,ename)',                                     indlvl, bindent, pbs.fidjl );
            end
            prevename  = f{2};                 %  Remember the SIF name of the current element
            ehas_eltyp = 0;                    %  Remember that the element type hasn't beeen defined yet
            ehas_elpar = 0;                    %  The element doesn't yet contain elemental parameters
         end
         
         switch ( act )

         %  Set the element's type.
         
         case 'T'

            tname = [ 'e', nlfname( f{3} ) ];  % the elftype name
            switch ( pbs.lang )
            case 'matlab'
               printmline( sprintf( 'pbm.elftype{ie} = ''%s'';',    tname ),                  indlvl, bindent, pbs.fidma );
               printmline( sprintf( 'ielftype(ie) = iet_(''%s'');', tname ),                  indlvl, bindent, pbs.fidma );
            case 'python'
               printpline( sprintf( 'self.elftype = arrset(self.elftype,ie,''%s'')', tname ), indlvl, bindent, pbs.fidpy );
               printpline( sprintf( 'ielftype = arrset(ielftype,ie,iet_["%s"])', tname ),     indlvl, bindent, pbs.fidpy );
            case 'julia'
               printjline( sprintf( 'arrset(pbm.elftype,ie,"%s")',    tname ),                indlvl, bindent, pbs.fidjl );
               printjline( sprintf( 'arrset(ielftype,ie,iet_["%s"])', tname ),                indlvl, bindent, pbs.fidjl );
            end
            ehas_eltyp = 1;                    %  The element type has now been defined
            
         %  Define the elemental variables from the problem's variables.
         
         case 'V'

            %  Has the element type been defined? If not, define it and identify its position in the list of element types
            %  and the associated list of elemental variables.

            if ( ~ehas_eltyp )
               if ( ~isempty( defelftype ) )
                  switch ( pbs.lang )
                  case 'matlab'
                     printmline( 'if(newelt)',                                                indlvl, bindent, pbs.fidma );
                     printmline( sprintf( 'pbm.elftype{ie} = ''%s'';',     defelftype ),  indlvl + 1, bindent, pbs.fidma );
                     printmline( sprintf( 'ielftype(ie) = iet_(''%s'');',  defelftype ),  indlvl + 1, bindent, pbs.fidma );
                     printmline( 'end',                                                       indlvl, bindent, pbs.fidma );
                  case 'python'
                     printpline( 'if newelt:',                                                indlvl, bindent, pbs.fidpy );
                     printpline( sprintf( 'self.elftype = arrset(self.elftype,ie,''%s'')', defelftype ), ...
                                                                                          indlvl + 1, bindent, pbs.fidpy );
                     printpline( sprintf( 'ielftype = arrset(ielftype,ie,iet_[''%s''])', defelftype ),   ...
                                                                                          indlvl + 1, bindent, pbs.fidpy );
                  case 'julia'
                     printjline( 'if newelt > 0',                                             indlvl, bindent, pbs.fidjl );
                     printjline( sprintf( 'arrset(pbm.elftype,ie,"%s")',    defelftype ), indlvl + 1, bindent, pbs.fidjl );
                     printjline( sprintf( 'arrset(ielftype,ie,iet_["%s"])', defelftype ), indlvl + 1, bindent, pbs.fidjl );
                     printjline( 'end',                                                       indlvl, bindent, pbs.fidjl );
                  end
                  ehas_eltyp = 1;               %  The element type has now been defined
               end
            end
             
            %  Assign the elemental variable. Note that vname may be the name of a 'nonlinear' variable not declared in
            %  the VARIABLES section, in which case vname must be added to the variables' dictionary ix_, with bounds,
            %  start point and types defined by their default settings (it is too late to define problem-specific values).
            %  This task is performed by s2mpj_nlx,( ...) and the different calls are needed to reflect the fact that these
            %  defaults may have been set or not (the default variable type is always defined to be 'r').
            
            if ( ~has_x0def )
                printpline( 'self.x0 = np.zeros((self.n,1))',                                 indlvl, bindent, pbs.fidpy );
                has_x0def = 1;
            end
            vname = s2mpjname('',f{5},pbs);
            printmline( sprintf( 'vname = %s;', vname ),                                      indlvl, bindent, pbs.fidma );
            printpline( sprintf( 'vname = %s',  vname ),                                      indlvl, bindent, pbs.fidpy );
            printjline( sprintf( 'vname = %s',  vname ),                                      indlvl, bindent, pbs.fidjl );
            if     (  isempty( xlowdef ) &&  isempty( xuppdef ) &&  isempty( x0def ) )
               printmline( sprintf( '[iv,ix_,pb] = s2mpjlib(''nlx'',vname,ix_,pb,%d,[],[],[]);',               ...
                           getxnames ),                                                       indlvl, bindent, pbs.fidma );
               printpline( sprintf( '[iv,ix_] = s2mpj_nlx(self,vname,ix_,%d,None,None,None)',                  ...
                           getxnames ),                                                       indlvl, bindent, pbs.fidpy );
               printjline( sprintf( 'iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,%d,nothing,nothing,nothing)',          ...
                           getxnames ),                                                       indlvl, bindent, pbs.fidjl );
            elseif ( ~isempty( xlowdef ) &&  isempty( xuppdef ) &&  isempty( x0def ) )
               printmline( sprintf( '[iv,ix_,pb] = s2mpjlib(''nlx'',vname,ix_,pb,%d,%s,[],[]);',               ...
                           getxnames, xlowdef ),                                              indlvl, bindent, pbs.fidma );
               printpline( sprintf( '[iv,ix_] = s2mpj_nlx(self, vname,ix_,%d,float(%s),None,None)',            ...
                           getxnames, xlowdef ),                                              indlvl, bindent, pbs.fidpy );
               printjline( sprintf( 'iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,%d,Float64(%s),nothing,nothing)',      ...
                           getxnames, xlowdef ),                                              indlvl, bindent, pbs.fidjl );
            elseif (  isempty( xlowdef ) && ~isempty( xuppdef ) &&  isempty( x0def ) )
               printmline( sprintf( '[iv,ix_,pb] = s2mpjlib(''nlx'',vname,ix_,pb,%d,[],%s,[]);',               ...
                           getxnames, xuppdef ),                                              indlvl, bindent, pbs.fidma );
               printpline( sprintf( '[iv,ix_] = s2mpj_nlx(self,vname,ix_,%d,None,float(%s),None)',             ...
                           getxnames, xuppdef ),                                              indlvl, bindent, pbs.fidpy );
               printjline( sprintf( 'iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,%d,nothing,Float64(%s),nothing)',      ...
                           getxnames, xuppdef ),                                              indlvl, bindent, pbs.fidjl );
            elseif (  isempty( xlowdef ) &&  isempty( xuppdef ) && ~isempty( x0def ) )
               printmline( sprintf( '[iv,ix_,pb] = s2mpjlib(''nlx'',vname,ix_,pb,%d,[],[],%s);',               ...
                           getxnames, x0def ),                                                indlvl, bindent, pbs.fidma );
               printpline( sprintf( '[iv,ix_] = s2mpj_nlx(self,vname,ix_,%d,None,None,float(%s))',             ...
                           getxnames, x0def ),                                                indlvl, bindent, pbs.fidpy );
               printjline( sprintf( 'iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,%d,nothing,nothing,Float64(%s))',      ...
                           getxnames, x0def ),                                                indlvl, bindent, pbs.fidjl );
            elseif ( ~isempty( xlowdef ) && ~isempty( xuppdef ) &&  isempty( x0def ) )
               printmline( sprintf( '[iv,ix_,pb] = s2mpjlib(''nlx'',vname,ix_,pb,%d,%s,%s,[]);',               ...
                           getxnames, xlowdef, xuppdef ),                                     indlvl, bindent, pbs.fidma );
               printpline( sprintf( '[iv,ix_] = s2mpj_nlx(self,vname,ix_,%d,float(%s),float(%s),None)',        ...
                           getxnames, xlowdef, xuppdef ),                                     indlvl, bindent, pbs.fidpy );
               printjline( sprintf( 'iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,%d,Float64(%s),Float64(%s),nothing)',           ...
                           getxnames, xlowdef, xuppdef ),                                     indlvl, bindent, pbs.fidjl );
            elseif ( ~isempty( xlowdef ) &&  isempty( xuppdef ) && ~isempty( x0def ) )
               printmline( sprintf( '[iv,ix_,pb] = s2mpjlib(''nlx'',vname,ix_,pb,%d,%s,[],%s);',               ...
                           getxnames, xlowdef, x0def ),                                       indlvl, bindent, pbs.fidma );
               printpline( sprintf( '[iv,ix_] = s2mpj_nlx(self,vname,ix_,%d,float(%s),None,float(%s))',            ...
                           getxnames, xlowdef, x0def ),                                       indlvl, bindent, pbs.fidpy );
               printjline( sprintf( 'iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,%d,Float64(%s),nothing,Float64(%s))',           ...
                           getxnames, xlowdef, x0def ),                                       indlvl, bindent, pbs.fidjl );
            elseif (  isempty( xlowdef ) && ~isempty( xuppdef ) && ~isempty( x0def ) )
               printmline( sprintf( '[iv,ix_,pb] = s2mpjlib(''nlx'',vname,ix_,pb,%d,[],%s,%s);',               ...
                           getxnames, xuppdef, x0def ),                                       indlvl, bindent, pbs.fidma );
               printpline( sprintf( '[iv,ix_] = s2mpj_nlx(self,vname,ix_,%d,None,float(%s),float(%s))',            ... 
                           getxnames, xuppdef, x0def ),                                       indlvl, bindent, pbs.fidpy );
               printjline( sprintf( 'iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,%d,nothing,Float64(%s),Float64(%s))',  ... 
                           getxnames, xuppdef, x0def ),                                       indlvl, bindent, pbs.fidjl );
            elseif ( ~isempty( xlowdef ) && ~isempty( xuppdef ) && ~isempty( x0def ) )
               printmline( sprintf( '[iv,ix_,pb] = s2mpjlib(''nlx'',vname,ix_,pb,%d,%s,%s,%s);', ...
                           getxnames, xlowdef, xuppdef, x0def ),                              indlvl, bindent, pbs.fidma );
               printpline( sprintf( '[iv,ix_] = s2mpj_nlx(self,vname,ix_,%d,float(%s),float(%s),float(%s))',   ...
                           getxnames, xlowdef, xuppdef, x0def ),                              indlvl, bindent, pbs.fidpy );
               printjline( sprintf( 'iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,%d,Float64(%s),Float64(%s),Float64(%s))',   ...
                           getxnames, xlowdef, xuppdef, x0def ),                              indlvl, bindent, pbs.fidjl );
            end

            %  Add the (new) variable at the correct position in pbm.elvar.
            
            switch ( pbs.lang )
            case 'matlab'
               printmline( sprintf( 'posev = find(strcmp(''%s'',elftv{ielftype(ie)}));',  nlfname( f{3} ) ),   ...
                                                                                              indlvl, bindent, pbs.fidma );
               printmline( 'pbm.elvar{ie}(posev) = iv;',                                      indlvl, bindent, pbs.fidma );
            case 'python'
               printpline( sprintf( 'posev = np.where(elftv[ielftype[ie]]==''%s'')[0]', nlfname( f{3} ) ), ...
                                                                                              indlvl, bindent, pbs.fidpy );
               printpline( 'self.elvar = loaset(self.elvar,ie,posev[0],iv)',                  indlvl, bindent, pbs.fidpy );
            case 'julia'
               printpline( sprintf( 'posev = findfirst(x->x=="%s",elftv[ielftype[ie]])', nlfname( f{3} ) ),    ...
                                                                                              indlvl, bindent, pbs.fidjl );
               printjline( 'loaset(pbm.elvar,ie,posev,iv)',                                   indlvl, bindent, pbs.fidjl );
            end

         %  Element's parameters
         
         case 'P'

            %  Has the element type been defined? If not, define it and identify its position in the list of element types.
            
            if ( ~ehas_eltyp )
               if ( ~isempty( defelftype ) )
                  switch ( pbs.lang )
                  case 'matlab'
                     printmline( 'if(newelt)',                                                indlvl, bindent, pbs.fidma );
                     printmline( sprintf( 'pbm.elftype{ie} = ''%s'';',    defelftype ),   indlvl + 1, bindent, pbs.fidma );
                     printmline( sprintf( 'ielftype(ie) = iet_(''%s'');', defelftype ),   indlvl + 1, bindent, pbs.fidma );
                     printmline( 'end',                                                       indlvl, bindent, pbs.fidma );
                  case 'python'
                     printpline( 'if newelt:',                                                indlvl, bindent, pbs.fidpy );
                     printpline( sprintf( 'self.elftype = arrset(self.elftype,ie,''%s'')',  defelftype ), ...
                                                                                          indlvl + 1, bindent, pbs.fidpy );
                     printpline( sprintf( 'ielftype = arrset(ielftype,ie,iet_[''%s''])', defelftype ),    ...
                                                                                          indlvl + 1, bindent, pbs.fidpy );
                  case 'julia'
                     printjline('if newelt > 0',                                              indlvl, bindent, pbs.fidjl );
                     printjline( sprintf( 'arrset( pbm.elftype,ie,"%s")',   defelftype ), indlvl + 1, bindent, pbs.fidjl );
                     printjline( sprintf( 'arrset(ielftype,ie,iet_["%s"])', defelftype ), indlvl + 1, bindent, pbs.fidjl );
                     printjline('end',                                                        indlvl, bindent, pbs.fidjl );
                  end
                  ehas_eltyp = 1;                  %  The element type is now defined
               end
            end

            %  Have elemental parameter(s) been defined for the current element yet? If not extract the list of elemental
            %  parameters corresponding to its type.

            %  Assign the elemental parameter(s).

            switch ( pbs.lang )
            case 'matlab'
               printmline( sprintf( '[~,posep] = ismember(''%s'',elftp{ielftype(ie)});', f{3} ),          ...
                                                                                              indlvl, bindent, pbs.fidma );
               printmline( sprintf( 'pbm.elpar{ie}(posep) = %s;', getv1r(f{1},f,pbs) ),       indlvl, bindent, pbs.fidma );
            case 'python'
               printpline( sprintf( 'posep = np.where(elftp[ielftype[ie]]==''%s'')[0]', f{3} ),       ...
                                                                                              indlvl, bindent, pbs.fidpy );
               printpline( sprintf( 'self.elpar = loaset(self.elpar,ie,posep[0],float(%s))', getv1r(f{1},f,pbs) ), ...
                                                                                              indlvl, bindent, pbs.fidpy );
            case 'julia'
               printjline( sprintf( 'posep = findfirst(x->x=="%s",elftp[ielftype[ie]])', f{3} ),          ...
                                                                                              indlvl, bindent, pbs.fidjl );
               printjline( sprintf( 'loaset(pbm.elpar,ie,posep,Float64(%s))', getv1r(f{1},f,pbs) ),       ...
                                                                                              indlvl, bindent, pbs.fidjl );
            end
            if ( nf > 5 )
               switch ( pbs.lang )
               case 'matlab'
                  printmline( sprintf( '[~,posep] = ismember(''%s'',elftp{ielftype(ie)});', f{5} ),     ...
                                                                                              indlvl, bindent, pbs.fidma );
                  printmline( sprintf( 'pbm.elpar{ie}(posep) = %s;', getv2r( f ) ),           indlvl, bindent, pbs.fidma );
               case 'python'
                  printpline( sprintf( 'posep = np.where(elftp[ielftype[ie]]==''%s'')[0]', f{5} )', ...
                                                                                              indlvl, bindent, pbs.fidpy );
                  printpline( sprintf( 'loaset(self.elpar,ie,posep[0],float(%s))', getv2r( f ) ),        ...
                                                                                              indlvl, bindent, pbs.fidpy );
               case 'julia'
                  printjline( sprintf( 'posep = findfirst(x->x=="%s",elftp[ielftype[ie]])', f{5} )',    ...
                                                                                              indlvl, bindent, pbs.fidjl );
                  printjline( sprintf( 'loaset(pbm.elpar,ie,posep,Float64(%s))', getv2r( f ) ),         ...
                                                                                              indlvl, bindent, pbs.fidjl );
               end
            end
            
         end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  GRFTYPE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      %  As for ELFTYPE, assigning group parameters correctly requires storing the
      %  names and associated list of parameters for the grtypes in the output file.
      %  Since each group only has a single variable, this is not necessary for group
      %  variables.
      
      case 'GRFTYPE'

         %  The default

         grftname = [ 'g', nlfname( f{2} ) ];
         pos = ismemstr( grftname, grftype, 'name' );
         if ( ~pos )
            pos = length( grftype)+1;
            grftype{pos}.name = grftname;
         end
         printmline( sprintf( '[it,igt_] = s2mpjlib(''ii'',''%s'',igt_);', grftname ),        indlvl, bindent, pbs.fidma );
         printpline( sprintf( '[it,igt_,_] = s2mpj_ii(''%s'',igt_)',       grftname ),        indlvl, bindent, pbs.fidpy );
         printjline( sprintf( 'it,igt_,_ = s2mpj_ii("%s",igt_)',           grftname ),        indlvl, bindent, pbs.fidjl );

         %   Specific values
         
         switch ( f{1} )

         %  The group's variable
         
         case 'GV'
            grftype{pos}.grvar = f{3};

         %  The group's parameters
         
         case 'GP'

               %  Have group parameters been defined for the current group?  If not, initialize their list.
               
               if ( ~has_grpar )
                  printpline( 'grftp = []',                                                   indlvl, bindent, pbs.fidpy );
                  printjline( 'grftp = Vector{Vector{String}}()',                             indlvl, bindent, pbs.fidjl );
                  has_grpar = 1;
               end
               
            if ( isfield( grftype{pos}, 'grpar' ) )
               posp = length( grftype{pos}.grpar ) + 1;
               grftype{pos}.grpar{posp} = f{3};
               printmline( sprintf( 'grftp{it}{%d} = ''%s'';',    posp,   f{3} ),             indlvl, bindent, pbs.fidma );
               printpline( sprintf( 'grftp = loaset(grftp,it,%d,''%s'')', posp-1, f{3} ),     indlvl, bindent, pbs.fidpy );
               printjline( sprintf( 'loaset(grftp,it,%d,"%s")', posp, f{3}  ),                indlvl, bindent, pbs.fidjl );
            else
               grftype{pos}.grpar{1} = f{3};
               printmline( sprintf( 'grftp{it}{1} = ''%s'';',    f{3} ),                      indlvl, bindent, pbs.fidma );
               printpline( sprintf( 'grftp = loaset(grftp,it,0,''%s'')', f{3} ),              indlvl, bindent, pbs.fidpy );
               printjline( sprintf( 'loaset(grftp,it,1,"%s")', f{3} ),                        indlvl, bindent, pbs.fidjl );
            end
            if ( nf > 4 )
               posp = length( grftype{pos}.grpar ) + 1;
               grftype{pos}.grpar{posp} = f{3};
               printmline( sprintf( 'grftp{it}{%d} = ''%s'';',    posp,   f{5} ),             indlvl, bindent, pbs.fidma );
               printpline( sprintf( 'grftp = loaset(grftp,it,%d,''%s'')', posp-1, f{5} ),     indlvl, bindent, pbs.fidpy );
               printjline( sprintf( 'loaset(grftp,it,%d,"%s")', posp, f{5} ),                 indlvl, bindent, pbs.fidjl );
            end
         end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  GRUSES  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      case 'GRUSES'

         [ ~, act ] = formact( f{1} );
         switch( f{2} )

         %  The defaults
         
         case '''DEFAULT'''
            switch ( act )  
            case 'T'
               defgrftype = [ 'g', nlfname(f{3}) ];
               if ( ~strcmp( defgrftype, 'TRIVIAL' ) )
                  switch ( pbs.lang )
                  case 'matlab'
                     printmline( 'for ig = 1:ngrp',                                       indlvl,     bindent, pbs.fidma );
                     printmline( sprintf( 'pbm.grftype{ig} = ''%s'';', defgrftype ),      indlvl + 1, bindent, pbs.fidma );
                     printmline( 'end',                                                   indlvl,     bindent, pbs.fidma );
                  case 'python' 
                     printpline( 'for ig in range(0,ngrp):',                              indlvl,     bindent, pbs.fidpy );
                     printpline( sprintf( 'self.grftype = arrset(self.grftype,ig,''%s'')', defgrftype ), ...
                                                                                          indlvl + 1, bindent, pbs.fidpy );
                  case 'julia'
                     printjline( 'for ig in 1:ngrp',                                      indlvl,     bindent, pbs.fidjl );
                     printjline( sprintf( 'arrset(pbm.grftype,ig,"%s")', defgrftype ),    indlvl + 1, bindent, pbs.fidjl );
                     printjline( 'end',                                                   indlvl,     bindent, pbs.fidjl );
                  end
               end
            end

         otherwise

            if ( nf > 2 )
               gname = s2mpjname('',f{2},pbs);    %  Note that the group must have been defined in the data GROUPS section.
               
               %  Is this the same group as the previous one?  If not, find its flat index.

               if ( ~strcmp( f{2}, prevgname ) )
                  printmline( sprintf( 'ig = ig_(%s);', gname ),                              indlvl, bindent, pbs.fidma );
                  printpline( sprintf( 'ig = ig_[%s]',  gname ),                              indlvl, bindent, pbs.fidpy );
                  printjline( sprintf( 'ig = ig_[%s]',  gname ),                              indlvl, bindent, pbs.fidjl );
                  prevgname = f{2};        %  Remember the SIF name of the current group.
               end
            
               switch( act )

               % The group's type
            
               case 'T'

                  tname = [ 'g', nlfname( f{3}) ];
                  printmline( sprintf( 'pbm.grftype{ig} = ''%s'';',tname),                    indlvl, bindent, pbs.fidma );
                  printpline( sprintf( 'self.grftype = arrset(self.grftype,ig,''%s'')',tname ), ...
                                                                                              indlvl, bindent, pbs.fidpy );
                  printjline( sprintf( 'arrset(pbm.grftype,ig,"%s")',tname ),                 indlvl, bindent, pbs.fidjl );

               %  The group's elements
                  
               case 'E'

                  %  Add an element to the group's list of elements.

                  ename = s2mpjname('',f{3},pbs);
                  switch( pbs.lang )
                  case 'matlab'
                     printmline( 'posel = length(pbm.grelt{ig})+1;',                          indlvl, bindent, pbs.fidma );
                     printmline( sprintf( 'pbm.grelt{ig}(posel) = ie_(%s);', ename ),         indlvl, bindent, pbs.fidma );
                  case 'python'
                     printpline( 'posel = len(self.grelt[ig])',                               indlvl, bindent, pbs.fidpy );
                     printpline( sprintf( 'self.grelt = loaset(self.grelt,ig,posel,ie_[%s])', ename ),     ...
                                                                                              indlvl, bindent, pbs.fidpy );
                  case 'julia'
                     printjline( 'posel = length(pbm.grelt[ig])+1',                           indlvl, bindent, pbs.fidjl );
                     printjline( sprintf( 'loaset(pbm.grelt,ig,posel,ie_[%s])', ename ),      indlvl, bindent, pbs.fidjl );
                  end
                  if ( has_constraints )
                     printmline( 'nlc = union(nlc,ig);',                                      indlvl, bindent, pbs.fidma );
                     printpline( 'nlc = np.union1d(nlc,np.array([ig]))',                      indlvl, bindent, pbs.fidpy );
                     printjline( 'arrset(nlc,length(nlc)+1,ig)',                              indlvl, bindent, pbs.fidjl );
                     has_nonlinc = 1;  % remember that nonlinear elements were found
                  end
                  
                  %  Define the associated element weight.

                  if ( nf > 3 )
                     if ( f{1}(1) == 'Z' )
                        printmline( sprintf( 'pbm.grelw{ig}(posel) = %s;', getv1(f{1},f,pbs) ),                         ...
                                                                                              indlvl, bindent, pbs.fidma );
                        printpline( sprintf( 'self.grelw = loaset(self.grelw,ig,posel,float(%s))', getv1(f{1},f,pbs) ), ...
                                                                                              indlvl, bindent, pbs.fidpy );
                        printjline( sprintf( 'loaset(pbm.grelw,ig,posel,Float64(%s))', getv1(f{1},f,pbs) ),             ...
                                                                                              indlvl, bindent, pbs.fidjl );
                     else
                        if( ~isempty( f{4} ) )
                           printmline( sprintf( 'pbm.grelw{ig}(posel) = %s;', getv1r(f{1},f,pbs) ),                     ...
                                                                                              indlvl, bindent, pbs.fidma );
                           printpline( sprintf( 'self.grelw = loaset(self.grelw,ig,posel,float(%s))',                   ...
                                                getv1r(f{1},f,pbs) ),                         indlvl, bindent, pbs.fidpy );
                           printjline( sprintf( 'loaset(pbm.grelw,ig,posel,Float64(%s))', getv1r(f{1},f,pbs) ),         ...
                                                                                              indlvl, bindent, pbs.fidjl );
                        else
                           printmline( 'pbm.grelw{ig}(posel) = 1.;',                          indlvl, bindent, pbs.fidma );
                           printpline( 'self.grelw = loaset(self.grelw,ig,posel,1.)',         indlvl, bindent, pbs.fidpy );
                           printjline( 'loaset(pbm.grelw,ig,posel,1.)',                       indlvl, bindent, pbs.fidjl );
                        end

                        %  Possibly add another element and its weight.
                        
                        if ( nf > 4 )
                           ename = s2mpjname('',f{5},pbs);
                           switch ( pbs.lang )
                           case 'matlab'
                              printmline( 'posel = posel+1;',                                 indlvl, bindent, pbs.fidma );
                              printmline( sprintf( 'pbm.grelt{ig}(posel) = ie_(%s);', ename ),                ...
                                                                                              indlvl, bindent, pbs.fidma );
                           case 'python'
                              printpline( 'posel = posel+1',                                  indlvl, bindent, pbs.fidpy ); 
                              printpline( sprintf( 'self.grelt = loaset(self.grelt,ig,posel,ie_[%s])', ename ), ...
                                                                                              indlvl, bindent, pbs.fidpy );
                           case 'julia'
                              printjline( 'posel = posel+1',                                  indlvl, bindent, pbs.fidjl );
                              printjline( sprintf( 'loaset(pbm.grelt,ig,posel,ie_[%s])', ename ),             ...
                                                                                              indlvl, bindent, pbs.fidjl );
                           end
                           if ( nf > 5 )
                              printmline( sprintf( 'pbm.grelw{ig}(posel) = %s;', getv2r(f) ), indlvl, bindent, pbs.fidma );
                              printpline( sprintf( 'self.grelw = loaset(self.grelw,ig,posel,float(%s))', getv2r(f) ),   ...
                                                                                              indlvl, bindent, pbs.fidpy );
                              printjline( sprintf( 'loaset(pbm.grelw,ig,posel,Float64(%s))', getv2r(f) ),             ...
                                                                                              indlvl, bindent, pbs.fidjl );
                           else
                              printmline( 'pbm.grelw{ig}(posel) = 1.;',                       indlvl, bindent, pbs.fidma );
                              printpline( 'self.grelw = loaset(self.grelw,ig,posel, 1.)',     indlvl, bindent, pbs.fidpy );
                              printjline( 'loaset(pbm.grelw,ig,posel, 1.)',                   indlvl, bindent, pbs.fidjl );
                           end
                        
                        end
                     end
                  else
                     printmline( 'pbm.grelw{ig}(posel) = 1.;',                                indlvl, bindent, pbs.fidma );
                     printpline( 'self.grelw = loaset(self.grelw,ig,posel,1.)',               indlvl, bindent, pbs.fidpy );
                     printjline( 'loaset(pbm.grelw,ig,posel,1.)',                             indlvl, bindent, pbs.fidjl );
                  end

               %  The group's parameters
            
               case 'P'

                  %  Assign values to 1 or 2 elemental parameters.
                  %
                  %  Note that the group ftype must have been defined because the TRIVIAL ftype has no parameter.

                  switch ( pbs.lang )
                  case 'matlab'
                     printmline( sprintf( '[~,posgp] = ismember(''%s'',grftp{igt_(pbm.grftype{ig})});',f{3})',        ...
                                                                                              indlvl, bindent, pbs.fidma );
                     printmline( sprintf( 'pbm.grpar{ig}(posgp) = %s;', getv1r( f{1}, f, pbs ) ),                     ...
                                                                                              indlvl, bindent, pbs.fidma );
                     if ( nf > 5 )
                        printmline( sprintf( '[~,posgp] = ismember(''%s'',grftp[igt_(pbm.grftype[ig])]);', f{5} ),    ...
                                                                                              indlvl, bindent, pbs.fidma );
                        printmline( sprintf( 'pbm.grpar{ig}(posgp) = %s;', getv2r( f ) ),     indlvl, bindent, pbs.fidma );
                     end
                  case 'python'
                     printpline( sprintf( 'posgp = np.where(grftp[igt_[self.grftype[ig]]]==''%s'')[0]', f{3})',    ...
                                                                                              indlvl, bindent, pbs.fidpy );
                     printpline( sprintf( 'self.grpar =loaset(self.grpar,ig,posgp[0],float(%s))',                       ...
                                           getv1r( f{1}, f, pbs ) ),                          indlvl, bindent, pbs.fidpy );
                     if ( nf > 5 )
                        printpline( sprintf( 'posgp = np.where(grftp[igt_[pbm.grftype[ig]]]==''%s'')[0]', f{5} ), ...
                                                                                              indlvl, bindent, pbs.fidpy );
                        printpline( sprintf( 'self.grpar = loaset(self.grpar,ig,posgp[0],float(%s))',getv2r( f ) ),     ...
                                                                                              indlvl, bindent, pbs.fidpy );
                     end
                  case 'julia'
                     printjline( sprintf( 'posgp = findfirst(x->x=="%s",grftp[igt_[pbm.grftype[ig]]])', f{3})',       ...
                                                                                              indlvl, bindent, pbs.fidjl );
                     printjline( sprintf( 'loaset(pbm.grpar,ig,posgp,Float64(%s))', getv1r( f{1}, f, pbs ) ),         ...
                                                                                              indlvl, bindent, pbs.fidjl );
                     if ( nf > 5 )
                        printjline( sprintf( 'posgp = findfirst(x=="%s",grftp[igt_[pbm.grftype[ig]]])', f{5} ),       ...
                                                                                              indlvl, bindent, pbs.fidjl );
                        printjline( sprintf( 'loaset(pbm.grpar,ig,posgp,Float64(%s))', getv2r( f ) ),                 ...
                                                                                              indlvl, bindent, pbs.fidjl );
                     end
                  end
               end
            end
         end   

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  OBOUND  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      case 'OBOUND'

         %  If requested, write the definitions for alternate sets of objective bounds as comments in the output file.
         
         if ( ~isempty( oname ) )
            if ( ~strcmp( f{2}, oname ) )
               if ( writealtsets ) 
                  printcmline( sprintf( '%s', line(2:end) ),                                          bindent, pbs.fidma );
                  if ( ~prevlineispass )
                     printpline( sprintf( 'pass # %s', line(2:end) ),                         indlvl, bindent, pbs.fidpy );
                     prevlineispass = 1;
                  end
                  printcjline( sprintf( '%s', line(2:end) ),                                          bindent, pbs.fidjl );
               else
                  if ( ~prevlineispass )
                     printpline( 'pass',                                                      indlvl, bindent, pbs.fidpy );
                     prevlineispass = 1;
                  end
               end
               continue
            end
         else
            oname = f{2};
         end

         %  Specific values
         
         switch( f{1} )
         case { 'LO', 'XL' }
            printmline( sprintf( 'pb.objlower = %s;',  d2e( f{4}) ),                          indlvl, bindent, pbs.fidma );
            printpline( sprintf( 'self.objlower = %s', d2e( f{4}) ),                          indlvl, bindent, pbs.fidpy );
            printjline( sprintf( 'pb.objlower = %s',   d2e( f{4}) ),                          indlvl, bindent, pbs.fidjl );
         case 'ZL'
            printmline( sprintf( 'pb.objlower = v_(%s);',  s2mpjname('',f{5},pbs) ),          indlvl, bindent, pbs.fidma );
            printpline( sprintf( 'self.objlower = v_[%s]', s2mpjname('',f{5},pbs) ),          indlvl, bindent, pbs.fidpy );
            printjline( sprintf( 'pb.objlower = v_[%s]',   s2mpjname('',f{5},pbs) ),          indlvl, bindent, pbs.fidjl );
         case { 'UP', 'XU' }
            printmline( sprintf( 'pb.objupper = %s;',  d2e( f{4}) ),                          indlvl, bindent, pbs.fidma );
            printpline( sprintf( 'self.objupper = %s', d2e( f{4}) ),                          indlvl, bindent, pbs.fidpy );
            printjline( sprintf( 'pb.objupper = %s',   d2e( f{4}) ),                          indlvl, bindent, pbs.fidjl );
         case 'ZU'
            printmline( sprintf( 'pb.objupper = v_(%s);',  s2mpjname('',f{5},pbs) ),          indlvl, bindent, pbs.fidma );
            printpline( sprintf( 'self.objupper = v_[%s]', s2mpjname('',f{5},pbs) ),          indlvl, bindent, pbs.fidpy );
            printjline( sprintf( 'pb.objupper = v_[%s]',   s2mpjname('',f{5},pbs) ),          indlvl, bindent, pbs.fidjl );
         end
         posdoll = strfind( line( 37:end ), '$' );
         if ( ~isempty( posdoll ) )
            printmline( sprintf( '                           %% %s;', line( posdoll+38:end ) ), ...
                                                                                              indlvl, bindent, pbs.fidma );
            printpline( sprintf( '                           # %s',   line( posdoll+38:end ) ), ...
                                                                                              indlvl, bindent, pbs.fidpy );
            printjline( sprintf( '                           # %s',   line( posdoll+38:end ) ), ...
                                                                                              indlvl, bindent, pbs.fidjl );
         end

      end
      
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Nonlinear elements definitions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   elseif ( inelements )

      [ f, nf ] = parselgroupline( line, lline );
      if ( showsiflines )
         printelgroupline( f, nf, pbs.nline );
      end

      switch ( section )

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  ETEMPS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      case 'ETEMPS'

         switch( f{1} )
         case 'F'
            pbs.errors{end+1} = ...
                        '*** UNSUPPORTED: calls to external functions in the element''s definitions are not supported ***';
            if ( pbs.disperrors)
               disp( pbs.errors{end} )
            end
            errors = pbs.errors;
            exitc  = length( errors );
%            return
         end
         
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EGLOBS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      case 'EGLOBS'

         switch( f{1} )
     
         case 'A'

            n_eglobs = n_eglobs + 1;
            evalstr  = from_fortran( f{4}, globdict, pbs );
            switch ( pbs.lang )
            case 'matlab'
               if ( n_eglobs == 1 )
                  printmline( 'pbm.efpar = [];',                                              indlvl, bindent, pbs.fidma );
               end
               pending{1}       = sprintf( 'pbm.efpar(%d) = %s;', n_eglobs, evalstr );
               globdict( f{2} ) = [ 'pbm.efpar(', int2str( n_eglobs ), ')' ];
            case 'python'
               if ( n_eglobs == 1 )
                  printpline( 'self.efpar = np.array([]);',                                   indlvl, bindent, pbs.fidpy );
               end
               pending{1}    = sprintf( 'self.efpar = arrset( self.efpar,%d,%s)', n_eglobs-1, evalstr );
               globdict( f{2} ) = [ 'self.efpar[', int2str( n_eglobs-1 ), ']' ];
            case 'julia'
               pending{1}    = sprintf( 'arrset(pbm.efpar,%d,%s)', n_eglobs, evalstr );
               globdict( f{2} ) = [ 'pbm.efpar[', int2str( n_eglobs ), ']' ];
            end
            pendingkey       = 'fortran';
            
         case { 'A+', 'I+', 'E+' }

            if ( nf >= 4 && ~isempty( f{4} ) )
               switch( pbs.lang )
               case 'matlab'
                  pending{1}  = [ pending{1}(1:end-1), from_fortran( f{4}, globdict, pbs ), ';' ];
               case { 'python', 'julia' }
                  pending{1}  = [ pending{1}, from_fortran( f{4}, globdict, pbs ) ];
               end
            end
            
         case 'I'

            if ( isKey( globdict, f{3} ) )
               parname    = globdict(f{3});
               [ ib, ie ] = regexp( parname, '[0_9]+' );
               i_eglobs   = round( str2double( parname(ib:ie) ) );  % the index of the glob parameter in pbm.efpar
            else
               n_eglobs   = n_eglobs + 1;
               i_eglobs   = n_eglobs;
               switch (pbs.lang )
               case 'matlab'
                  globdict( f{3} ) =  [ 'pbm.efpar(', int2str( n_eglobs ), ')' ];
               case 'python'
                  globdict( f{3} ) =  [ 'self.efpar[', int2str( n_eglobs ), ']' ];
               case 'julia'
                  globdict( f{3} ) =  [ 'pbm.efpar[', int2str( n_eglobs ), ']' ];
               end
            end
            evalstr = from_fortran( f{4}, globdict, pbs );
            switch( pbs.lang )
            case 'matlab'
               if ( n_eglobs == 1 )
                  printmline( 'pbm.efpar = [];',                                              indlvl, bindent, pbs.fidma );
               end
              printmline( sprintf( 'if(%s)', globdict( f{2} ) ),                              indlvl, bindent, pbs.fidma );
              pending{1}       = sprintf( '%spbm.efpar(%d) = %s;    %% this is %s', bindent, i_eglobs, evalstr, f{3} );
              pending{2}       = 'end';
            case 'python'
               if ( n_eglobs == 1 )
                  printpline( 'self.efpar = np.array([]);', indlvl, bindent, pbs.fidpy );
               end
              printpline( sprintf( 'if %s!=0:', globdict( f{2} ) ),                           indlvl, bindent, pbs.fidpy );
              pending{1}     = sprintf( '%sself.efpar = arrset(self.efpar,%d,%s)     # this is %s',  ...
                                        bindent, i_eglobs-1, evalstr, f{3} );
              if ( i_eglobs == n_eglobs )
                 globdict( f{3} ) =  [ 'self.efpar[', int2str( n_eglobs-1 ), ']' ];
              end
            case 'julia'
              printjline( sprintf( 'if %s!=0', globdict( f{2} ) ),                            indlvl, bindent, pbs.fidjl );
              pending{1}     = sprintf( '%sarrset(pbm.efpar,%d,%s)    # this is %s',  bindent, i_eglobs, evalstr, f{3} );
              pending{2}     = 'end';
            end
            pendingkey = 'fortran';

         case 'E'

            if ( isKey( globdict, f{3} ) )
               parname    = globdict( f{3} );
               [ ib, ie ] = regexp( parname, '[0_9]+' ) 
               i_eglobs   = round( str2double( parname(ib:ie) ) );  % the index of the glob parameter in pbm.efpar
            else
               n_eglobs   = n_eglobs + 1;
               i_eglobs   = n_eglobs;
               switch ( pbs.lang )
               case 'matlab'
                  globdict( f{3} ) =  [ 'pbm.efpar(', int2str( n_eglobs ), ')' ];
               case 'python'
                  globdict( f{3} ) =  [ 'self.efpar[', int2str( n_eglobs ), ']' ];
               case 'julia'
                  globdict( f{3} ) =  [ 'pbm.efpar[', int2str( n_eglobs ), ']' ];
               end
            end
            evalstr = from_fortran( f{4}, globdict, pbs );
            switch( pbs.lang )
            case 'matlab'
               if ( n_eglobs == 1 )
                  printmline( 'pbm.efpar = [];',                                              indlvl, bindent, pbs.fidma );
               end
              printmline( sprintf( 'if(~%s)',  globdict( f{2} ) ),                            indlvl, bindent, pbs.fidma );
              pending{1}       = sprintf( '%spbm.efpar(%d) = %s;    %% this is %s', bindent, i_eglobs, evalstr, f{3} );
              pending{2}       = 'end';
            case 'python'
               if ( n_eglobs == 1 )
                  printpline( 'self.efpar = np.array([]);',                                   indlvl, bindent, pbs.fidpy );
               end
              printpline( sprintf( 'if %s==0:',  globdict( f{2} ) ),                          indlvl, bindent, pbs.fidpy );
              pending{1}     = sprintf( '%sself.efpar = arrset(self.efpar,%d,%s)     # this is %s',  ...
                                        bindent, i_eglobs-1, evalstr, f{3} );
            case 'julia'
              printjline( sprintf( 'if !%s',  globdict( f{2} ) ),                             indlvl, bindent, pbs.fidjl );
              pending{1}     = sprintf( '%sarrset(pbm.efpar,%d,%s)    # this is %s',  bindent, i_eglobs, evalstr, f{3} );
              pending{2}     = 'end';
            end
            pendingkey = 'fortran';

         end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  EINDIVS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      case 'EINDIVS'

         switch ( f{1} )
         case 'T'

            %  Conclude the code for the previous element, if any.
            
            conclude_nlf( codeg, codeH, isempty( todefine ), indlvl, bindent, pbs );

            ename = [ 'e', nlfname(f{2}) ];
            switch ( pbs.lang )
            case 'matlab'
               printmline( ' ',                                                               0,      bindent, pbs.fidma );
               printmline( sprintf( 'case ''%s''', ename ),                                   1,      bindent, pbs.fidma );
               printmline( ' ',                                                               0,      bindent, pbs.fidma );
            case 'python'
               printpline( ' ',                                                               0,      bindent, pbs.fidpy );
               printpline( '@staticmethod',                                                   1,      bindent, pbs.fidpy );
               printpline( sprintf( 'def %s(self, nargout,*args):', ename ),                  1,      bindent, pbs.fidpy );
               printpline( ' ',                                                               0,      bindent, pbs.fidpy );
               printpline( 'import numpy as np',                                              2,      bindent, pbs.fidpy );
            case 'julia'
               printjline( ' ',                                                               0,      bindent, pbs.fidjl );
               printjline( sprintf( 'elseif action == "%s"', ename ),                         1,      bindent, pbs.fidjl );
               printjline( ' ',                                                               0,      bindent, pbs.fidjl );
            end
            
            %  Get the elemental variables and the element index.
            
            switch( pbs.lang )
            case 'matlab'
               printmline( 'EV_  = varargin{1};',                                             indlvl, bindent, pbs.fidma );
               printmline( 'iel_ = varargin{2};',                                             indlvl, bindent, pbs.fidma ); 
            case 'python'
               printpline( 'EV_  = args[0]',                                                  indlvl, bindent, pbs.fidpy );
               printpline( 'iel_ = args[1]',                                                  indlvl, bindent, pbs.fidpy );
            case 'julia'
               printjline( 'EV_     = args[1]',                                               indlvl, bindent, pbs.fidjl );
               printjline( 'iel_    = args[2]',                                               indlvl, bindent, pbs.fidjl );
               printjline( 'nargout = args[3]',                                               indlvl, bindent, pbs.fidjl );
               printjline( 'pbm     = args[4]',                                               indlvl, bindent, pbs.fidjl );
            end

            %  Find the index of the element type specified by f{2}.

            post = ismemstr( ename, elftype, 'name' );
            if ( ~post )
               pbs.errors{end+1} = sprintf( 'ERROR in line %d: element type %s not found', pbs.nline, f{2} );
               if ( pbs.disperrors )
                  disp( pbs.errors{end} )
               end
            end

            %  Build the element-dependent dictionary of reserved names for this type of element:
            %  search for elemental/internal variable's names and elemental parameters' names

            locdict = globdict;

            nel = length( elftype{post}.elvar );
            for i = 1:nel
               switch ( pbs.lang )
               case 'matlab'
                  locdict( elftype{post}.elvar{i} ) = [ 'EV_(', int2str(i), ')' ];
               case 'python'
                  locdict( elftype{post}.elvar{i} ) = [ 'EV_[', int2str(i-1), ']' ];
               case 'julia'
                  locdict( elftype{post}.elvar{i} ) = [ 'EV_[', int2str(i), ']' ];
               end
            end
            if ( isfield( elftype{post}, 'invar' ) )
               for i = 1:length( elftype{post}.invar )
                  switch ( pbs.lang )
                  case 'matlab'
                     locdict( elftype{post}.invar{i} ) = [ 'IV_(', int2str(i), ')' ];
                  case 'python'
                     locdict( elftype{post}.invar{i} ) = [ 'IV_[', int2str(i-1), ']' ];
                  case 'julia'
                     locdict( elftype{post}.invar{i} ) = [ 'IV_[', int2str(i), ']' ];
                  end
               end
            end
            if ( isfield( elftype{post}, 'elpar' ) )
               for i = 1:length( elftype{post}.elpar )
                  switch( pbs.lang )
                  case 'matlab'
                     locdict( elftype{post}.elpar{i} ) = sprintf( 'pbm.elpar{iel_}(%d)', i );
                  case 'python'
                     locdict( elftype{post}.elpar{i} ) = sprintf( 'self.elpar[iel_][%d]', i-1 );
                  case 'julia'
                     locdict( elftype{post}.elpar{i} ) = sprintf( 'pbm.elpar[iel_][%d]', i );
                  end                  
               end
            end

            %  Reset the counter of local Fortran parameters to the number of global Fortran parameters for the nonlinear
            %  elements.

            n_fpar = n_eglobs;
            
            %  Initialize the holders for the code describing the derivatives of the current group.
            
            ing = 0;
            inH = 0;
            codeg = {};
            codeH = {};

            %  Manage internal variables. Because the SIF standard does not specify that internal variables must be defined
            %  before assignments, one needs to check the rhs of every assignment for the presence of internal variables
            %  and (assuming the SIF file is correct) define them before their use.
            
            todefine = [];  % if = i, the i-th internal variable must still be defined,
                            % if 0, it has been defined. Empty if no internal variables.

         %  Assign a parameter.
            
        case 'A'

            evalstr         = from_fortran( f{4}, locdict, pbs );
            todefine        = initIV( todefine, pbs, indlvl, bindent );
            parname         = s2mpjname( 'u', f{2}, pbs) ;
            switch( pbs.lang )
            case 'matlab'
               pending{1}   = sprintf( '%s = %s;', parname, evalstr );
            case { 'python', 'julia' }
               pending{1}   = sprintf( '%s = %s',  parname, evalstr );
            end
            pendingkey      = 'fortran';
            locdict( f{2} ) = parname;

         %  Assign a parameter conditionally (if)
         
         case 'I'
     
            evalstr         = from_fortran( f{4}, locdict, pbs );
            todefine        = initIV( todefine, pbs, indlvl, bindent );
            parname         = s2mpjname( 'u', f{3}, pbs );
            switch ( pbs.lang )
            case 'matlab'
               printmline( sprintf( 'if(%s)', locdict( f{2} ) ),                              indlvl, bindent, pbs.fidma );
               pending{1}   = sprintf( '%s%s = %s;', bindent, parname, evalstr );
               pending{2}   = 'end';
            case 'python'
               printpline( sprintf( 'if %s!=0:', locdict( f{2} ) ),                           indlvl, bindent, pbs.fidpy );
               pending{1}   = sprintf( '%s%s = %s', bindent, parname, evalstr );
            case 'julia'
               printjline( sprintf( 'if %s', locdict( f{2} ) ),                               indlvl, bindent, pbs.fidjl );
               pending{1}   = sprintf( '%s%s = %s', bindent, parname, evalstr );
               pending{2}   = 'end';
            end
            pendingkey      = 'fortran';
            locdict( f{3} ) = parname;

         %  Assign a parameter conditionally (else)
         
         case 'E'

            evalstr       = from_fortran( f{4}, locdict, pbs );
            todefine      = initIV( todefine, pbs, indlvl, bindent );
            parname       = s2mpjname( 'u', f{3}, pbs );
            switch ( pbs.lang )
            case 'matlab'
               printmline( sprintf( 'if(~%s)', locdict( f{2} ) ),                             indlvl, bindent, pbs.fidma );
               pending{1} = sprintf( '%s%s = %s;', bindent, parname, evalstr );
               pending{2} = 'end';
            case 'python'
               printpline( sprintf( 'if %s==0:', locdict( f{2} ) ),                           indlvl, bindent, pbs.fidpy );
               pending{1} = sprintf( '%s%s = %s', bindent, parname, evalstr );
            case 'julia'
               printjline( sprintf( 'if !%s',  locdict( f{2} ) ),                             indlvl, bindent, pbs.fidjl );
               pending{1} = sprintf( '%s%s = %s', bindent, parname, evalstr );
               pending{2} = 'end';
            end
            pendingkey      = 'fortran';
            locdict( f{3} ) = parname;

         %  Define the internal range of the element.

         case 'R'

            [ f, nf ]  = parsedataline( line );                   % reinterpret the line as a data line
            [ ~, pos ] = ismember( f{2}, elftype{post}.invar );   %  --> IV(pos)
            [ ~, ev1 ] = ismember( f{3}, elftype{post}.elvar );
            if ( isempty( todefine ) )
                switch( pbs.lang )
                case 'matlab'
                   printmline( sprintf( 'U_ = zeros(%d,%d);', length( elftype{post}.invar ), ...
                                        length( elftype{post}.elvar ) ),                      indlvl, bindent, pbs.fidma );
                case 'python'
                   printpline( sprintf( 'U_ = np.zeros((%d,%d))', length( elftype{post}.invar ), ...
                                        length( elftype{post}.elvar ) ),                      indlvl, bindent, pbs.fidpy );
                   printpline( sprintf( 'IV_ = np.zeros(%d)',length( elftype{post}.invar ) ), indlvl, bindent, pbs.fidpy );
                case 'julia'
                   printjline( sprintf( 'U_ = zeros(Float64,%d,%d)', length( elftype{post}.invar ),   ... 
                                        length( elftype{post}.elvar ) ),                      indlvl, bindent, pbs.fidjl );
                   printjline( sprintf( 'IV_ =  zeros(Float64,%d)',     length( elftype{post}.invar ) ), ...
                                                                                              indlvl, bindent, pbs.fidjl );
                end                
            end
            switch (pbs.lang )
            case 'matlab'
               updU = replace( sprintf( 'U_(%d,%d) = U_(%d,%d)+%d;', pos, ev1, pos, ev1, str2double( d2e(f{4}) ) ), ...
                           '+-', '-' );
               printmline( updU,                                                              indlvl, bindent, pbs.fidma );
            case 'python'
               updU = replace( sprintf( 'U_[%d,%d] = U_[%d,%d]+%d',  pos-1, ev1-1, pos-1, ev1-1, str2double( d2e(f{4}) ) ), ...
                           '+-', '-' );
               printpline( updU,                                                              indlvl, bindent, pbs.fidpy );
            case 'julia'
               updU = replace( sprintf( 'U_[%d,%d] = U_[%d,%d]+%d',  pos, ev1, pos, ev1, str2double( d2e(f{4}) ) ), ...
                           '+-', '-' );
               printjline( updU,                                                              indlvl, bindent, pbs.fidjl );
            end
            if ( nf > 4 )
               [ ~, ev2 ] = ismember( f{5}, elftype{post}.elvar);
               switch ( pbs.lang )
               case 'matlab'
                  updU = replace( sprintf( 'U_(%d,%d) = U_(%d,%d)+%d;', pos, ev2, pos, ev2, str2double( d2e(f{6}) ) ), ...
                           '+-', '-' );
                  printmline( updU,                                                           indlvl, bindent, pbs.fidma );
               case 'python'
                  updU = replace( sprintf( 'U_[%d,%d] = U_[%d,%d]+%d', pos-1, ev2-1, pos-1, ev2-1, ...
                                  str2double( d2e(f{6}) ) ), '+-', '-' );
                  printpline( updU,                                                           indlvl, bindent, pbs.fidpy );
               case 'julia'
                  updU = replace( sprintf( 'U_[%d,%d] = U_[%d,%d]+%d', pos, ev2, pos, ev2, str2double( d2e(f{6}) ) ), ...
                           '+-', '-' );
                  printjline( updU,                                                           indlvl, bindent, pbs.fidjl );
               end
            end
            todefine( pos ) = 1;
         %  Compute the element's value.
            
         case 'F'

            %  Assign the evaluation string for the element function value.

            evalstrf  = from_fortran( f{4}, locdict, pbs );
            todefine  = initIV( todefine, pbs, indlvl, bindent );
            switch ( pbs.lang )
            case 'matlab' 
               pending{1} = sprintf( 'varargout{1} = %s;', evalstrf );
            case 'python'
               pending{1} = sprintf( 'f_   = %s', evalstrf );
               pending{2} = 'if not isinstance( f_, float ):';
               pending{3} = sprintf( '    f_   = f_.item();');  % force f_ to be a scalar
            case 'julia'
               pending{1} = sprintf( 'f_   = %s', evalstrf );
            end
            pendingkey    = 'fortran';

         %  Continuation lines
            
         case { 'A+', 'I+', 'E+', 'F+' }

            if ( nf >= 4 && ~isempty( f{4} ) )
               switch ( pbs.lang )
               case 'matlab'
                  pending{1} = [ pending{1}(1:end-1), from_fortran( f{4}, locdict, pbs ), ';' ];
               case { 'python', 'julia' }
                  pending{1} = [ pending{1}, from_fortran( f{4}, locdict, pbs ) ];
               end
             end

         %  Compute the element's gradient.
         %
         %  Note that, in Python, the internal gradient vector must be initilaized with a dimension which is the smallest
         %  of the length of the internal variables vector (IV_) and the length of the elemental variables vector (EV_).
         %  Since the same dimension is also used to initialize the internal Hessian (in Python), one also remembers the
         %  index pbs.gstrt in codeg at which the effective gradient code starts, and the internal dimension is then given
         %  as the length of codeg minus pb.gstrt.
             
         case 'G'

            if ( ing == 0 )
               switch ( pbs.lang )
               case 'matlab'
                  codeg{1}  = 'if(nargout>1)';
               case 'python'
                  codeg{1}  = 'if nargout>1:';
                  codeg{2}  = 'try:';               %  Detect the internal dimension, if any
                  codeg{3}  = '    dim = len(IV_)';
                  codeg{4}  = 'except:';            %  No internal variables
                  codeg{5}  = '    dim = len(EV_)';
                  codeg{6}  = 'g_ = np.zeros(dim)';
               case 'julia'
                  codeg{1}  = 'if nargout>1';
                  codeg{2}  = 'dim = try length(IV_) catch; length(EV_) end';   %  Detect the internal dimension, if any
                  codeg{3}  = 'g_  = zeros(Float64,dim)';
               end
               ing = length( codeg );
               pbs.gstrt = ing;
            end
            ing = ing + 1;
            if ( isempty( todefine ) )
               [ vf, pos ] = ismember( f{2}, elftype{post}.elvar );
            else
               [ vf, pos ] = ismember( f{2}, elftype{post}.invar );
            end
            if ( ~vf )
               pbs.errors{end+1} = sprintf(' ERROR in line %d: can''t find the variable %s in element type %s', ...
                                        pbs.nline, f{2}, elftype{post}.name );
               if ( pbs.disperrors )
                  disp( pbs.errors{end} )
               end
            end
            switch ( pbs.lang )
            case 'matlab'
               codeg{ing} = sprintf( 'g_(%d,1) = %s;', pos, from_fortran( f{4}, locdict, pbs ) );
            case 'python'
               codeg{ing} = sprintf( 'g_[%d] = %s',  pos-1, from_fortran( f{4}, locdict, pbs ) );
            case 'julia'
               codeg{ing} = sprintf( 'g_[%d] = %s',  pos,   from_fortran( f{4},locdict, pbs ) );
            end            
            
         case 'G+'

            if ( nf >= 4 && ~isempty( f{4} ) )
               switch( pbs.lang )
               case 'matlab'
                  codeg{ing} = [ codeg{ing}(1:end-1), from_fortran( f{4}, locdict, pbs ), ';' ];
               case { 'python', 'julia' }
                  codeg{ing} = [ codeg{ing},  from_fortran( f{4}, locdict, pbs ) ];
               end               
            end

         %  Compute the element's Hessian.
         
         case 'H'

            if ( inH == 0 )
               switch ( pbs.lang )
               case 'matlab'
                  codeH{1} = sprintf( 'if(nargout>2)' );
               case 'python'
                  codeH{1} = sprintf( 'if nargout>2:' );
               case 'julia'
                  codeH{1} = sprintf( 'if nargout>2' );
               end
               inH = 1;
            end
            inH = inH + 1;
            if ( isempty( todefine ) )
               [ vf1, pos1 ] = ismember( f{2}, elftype{post}.elvar );
               [ vf2, pos2 ] = ismember( f{3}, elftype{post}.elvar );
            else
               [ vf1, pos1 ] = ismember( f{2}, elftype{post}.invar );
               [ vf2, pos2 ] = ismember( f{3}, elftype{post}.invar );
            end            
            if ( ~vf1 )
               pbs.errors{end+1} = sprintf( 'ERROR in line %d: can''t find the variable %s in element type %s',...
                                         pbs.nline, f{2}, elftype{post}.name );
               if ( pbs.disperrors )
                  disp( pbs.errors{end} )
               end
            end
            if ( ~vf2 )
               pbs.errors{end+1} = sprintf( 'ERROR in line %d: can''t find the variable %s in element type %s',...
                                         pbs.nline, f{3}, elftype{post}.name );
               if ( pbs.disperrors )
                  disp( pbs.errors{end} )
               end
            end
            switch( pbs.lang )
            case 'matlab'
               codeH{inH} = sprintf( 'H_(%d,%d) = %s;', pos1  , pos2  , from_fortran( f{4}, locdict, pbs) );
            case 'python'
               codeH{inH} = sprintf( 'H_[%d,%d] = %s',  pos1-1, pos2-1, from_fortran( f{4}, locdict, pbs) );
            case { 'julia' }
               codeH{inH} = sprintf( 'H_[%d,%d] = %s',  pos1  , pos2,   from_fortran( f{4}, locdict, pbs) );
            end

         case 'H+'

            if ( nf >= 4 && ~isempty( f{4} ) )
               switch ( pbs.lang )
               case 'matlab'
                  codeH{inH} = [ codeH{inH}(1:end-1), from_fortran( f{4}, locdict, pbs ), ';' ];
               case { 'python', 'julia' }
                  codeH{inH} = [ codeH{inH}, from_fortran( f{4}, locdict, pbs ) ];
               end
            end
             
         end
      end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Nonlinear group's definitions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   elseif( ingroups )
   
      [ f, nf ] = parselgroupline( line, lline );
      if ( showsiflines )
         printelgroupline( f, nf, pbs.nline );
      end
      
      switch ( section )

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  GTEMPS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      case 'GTEMPS'
      
         switch( f{1} )
         case 'F'
            pbs.errors{end+1} = ...
                        ' *** UNSUPPORTED: calls to external functions in the group''s definitions are not supported ***';
            if ( pbs.disperrors )
               disp( pbs.errors{end} )
            end
            switch( pbs.lang )
            case { 'matlab', 'julia' }
               errors = pb.errors;
            case 'python'
               errors = self.errors;
            end
            exitc  = length( errors );
            return
         end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  GGLOBS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      case 'GGLOBS'

         switch( f{1} )
     
         case 'A'

            evalstr = from_fortran( f{4}, globdict, pbs );
            n_gglobs = n_gglobs + 1;
            switch ( pbs.lang )
            case 'matlab'
               pending{1}       = sprintf( 'pbm.gfpar(%d) = %s;    %% this is %s', n_gglobs, evalstr, f{2} );
               globdict( f{2} ) = [ 'pbm.gfpar(', int2str( n_gglobs ), ')' ];
            case 'python'
               if ( n_gglobs == 1 )
                  printpline( 'self.gfpar = np.array([]);',                                   indlvl, bindent, pbs.fidpy );
               end
               pending{1}    = sprintf( 'self.gfpar = arrset(self.gfpar,%d,%s)     # this is %s', ...
                                         n_gglobs-1, evalstr, f{2} );
               globdict( f{2} ) = [ 'self.gfpar[', int2str( n_gglobs-1 ), ']' ];
            case 'julia'
               pending{1}    = sprintf( 'arrset(pbm.gfpar,%d,%s)    # this is  %s', n_gglobs, evalstr, f{2} );
               globdict( f{2} ) = [ 'pbm.gfpar[', int2str( n_gglobs ), ']' ];
            end
            pendingkey       = 'fortran';
            
         case { 'A+', 'I+', 'E+' }

            if ( nf >= 4 && ~isempty( f{4} ) )
               switch ( pbs.lang )
               case 'matlab'
                  pending{1}  = [ pending{1}(1:end-1), from_fortran( f{4}, globdict, pbs ), ';' ];
               case { 'python', julia' }
                  pending{1}  = [ pending{1}, from_fortran( f{4}, globdict, pbs ) ];
               end
            end

         case 'I'
     
            if ( isKey( globdict, f{3} ) )
               parname    = globdict( f{3} );
               [ ib, ie ] = regexp( parname, '[0_9]+' ) 
               i_eglobs   = round( str2double( parname(ib:ie) ) );  % the index of the glob parameter in pbm.gfpar
            else
               n_gglobs   = n_gglobs + 1;
               i_gglobs   = n_gglobs;
               switch (pbs.lang )
               case 'matlab'
                  globdict( f{3} ) =  [ 'pbm.gfpar(', int2str( n_gglobs ), ')' ];
               case 'python'
                  globdict( f{3} ) =  [ 'self.gfpar[', int2str( n_gglobs ), ']' ];
               case 'julia'
                  globdict( f{3} ) =  [ 'pbm.gfpar[', int2str( n_gglobs ), ']' ];
               end
            end
            evalstr = from_fortran( f{4}, globdict, pbs );
            switch( pbs.lang )
            case 'matlab'
               if ( n_gglobs == 1 )
                  printmline( 'pbm.gfpar = [];',                                              indlvl, bindent, pbs.fidma );
               end
              printmline( sprintf( 'if(%s)', globdict( f{2} ) ),                              indlvl, bindent, pbs.fidma );
              pending{1} = sprintf( '%spbm.gfpar(%d) = %s;    %% this is %s', bindent, i_gglobs, evalstr, f{3} );
              pending{2} = 'end';
            case 'python'
               if ( n_gglobs == 1 )
                  printpline( 'self.gfpar = np.array([]);',                                   indlvl, bindent, pbs.fidpy );
               end
              printpline( sprintf( 'if %s!=0:', globdict( f{2} ) ),                           indlvl, bindent, pbs.fidpy );
              pending{1} = sprintf( '%sself.gfpar = arrset(self.gfpar,%d,%s)     # this is %s', ...
                                     bindent, i_gglobs-1, evalstr, f{3} );
            case 'julia'
              printjline( sprintf( 'if %s!=0', globdict( f{2} ) ),                            indlvl, bindent, pbs.fidjl );
              pending{1}  = sprintf( '%sarrset(pbm.gfpar,%d,%s)', bindent, i_gglobs, evalstr, f{3} );
              pending{2}  = 'end';
            end
            pendingkey = 'fortran';

         case 'E'

            if ( isKey( globdict, f{3} ) )
               parname    =  globdict( f{3} );
               [ ib, ie ] = regexp( parname, '[0_9]+' ) 
               i_eglobs   = round( str2double( parname(ib:ie) ) );  % the index of the glob parameter in pbm.gfpar
            else
               n_gglobs   = n_gglobs + 1;
               i_gglobs   = n_gglobs;
               switch( pbs.lang )
               case 'matlab'
                  globdict( f{3} ) =  [ 'pbm.gfpar(', int2str( n_gglobs ), ')' ];
               case 'python'
                  globdict( f{3} ) =  [ 'self.gfpar[', int2str( n_gglobs ), ']' ];
               case 'matlab'
                  globdict( f{3} ) =  [ 'pbm.gfpar[', int2str( n_gglobs ), ']' ];
               end
            end
            evalstr =  from_fortran( f{4}, globdict, pbs );
            switch( pbs.lang )
            case 'matlab'
               if ( n_gglobs == 1 )
                  printmline( 'pbm.gfpar = {};',                                              indlvl, bindent, pbs.fidma );
               end
               printmline( sprintf( 'if(2%s)', globdict( f{2} ) ),                            indlvl, bindent, pbs.fidma );
               pending{1} = sprintf( '%spbm.gfpar(%d) = %s;    %% this is %s', bindent, i_gglobs, evalstr, f{3} );
               pending{2} = 'end';
            case 'python'
               if ( n_gglobs == 1 )
                  printpline( 'self.gfpar = np.array([]);',                                   indlvl, bindent, pbs.fidpy );
               end
              printpline( sprintf( 'if %s==0:', globdict( f{2} ) ),                           indlvl, bindent, pbs.fidpy );
              pending{1} = sprintf( '%sself.gfpar = arrset(self.gfpar,%d,%s)     # this is  %s',  ...
                                     bindent, i_gglobs-1, evalstr, f{3} );
            case 'julia'
              printjline( sprintf( 'if !%s', globdict( f{2} ) ),                              indlvl, bindent, pbs.fidjl );
              pending{1} = sprintf( '%sarrset(pbm.gfpar,%d,%s)    # this is %s',  bindent, i_gglobs, evalstr, f{3} );
              pending{2} = 'end';
            end
            pendingkey = 'fortran';

         end
         
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  GINDIVS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      case 'GINDIVS'

         switch ( f{1} )
         case 'T'

            %  Conclude the code for the previous group, if any.

            conclude_nlf( codeg, codeH, 1, indlvl, bindent, pbs );

            %  New group

            gname = [ 'g', nlfname(f{2}) ];
            switch ( pbs.lang )
            case 'matlab'
               printmline( ' ',                                                               0,      bindent, pbs.fidma );
               printmline( sprintf( 'case ''%s''', gname ),                                   1,      bindent, pbs.fidma ); 
               printmline( ' ',                                                               0,      bindent, pbs.fidma );
            ;case 'python'
               printpline( ' ',                                                               0,      bindent, pbs.fidpy );
               printpline( '@staticmethod',                                                   1,      bindent, pbs.fidpy );
               printpline( sprintf( 'def %s(self,nargout,*args):', gname ),                   1,      bindent, pbs.fidpy );
               printpline( ' ',                                                               0,      bindent, pbs.fidpy );
            case 'julia'
               printjline( ' ',                                                               0,      bindent, pbs.fidjl );
               printjline( sprintf( 'elseif action == "%s"', gname ),                         1,      bindent, pbs.fidjl ); 
               printjline( ' ',                                                               0,      bindent, pbs.fidjl );
            end
            
            %  Get the group variable and the group index

            switch ( pbs.lang )
            case 'matlab'
               printmline( 'GVAR_ = varargin{1};',                                            indlvl, bindent, pbs.fidma );
               printmline( 'igr_  = varargin{2};',                                            indlvl, bindent, pbs.fidma );
            case 'python'
               printpline( 'GVAR_ = args[0]',                                                 indlvl, bindent, pbs.fidpy );
               printpline( 'igr_  = args[1]',                                                 indlvl, bindent, pbs.fidpy );
            case 'julia'
               printjline( 'GVAR_   = args[1]',                                               indlvl, bindent, pbs.fidjl );
               printjline( 'igr_    = args[2]',                                               indlvl, bindent, pbs.fidjl );
               printjline( 'nargout = args[3]',                                               indlvl, bindent, pbs.fidjl );
               printjline( 'pbm     = args[4]',                                               indlvl, bindent, pbs.fidjl );
            end
            
            %  Find the index of the group type specified by f{2}.

            post = ismemstr( gname, grftype, 'name' );
            if ( ~post )
               pbs.errors{end+1} = sprintf( 'ERROR in line %d: group type %s not found', pbs.nline, f{2} );
               disp( pbs.errors{end} )
            end

            %  Build the group-dependent dictionary of reserved names for this type of group:
            %  search for group variable's name and group parameters' names

            locdict = globdict;
            locdict( grftype{post}.grvar ) = 'GVAR_';
            
            if ( isfield( grftype{post}, 'grpar' ) )
               for i = 1:length( grftype{post}.grpar )
                  switch( pbs.lang )
                  case 'matlab'
                     locdict( grftype{post}.grpar{i} ) = sprintf( 'pbm.grpar{igr_}(%d)', i );
                  case 'python'
                     locdict( grftype{post}.grpar{i} ) = sprintf( 'self.grpar[igr_][%d]', i-1 );
                  case 'julia'
                     locdict( grftype{post}.grpar{i} ) = sprintf( 'pbm.grpar[igr_][%d]', i );
                  end
               end
            end


            %  Initialize the holders for the code describing the derivatives of the current group.
            
            ing     = 0;
            inH     = 0;
            codeg = {};
            codeH = {};
            
         case 'A'

            switch( pbs.lang )
            case 'matlab'
               pending{1}   = sprintf( '%s = %s;', s2mpjname( 'u', f{2}, pbs ), from_fortran( f{4}, locdict, pbs ) );
            case { 'python', 'julia' }
               pending{1}   = sprintf( '%s = %s',  s2mpjname( 'u', f{2}, pbs ), from_fortran( f{4}, locdict, pbs ) );
            end
            pendingkey      = 'fortran';
            locdict( f{2} ) = s2mpjname( 'u', f{2}, pbs );
            
         case 'I'

            parname = s2mpjname( 'u', f{3}, pbs );
            switch( pbs.lang )
            case 'matlab'
               printmline( sprintf( 'if(%s)', locdict( f{2} ) ),                              indlvl, bindent, pbs.fidma );
               pending{1}   = sprintf( '%s%s = %s;', bindent, parname, from_fortran( f{4}, locdict, pbs ) );
               pending{2}   = 'end';
            case 'python'
               printpline( sprintf( 'if %s!=0:', locdict( f{2} ) ),                           indlvl, bindent, pbs.fidpy );
               pending{1}   = sprintf( '%s%s = %s', bindent, parname, from_fortran( f{4}, locdict, pbs ) );
            case 'julia'
               printjline( sprintf( 'if %s', locdict( f{2} ) ),                               indlvl, bindent, pbs.fidjl );
               pending{1}   = sprintf( '%s%s = %s', bindent, parname, from_fortran( f{4}, locdict, pbs ) );
               pending{2}   = 'end';
            end
            pendingkey      = 'fortran';
            locdict( f{3} ) = parname;
            
         case 'E'

            parname = s2mpjname( 'u', f{3}, pbs );
            switch ( pbs.lang )
            case 'matlab'
               printmline( sprintf( 'if(~%s)', locdict( f{2} ) ),                             indlvl, bindent, pbs.fidma );
               pending{1}   = sprintf( '%s%s = %s;', bindent, parname, from_fortran( f{4}, locdict, pbs ) );
               pending{2}   = 'end';
            case 'python'
               printpline( sprintf( 'if %s==0:', locdict( f{2} ) ),                           indlvl, bindent, pbs.fidpy );
               pending{1}   = sprintf( '%s%s = %s', bindent, parname, from_fortran( f{4}, locdict, pbs ) );
            case 'julia'
               printjline( sprintf( 'if !%s', locdict( f{2} ) ),                              indlvl, bindent, pbs.fidjl );
               pending{1}   = sprintf( '%s%s = %s', bindent, parname, from_fortran( f{4}, locdict, pbs ) );
               pending{2}   = 'end';
            end
            pendingkey      = 'fortran';
            locdict( f{3} ) =  parname;
            
         case 'F'

            %  Assign the evaluation string for the element function value.

            switch ( pbs.lang )
            case 'matlab'
               pending{1} = sprintf( 'varargout{1} = %s;', from_fortran( f{4}, locdict, pbs ) );
            case { 'python', 'julia' }
               pending{1} = sprintf( 'f_= %s', from_fortran( f{4}, locdict, pbs ) );
            end
            pendingkey    = 'fortran';
            
         case { 'A+', 'I+', 'E+', 'F+' }

            if ( nf >= 4 && ~isempty( f{4} ) )
               switch( pbs.lang )
               case 'matlab'
                  pending{1} = [ pending{1}(1:end-1), from_fortran( f{4}, locdict, pbs ), ';' ];
               case { 'python', 'julia' }
                  pending{1} = [ pending{1}, from_fortran( f{4}, locdict, pbs ) ];
               end
            end
            
         case 'G'

            if ( ing == 0 )
               switch ( pbs.lang )
               case 'matlab'
                  codeg{1} = 'if(nargout>1)';
               case 'python'
                  codeg{1} = 'if nargout>1:';
               case 'julia'
                  codeg{1} = 'if nargout>1';
               end
               ing       = 1;
               pbs.gstrt = ing;
            end
            ing = ing + 1;
            switch ( pbs.lang )
            case 'matlab'
               codeg{ing} = sprintf( 'g_ = %s;', from_fortran( f{4}, locdict, pbs ) );
            case { 'python', 'julia' }
               codeg{ing} = sprintf( 'g_ = %s',  from_fortran( f{4}, locdict, pbs ) );
            end            
         case 'G+'

            if ( nf >= 4 && ~isempty( f{4} ) )
               switch( pbs.lang )
               case 'matlab'
                  codeg{ing} = [ codeg{ing}(1:end-1), from_fortran( f{4}, locdict, pbs ), ';' ];
               case { 'python', 'julia' }
                  codeg{ing} = [ codeg{ing}, from_fortran( f{4}, locdict, pbs ) ];
               end
            end

         case 'H'
        
            if ( inH == 0 )
               switch ( pbs.lang )
               case 'matlab'
                  codeH{1} = sprintf( 'if(nargout>2)' );
               case 'python'
                  codeH{1} = sprintf( 'if nargout>2:' );
               case 'julia'
                  codeH{1} = sprintf( 'if nargout>2' );
               end
               inH = 1;
            end
            inH = inH + 1;
            switch ( pbs.lang )
            case 'matlab'
               codeH{inH} = sprintf( 'H_ = %s;', from_fortran( f{4}, locdict, pbs ) );
            case { 'python', 'julia' }
               codeH{inH} = sprintf( 'H_ = %s',  from_fortran( f{4}, locdict, pbs ) );
            end
            
         case 'H+'

            if ( nf >= 4 && ~isempty( f{4} ) )
               switch ( pbs.lang )
               case 'matlab'
                  codeH{inH} = [ codeH{inH}(1:end-1),from_fortran( f{4}, locdict, pbs ), ';' ];
               case { 'python', 'julia' }
                  codeH{inH} = [codeH{inH}, from_fortran( f{4}, locdict, pbs ) ];
               end
            end
            
         end
      end
      
   end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  Write the instructions to call the high-level actions. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch( pbs.lang )
case 'matlab'
   printmline( ' ',                                                                           0,      bindent, pbs.fidma );
   printmline( '%%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%',                           1,      bindent, pbs.fidma );
   printmline( ' ',                                                                           0,      bindent, pbs.fidma );
   printmline(                                                                                                     ...
      'case {''fx'',''fgx'',''fgHx'',''cx'',''cJx'',''cJHx'',''cIx'',''cIJx'',''cIJHx'',''cIJxv'',''fHxv'',...' ,  ...
                                                                                              1,      bindent, pbs.fidma );
   printmline( '      ''cJxv'',''cJtxv'',''cIJtxv'',''Lxy'',''Lgxy'',''LgHxy'',''LIxy'',''LIgxy'',''LIgHxy'',...', ...
                                                                                              1,      bindent, pbs.fidma );
   printmline( '      ''LHxyv'',''LIHxyv''}',                                                 1,      bindent, pbs.fidma );
   printmline( ' ',                                                                           0,      bindent, pbs.fidma );
   printmline( 'if(isfield(pbm,''name'')&&strcmp(pbm.name,name))',                            2,      bindent, pbs.fidma );
   printmline( sprintf( 'pbm.has_globs = [%d,%d];', n_eglobs, n_gglobs ),                     3,      bindent, pbs.fidma );
   printmline( '[varargout{1:max(1,nargout)}] = s2mpjlib(action,pbm,varargin{:});',           3,      bindent, pbs.fidma );
   printmline( 'else',                                                                        2,      bindent, pbs.fidma );
   printmline( 'disp([''ERROR: please run '',name,'' with action = setup''])',                3,      bindent, pbs.fidma );
   printmline( '[varargout{1:nargout}] = deal(NaN);',                                         3,      bindent, pbs.fidma );
   printmline( 'end',                                                                         2,      bindent, pbs.fidma );
case 'python'
   %  Nothing to do here: everything is managed within s2mpjlib.py.
case 'julia'
   printjline( ' ',                                                                           0,      bindent, pbs.fidjl );
   printjline( '#%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%',                           1,      bindent, pbs.fidjl );
   printjline( ' ',                                                                           0,      bindent, pbs.fidjl );
   printjline( 'elseif action in  ["fx","fgx","fgHx","cx","cJx","cJHx","cIx","cIJx","cIJHx","cIJxv","fHxv",', ...
                                                                                              1,      bindent, pbs.fidjl ); 
   printjline( '                   "cJxv","cJtxv","cIJtxv","Lxy","Lgxy","LgHxy","LIxy","LIgxy","LIgHxy",',    ...
                                                                                              1,      bindent, pbs.fidjl ); 
   printjline( '                   "LHxyv","LIHxyv"]',                                        1,      bindent, pbs.fidjl ); 
   printjline( ' ',                                                                           0,      bindent, pbs.fidjl );
   printjline( 'pbm = args[1]',                                                               2,      bindent, pbs.fidjl );
   printjline( 'if pbm.name == name',                                                         2,      bindent, pbs.fidjl );
   printjline( sprintf( 'pbm.has_globs = [%d,%d]', n_eglobs, n_gglobs ),                      3,      bindent, pbs.fidjl );
   printjline( 'return s2mpj_eval(action,args...)',                                           3,      bindent, pbs.fidjl );
   printjline( 'else',                                                                        2,      bindent, pbs.fidjl );
   printjline( 'println("ERROR: please run "*name*" with action = setup")',                   3,      bindent, pbs.fidjl );
   printjline( 'return ntuple(i->undef,args[end])',                                           3,      bindent, pbs.fidjl );
   printjline( 'end',                                                                         2,      bindent, pbs.fidjl );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Conclude %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Close the problem's SIF file.

fclose( fidSIF );

%  Terminate and close the problem's output file.

%  1) the list of actions

switch( pbs.lang )
case 'matlab'
   printmline( ' ',                                                                           0,      bindent, pbs.fidma );
   printmline( 'otherwise',                                                                   1,      bindent, pbs.fidma );
   printmline( 'disp(['' ERROR: action '',action,'' unavailable for problem '',name,''.m''])', 2,     bindent, pbs.fidma );
   printmline( 'end',                                                                         1,      bindent, pbs.fidma );
case 'python'
   %  The python output files are pretty terse...   
case 'julia'
   printjline( ' ',                                                                           0,      bindent, pbs.fidjl );
   printjline( 'else',                                                                        1,      bindent, pbs.fidjl );
   printjline( 'println("ERROR: action "*action*" unavailable for problem "*name*".jl")',     2,      bindent, pbs.fidjl );
   printjline( 'return ntuple(i->undef,args[end])',                                           2,      bindent, pbs.fidjl );
   printjline( 'end',                                                                         1,      bindent, pbs.fidjl );
end

%  2) the file itself

switch( pbs.lang )
case 'matlab'
   printmline( ' ',                                                                           0,      bindent, pbs.fidma );
   printmline( 'return',                                                                      0,      bindent, pbs.fidma );
   printmline( ' ',                                                                           0,      bindent, pbs.fidma );
   printmline( 'end',                                                                         0,      bindent, pbs.fidma );
   printmline( ' ',                                                                           0,      bindent, pbs.fidma );
   printmline( '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',                        0,      bindent, pbs.fidma );
   printmline( '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',                        0,      bindent, pbs.fidma );
   printmline( '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',                        0,      bindent, pbs.fidma );
   printmline( ' ',                                                                           0,      bindent, pbs.fidma );
case 'python'
   printpline( ' ',                                                                           0,      bindent, pbs.fidpy );
   printpline( '#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',                        0,      bindent, pbs.fidpy );
   printpline( '#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',                        0,      bindent, pbs.fidpy );
   printpline( '#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',                        0,      bindent, pbs.fidpy );
   printpline( ' ',                                                                           0,      bindent, pbs.fidpy );
case 'julia'
   printjline( ' ',                                                                           0,      bindent, pbs.fidjl );
   printjline( 'end',                                                                         0,      bindent, pbs.fidjl );
   printjline( ' ',                                                                           0,      bindent, pbs.fidjl );
   printjline( '#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',                        0,      bindent, pbs.fidjl );
   printjline( '#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',                        0,      bindent, pbs.fidjl );
   printjline( '#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',                        0,      bindent, pbs.fidjl );
   printjline( ' ',                                                                           0,      bindent, pbs.fidjl );
end

if ( pbs.fidma && pbs.fidma ~= 1 )
   fclose( pbs.fidma );
end
if ( pbs.fidpy && pbs.fidpy ~= 1 )
   fclose( pbs.fidpy );
end
if ( pbs.fidjl && pbs.fidjl ~= 1 )
   fclose( pbs.fidjl );
end

%  Set the ouputs.

errors = pbs.errors;
exitc  = length( errors );

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ f, nf ] = parsedataline( line )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Parses a data line into its up to 6 different fields.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = {};

line  = regexprep( line, '\s*$','');       %  remove trailing blanks
lline = length( line );
if ( lline < 3 )
   return
end
f{1} = strtrim( line(2:3) );
if ( lline > 4 )
   f{2} = strtrim( line( 5:min( lline, 14 ) ) );
end
if ( lline > 14 )
   f{3} = strtrim( line(15:min( lline, 24 ) ) );
   if( ~isempty( f{3} ) && f{3}(1) == '$' && ( lline < 25 || ~strcmp(line(15:15),'$-PARAMETER' ) ) )
     f{3} = '';
     nf = 2;
     return
   end
end
if( lline > 24 )
   f{4} = strtrim( line(25:min( lline, 36 ) ) );
   if( ~isempty( f{4} ) && f{4}(1) == '$' && ( lline < 35 || ~strcmp(line(25:35),'$-PARAMETER' ) ) )
     f{4} = '';
     nf = 3;
     return
   end
end
if( lline > 39 )
   f{5} = strtrim( line(40:min( lline, 49 ) ) );
   if( ~isempty( f{5} ) && f{5}(1) == '$' && ( lline < 50 || ~strcmp(line(40:50),'$-PARAMETER' ) ) )
     f{5} = '';
     nf = 4;
     return
   end
end
if( lline > 49 )
   f{6} = strtrim( line(50:min( lline, 61 ) ) );
   if( ~isempty( f{6} ) && f{6}(1) == '$' && (lline < 60 || ~strcmp(line(50:60),'$-PARAMETER' ) ) )
     f{6} = '';
     nf = 5;
     return
   end
end
nf = length( f );

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ f, nf ] = parselgroupline( line, lline )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Parse a line in the nonlinear ELEMENTS and GROUPS sections into its up to 4 different fields.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = {};
if ( lline < 3 )
   return
end

f{1} = strtrim( line(2:3) );
if ( lline > 4 )
   f{2} = strtrim( line( 5:min( lline, 14 ) ) );
end
if ( lline > 14 )
   f{3} = strtrim( line(15:min( lline, 24 ) ) );
end
if( lline > 24 )
   f{4} = strtrim( line(25:min( lline, 65 ) ) );
end
nf = length( f );

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = printelgroupline( f, nf, nline )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Print a line in the nonlinear elements and groups sections
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch( nf )
case 1
   fprintf( ' line %d:   %-2s\n', nline, f{1} );
case 2
   fprintf( ' line %d:   %-2s %-10s\n', nline, f{1}, f{2} );
case 3
   fprintf( ' line %d:   %-2s %-10s%-10s\n', nline, f{1}, f{2}, f{3} );
case 4
   fprintf( ' line %d:   %-2s %-10s%-10s%-40s\n', nline, f{1}, f{2}, f{3}, f{4} );
end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ form, act ] = formact( f1 )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Find the type and type of action for a line, depending on its first field f1.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lf1 = length( f1 );
if ( lf1 )
   if ( f1(1) == 'X' || f1(1) == 'Z' )
      form = 'a';
      if( length( f1) > 1 )
         act  = f1(2);
      else
         act = '';
      end
   else
      form = '';
      act  = f1;
   end
else  
   form = '';
   act  = '';
end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function value = getv1( form, f, pbs )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Get the first data value, irrespective of form.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( ~isempty( form ) && form(1) == 'Z' )
   value = s2mpjvalue( '', f{5}, pbs );
else
   value = d2e( f{4} );
end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function value = getv1r( form, f, pbs )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Get the first data value, irrespective of form, making sure the output string can be interpreterd as a real.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( ~isempty( form ) && form(1) == 'Z' )
   value = s2mpjvalue( '', f{5}, pbs );
else
   value = d2e( f{4} );
   if ( isempty( strfind( value, '.' ) ) )
      value = [value, '.' ];
   end
end
return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function value = getv2( f )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Get the second data value from field 6: it must be numerical.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

value = d2e( f{6} );

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function value = getv2r( f )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Get the second data value from field 6: it must be numerical. Make sure the output string can be interpreterd as a real.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

value = d2e( f{6} );
if ( isempty( strfind( value, '.' ) ) )
   value = [value, '.' ];
end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function string = from_fortran( string, dict, pbs )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Transform a Fortran expression to one (hopefully) acceptable for S2MPJ, depending on the output file language.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

string = regexprep( string, '\s*','' );           % remove blanks
string = d2ef( string );                          % transform exponent notation in Fortran
string = replace( string, '+-', '-'   );
string = replace( string, '= +', '= ' );

switch ( pbs.lang )

case { 'matlab', 'julia' }

   string = replace( string, '**', '^'  );           % replace ** by ^
   string = replace( string, 'MOD('   , 'mod('   ); 
   string = replace( string, 'DABS('  , 'abs('   ); 
   string = replace( string, 'ABS('   , 'abs('   );  % lower case functions
   string = replace( string, 'DSQRT(' , 'sqrt('  );
   string = replace( string, 'SQRT('  , 'sqrt('  );
   string = replace( string, 'DEXP('  , 'exp('   );
   string = replace( string, 'EXP('   , 'exp('   );
   string = replace( string, 'DLOG10(', 'log10(' );
   string = replace( string, 'LOG10(' , 'log10(' );
   string = replace( string, 'DLOG('  , 'log('   );
   string = replace( string, 'LOG('   , 'log('   );
   string = replace( string, 'DASIN(' , 'asin('  );
   string = replace( string, 'ASIN('  , 'asin('  );
   string = replace( string, 'DACOS(' , 'acos('  );
   string = replace( string, 'ACOS('  , 'acos('  );
   string = replace( string, 'DATAN(' , 'atan('  );
   string = replace( string, 'ATAN('  , 'atan('  );
   string = replace( string, 'DSINH(' , 'sinh('  );
   string = replace( string, 'SINH('  , 'sinh('  );
   string = replace( string, 'DCOSH(' , 'cosh('  );
   string = replace( string, 'COSH('  , 'cosh('  );
   string = replace( string, 'TANH('  , 'tanh('  );
   string = replace( string, 'DTANH(' , 'tanh('  );
   string = replace( string, 'DSIN('  , 'sin('   );
   string = replace( string, 'SIN('   , 'sin('   );
   string = replace( string, 'DCOS('  , 'cos('   );
   string = replace( string, 'COS('   , 'cos('   );
   string = replace( string, 'DTAN('  , 'tan('   );
   string = replace( string, 'TAN('   , 'tan('   );
   string = replace( string, 'DMAX('  , 'max('   );
   string = replace( string, 'AMAX('  , 'max('   );
   string = replace( string, 'MAX('   , 'max('   );
   string = replace( string, 'DMIN('  , 'min('   );
   string = replace( string, 'AMIN('  , 'min('   );
   string = replace( string, 'MIN('   , 'min('   );
   string = replace( string, '.LE.'   , '<='     );
   string = replace( string, '.LT.'   , '<'      );
   string = replace( string, '.EQ.'   , '=='     );
   string = replace( string, '.GE.'   , '>='     );
   string = replace( string, '.GT.'   , '>'      );
   string = replace( string, '.AND.'  , '&&'     );
   string = replace( string, '.OR.'   , '||'     );
   switch ( pbs.lang )
   case 'matlab' 
      string = replace( string, 'DATAN2(', 'atan2(' );
      string = replace( string, 'ATAN2(' , 'atan2(' );
      string = replace( string, '.NOT.'  , '~'      );
   case 'julia'
      string = replace( string, 'DATAN2(', 'atan(' );
      string = replace( string, 'ATAN2(' , 'atan(' );
      string = replace( string, '.NOT.'  , '!'     );
   end
   
case 'python'

   string = replace( string, 'mod('  ,  'np.mod(' );
   string = replace( string, 'MOD('  ,  'np.mod(' );

   string = replace( string, 'abs('  ,  'np.absolute#(' );
   string = replace( string, 'DABS(' ,  'np.absolute#(' );
   string = replace( string, 'ABS('  ,  'np.absolute#(' );
   
   string = replace( string, 'sqrt(' ,  'np.sqrt#('     );
   string = replace( string, 'DSQRT(',  'np.sqrt#('     );
   string = replace( string, 'SQRT(' ,  'np.sqrt#('     );
   
   string = replace( string, 'exp('  ,  'np.exp#('      );
   string = replace( string, 'DEXP(' ,  'np.exp#('      );
   string = replace( string, 'EXP('  ,  'np.exp#('      );
   
   string = replace( string, 'log10(',  'np.log10#('    );
   string = replace( string, 'DLOG10(', 'np.log10#('    );
   string = replace( string, 'LOG10(',  'np.log10#('    );

   string = replace( string, 'log('  ,  'np.log#('      );
   string = replace( string, 'DLOG(' ,  'np.log#('      );
   string = replace( string, 'LOG('  ,  'np.log#('      );

   string = replace( string, 'asin(' ,  'np.arcsin#('   );
   string = replace( string, 'DASIN(',  'np.arcsin#('   );
   string = replace( string, 'ASIN(' ,  'np.arcsin#('   );
   
   string = replace( string, 'acos(' ,  'np.arccos#('   );
   string = replace( string, 'DACOS(',  'np.arccos#('   );
   string = replace( string, 'ACOS(' ,  'np.arccos#('   );
   
   string = replace( string, 'atan(' ,  'np.arctan#('   );
   string = replace( string, 'DATAN(',  'np.arctan#('   );
   string = replace( string, 'ATAN(' ,  'np.arctan#('   );
   
   string = replace( string, 'atan2(',  'np.arctan2#('  );
   string = replace( string, 'DATAN2(', 'np.arctan2#('  );
   string = replace( string, 'ATAN2(',  'np.arctan2#('  );
   
   string = replace( string, 'sinh(' ,  'np.sinh#('     );
   string = replace( string, 'DSINH(',  'np.sinh#('     );
   string = replace( string, 'SINH(' ,  'np.sinh#('     );

   string = replace( string, 'cosh(' ,  'np.cosh#('     );
   string = replace( string, 'DCOSH(',  'np.cosh#('     );
   string = replace( string, 'COSH(' ,  'np.cosh#('     );

   string = replace( string, 'tanh(' ,  'np.tanh#('     );
   string = replace( string, 'DTANH(',  'np.tanh#('     );
   string = replace( string, 'TANH(' ,  'np.tanh#('     );
   
   string = replace( string, 'sin('  ,  'np.sin#('      );
   string = replace( string, 'DSIN('  , 'np.sin#('      );
   string = replace( string, 'SIN('  ,  'np.sin#('      );
   
   string = replace( string, 'cos('  ,  'np.cos#('      );
   string = replace( string, 'DCOS('  , 'np.cos#('      );
   string = replace( string, 'COS('  ,  'np.cos#('      );
   
   string = replace( string, 'tan('  ,  'np.tan#('      );
   string = replace( string, 'DTAN(' ,  'np.tan#('      );
   string = replace( string, 'TAN('  ,  'np.tan#('      );
   
   string = replace( string, 'DMAX(' ,  'max#('         );
   string = replace( string, 'AMAX(' ,  'max#('         );
   string = replace( string, 'MAX('  ,  'max#('         );
   string = replace( string, 'DMIN(' ,  'min#('         );
   string = replace( string, 'AMIN(' ,  'min#('         );
   string = replace( string, 'MIN('  ,  'min#('         );
   string = replace( string, '.LE.'  ,  '<='            );
   string = replace( string, '.LT.'  ,  '<'             );
   string = replace( string, '.EQ.'  ,  '=='            );
   string = replace( string, '.GE.'  ,  '>='            );
   string = replace( string, '.GT.'  ,  '>'             );
   string = replace( string, '.NOT.' ,  'not'           );
   string = replace( string, '.AND.' ,  ' and '         );
   string = replace( string, '.OR.'   , ' or '          );
   string = replace( string, '#','' );

end

%  Check for the SIGN FORTRAN function with two arguments

string = replace( string, 'DSIGN(', 'SIGN(' );
ia = strfind( string, 'SIGN(' );
if ( ~isempty( ia ) )
   string = decode_sign( string, ia, pbs );
end

%  Substitute the internal references to reserved names
%  according to the supplied dictionary
%
%  First sort the dictionary by order of decreasing key length.

switch ( pbs.dicttype )
case 'native'
   ldict = numEntries( dict );
case 'custom'
   ldict = dict.Count ;
end
lend    = zeros( ldict, 1 );
thekeys = cellstr( keys( dict ) );
for i = 1:ldict
   lend(i)  = length( thekeys{i} );
end
[ ~, perm ] = sort( lend, 'descend' );

%  Perform the substitution(s) in two passes in order to avoid recursive loops.

for i = 1:ldict
   ii     = perm(i);
   string = replace( string, thekeys(ii), ['#{',int2str( ii ),'}']);
end
for i = 1:ldict
   string = replace( string, ['#{',int2str( i ),'}'], dict( thekeys{ i } ) );
end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ string, poss, posb ] = d2e( string )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Reformats numbers in string,.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Change D exponent notation to e.

[ poss, posb ] = regexp( string, '[0-9]*\.[0-9]*[Dd][-|\+]?[0-9]{1,2}' );
if ( ~isempty( poss ) )
   posD = regexp( string( poss(1):posb(1)), '(D|d)' );
   string(poss(1)+posD(1)-1) = 'e';
   exps = regexp( string( poss(1):posb(1)), 'e[-|\+]?00' );
   if ( ~isempty( exps ) )
      string = [ string(1:poss(1)+exps-2), string(posb(1)+1:end) ];
   end
   string = d2e( string );
else

   %  Add a trailing 0 to reals

   if ( string(end) == "." )
      string = [string, '0' ];
   end

   %  Avoid numbers like 1234+04
  
   poss = regexp( string, '[0-9][-|\+][0-9]{1,2}' );
   if ( ~isempty( poss ) )
      string = [ string(1:poss), 'e', string(poss+1:end) ];
   end
end
string = replace( string, '.e0', '.0' );

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ string, poss, posb ] = d2ef( string )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Reformats numbers in string in Fortran
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Change D exponent notation to e.

[ poss, posb ] = regexp( string, '[0-9]*\.[0-9]*[Dd][-|\+]?[0-9]{1,2}' );
if ( ~isempty( poss ) )
   posD = regexp( string( poss(1):posb(1)), '(D|d)' );
   string(poss(1)+posD(1)-1) = 'e';
   exps = regexp( string( poss(1):posb(1)), 'e[-|\+]?00' );
   if ( ~isempty( exps ) )
      string = [ string(1:poss(1)+exps-2), string(posb(1)+1:end) ];
   end
   string = d2ef( string );
else

   %  Add a trailing 0 to reals
   
   if ( string(end) == "." )
      string = [string, '0' ];
   end
end
string = replace( string, '.e0', '.0' );

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function res = decode_sign( string, ia, pbs )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Decodes the old FORTRAN 77 SIGN function with two arguments.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opar = 1;
coma = 0;
iend = 0;
ss   = string(ia+5:end);
for i = 1:length( ss )
   switch( ss(i) )
   case ' '
      continue
   case '('
      opar = opar + 1;
   case ')'
      opar = opar - 1;
   case ','
      if ( opar == 1 )
         coma = i;
      end
   end
   if ( opar == 0 )
      iend = i;
      break
   end
end
arg1 = ss(1:coma-1);
ia1  = strfind( arg1, 'SIGN(');
if ( ~isempty( ia1 ) )
   arg1 = decode_sign( arg1, ia1, pbs );
end
arg2 =  ss(coma+1:iend-1);
ia2  = strfind( arg2, 'SIGN(');
if ( ~isempty( ia2 ) )
   arg2 = decode_sign( arg2, ia2, pbs );
end
switch( pbs.lang )
case 'matlab'
   res = [ string(1:ia-1), 'sign(', arg2, ')*(', arg1, ')',    string(ia+iend+5:end) ];
case 'python'
   res = [ string(1:ia-1), 'np.sign(', arg2, ')*(', arg1, ')', string(ia+iend+5:end) ];
 case 'julia'
   res = [ string(1:ia-1), 'sign(', arg2, ')*(', arg1, ')',    string(ia+iend+5:end) ];
end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function conclude_nlf( codeg, codeH, no_internals, indlvl, bindent, pbs )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Terminate writing the pending statements for the evaluation of nonlinear functions in the ELEMENTS and GROUPS sections.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  The gradient

if( ~isempty( codeg ) )

  % the test on nargout
  
  printmline( codeg{1},                                                                   indlvl,     bindent, pbs.fidma );
  printpline( codeg{1},                                                                   indlvl,     bindent, pbs.fidpy );
  printjline( codeg{1},                                                                   indlvl,     bindent, pbs.fidjl );
  for ig = 2:length( codeg )
     printmline( codeg{ig},                                                               indlvl + 1, bindent, pbs.fidma );
     printpline( codeg{ig},                                                               indlvl + 1, bindent, pbs.fidpy );
     printjline( codeg{ig},                                                               indlvl + 1, bindent, pbs.fidjl );
  end
  if ( no_internals )
     printmline( 'varargout{2} = g_;',                                                    indlvl + 1, bindent, pbs.fidma );
  else
     printmline( 'varargout{2} = U_.''*g_;',                                              indlvl + 1, bindent, pbs.fidma );
     printpline( 'g_ =  U_.T.dot(g_)',                                                    indlvl + 1, bindent, pbs.fidpy );
     printjline( 'g_ =  U_''*g_',                                                         indlvl + 1, bindent, pbs.fidjl );
 end

  %  The Hessian
  %
  %  Note that the initialization of the Hessian is made using the internal dimension as given by the number of internal
  %  gradient components, which is the length of codeg  minus pbs.gstrt.  This initailization is unnecessary in Matlab
  %  if the internal dimension is 1.

  if( ~isempty( codeH ) )
     dim = length( codeg ) - pbs.gstrt;                               % the length of the gradient
     switch( pbs.lang )
     case 'matlab'

        % the test on nargout
           
        printmline( codeH{1},                                                             indlvl + 1, bindent, pbs.fidma );
        if ( dim > 1 )
           printmline( sprintf( 'H_ = sparse(%d,%d);', dim, dim )',                       indlvl + 2, bindent, pbs.fidma );
        end
     case 'python'
        printpline( codeH{1},                                                             indlvl + 1, bindent, pbs.fidpy );
        printpline( sprintf( 'H_ = np.zeros((%d,%d))', dim, dim )',                       indlvl + 2, bindent, pbs.fidpy );
     case 'julia'
        printjline( codeH{1},                                                             indlvl + 1, bindent, pbs.fidjl );
        printjline( sprintf( 'H_ = zeros(Float64,%d,%d)', dim, dim )',                    indlvl + 2, bindent, pbs.fidjl );
     end

     for ih = 2:length( codeH )

        %  the specified Hessian entry
        
        printmline( codeH{ih},                                                            indlvl + 2, bindent, pbs.fidma );
        printpline( codeH{ih},                                                            indlvl + 2, bindent, pbs.fidpy );
        printjline( codeH{ih},                                                            indlvl + 2, bindent, pbs.fidjl );

        %  its symmetric (if non diagonal)

        switch( pbs.lang )
        case 'matlab'
           [ opar, cpar ] = regexp( codeH{ih}, '\(\S+\)' );
        case { 'python', 'julia' }
           [ opar, cpar ] = regexp( codeH{ih}, '\[\S+\]' );
        end
        if ( ~isempty( cpar ) && cpar(1) < strfind( codeH{ih}, '=' ) )  % no symmetric for H(1,1) in groups
           inds = split( codeH{ih}(opar(1)+1:cpar(1)-1), ',' );
           if ( ~strcmp( inds{1}, inds{2} ) )
              printmline( sprintf( 'H_(%s,%s) = H_(%s,%s);', inds{2}, inds{1}, inds{1}, inds{2} ), ...
                                                                                          indlvl + 2, bindent, pbs.fidma );
              printpline( sprintf( 'H_[%s,%s] = H_[%s,%s]',  inds{2}, inds{1}, inds{1}, inds{2} ), ...
                                                                                          indlvl + 2, bindent, pbs.fidpy );
              printjline( sprintf( 'H_[%s,%s] = H_[%s,%s]',  inds{2}, inds{1}, inds{1}, inds{2} ), ...
                                                                                          indlvl + 2, bindent, pbs.fidjl );
           end                          
        end
     end
     if ( no_internals )
        printmline( 'varargout{3} = H_;',                                                 indlvl + 2, bindent, pbs.fidma );
     else
        printmline( 'varargout{3} = U_.''*H_*U_;',                                        indlvl + 2, bindent, pbs.fidma );
        printpline( 'H_ = U_.T.dot(H_).dot(U_)',                                          indlvl + 2, bindent, pbs.fidpy );
        printjline( 'H_ = U_''*H_*U_',                                                    indlvl + 2, bindent, pbs.fidjl );
     end
     printmline( 'end',                                                                   indlvl + 1, bindent, pbs.fidma );
     printjline( 'end',                                                                   indlvl + 1, bindent, pbs.fidjl );
  end
  switch ( pbs.lang )
  case 'matlab'
     printmline( 'end',                                                                   indlvl,     bindent, pbs.fidma );
  case 'python'
     printpline( 'if nargout == 1:',                                                      indlvl,     bindent, pbs.fidpy );
     printpline( 'return f_',                                                             indlvl + 1, bindent, pbs.fidpy );
     printpline( 'elif nargout == 2:',                                                    indlvl,     bindent, pbs.fidpy );
     printpline( 'return f_,g_',                                                          indlvl + 1, bindent, pbs.fidpy );
     printpline( 'elif nargout == 3:',                                                    indlvl,     bindent, pbs.fidpy );
     printpline( 'return f_,g_,H_',                                                       indlvl + 1, bindent, pbs.fidpy );
  case 'julia'
     printjline( 'end',                                                                   indlvl,     bindent, pbs.fidjl );
     printjline( 'if nargout == 1',                                                       indlvl,     bindent, pbs.fidjl );
     printjline( 'return f_',                                                             indlvl + 1, bindent, pbs.fidjl );
     printjline( 'elseif nargout == 2',                                                   indlvl,     bindent, pbs.fidjl );
     printjline( 'return f_,g_',                                                          indlvl + 1, bindent, pbs.fidjl );
     printjline( 'elseif nargout == 3',                                                   indlvl,     bindent, pbs.fidjl );
     printjline( 'return f_,g_,H_',                                                       indlvl + 1, bindent, pbs.fidjl );
     printjline( 'end',                                                                   indlvl,     bindent, pbs.fidjl );
  end
end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ posact, posall ] = isactloop( ind, pbs )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  If ind is the index of an active loop, returns it position in pbs.actloop (in posact) and in pbs.loop (in posall).
%  Returns 0 otherwise.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

posact = 0;
posall = 0;
if ( isfield( pbs, 'loop' ) && ~isempty( pbs.loop ) )
   for i = 1:length( pbs.actloop )
       ii = pbs.actloop(i);
       if ( strcmp( ind, pbs.loop{ii}.index ) )
          posact = i;
          posall = ii;
       end
   end
end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pos = ismemstr( target, structarray, field )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Find the position, in a struct array, for the struct whose field 'field' is the string 'target'.  Returns 0 if none
%  is found.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pos = 0;
for i = 1:length( structarray )
   if ( isfield( structarray{i}, field ) )
      thefield = getfield( structarray{i}, field );
      if ( ischar( thefield ) )
         if ( strcmp( target, thefield ) )
            pos = i;
         end
      else
         [ ~, pos ] = ismember( target, thefield );
      end
      if ( pos )
         return
      end
   end
end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function printcpline( line, bindent, fidpy )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Prints a Python comment line, breaking the lines when possible to make it readable.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Formatting
if ( ~fidpy )
   return
end

width   = 100;                %  the target width of the Matlab file
defind  = '     ';            %  the additional indentation space on 
                              %  continuation after a break in a line
                              
%  If the line is short enough, print it as such.

if ( length( line ) <= width )
   fprintf( fidpy, '#%s\n', line );
   
%  The line is too long.: break at blanks.

else
   iblanks = find( line == ' ' );
   if ( length( iblanks ) > 1 && iblanks(end) > 80 )
      for i= 1:length( iblanks )
         if ( iblanks( i ) > width )
            breaks = iblanks(i-1);
            break
         end
      end
      fprintf( fidpy, '#%s\n', line(1:breaks-1) );
      printcpline( [ defind, line(breaks+1:end) ], bindent, fidpy );
   else  %  no blanks: print as such.
      fprintf( fidpy, '#%s\n', line );
   end
end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function printpline( line, indlvl, bindent, fidpy )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Prints a Python line, breaking the lines when possible to make it readable.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( ~fidpy )
   return
end
%  Formatting

width   = 80;                           %  the target width of the Matlab file
defind  = '     ';                      %  the additional indentation space on 
                                        %  continuation after a break in a line
indent  = repmat( bindent, indlvl, 1 ); %  the current indentation string
lindent = length( indent );             %  the amount of indentation
lenline = length( line );

%  Avoid the final semicolon on non-assignment Matlab statements.

sline   = strtrim(line);     %  remove leading blanks
if( isempty( sline) || length( sline ) < 3 )
   fprintf( fidpy, '\n' );
   return
end

%  If the line is short enough, print it without processing
%  (but with proper indentation).

if ( lindent+lenline <= width )

   fprintf( fidpy, '%s%s\n',  indent, line );
   
%  The line is too long.

else

   %  Search for a possible break point.  Break points are either at the =  sign in an assignment or at a position in the
   %  rhs which is a 'level 0' of bracketing.

   nopenpar = 0;       % the number of open paranthesis
   breaks   = 0;       % the position(s) of potential breakpoints
   poseq    = regexp( line, '[^<=>]=[^=]' ) + 1; % the first position after an = sign
   nopenstr = 0;       % the number of open strings
   for i = 1:lenline

      %  Track opening and closing paranthesis.
      
      switch( line(i) )
      case '('
         nopenpar = nopenpar + 1;
      case ')'
         nopenpar = nopenpar - 1;
      case { '"', '''' }
         if ( nopenstr )
            nopenstr = nopenstr - 1;
         else
            nopenstr = nopenstr + 1;
         end
      end

      %  Collect potential break points.

      if ( nopenstr == 0                             && ...  % no string open
           ~isempty( poseq )                         && ...  % there is an = sign
           i > poseq(1)                              && ...  % i is beyond it
           nopenpar == 0                             && ...  % there is no open paranthesis
           ( ismember( line(i), { '+', '-', '/' } )  && line(i-1) ~= 'e' ) && ...  % pos i is a suitable math break
           i-breaks(end) > 0.8*width                 && ...  % i is close enough to max length
           i < width                  )                      % but before it
         breaks(end+1) = i;
      end
   end

   %  Print the line using continuation paranthesis.
   
   if ( length( breaks ) > 1 )
      fprintf( fidpy, '%s%s = (%s\n', indent, line(1:poseq-2), line(poseq+2:breaks(2) ) );
      for j = 2:length( breaks ) - 1
         fprintf( fidpy, '%s%s%s\n', indent, defind, line(breaks(j)+1:breaks(j+1)) );
      end
      fprintf( fidpy, '%s%s%s)\n', indent, defind, line(breaks(end)+1:end ) );
      
   %  No break found: print the long line by breaking at the = sign, if possible.
   
   elseif ( poseq < 77 )
      fprintf( fidpy, '%s%s = (\n', indent, line(1:poseq-1) );
      fprintf( fidpy, '%s%s%s)\n', indent, defind, line(poseq+1:end) );
   else
      fprintf( fidpy, '%s%s\n',  indent, line );
   end
end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function printcjline( line, bindent, fidjl )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Prints a Julia comment line, breaking the lines when possible to make it readable.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Formatting

if ( ~fidjl )
   return
end

width   = 100;                %  the target width of the Matlab file
defind  = '     ';            %  the additional indentation space on 
                              %  continuation after a break in a line
                              
%  If the line is short enough, print it as such.

if ( length( line ) <= width )
   fprintf( fidjl, '#%s\n', line );
   
%  The line is too long.: break at blanks.

else
   iblanks = find( line == ' ' );
   if ( length( iblanks ) > 1 && iblanks(end) > 80 )
      for i= 1:length( iblanks )
         if ( iblanks( i ) > width )
            breaks = iblanks(i-1);
            break
         end
      end
      fprintf( fidjl, '#%s\n', line(1:breaks-1) );
      printcpline( [ defind, line(breaks+1:end) ], bindent, fidjl );
   else  %  no blanks: print as such.
      fprintf( fidjl, '#%s\n', line );
   end
end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function printjline( line, indlvl, bindent, fidjl )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Prints a Julia line, breaking the lines when possible to make it readable.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( ~fidjl )
   return
end

%  Formatting

width   = 80;                           %  the target width of the Matlab file
defind  = '     ';                      %  the additional indentation space on 
                                        %  continuation after a break in a line
indent  = repmat( bindent, indlvl, 1 ); %  the current indentation string
lindent = length( indent );             %  the amount of indentation
lenline = length( line );

%  Avoid the final semicolon on non-assignment Matlab statements.

sline   = strtrim(line);     %  remove leading blanks
if( isempty( sline) || length( sline ) < 3 )
   fprintf( fidjl, '\n' );
   return
end

%  If the line is short enough, print it without processing
%  (but with proper indentation).

if ( lindent+lenline <= width )

   fprintf( fidjl, '%s%s\n',  indent, line );
   
%  The line is too long.

else

   %  Search for a possible break point.  Break points are either at the =  sign in an assignment or at a position in the
   %  rhs which is a 'level 0' of bracketing.

   nopenpar = 0;       % the number of open paranthesis
   breaks   = 0;       % the position(s) of potential breakpoints
   poseq    = regexp( line, '[^<=>!]=[^=]' ) + 1; % the first position after an = sign
   nopenstr = 0;       % the number of open strings
   for i = 1:lenline

      %  Track opening and closing paranthesis.
      
      switch( line(i) )
      case '('
         nopenpar = nopenpar + 1;
      case ')'
         nopenpar = nopenpar - 1;
      case '"'
         if ( nopenstr )
            nopenstr = nopenstr - 1;
         else
            nopenstr = nopenstr + 1;
         end
      end

      %  Collect potential break points.

      if ( nopenstr == 0                             && ...  % no string open
           ~isempty( poseq )                         && ...  % there is an = sign
           i > poseq(1)                              && ...  % i is beyond it
           nopenpar == 0                             && ...  % there is no open paranthesis
           ( ismember( line(i), { '+', '-', '/' } )  && line(i-1) ~= 'e' ) && ...  % pos i is a suitable math break
           i-breaks(end) > 0.8*width                 && ...  % i is close enough to max length
           i < width                  )                      % but before it
         breaks(end+1) = i;
      end
   end

   %  Print the line using continuation paranthesis.
   
   if ( length( breaks ) > 1 )
      fprintf( fidjl, '%s%s = (%s\n', indent, line(1:poseq-2), line(poseq+2:breaks(2) ) );
      for j = 2:length( breaks ) - 1
         fprintf( fidjl, '%s%s%s\n', indent, defind, line(breaks(j)+1:breaks(j+1)) );
      end
      fprintf( fidjl, '%s%s%s)\n', indent, defind, line(breaks(end)+1:end ) );
      
   %  No break found: print the long line by breaking at the = sign, if possible.
   
   elseif ( poseq < 77 )
      fprintf( fidjl, '%s%s = (\n', indent, line(1:poseq-1) );
      fprintf( fidjl, '%s%s%s)\n',  indent, defind, line(poseq+1:end) );
   else
      fprintf( fidjl, '%s%s\n',  indent, line );
   end
end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function printcmline( line, bindent, fidma )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Prints a Matlab comment line, breaking the lines when possible to make it readable.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Do nothing if pbs.fidma is zero.

if ( ~fidma )
   return
end

%  Formatting

width   = 80;                 %  the target width of the Matlab file
defind  = '     ';            %  the additional indentation space on 
                              %  continuation after a break in a line
                              
%  If the line is short enough, print it as such.

if ( length( line ) <= width )
   fprintf( fidma, '%%%s\n', line );
   
%  The line is too long.: break at blanks, if any.

else
   iblanks = find( line == ' ' );
   if ( length( iblanks ) > 1 && iblanks(end) > 80 )
      for i= 1:length( iblanks )
         if ( iblanks( i ) > width )
            breaks = iblanks(i-1);
            break
         end
      end
      fprintf( fidma, '%%%s\n', line(1:breaks-1) );
      printcmline( [ defind, line(breaks+1:end) ], bindent, fidma );
   else  %  no blanks: print as such.
      fprintf( fidma, '%%%s\n', line );
   end
end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function printmline( line, indlvl, bindent, fidma )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Prints a Matlab line, breaking the lines when possible to make it readable.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Do nothing if pbs.fidma is zero.

if ( ~fidma )
   return
end

%  Formatting

width   = 80;                 %  the target width of the Matlab file
defind  = '     ';            %  the additional indentation space on 
                              %  continuation after a break in a line
indent  = repmat( bindent, indlvl, 1 ); %  the current indentation string
lindent = length( indent );   %  the amount of indentation
lenline = length( line );

%  Avoid the final semicolon on non-assignment Matlab statements.

sline   = strtrim(line);     %  remove leading blanks
if( isempty( sline) || length( sline ) < 3 )
   fprintf( fidma, '\n' );
   return
end

%  If the line is short enough, print it without processing
%  (but with proper indentation).

if ( lindent+lenline <= width )
   fprintf( fidma, '%s%s\n',  indent, line );
   
%  The line is too long.

else

   %  Search for a possible break point.  Break points are either at the = sign in an assignment or at a position in the
   %  rhs which is a 'level 0' of bracketing.

   nopenpar = 0;       % the number of open paranthesis
   breaks   = 0;       % the position(s) of potential breakpoints
   poseq    = regexp( line, '[^<=>]=[^=]' ) + 1; % the first position after an = sign
   for i = 1:lenline

      %  Track opening and closing paranthesis.
      
      switch( line(i) )
      case '('
         nopenpar = nopenpar + 1;
      case ')'
         nopenpar = nopenpar - 1;
      end

      %  Collect potential break points.
      
      if ( ~isempty( poseq )                             && ...  % there is an == sign
           i > poseq(1)                                  && ...  % i is beyond it
           nopenpar == 0                                 && ...  % there is no open paranthesis
           ( ismember( line(i), { '+', '-', '*', '/' } ) &&  line(i-1)~='e') && ...  % pos i is a suitable math break
           i-breaks(end) > 0.8*width                     && ...  % i is close enough to max length
           i < width                  )                          % but before it
         breaks(end+1) = i;
      end
   end

   %  Print the line using continuation lines.
   
   if ( length( breaks ) > 1 )
      fprintf( fidma, '%s%s...\n', indent, line(breaks(1)+1:breaks(2)) );
      for j = 2:length( breaks ) - 1
         fprintf( fidma, '%s%s%s...\n', indent, defind, line(breaks(j)+1:breaks(j+1)) );
      end
      fprintf( fidma, '%s%s%s\n', indent, defind, line(breaks(end)+1:end ) );
        
   %  No break found: print the long line by breaking at the = sign, if possible.
     
   elseif ( poseq < 77 )
      fprintf( fidma, '%s%s...\n', indent, line(1:poseq) );
      fprintf( fidma, '%s%s%s\n', indent, defind, line(poseq+1:end) );
   else
      fprintf( fidma, '%s%s\n',  indent, line );
   end

end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ thename, isidx, explindex ] = s2mpjname( nametype, sifname, pbs )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Return the name of a SIF entity specified by sifname as a string which can be evaluated unambiguously at runtime of the
%  output file (using the effective values of parameters and loop indeces). The flag nametype is used to distinguish names
%  of scalar integer or real parameters, as well as names occuring in nonlinear ELEMENTS and GROUPS section (in
%  Fortran-based statements), for which an unquoted naming convention is necessary to pass their values to s2mpjlib. The
%  s2mpj internal pbs structure is needed in input because knowledge of loop indeces is necessary to verify if sifname is
%  one of them. The flag isidx is returned true iff sifname is the name of an active loop' index, in which case its value
%  is well-defined in the output file and can be used in the s2mpjvalue function without recourse to the dictionary
%  of parameters. Finally, the flag explindex indicates if the name explictly depends on an index which is not that of
%  active loop.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Avoid some silly (but SIF compliant) confusion-inducing naming.

sifname = replace( sifname, ',(),', ',' );
sifname = replace( sifname, '()', '' );
sifname = replace( sifname, ',,', ',' );
sifname = replace( sifname, '_', 'u' );
sifname = replace( sifname, '''','p' );
switch( sifname )
case 'PB'
   sifname = 'loc_PB';
case 'PBM'
   sifname = 'loc_PBM';
end
thename = sifname;

%  Initialization

isidx     = 0;
explindex = 0;
opar      = [];

%  Identify scalar parameters.

scalar = ~isempty( nametype ) && ( nametype == 'R' || nametype == 'I' );

%  Check for opening and closing paranthesis.

if ( isempty( nametype ) || ~scalar )
   [ opar, cpar ] = regexp( sifname, '\(\S+\)' );
   if ( ~isempty( opar )  &&  length( opar ) ~= 1 )  %  more than one opening bracket
      pbs.errors{end+1} = sprintf( 'ERROR in line %d: more than 1 opening paranthesis in %s', pbs.nline, sifname );
      if ( pbs.disperrors )
         disp( pbs.errors{end} )
      end
      return
   end
   if ( ~isempty( cpar) &&  length( cpar ) ~= 1 )  %  more than one closing bracket
      pbs.errors{end+1} = sprintf( 'ERROR in line %d: more than 1 closing paranthesis in %s', pbs.nline, sifname );
      if ( pbs.disperrors )
         disp( pbs.errors{end} )
      end
      return
   end
end

%  The name has no indeces.

if (  scalar || isempty( opar ) || opar(1) == 1 )

   %  Check whether sifname refers to an active loop index.
   
   [ posact, posall ] = isactloop( sifname, pbs );
   if ( posact )
      thename = pbs.loop{posall}.index;
      isidx   = 1;
      return
   end

   %  A scalar name is returned. If this name is checked for in the nonlinear ELEMENTS or GROUPS sections, no surrounding
   %  quotes are needed because it occurs in FORTRAN-based statements. 

   switch( nametype )
   case 'u'
      thename = sifname;
   otherwise
      switch ( pbs.lang )
      case { 'matlab', 'python' }
         thename = [ '''', sifname, '''' ];    % first time this sifname is seen
      case 'julia'
         thename = [ '"', sifname, '"' ]; 
      end
   end

%  The name has indeces. 

else

   thename = sifname(1:opar-1);     % the name before the opening paranthesis

   %  Find the individual indeces. These must be either active loop indeces or integer parameters.
   
   theinds = split( sifname(opar+1:cpar-1), ',' );

   %  Check if any index is that of an active loop. If yes, prepare a variant of the name which can be used at runtime
   %  of the output file.
   
   for ii = 1:length( theinds )
      foundii = 0;
      [ posact, posall ] = isactloop( theinds{ii}, pbs );
      if ( posact )
         switch( pbs.lang )
         case 'matlab'
            pinds{ii} = sprintf( 'int2str(%s)', pbs.loop{posall}.index );
         case 'python'
            pinds{ii} = sprintf( 'str(%s)'    , pbs.loop{posall}.index );
         case 'julia'
            pinds{ii} = sprintf( 'string(%s)' , pbs.loop{posall}.index );
         end
         foundii = 1;
      end

      %  Index ii is not the index of an active loop.  Its (string) value must thus be available in the dictionary of
      %  integer parameters pbs.irpdict. Note that, because active loop indeces are checked for first, this allows
      %  an inactive loop index or standard integer parameter to be used as an active loop index.

      if ( ~foundii )
         switch ( pbs.lang )
         case 'matlab'
             pinds{ii} = sprintf( 'int2str(round(v_(%s)))',        pbs.irpdict( theinds{ii} ) );
         case 'python'
             pinds{ii} = sprintf( 'str(int(v_[%s]))'    ,        pbs.irpdict( theinds{ii} ) );
         case 'julia'
             pinds{ii} = sprintf( 'string(Int64(v_[%s]))' , pbs.irpdict( theinds{ii} ) );
         end
         explindex = 1;
      end
   end

   switch ( pbs.lang )
   case 'matlab'
      switch ( length( theinds ) )
      case 0
         thename = [ '''', thename, '''' ];
      case 1
         thename = [ '[',  '''', thename, '''',',', pinds{1},  ']' ];
      case 2 
         thename = [ '[',  '''', thename, '''',',', pinds{1}, ','','',', pinds{2}, ']' ];
      case 3
         thename = [ '[',  '''', thename, '''',',', pinds{1}, ','','',', pinds{2}, ','','',', pinds{3}, ']' ];
      end
   case 'python'
      switch ( length( theinds ) )
      case 0
         thename = [ '''', thename, '''' ];
      case 1
         thename = [ '''', thename, '''', '+', pinds{1}];
      case 2
         thename = [ '''', thename, '''', '+', pinds{1}, '+', ''',''', '+',pinds{2}  ];
      case 3
         thename = [ '''', thename, '''', '+', pinds{1}, '+', ''',''', '+', pinds{2}, '+', ''',''', '+', pinds{3} ];
      end
   case 'julia'
      switch ( length( theinds ) )
      case 0
         thename = [ '"', thename, '"' ];
      case 1
         thename = [ '"', thename, '"', '*', pinds{1} ];
      case 2
         thename = [ '"', thename, '"', '*', pinds{1}, '*', '","', '*', pinds{2} ];
      case 3
         thename = [ '"', thename, '"', '*', pinds{1}, '*', '","', '*', pinds{2}, '*', '","', '*', pinds{3} ];
      end
   end
end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pvalue = s2mpjvalue( nametype, sifname, pbs )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Returns a string defining the value of the sifname, which is either an active loop's index name, or the value of a
%  parameter extracted from the parameters' dictionary.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ thename, isidx ] = s2mpjname( nametype, sifname, pbs );
if ( isidx  )
   pvalue = thename;
else
   if ( isKey( pbs.irpdict, thename(2:end-1) ) )
      thename = char( pbs.irpdict( thename(2:end-1) ) );
   end
   switch ( pbs.lang )
   case 'matlab'
      pvalue = sprintf( 'v_(%s)', thename );
   case { 'python', 'julia' }
      pvalue = sprintf( "v_[%s]", thename );
   end
end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function initA( indlvl, bindent, pbs )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Initialize the three temporary vectors used to store the row and column indeces and values for A before
%  its sonstruction as a sparse matrix.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch ( pbs.lang )
case 'matlab'
   printmline( 'irA  = [];',                                                                  indlvl, bindent, pbs.fidma );
   printmline( 'icA  = [];',                                                                  indlvl, bindent, pbs.fidma );
   printmline( 'valA = [];',                                                                  indlvl, bindent, pbs.fidma );
case 'python'
   printpline( 'irA          = np.array([],dtype=int)',                                       indlvl, bindent, pbs.fidpy );
   printpline( 'icA          = np.array([],dtype=int)',                                       indlvl, bindent, pbs.fidpy );
   printpline( 'valA         = np.array([],dtype=float)',                                     indlvl, bindent, pbs.fidpy );
case 'julia'
   printjline( 'irA   = Int64[]',                                                             indlvl, bindent, pbs.fidjl );
   printjline( 'icA   = Int64[]',                                                             indlvl, bindent, pbs.fidjl );
   printjline( 'valA  = Float64[]',                                                           indlvl, bindent, pbs.fidjl );
end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function add2A( rowcol, idx, val, indlvl, bindent, pbs )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Write the code to add an entry of value val to the linear term, summing with an existing entry if present.
%  This can be done row or column-wise, depending on the value of rowcol. idx is the name of the considered
%  column/variable (rowcol = 'col') or row/group (rowcol = 'row').
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch( rowcol )
case 'row'
   printmline( 'icA(end+1) = iv;',                                                            indlvl, bindent, pbs.fidma );
   printmline( sprintf( 'irA(end+1) = ig_(%s);', idx ),                                       indlvl, bindent, pbs.fidma );
   printpline( 'icA  = np.append(icA,[iv])',                                                  indlvl, bindent, pbs.fidpy );
   printpline( sprintf( 'irA  = np.append(irA,[ig_[%s]])', idx),                              indlvl, bindent, pbs.fidpy );
   printjline( 'push!(icA,iv)',                                                               indlvl, bindent, pbs.fidjl );
   printjline( sprintf( 'push!(irA,ig_[%s])', idx ),                                          indlvl, bindent, pbs.fidjl );
case 'col'
   
   printmline( 'irA(end+1)  = ig;',                                                           indlvl, bindent, pbs.fidma );
   printmline( sprintf( 'icA(end+1)  = ix_(%s);', idx ),                                      indlvl, bindent, pbs.fidma );
   printpline( 'irA  = np.append(irA,[ig])',                                                  indlvl, bindent, pbs.fidpy );
   printpline( sprintf( 'icA  = np.append(icA,[ix_[%s]])', idx),                              indlvl, bindent, pbs.fidpy );
   printjline( 'push!(irA,ig)',                                                               indlvl, bindent, pbs.fidjl );
   printjline( sprintf( 'push!(icA,ix_[%s])',  idx ),                                         indlvl, bindent, pbs.fidjl );
end
printmline( sprintf( 'valA(end+1) = %s;' ,val ),                                              indlvl, bindent, pbs.fidma );
printpline( sprintf( 'valA = np.append(valA,float(%s))' ,val),                                indlvl, bindent, pbs.fidpy );
printjline( sprintf( 'push!(valA,Float64(%s))', val ),                                        indlvl, bindent, pbs.fidjl );

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function add2H( vname1, vname2, Hval, indlvl, bindent, pbs )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Write the code to add an entry of value val to the Hessian term, summing with an existing entry if present.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch( pbs.lang )
case 'matlab'
   printmline( sprintf( 'irH(end+1)  =  ix_(%s);',vname1 ),                                   indlvl, bindent, pbs.fidma );
   printmline( sprintf( 'icH(end+1)  =  ix_(%s);',vname2 ),                                   indlvl, bindent, pbs.fidma );
   printmline( sprintf( 'valH(end+1) =  %s;', Hval ),                                         indlvl, bindent, pbs.fidma );
   if ( ~strcmp( vname1, vname2 ) )
      printmline( sprintf( 'irH(end+1)  =  ix_(%s);',vname2 ),                                indlvl, bindent, pbs.fidma );
      printmline( sprintf( 'icH(end+1)  =  ix_(%s);',vname1 ),                                indlvl, bindent, pbs.fidma );
      printmline( sprintf( 'valH(end+1) =  %s;', Hval ),                                      indlvl, bindent, pbs.fidma );
   end
case 'python'
   printpline( sprintf( 'irH  = np.append(irH,[ix_[%s]])', vname1 ),                          indlvl, bindent, pbs.fidpy );
   printpline( sprintf( 'icH  = np.append(icH,[ix_[%s]])', vname2 ),                          indlvl, bindent, pbs.fidpy );
   printpline( sprintf( 'valH = np.append(valH,float(%s))', Hval ),                           indlvl, bindent, pbs.fidpy );
   if ( ~strcmp( vname1, vname2 ) )
      printpline( sprintf( 'irH  = np.append(irH,[ix_[%s]])', vname2 ),                       indlvl, bindent, pbs.fidpy );
      printpline( sprintf( 'icH  = np.append(icH,[ix_[%s]])', vname1 ),                       indlvl, bindent, pbs.fidpy );
      printpline( sprintf( 'valH = np.append(valH,float(%s))', Hval ),                        indlvl, bindent, pbs.fidpy );
   end
case 'julia'
   printjline( sprintf( 'push!(irH,ix_[%s])', vname1 ),                                       indlvl, bindent, pbs.fidjl );
   printjline( sprintf( 'push!(icH,ix_[%s])', vname2 ),                                       indlvl, bindent, pbs.fidjl );
   printjline( sprintf( 'push!(valH,Float64(%s))', Hval ),                                    indlvl, bindent, pbs.fidjl );
   if ( ~strcmp( vname1, vname2 ) )
      printjline( sprintf( 'push!(irH,ix_[%s])', vname2 ),                                    indlvl, bindent, pbs.fidjl );
      printjline( sprintf( 'push!(icH,ix_[%s])', vname1 ),                                    indlvl, bindent, pbs.fidjl );
      printjline( sprintf( 'push!(valH,Float64(%s))', Hval ),                                 indlvl, bindent, pbs.fidjl );
   end
end



return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function todefine = initIV( todefine, pbs, indlvl, bindent )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Defines yet undefined internal variables from the U_ matrix and the elemental variables, and unset the corresponding
%  entry of todefine to avoid doing this more than once.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for indiv = 1:length( todefine )
   if ( todefine( indiv ) )
      printmline( sprintf('IV_(%d) = U_(%d,:)*EV_;', indiv, indiv ),                          indlvl, bindent, pbs.fidma );
      printpline( sprintf('IV_[%d] = U_[%d:%d,:].dot(EV_)', indiv-1, indiv-1, indiv ),        indlvl, bindent, pbs.fidpy );
      printjline( sprintf('IV_[%d] = dot(U_[%d,:],EV_)', indiv, indiv ),                      indlvl, bindent, pbs.fidjl );
      todefine( indiv ) = 0;
   end
end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function okname = nlfname( sifname )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Find a name which can be used by Matlab and Python as a function name.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

okname = sifname;

%  Avoid names starting with a digit.

[ ib,ie ] = regexp( okname, '^[0-9]+' );
if ( ib )
   okname = [ 'n', okname ];
end

%  Avoid names with +,-,*,/

okname = replace( okname, '+', 'p' );
okname = replace( okname, '-', 'm' );
okname = replace( okname, '*', 't' );
okname = replace( okname, '/', 'd' );

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End of S2MPJ  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
