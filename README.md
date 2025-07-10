# S2MPJ

Thank you for your interest in S2MPJ, a collection of test problems for benchmarking
optimization software, available in Matlab, Python and Julia. The collection contains
most of the problems in the CUTEst collection, making them accessible without any use
of Fortran or translation from Fortran.

The collection provide simple tools to set up a problem's data structure and to evaluate
its objective function and/or constraints as well as their derivatives.

A full description of the collection and associated tools is available in the paper
's2mpj.pdf'.

## NOTE:

Due to a limitation on the length of displayed directory content in Github, *only a
fraction of the test problems are listed* when showing the content of the matlab_problems,
python_problems and julia_problems directories on the Github web site. The lists of all
available problems can be accessed in `list_of_matlab_problems`, `list_of_python_problems` and
`list_of_julia_problems`, respectively.  The complete set of test problems can be obtained
by downloading `all_matlab_problems.tar.gz`, `all_python_problems.tar.gz` and/or
`all_julia_problems.tar.gz`.

## Using the collection in Matlab

To use S2MPJ in Matlab, download the contents of this repository,
open Matlab and make sure the local directory containing the 's2mpjlib.m' file  and its
subdirectory 'matlab_problems' are in the Matlab path.

Suppose now  we have an optimization code MyOptimizer.m which we would like to apply to
one on the S2MPJ problems, let's say `PROBLEM.m`.  We
first setup the problem data structure, along with the
starting point and bounds on the variables by issuing the command
``` 
[ pb, pbm ] = PROBLEM( 'setup', args{:} );
```      
at the beginning of MyOptimizer.m, before any call to function/constraint
values or derivatives.  In this call, `args{:}` is an optional comma-separated
list of  problem-dependent arguments (such as problem's dimension, for
instance) identified by the string "$-PARAMETER" in the SIF file. If `args{:}`
is missing, setup is performed using the SIF file defaults. The starting point
is now given by pb.x0, the lower bounds by pb.xlower and the upper bounds by pb.xupper.
   
Every time values of the objective function, and possibly of its derivatives,
at the vector `x` are needed in MyOptimizer.m, they are obtained by
issuing one of the commands
```  
fx = PROBLEM( 'fx', x );
[ fx, gx ] = PROBLEM( 'fgx', x );
[ fx, gx, Hx ]  = PROBLEM( 'fgHx', x );
```      
The product of the objective function's Hessian (at `x`) times a
user-supplied vector `v` can be obtained by the command
```  
Hxv = PROBLEM( 'fHxv', x, v );
```      
If the problem has constraints (other than bounds on the variables),
their values (and that of derivatives, if desired) are obtained by
one of the commands
```
cx = PROBLEM( 'cx', x );
[ cx, Jx ] = PROBLEM( 'cJx', x );
[ cx, Jx, Hx ]  = PROBLEM( 'cJHx', x );
```      
The product of the constraints' Jacobian (at `x`) times `v` is computed
by issuing the  command
```
Jxv  = PROBLEM( 'cJxv', x, v );
```      
and the product of the Jacobian transposed (at `x`) times `v` is obtained by
```   
Jxv  = PROBLEM( 'cJtxv', x, v );
```
The value (and derivatives) of the Lagrangian function L(x,y) at `(x,y)`
is obtained by one of the commands
```   
Lxy = PROBLEM( 'Lxy', x, y );}
[ Lxy, Lgxy ] = PROBLEM( 'Lgxy', x, y );
[ Lxy, Lgxy, LHxy ] = PROBLEM( 'LgHxy', x, y );
```      
while the product of the Lagrangian's Hessian at `(x,y)` times a vector `v`
can be computed with the command
```   
HLxyv = PROBLEM( 'HLxyv', x, y, v );
```
The list of Matlab test problems is available in the "list_of_matlab_problems" file.

## Using the collection in Python

The use of the collection in Python is very similar to that described for Matlab.
After downloading the content of the repository, make sure the directory containing
the file "s2mpjlib.py" and the subdirectory 'python_problems" are in sys.path.
Consider now using the problem `PROBLEM.py` in an optimization code.
After importing the functions of the PROBLEM module by
```
from PROBLEM import *
```
the setup of the problem is now performed by a call of the form
```
PROB = PROBLEM( args[:] )
```
where `args[:]` is an optional comma separated list of problem-dependent
arguments. The subsequent evaluation tasks are then obtained by
calling one (or more) of
```
fx          = PROB.fx( x )
fx, gx      = PROB.fgx( x )
x, gx, Hx   = PROB.fgHx( x )
Hxv         = PROB.fHxv( x, v )
cx          = PROB.cx( x )
cx, Jx      = PROB.cJx( x )
cx, Jx, HXs = PROB.cJHx( x )
cJxv        = PROB.cJxv( x, v )
cJtxv       = PROB.cJtxv( x, v )
cIx         = PROB.cIx( x, I )
cIx, cIJx   = PROB.cIJx( x, I )
cIx, cIJx, cIJHx = PROB.cIJHx( x, I )
cIJv        = PROB.cIJxv( x, v, I )
cIJtv       = PROB.cIJtxv( x, v, I )
Lxy         = PROB.Lxy( x, y )
Lxy, Lgxy   = PROB.Lgxy( x, y )
Lxy, gLxy, HLxy = PROB.LgHxy( x, y )
LHxyv       = PROB.LHxyv( x, y, v )
LIxy        = PROB.LIxy( x, y, I )
LIxy, gLIxy = PROB.LIgxy( x, y, I )
LIxy, gLIxy, HLIxy = PROB.LIgHxy( x, y, I )
LIHxyv      = PROB.LIHxyv( x, y, v, I )
```
The list of Python test problems is available in the "list_of_python_problems" file.

## Using the collection in Julia

To use of the problem `PROBLEM.jl` in an optimization code (after dowloading the
content of the repository), one first includes the
s2mpj Julia library and problem as in
```
include( "s2mpjlib.jl" )
include( "PROBLEM.jl" )
```
(possibly giving the full path to these files if they are not in the current directory).

NOTE: this differs from the obsolete procedure described in the published paper!

The setup of the problem is then performed by a call 
```
pb, pbm = PROBLEM( "setup", args[:] )
```
where `args[:]` is an optional comma separated list of
problem-dependent arguments. The subsequent evaluation tasks are
then obtained by calling one (or more) of 
```
fx          = PROBLEM( "fx", pbm, x )
fx, gx      = PROBLEM( "fgx", pbm, x )
fx, gx, Hx  = PROBLEM( "fgHx", pbm, x )
Hxv         = PROBLEM( "fHxv", pbm x, v )
cx          = PROBLEM( "cx", pbm, x )
cx, Jx      = PROBLEM( "cJx", pbm, x )
cx, Jx, HXs = PROBLEM( "cJHx", pbm,  x )
cJxv        = PROBLEM( "cJxv", pbm, x, v )
cJtxv       = PROBLEM( "cJtxv", pbm, x, v )
cIx         = PROBLEM( "cIx", pbm, x, I )
cIx, cIJx   = PROBLEM( "cIJx", pbm, x, I )
cIx, cIJx, cIJHx = PROBLEM( "cIJHx", pbm,  x, I )
cIJv        = PROBLEM( "cIJxv", pbm, x, v, I )
cIJtv       = PROBLEM( "cIJtxv", pbm, x, v, I )
Lxy         = PROBLEM( "Lxy", pbm, x, y )
Lxy, Lgxy   = PROBLEM( "Lgxy", pbm, x, y )
Lxy, Lgxy, LgHxy = PROBLEM( "LgHxy", pbm, x, y )
LHxyv       = PROBLEM( "LHxyv", pbm, x, y, v )
LIxy        = PROBLEM( "LIxy", pbm, x, y, I )
LIxy, LIgxy = PROBLEM( "LIgxy", pbm, x, y, I )
LIxy, LIgxy, LIgHxy = PROBLEM( "LIgHxy", pbm, x, y, I )
LIHxyv      = PROBLEM( "LIHxyv", pbm, x, y, v, I )
```
Note that `pbm` always occurs as the second input argument. 

The list of Julia test problems is available in the "list_of_julia_problems" file.

## Authors

We hope that S2MPJ will  prove useful in your benchmarking.  If you have
any bug reports or comments, please feel free to email one of the
toolbox authors:
```
  Serge Gratton <serge.gratton@enseeiht.fr>
  Philippe Toint <philippe.toint@unamur.be>
```

## Thanks

Many thanks to
```
Nick Gould, Zaikun Zhang, Cunxin Huang, Filippo Marini, Greta Malaspina, Mihai Costin
```
