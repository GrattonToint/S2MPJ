function [ err_gf, err_Hf, err_gc, err_Hc ] = test_fd ( probname, verbose, pb, pbm, atx0 )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Test the derivatives of the functions for problem probname in S2X format, using finite differences
%   in continuous  variables at a random point feasible for the problems bounds or at the problem's
%   starting point. If the problem's dimension is too large (> 15), select a random set of components
%   for testing gradients and Hessians.

%   INPUT:
%
%   probname : the name of the considered optimization problem
%   verbose  : 1 if output and the standard output is desired, 0 for silent execution
%   n        : the dimension
%   pb       : the pb struct as output from probname('setup')
%   pbm      : the pbm struct as output from probname('setup')
%   atx0     : 1 = performs the test at pb.x0, 0 = at a random point
%
%   OUTPUT:
%
%   err_gf : the max absolute value of the error on the objective function's gradient components
%   err_Hf : the max absolute value of the error on the objective function's Hessian components
%   err_gc : the max absolute value of the error on the constraints' function's gradient components
%   err_Hc : the max absolute value of the error on the constraints' function's Hessian components
%
%   Programming : Ph. Toint and S. Gratton, July 2018
%   This version: 21 II 2024
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxsize = 15;
maxcons = 20;

has_obj = ( isfield( pbm, 'objgrps') && ~isempty( pbm.objgrps ) ) || isfield( pbm, H );
has_c   = pb.m;

version = 'complex';
%version = 'real';

switch ( version )
case 'complex'
   ep = 1e-12;
case 'real'
   ep = 1e-7;
end

%  The problem's dimension

n = pb.n;

%  If n exceeds maxsize, select maxsize components at random between 1 and n.

if ( n > maxsize )
  tc = sort( randi( n, maxsize, 1 ) );
else
  tc = [1:n];
end
dim  = length( tc );

%  Choose a random bound-feasible evaluation point, if requested, or the problem's
%  starting point otherwise.

if ( atx0 )
   x = pb.x0;
else
   x = -ones(n,1)+2*rand( n, 1 );
end

%  Project the evaluation point onto the feasible set for bound constraints.

if ( ~isempty( pb.xlower ) )
   x  = max( [ x'; pb.xlower' ] )';
end
if ( ~isempty( pb.xupper ) )
   x  = min( [ pb.xupper';  x' ] )';
end

v = zeros( n, 1 );

%%%%%%%%%%%%%%%%%%%%%%%%% Objective function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Get the values of the objective function at x.

if ( has_obj )
   [ f, g, H ] = feval( probname, 'fx', x );

   %   Test its gradient.

   err_gf  = 0;
   err_gfi = zeros( n, 1 );
   for ii   = 1:dim
      i = tc( ii );
      if ( length( pb.xlower ) == 0 || abs( pb.xupper( i ) - pb.xlower( i ) ) > 1e-15 )
         v( i ) = 1;
         switch( version )
         case 'real'
            fp = feval( probname, 'fx', x + ep*v );
            err_gfi( i )   = abs( (1/ep)*(fp-f) - g(i) );
         case 'complex'
            fp = feval( probname, 'fx', x + ep*sqrt(-1)*v );
            err_gfi( i ) = abs( (1/ep)*imag( fp ) - g(i) );
         end
         v( i ) = 0;
      end
   end
   %[ errgmax, ierrgmax ] = max( err_gfi );

   if ( verbose )
      disp(    ' ' )
      disp(    '*** FD testing problem  ***' )
      disp(    ' ' )
      disp(    '    objective function : errors on the gradient' )
      disp(    ' ' )
      fprintf( '                   ' ); fprintf( '%0.1e  ', err_gfi(tc) )
      disp(    ' ')
   %   disp(    ' ' )
   %   fprintf( '    max = %0.1e for component %2d\n', errgmax, ierrgmax )
   %   disp(    ' ')
   end
   err_gf = norm( err_gfi, Inf );

   %   Test its Hessian.

   err_Hf = zeros( n, n );
   for ii = 1:dim
      i = tc( ii );
      if ( length( pb.xlower ) == 0 || abs( pb.xupper( i ) - pb.xlower( i ) ) > 1e-15 )
         for jj = 1:dim
            j = tc( jj );
            if ( length( pb.xlower ) == 0 || abs( pb.xupper( j ) - pb.xlower( j ) ) > 1e-15 )
               v( j ) = 1;
   	    switch( version )
   	    case 'real'
                  [ ~, gp ]      = feval( probname, 'fx', x + ep*v );
                  err_Hf( i, j ) = abs( (1/ep)*( gp(i) -g(i) ) - H(i,j) );
   	    case 'complex'
                  [ ~, gp ]      = feval( probname, 'fx', x + ep*sqrt(-1)*v );
                  err_Hf( i, j ) = abs( (1/ep)*imag( gp(i) ) - H(i,j) );
   	    end
               v( j ) = 0;
   	 end
         end
      end
   end
   %[ errHmax, ierrHmax, jerrHmax ] = max( err_Hf );

   if ( verbose )
      disp(    ' ' )
      disp(    '    objective function : errors on the Hessian' )
      disp(    ' ' )
      fprintf( '      cols ----->  ' ); fprintf( '%7d  ', tc ); fprintf( '\n' );
      for ii = 1:dim
         i = tc(  ii );
         fprintf( '      row %4d  :  ', i ); fprintf( '%0.1e  ', err_Hf(i,tc) ); fprintf( '\n' )
      end
   %   fprintf( '    max = %0.1e for component %2d, %2d\n', errHmax, ierrHmax, jerrHmax )
   %   disp(    ' ')
   end
   err_Hf = norm( err_Hf, Inf );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Constraints  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( ~has_c )
   err_gc = NaN;
   err_Hc = NaN;
   return
end

err_gc  = 0;
err_Hc  = 0;

[ cx, Jx, Hcx ]  = feval( probname, 'cx',  x );
ncon = length( cx );

%  If ncon exceeds maxcons, select maxcons constraints at random between 1 and ncon.

if ( ncon > maxcons )
  mc = sort( randi( ncon, maxcons, 1 ) );
else
  mc = [1:ncon];
end

for iicon = 1:length(mc)

   icon = mc(iicon );
   c    = cx(icon);
   gc   = Jx(icon,:);
   Hc   = Hcx{icon};

%  Test the gradient of the icon-th constraint

   err_gci = zeros( n, 1 );
   for ii  = 1:dim
      i = tc( ii );
      if ( abs( pb.xupper(i) - pb.xlower(i) ) > 1e-15 )
         v( i ) = 1;
	     switch( version )
	     case 'real'
            cp = feval( probname, 'cIx', x + ep*v, icon);
            err_gci( i ) = abs( (1/ep)*(cp-c) - gc(i) );
	     case 'complex'
            cp = feval( probname, 'cIx', x + ep*sqrt(-1)*v, icon );
            err_gci( i ) = abs( (1/ep)*imag( cp ) - gc(i) );
	     end
         v( i ) = 0;
      end
   end
   err_gc = max( err_gc, err_gci );

   if ( verbose )
      disp(    ' ' )
      disp(  [ '    constraint ', int2str( icon ),': errors on the gradient' ] )
      disp(    ' ' )
      fprintf( '                   ' ); fprintf( '%0.1e  ', err_gci(tc) )
      disp(    ' ' )
   end

%  Test its Hessian.

   err_Hci = zeros( n, n );
   for ii = 1:dim
      i = tc( ii );
      if ( abs( pb.xupper( i ) - pb.xlower( i ) ) > 1e-15 )
         for jj = 1:dim
            j = tc( jj );
            if ( abs( pb.xupper( j ) - pb.xlower( j ) ) > 1e-15 )
               v( j ) = 1;
	           switch( version )
	           case 'real'
                  [ ~, gcp ]      = feval( probname, 'cIx', x + ep*v, icon );
                  err_Hci( i, j ) = abs( (1/ep)*( gcp(i)-gc(i) ) - Hc(i,j) );
	           case 'complex'
                  [ ~, gcp ]      = feval( probname, 'cIx', x + ep*sqrt(-1)*v, icon );
                  err_Hci( i, j ) = abs( (1/ep)*imag( gcp(i) ) - Hc(i,j) );
	           end
               v( j ) = 0;
	        end
         end
      end
   end
   if ( verbose )
      disp(    ' ' )
      disp(  [ '    constraint ', int2str(icon),': errors on the Hessian' ] )
      disp(    ' ' )
      fprintf( '      cols ----->  ' ); fprintf( '%7d  ', tc ); fprintf( '\n' );
      for ii = 1:dim
         i = tc(ii );
         fprintf( '      row %4d  :  ', i ); fprintf( '%0.1e  ', err_Hci(i,tc) ); fprintf( '\n' )
      end
   end
   err_Hc = max( err_Hc, norm( err_Hci, Inf ) );

end

return

end
