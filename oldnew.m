%
%  Compares the old and new versions of the
%  Matlab problems
%
%  Programming: Ph. Toint
%  This version :26 V 2024
%
addpath ./sif

start = 1070;
stop  =  Inf;

second_evaluation = 1;
eval_matvec       = 1;
printres          = 1;

filename = 'fullproblist';
fid = fopen(filename, "r");

nline = 0;
while( ~feof( fid ) )

   line = fgetl( fid );
   nline = nline + 1;
   if ( line(1) == '#' || nline < start || nline > stop )
      continue
   end
   fields = split( line );
   problem = fields{1};
   invalue = fields{2};
   args    = split( invalue, ',' );
   for i = 1:length(args)
      args{i} = str2double( args{i} );
   end
   disp( ['=== Checking problem ', int2str( nline ), ': ', problem ] )
   prob = str2func( problem );
   
   %%%%%%%%%%%%%%  new version %%%%%%%%%%%%%%%

   options.dicttype = 'custom';
   tic
   s2mpj( problem, options );
   ctime0 = toc;

   prob = str2func( problem );

   tic
   [ pb, pbm ] = prob( 'setup', args{:} );
   ctime1 = toc;

   ifixed = find( pb.xlower == pb.xupper );
   pb.x0( ifixed ) = pb.xlower( ifixed );
   y= ones(pb.m,1);
   
   tic
   [ Lnew, Lgnew, LHnew ] = prob( 'LgHxy', pb.x0, y );
   ctime2 = toc;
   if ( second_evaluation )
      tic
      [ f, g, H ] = prob( 'LgHxy', pb.x0, y);
      ctime3 = toc;
   end

   v= ones(pb.n,1);
   v(ifixed) = zeros( length(ifixed), 1 );

   if ( eval_matvec )
      tic
      Hvnew = prob( 'LHxyv', pb.x0, y, v );
      ctime4 = toc;
   end
   
   %%%%%%%%%%%%%%  old version %%%%%%%%%%%%%%%
   
   system( [ 'rm ', problem, '.m' ] );
   options.dicttype = 'native';
   tic
   s2mpj( problem, options );
   ntime0 = toc;

   prob = str2func( problem );

   tic
   [ pb, pbm ] = prob( 'setup', args{:} );
   ntime1 = toc;

   ifixed = find( pb.xlower == pb.xupper );
   pb.x0( ifixed ) = pb.xlower( ifixed );
   y= ones(pb.m,1);
   
   tic
   [ Lold, Lgold, LHold ] = prob( 'LgHxy', pb.x0, y );
   ntime2 = toc;
   if ( second_evaluation )
      tic
      [ f, g, H ] = prob( 'LgHxy', pb.x0, y);
      ntime3 = toc;
   end

   v= ones(pb.n,1);
   v(ifixed) = zeros( length(ifixed), 1 );

   if ( eval_matvec )
      tic
      Hvold = prob( 'LHxyv', pb.x0, y, v );
      ntime4 = toc;
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   fprintf( 1, '                 Custom      Native            Cust/Nat\n' );  
   fprintf( 1, '   Decode      : %f   %f seconds   %f\n', ctime0, ntime0, ctime0/ntime0 );
   fprintf( 1, '   Setup       : %f   %f seconds   %f\n', ctime1, ntime1, ctime1/ntime1 );
   fprintf( 1, '   Evaluation 1: %f   %f seconds   %f\n', ctime2, ntime2, ctime2/ntime2 );
   fprintf( 1, '   Evaluation 2: %f   %f seconds   %f\n', ctime3, ntime3, ctime3/ntime3 );
   fprintf( 1, '   H Mtavec    : %f   %f seconds   %f\n', ctime4, ntime4, ctime4/ntime4 );

   errf = abs( Lnew - Lold );
   errg = norm( Lgnew - Lgold );
   errH = norm( LHnew - LHold, 'fro' );
   errHv= norm( Hvnew - Hvold );
   fprintf( 1, '   Errors on f, g, H and Hv: %+8.3e %+8.3e %+8.3e %+8.3e\n', errf, errg, errH, errHv )
   system( [ 'rm ', problem, '.m' ] );
end

fclose(fid);

