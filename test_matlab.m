%
%  Test the setup and evaluation of Matlab versions of problems
%
%  Programming: Ph. Toint
%  This version : 19 V 2024
%
addpath ./matlab_problems

start = 1;
stop =  Inf;

second_evaluation = 1;
eval_matvec       = 1;
printres          = 1;

filename = 'fullproblist';
fid = fopen(filename, "r");

if ( printres )
   fidres = fopen('test_matlab.res','a');
%    fidres = 1:  % for STD output
end

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
   if ( printres )
   end
   prob = str2func( problem );

   tic
   [ pb, pbm ] = prob( 'setup', args{:} );
   time1 = toc;
   if ( printres )
      fprintf(fidres, '=== Checking problem  %s :  %s ( n = %d  m = %d )\n', int2str(nline), problem, pb.n, pb.m );
      for j = 1:5:pb.n
         fprintf(fidres, 'x0 = %+18.15e %+18.15e %+18.15e %+18.15e %+18.15e', pb.x0(j:min(j+4,pb.n)) );
         fprintf(fidres, '\n' );
      end
   end
   y = ones( pb.m, 1 );
   fprintf( 1, '   Matlab setup      : %f seconds\n', time1 );

   ifixed = find( pb.xlower == pb.xupper );
   pb.x0( ifixed ) = pb.xlower( ifixed );
   
   tic
   [ f, g, H ] = prob( 'LgHxy', pb.x0, y );
   time2 = toc;
   if ( printres )
      fprintf(fidres,'L0 = %+18.15e\n', f );
      for j = 1:5:pb.n
         fprintf(fidres, 'g0 = %+18.15e %+18.15e %+18.15e %+18.15e %+18.15e', g(j:min(j+4,pb.n)) );
         fprintf(fidres, '\n' );
      end
   end
   fprintf( 1, '   Matlab evaluation : %f seconds\n', time2 );
   if ( second_evaluation )
      tic
      [ f, g, H ] = prob( 'LgHxy', pb.x0, y);
      time3 = toc;
      fprintf( 1, '   Matlab evaluation : %f seconds\n', time3 );
   end

   v= ones(pb.n,1);
   v(ifixed) = zeros( length(ifixed), 1 );

   if ( eval_matvec )
      tic
      Hv = prob( 'LHxyv', pb.x0, y, v );
      time4 = toc;
      if ( printres )
         for j = 1:5:pb.n
            fprintf(fidres, 'Hv = %+18.15e %+18.15e %+18.15e %+18.15e %+18.15e', Hv(j:min(j+4,pb.n)) );
            fprintf(fidres, '\n' );
         end
      end
      fprintf( 1, '   Matlab prod H*v   : %f seconds\n', time4 );
%      fprintf( 1, '   Verif matvec    = %f\n', norm(H*v-Hv) );
%      fprintf( 1, '   Matlab err on prod: %f \n', norm( Hv - H*v) );
   end
   
end
if ( printres )
   fclose(fidres);
end
fclose(fid);

