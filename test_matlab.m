%
%  Test the setup and evaluation of Matlab versions of problems
%
%  Programming  : Ph. Toint
%  This version : 2 IX 2024
%
addpath ./matlab_problems

start = 1;
stop =  Inf;

second_evaluation = 1;
eval_matvec       = 1;
printtimes        = 1;   % print timings in test_matlab.times
printres          = 1;   % print values of x0, L0, g0, and Hv in test_matlab.data

filename = 'fullproblist';
fid = fopen(filename, "r");

if ( printres )
   fidres = fopen('test_matlab.data','a');
%    fidres = 1:  % for STD output
end
if ( printtimes )
   fidtimes = fopen('test_matlab.times','a');
%    fidtimes = 1:  % for STD output
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
   if ( printtimes )
      fprintf(fidtimes, '=== Checking problem  %s : %s \n', int2str(nline), problem);
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
   if ( printtimes )
      fprintf( fidtimes, '   Matlab setup      : %f seconds\n', time1 );
   end
   
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
   if ( printtimes )
      fprintf( fidtimes, '   Matlab evaluation : %f seconds\n', time2 );
   end
   if ( second_evaluation )
      tic
      [ f, g, H ] = prob( 'LgHxy', pb.x0, y);
      time3 = toc;
      fprintf( 1, '   Matlab evaluation : %f seconds\n', time3 );
      if ( printtimes )
         fprintf( fidtimes, '   Matlab evaluation : %f seconds\n', time2 );
      end
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
      if ( printtimes )
         fprintf( fidtimes, '   Matlab prod H*v   : %f seconds\n', time4 );
      end
   end
   
end
if ( printres )
   fclose(fidres);
end
if ( printtimes )
   fclose(fidtimes)
end

fclose(fid);

