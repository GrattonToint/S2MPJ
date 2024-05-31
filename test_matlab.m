%
%  Test the setup and evaluation of Matlab versions of problems
%
%  Programming: Ph. Toint
%  This version : 15 V 2024
%

start = 1;
stop =  10000;

second_evaluation = 1;
eval_matvec       = 1;

%filename = 'problist';
%filename = 'problist2';
%filename = 'problist3';
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

   tic
   [ pb, pbm ] = prob( 'setup', args{:} );
   y = ones( pb.m, 1 );
   time1 = toc;
   fprintf( '   Matlab setup      : %f seconds\n', time1 );
   tic
   [ f, g, H ] = prob( 'LgHxy', pb.x0, y );
   time2 = toc;
   fprintf( '   Matlab evaluation : %f seconds\n', time2 );
   if ( second_evaluation )
      tic
      [ f, g, H ] = prob( 'LgHxy', pb.x0, y);
      time3 = toc;
      fprintf( '   Matlab evaluation : %f seconds\n', time3 );
   end

   v= ones(pb.n,1);

   if ( eval_matvec )
      tic
      Hv = prob( 'LHxyv', pb.x0, y, v );
      time4 = toc;
      fprintf( '   Matlab prod H*v   : %f seconds\n', time4 );
%      fprintf( '   Verif matvacec    = %f\n', norm(H*v-Hv) );
%      fprintf( '   Matlab err on prod: %f \n', norm( Hv - H*v) );
   end
   
end
