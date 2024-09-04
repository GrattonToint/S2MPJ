%
%  list the indeces of unfixed variables for all problems
%
%  Programming: Ph. Toint
%  This version : 18 VI 2024
%
addpath ./matlab_problems

start = 221
stop =  10000;

filename = 'fullproblist';
fid = fopen(filename, "r");

fidres = fopen('fixed.lst','a');

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
   pb     = prob( 'setup', args{:} );
   fixed  = find( pb.xlower == pb.xupper );
   nfixed = length( fixed );
   fprintf( fidres, '=== Checking problem  %s :  %s ( n = %d  m = %d  nfixed = %d )\n', ...
                    int2str(nline), problem, pb.n, pb.m, nfixed );
   if ( nfixed )                    
      for j = 1:20:nfixed
         fprintf( fidres, 'fix = %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d', ...
         fixed(j:min(j+19,nfixed)) );
         fprintf( fidres, '\n' );
      end
   else
      fprintf( fidres, 'fix =\n' );
   end
end
fclose(fidres);
fclose(fid);

