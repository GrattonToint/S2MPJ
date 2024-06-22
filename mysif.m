%
%   Depending on the action argument:
%     action = "build' : copies the official SIF files as determined by fullproblsit
%                        from ~/problems/sif to the "private" version in ~/problems/mysif
%     action = 'check' : checks if official and private SIF files (in fullproblist) differ
%
%   Programming: Ph. Toint
%   This version 21 VI 2024
%

function mysif( action )

start = 1;

theproblems = 'fullproblist';

sifdir   = '/home/philippe/problems/sif/';
mysifdir = '/home/philippe/problems/mysif/';

fid = fopen( theproblems, 'r' );
while( ~feof(fid ) )
   line   = fgetl( fid );
   fields = split( line );
   if ( fields{1}(1) == '#' )
      continue
   end
   pname  = fields{1};
   sifname = replace( pname, 'n','' );
   sifname = replace( sifname,'p','+' );
   sifname = replace( sifname,'m','-' );
   sifname = replace( sifname,'t','*' );
   sifname = replace( sifname,'d','/' );
   ofile   = [ sifdir, sifname, '.SIF'];
   myfile  = [ mysifdir, sifname, '.SIF'];
   if ( exist( myfile ) )
      is_diff = system([ 'diff -q ', ofile, ' ', myfile, ' > /dev/null' ] );
      if ( is_diff )
         switch( action )
         case 'build'
            disp( [ 'Copying ', ofile, ' to ', myfile ] )
            system( [ 'cp ', ofile, ' ', myfile ] );
         case 'check'
            disp( [ 'files ', ofile, ' and ', myfile,' are different.'] )
         end
      end
   else
      switch ( action )
      case { 'build', 'check' }
         disp( [ 'Copying ', ofile, ' to ', myfile ] )
         system( [ 'cp ', ofile, ' ', myfile ] );
      end
   end
end
fclose(fid);
disp( 'Done.' )

return

end
