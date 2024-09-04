%
%   Update the docker files
%
%   Programming: Ph. Toint
%   This version 21 VI 2024
%

function update_docker()

   system( [ 'cp ./fullproblist /home/philippe/docker/mycute_2/newpy' ] );
%   system( [ 'cp ./s2mpjlib.py /home/philippe/docker/mycute_2/newpy' ] );
%   system( [ 'cp ./test_fortran_python.py /home/philippe/docker/mycute_2/newpy' ] );
   system( [ 'cp ./test_fortran.py /home/philippe/docker/mycute_2/newpy' ] );

   theproblems = 'fullproblist';
   
   fid = fopen(theproblems, 'r' );
   while( ~feof(fid ) )
       line   = fgetl( fid );
       fields = split( line );
       if ( fields{1}(1) ~= '#' )
          problem = fields{1};
          system( [ 'cp ./python_problems/', problem, '.py /home/philippe/docker/mycute_2/newpy' ] );
          if ( ismember( problem, { 'n10FOLDTR', 'n10FOLDTRLS','n3PK' } ) )
             system( [ 'cp ./sif/',problem,'.SIF /home/philippe/docker/mycute_2/newpy/',problem(2:end), '.SIF' ] );
          else
             system( [ 'cp ./sif/',problem,'.SIF /home/philippe/docker/mycute_2/newpy/' ] );
          end
       end
   end
   fclose(fid);

   disp( 'Done.' )

return

end
