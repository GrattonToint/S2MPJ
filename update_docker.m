%
%   Update the docker files
%
%   Programming: Ph. Toint
%   This version 21 V 2024
%

function update_docker()

   system( [ 'cp ./fullproblist /home/philippe/docker/mycute_2/newpy' ] );
   system( [ 'cp ./s2mpjlib.py /home/philippe/docker/mycute_2/newpy' ] );
   system( [ 'cp ./test_fortran_python.py /home/philippe/docker/mycute_2/newpy' ] );

   theproblems = 'fullproblist';
   
   fid = fopen(theproblems, 'r' );x
   while( ~feof(fid ) )
       line   = fgetl( fid );
       fields = split( line );
       if ( fields{1}(1) ~= '#' )
          problem = fields{1};
       end
       system( [ 'cp ./sif/',problem,'.SIF /home/philippe/docker/mycute_2/newpy/' ] );
       system( [ 'cp /python_problems/', problem, '.py /home/philippe/docker/mycute_2/newpy' ] );
   end
   fclose(fid);

   disp( 'Done.' )

return

end
