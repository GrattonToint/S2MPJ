%
%   Regenerate .m, .jl and .py files and copy the latter to the docker.
%
%   Programming: Ph. Toint
%   This version 28 V 2024
%

function regenerate( varargin )

addpath /home/philippe/s2x/sif
inpy.language = 'python';
injl.language = 'julia';

system( [ 'cp fullproblist /home/philippe/docker/mycute_2/newpy' ] );
system( [ 'cp s2xlib.py /home/philippe/docker/mycute_2/newpy' ] );
system( [ 'cp test_fortran_python.py /home/philippe/docker/mycute_2/newpy' ] );

if ( nargin )

   for i= 1:nargin
      pname = varargin{i};
      fprintf( ' Decoding problem %s\n', pname );
      thename = s2x( pname );
      thename = s2x( pname, inpy );
      thename = s2x( pname, injl );
      system( [ 'cp /home/philippe/s2x/sif/',pname,'.SIF /home/philippe/docker/mycute_2/newpy' ] );
      system( [ 'cp ', thename, '.py /home/philippe/docker/mycute_2/newpy' ] );
   end


else

   start = 1;

   theproblems = 'fullproblist';
   
   fid = fopen(theproblems, 'r' );
   problems = {};
   iprob = 0;
   while( ~feof(fid ) )
       line   = fgetl( fid );
       fields = split(line );
       iprob  = iprob + 1;
       if ( fields{1}(1) == '#' )
          if ( fields{1}(2) ~= '#' )
             problems{iprob} = fields{1}(2:end);
          end
       else
          problems{iprob} = fields{1};
       end
   end
   fclose(fid);
   
   for i = start:length(problems)
      pname = problems{i};
      fprintf( ' Decoding problem %4d: %s\n', i , pname );
      if ( ~ismember( pname, {'launchPB.py', 's2xlib.py', 'sifdims.py', 'test_fortran_python.py' } ) )
         thename = s2x( pname );
         thename = s2x( pname, injl );
         thename = s2x( pname, inpy );
         system( [ 'cp /home/philippe/s2x/sif/',pname,'.SIF /home/philippe/docker/mycute_2/newpy/' ] );
         system( [ 'cp ', thename, '.py /home/philippe/docker/mycute_2/newpy' ] );
      end
   end
   
end

