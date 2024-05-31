%
%   Regenerate .m, .jl and .py files and copy the latter to the docker.
%
%   Programming: Ph. Toint
%   This version 28 V 2024
%

function regenerate( varargin )

addpath ./sif
inma.language = 'matlab';
inma.outdir   = './matlab_problems';
inpy.language = 'python';
inpy.outdir   = './python_problems';
injl.language = 'julia';
injl.outdir   = './julia_problems';

if ( nargin )

   for i= 1:nargin
      pname = varargin{i};
      fprintf( ' Decoding problem %s\n', pname );
      thename = s2x( pname, inma );
      thename = s2x( pname, inpy );
      thename = s2x( pname, injl );
      system( [ 'cp ./sif/',pname,'.SIF /home/philippe/docker/mycute_2/newpy' ] );
      system( [ 'cp ./python_problems/', thename, '.py /home/philippe/docker/mycute_2/newpy' ] );
   end


else

   start = 1;

   system( [ 'cp fullproblist /home/philippe/docker/mycute_2/newpy' ] );
   system( [ 'cp s2xlib.py /home/philippe/docker/mycute_2/newpy' ] );
   system( [ 'cp test_fortran_python.py /home/philippe/docker/mycute_2/newpy' ] );

   theproblems = 'fullproblist';
   
   fid = fopen(theproblems, 'r' );
   problems = {};
   iprob = 0;
   while( ~feof(fid ) )
       line   = fgetl( fid );
       fields = split(line );
       if ( fields{1}(1) ~= '#' )
          iprob  = iprob + 1;
          problems{iprob} = fields{1};
       end
   end
   fclose(fid);
   
   for i = start:length(problems)
      pname = problems{i};
      fprintf( ' Decoding problem %4d: %s\n', i , pname );
      if ( ~ismember( pname, { 's2xlib.py', 'test_fortran_python.py' } ) )
         thename = s2x( pname, inma );
         thename = s2x( pname, injl );
         thename = s2x( pname, inpy );
         system( [ 'cp ./sif/',pname,'.SIF /home/philippe/docker/mycute_2/newpy/' ] );
         system( [ 'cp ./python_problems/', thename, '.py /home/philippe/docker/mycute_2/newpy' ] );
      end
   end
   
end

