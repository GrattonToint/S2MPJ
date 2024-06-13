function regenerate( varargin )
%
%   Regenerate .m, .jl and .py files from the SIF files.
%
%   If varargin contains a list of problem names, only those are decoded.
%   If regenerate is used without arguments, the list of problems to regenerate
%   is given by the fullproblist file, which also specifies the problems
%   arguments The values of start and stop may be modified to select a sublist.
%
%   Programming: Ph. Toint
%   This version 7 V 2024
%

addpath ./sif  % the directory containing the problem's SIF files

inma.language = 'matlab';
inpy.language = 'python';
injl.language = 'julia';

if ( nargin )  % Decode only a few problems

   %  Options for decoding a small set of problems
   
   inma.outdir   = './matlab_problems';
   inpy.outdir   = './python_problems';
   injl.outdir   = './julia_problems';

   for i= 1:nargin
      pname = varargin{i};
      fprintf( ' Decoding problem %s\n', pname );
      thename = s2mpj( pname, inma );
      thename = s2mpj( pname, inpy );
      thename = s2mpj( pname, injl );
   end

else       % Decode the full problem list

   theproblems = './fullproblist';
   start = 1;
   stop  = 10000;

   %  Provide a backup, just in case...
   
   if ( start == 1 )
      if( ~exist( 'matlab_problems_bak' ) )
         system( 'mkdir matlab_problems_bak' );
      end
      system( 'cp -r matlab_problems matlab_problems_bak' );
      if( ~exist( 'python_problems_bak' ) )
         system( 'mkdir python_problems_bak' );
      end
      system( 'cp -r python_problems python_problems_bak' );
      if( ~exist( 'julia_problems_bak' ) )
         system( 'mkdir julia_problems_bak' );
      end
      system( 'cp -r julia_problems julia_problems_bak' );
   end

   %  Create a temporary working directory.
   
   if ( ~exist( './temporaries' ) )
      system( 'mkdir ./temporaries' );
   end
   inma.outdir   = './temporaries';
   inpy.outdir   = './temporaries';
   injl.outdir   = './temporaries';

   % Loop on the problems, skipping entries starting with #

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

   for i = start:min(stop,length(problems))
      pname = problems{i};
      fprintf( ' Decoding problem %4d: %s\n', i , pname );

      %  Decode to a Matlab file

      thename = s2mpj( pname, inma );
      if ( exist( [ './matlab_problems/', thename, '.m' ] ) )
         has_changed = system([ 'diff -q ./matlab_problems/', thename, '.m  ./temporaries/', pname, '.m > /dev/null' ] );
      else
         has_changed = 1;
      end
      if ( has_changed )
         system( [ 'mv ./temporaries/', thename, '.m ./matlab_problems/', thename, '.m' ] );
      end

      %  Decode to a Python file
      
      thename = s2mpj( pname, inpy );
      if ( exist( [ './python_problems/', thename, '.py' ] ) )
         has_changed = system([ 'diff -q ./python_problems/', thename, '.py  ./temporaries/', pname, '.py > /dev/null' ] );
      else
         has_changed = 1;
      end
      if ( has_changed )
         system( [ 'mv ./temporaries/', thename, '.py ./python_problems/', thename, '.py' ] );
      end

      %  Decode to a Julia file

      thename = s2mpj( pname, injl );
      if ( exist( [ './julia_problems/', thename, '.jl' ] ) )
         has_changed = system([ 'diff -q ./julia_problems/', thename, '.jl  ./temporaries/', pname, '.jl > /dev/null' ] );
      else
         has_changed = 1;
      end
      if ( has_changed )
         system( [ 'mv ./temporaries/', thename, '.jl ./julia_problems/', thename, '.jl' ] );
      end
   end

   %  Clean up.

   system( 'rm -r ./temporaries' );

   %  Update the lists of available problems.
   
   system( 'ls matlab_problems > list_of_matlab_problems' );
   system( 'ls python_problems > list_of_python_problems' );
   system( 'ls julia_problems  > list_of_julia_problems'  );

end
disp( 'Done.' )
   

