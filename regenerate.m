function regenerate( varargin )
%
%   Regenerate .m, .jl and .py files from the SIF files.
%
%   If varargin contains one or more of 'matlab', 'python' or 'julia', problems
%   are decoded to the corresponding language(s).  Otherwise problems are
%   decoded to the three languages.
%
%   If varargin contains the name(s) of one or more problems (e.g. 'ROSENBR'), the
%   corresponding problems are decoded. Otherwise, all problems specified in the
%   ./fullproblist file are decoded.
%
%   If varargin contains a string of the form './problist', the list of problems yo
%   decode is read from the ./problist file (the ./ are mandatory) and this list
%   is used only  if no problem is explictly specified. (Default: ./fullproblist)
%   Only the last occurence of a ./problist argument is taken into account.
%   The ./problist file must contain a problem name as the first (blank separated)
%   field on each line.
%
%   Examples:
%     regenerate()           decodes all problems for the three languages
%     regenerate('python')   decodes the python version of all problems
%     regenerate('ROSENBR')  decodes the ROSENBR.SIF file for the three languages
%     regenerate('ROSENBR','JUNKTURN','matlab','julia') decodes the ROSENBR.SIF
%                            and JUNKTURN.SIF files for matlab and julia
%     regenerate('./problist','matlab','python') decodes the problems specified in the
%                            ./problist file for matlab and python
%
%   Programming: Ph. Toint
%   This version 25 VI 2024
%

addpath ./sif   % the directory containing the problem's SIF files

start  = 1;     % the index of the first problem (in the specified list of problems)
                % to be decoded
stop   = 10000; % the index of the lasst problem (in the specified list of problems)
                % to be decoded
backup = 1;     % true if a backup problems directories should be created
                % whenever start = 1

%  Process the list of arguments, if any.

problist = './fullproblist';
versions = {};
problems = {};

for i = 1:nargin
   if ( strcmp( varargin{i}(1:2), './' ) )
      problist = varargin{i};
      continue
   end
   switch( varargin{i} )
   case 'all'
      versions = { 'matlab', 'python', 'julia' };
   case 'matlab'
      versions = union( versions, 'matlab' );
   case 'python'
      versions = union( versions, 'python' );
   case 'julia'
      versions = union( versions, 'julia'  );
   otherwise
      problems = union( problems, varargin{i} );
   end
end
if ( isempty( versions ) )
   versions = { 'matlab', 'python', 'julia' };
end

%  Provide a backup, just in case...

if ( backup && start == 1 )
   if( ismember( 'matlab', versions ) && ~exist( 'matlab_problems_bak' ) )
      system( 'mkdir matlab_problems_bak' );
   end
   system( 'cp -r matlab_problems matlab_problems_bak' );
   if( ismember( 'python', versions ) && ~exist( 'python_problems_bak' ) )
      system( 'mkdir python_problems_bak' );
   end
   system( 'cp -r python_problems python_problems_bak' );
   if( ismember( 'julia', versions ) && ~exist( 'julia_problems_bak' ) )
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

%  Read the list of problems from ./problist, if none was explicitly specified.

if ( isempty( problems ) )
   fid = fopen(problist, 'r' );
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
end

%  Create s2mpj options to specify languages.

inma.language = 'matlab';
inpy.language = 'python';
injl.language = 'julia';

% Loop on the problems

for i = start:min(stop,length(problems))
   pname = problems{i};
   fprintf( ' Decoding problem %4d: %s\n', i , pname );

   %  Decode to a Matlab file

   if ( ismember( 'matlab', versions ) )
      thename = s2mpj( pname, inma );
      if ( exist( [ './matlab_problems/', thename, '.m' ] ) )
         has_changed = system([ 'diff -q ./matlab_problems/', thename, '.m  ./temporaries/', pname, '.m > /dev/null' ] );
      else
         has_changed = 1;
      end
      if ( has_changed )
         system( [ 'mv ./temporaries/', thename, '.m ./matlab_problems/', thename, '.m' ] );
     end
   end
   
   %  Decode to a Python file
   
   if ( ismember( 'python', versions ) )
      thename = s2mpj( pname, inpy );
      if ( exist( [ './python_problems/', thename, '.py' ] ) )
         has_changed = system([ 'diff -q ./python_problems/', thename, '.py  ./temporaries/', pname, '.py > /dev/null' ] );
      else
         has_changed = 1;
      end
      if ( has_changed )
         system( [ 'mv ./temporaries/', thename, '.py ./python_problems/', thename, '.py' ] );
      end
   end
   
   %  Decode to a Julia file

   if ( ismember( 'julia', versions ) )
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

end

%  Clean up.

system( 'rm -r ./temporaries' );

%  Update the lists of available problems.

system( 'ls matlab_problems > list_of_matlab_problems' );
system( 'ls python_problems > list_of_python_problems' );
system( 'ls julia_problems  > list_of_julia_problems'  );
disp( 'Done.' )
   
return

end

