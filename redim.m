%
%   Reduce dimension by setting problems arguments to coincide with these
%   specified in fullproblist
%
%   Programming: Ph. Toint
%   This version 20 VI 2024
%

function redim( varargin )

sifdir = '/home/philippe/problems/mysif/';

if ( nargin )

   myprobs = {};
   for i= 1:nargin
      myprobs{end+1} = varargin{i};
   end
   
end

start = 1;

theproblems = 'fullproblist';

fid = fopen( theproblems, 'r' );
iprob = 0;
while( ~feof(fid ) )
   line   = fgetl( fid );
   fields = split( line );
   if ( fields{1}(1) == '#' )
      continue
   end
   pname = fields{1};
   if ( nargin == 0 || ismember( pname, myprobs ) )
      iprob = iprob + 1;
      fprintf( ' Redim problem %4d: %s\n', iprob , pname );
      args = split( fields{2}, ',' );
      redimprob( sifdir, pname, args )
   end
end
fclose(fid);
   

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function redimprob( sifdir, pname, args )

okname  = pname;
sifname = replace( pname, 'n','' );
sifname = replace( sifname,'p','+' );
sifname = replace( sifname,'m','-' );
sifname = replace( sifname,'t','*' );
sifname = replace( sifname,'d','/' );

fidin    = fopen( [ sifdir, sifname, '.SIF' ], 'r' );
fidout   = fopen( ['./sif/', okname, '.SIF' ], 'w' );
nargs    = length( args );
iarg     = 0;
assigned = {};
while( ~feof( fidin ) )
   line   = fgetl( fidin );
   lline  = length( line );
   if ( lline &&  line(1) ~= '*' && contains( line, '$-PARAMETER' ) )
      parname = line(5:14);
      if ( ismember( parname, assigned ) )
         line(1) = '*';
      elseif ( iarg < nargs )
         iarg  = iarg + 1;
         argi  = args{ iarg };
         largi = length( argi );
         argi(largi+1:max(largi,12) ) = " ";
         cline    = line;
         cline(1) = '*';
         fprintf( fidout, '%s\n', cline );
         assigned{end+1} = parname;
         line(25:36)     = argi;
         line = [ line(1:strfind( line, '$-PARAMETER' )+10), '     modified for S2MPJ tests' ];
      end
   elseif ( ( lline >4 && strcmp(line(1:4),'NAME') ) ||     ...
            ( lline >8 && strcmp(line(1:8),'ELEMENTS') ) || ...
            ( lline >6 && strcmp(line(1:6),'GROUPS') )      )
      line = replace( line, sifname, okname );
   elseif ( lline >= 15 && strcmp( line(5:length(sifname)), sifname ) )
      line(5:length(okname)) = okname;
   end
   fprintf( fidout, '%s\n', line );
end
fclose( fidin );
fclose( fidout );

return

end
