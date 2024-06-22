%
%   Compares the results obtained by running test_matlab, test_python, test_julia and test_fortran
%
%   Poorgramming: Ph Toint, this version 21 VI 2024
%
%   Example: compare_res( 'python', 'matlab', 'julia' )
%

function compare_res( v0, varargin )

%start = 1;   % not implemented
stop  = Inf;

unfixed = 1;     % Compare only components that are not fixed

%   Define the base version, open the corresponding file and parse its first line.

basefile  = [ 'test_', v0, '.res' ];
vs{1}     = v0;
fidv0     = fopen( basefile, 'r' );
line      = fgetl( fidv0 );
fields    = split( line );

%   Do the same for the other versions.

if ( nargin == 1 )
   disp( 'Cannot compare a single version. Aborting.' )
   return
else
   nv = nargin - 1;
   if ( nargin > 1 )
      v1        = varargin{1};
      v1file    = [ 'test_',v1, '.res' ];
      fidv1     = fopen( v1file, 'r' );
      linev1    = fgetl( fidv1 );
      fieldsv1  = split( linev1 );
      vs{end+1} = v1;
      if ( nargin > 2 )
         v2        = varargin{2};
         v2file    = [ 'test_',v2, '.res' ];
         fidv2     = fopen( v2file, 'r' );
         linev2    = fgetl( fidv2 );
         fieldsv2  = split( linev2 );
         vs{end+1} = v2;
         if ( nargin > 3 )
            v3        = varargin{3};
            v3file    = [ 'test_',v3, '.res' ];
            fidv3     = fopen( v3file, 'r' );
            linev3    = fgetl( fidv3 );
            fieldsv3  = split( linev3 );
            vs{end+1} = v3;
         end
      end
   end
end

if ( ~unfixed && ismember( 'fortran', vs ) )
   disp( 'Cannot compare values for fixed variables in Fortran. Aborting.' )
   return
end

if ( unfixed)
   fidu    = fopen('fixed.lst','r');
   lineu   = fgetl( fidu );
   fieldsu = split( lineu );
end

%  Print header.

v0r      = sprintf( '%7s', v0 );
v1l      = sprintf( '%-7s', v1 );
subtitle  = '     x0          L0            g0          Hv      ';
subt2     = ' abs   rel    abs   rel    abs   rel    abs   rel  ';
fmt      = '%4.0e(%4.0e) %4.0e(%4.0e) %4.0e(%4.0e) %4.0e(%4.0e)';

blanks   = '                ';
subfmt   = '----------------- %s-%s -----------------';
switch( nv )
case 1
   fprintf(1,[ '%s ', subfmt, '\n'], blanks, v0r, v1l )  
   fprintf(1,'%s %s\n', blanks, subtitle );
   fprintf(1,'%s %s\n', blanks, subt2);
case 2
   v2l = sprintf( '%-7s', v2 );
   fprintf(1,[ blanks, ' ', subfmt, ' | ', subfmt, '\n'], v0r, v1l, v0r, v2l );
   fprintf(1,'%s %s | %s \n', blanks, subtitle, subtitle );
   fprintf(1,'%s %s | %s \n', blanks, subt2, subt2 );
case 3
   v2l = sprintf( '%-7s', v2 );
   v3l = sprintf( '%-7s', v3 );
   fprintf(1,[ blanks, ' ', subfmt, ' | ', subfmt, ' | ', subfmt, '\n'], v0r, v1l, v0r, v2l, v0r, v3l );
   fprintf(1,'%s %s | %s | %s \n', blanks, subtitle, subtitle, subtitle );
   fprintf(1,'%s %s | %s | %s \n', blanks, subt2, subt2, subt2 );
end
fprintf(1,'\n' );

%   Loop over the file (synchronously) to compare for each problem in turn.

while( ~feof( fidv0 ) )

   %  Read the file describing the active variables.

   if ( unfixed )
      n       = fieldsu{10};
      n_f     = fieldsu{16};
      lineu   = fgetl( fidu );
      fieldsu = split( strtrim( lineu ) );
      [ fixed, fieldsu ] = getvect( fieldsu, 'fix', fidu );
      active  = setdiff( [1:round(str2double(n))], fixed );
   else
      active  = [];
   end

   %  Find the next problem in the base version (v0) and get its data.

   if ( strcmp( fields{1},'===' ) )
      pbind = round( str2double(fields{4}));
      if ( pbind > stop )
         return
      end
      pname   = fields{6};
      n       = str2double( fields{10} );
      m       = str2double( fields{13} );
      line    = fgetl( fidv0 );
      fields  = split( strtrim( line ) );
      [ x0_v0, fields ] = getuvect( fields, 'x0', fidv0, active, strcmp( v0, 'fortran' ), unfixed );
     L0_v0 = str2double( fields{3} );
      line    = fgetl( fidv0 );
      fields  = split( strtrim( line ) );
      [ g0_v0, fields ] = getuvect( fields, 'g0', fidv0, active, strcmp( v0, 'fortran' ), unfixed );
      [ Hv_v0, fields ] = getuvect( fields, 'Hv', fidv0, active, strcmp( v0, 'fortran' ), unfixed );
   end

%  Find the same problem in all other versions and get their data.

   namev1      = fieldsv1{6};
   if ( strcmp( namev1, pname ) )
      n_v1     = str2double( fieldsv1{10} );
      m_v1     = str2double( fieldsv1{13} );
      line     = fgetl( fidv1 );
      fieldsv1 = split( strtrim( line ) );
      [ x0_v1, fieldsv1 ] = getuvect( fieldsv1, 'x0', fidv1, active, strcmp( v1, 'fortran' ), unfixed );
      L0_v1    = str2double( fieldsv1{3} );
      line     = fgetl( fidv1 );
      fieldsv1 = split( strtrim( line ) );
      [ g0_v1, fieldsv1 ] = getuvect( fieldsv1, 'g0', fidv1, active, strcmp( v1, 'fortran' ), unfixed );
      [ Hv_v1, fieldsv1 ] = getuvect( fieldsv1, 'Hv', fidv1, active, strcmp( v1, 'fortran' ), unfixed );
   else
      disp( [ 'Out of sync: can''t find problem ', pname, ' in ', v1file, ' !' ] )
      return
   end

   if ( nv > 1 )
      namev2   = fieldsv2{6};
      if ( strcmp( namev2, pname ) )
         n_v2     = str2double( fieldsv2{10} );
         m_v2     = str2double( fieldsv2{13} );
         line     = fgetl( fidv2 );
         fieldsv2 = split( strtrim( line ) );
         [ x0_v2, fieldsv2 ] = getuvect( fieldsv2, 'x0', fidv2, active, strcmp( v2, 'fortran' ), unfixed );
         L0_v2    = str2double( fieldsv2{3} );
         line     = fgetl( fidv2 );
         fieldsv2 = split( strtrim( line ) );
         [ g0_v2, fieldsv2 ] = getuvect( fieldsv2, 'g0', fidv2, active, strcmp( v2, 'fortran' ), unfixed );
         [ Hv_v2, fieldsv2 ] = getuvect( fieldsv2, 'Hv', fidv2, active, strcmp( v2, 'fortran' ), unfixed );
      else
         disp( [ 'Out of sync: can''t find problem ', pname, ' in ', v2file,' !' ] )
         return
      end
      if ( nv > 2 )
         namev3      = fieldsv3{6};
         if ( strcmp( namev3, pname ) )
            n_v3     = str2double( fieldsv3{10} );
            m_v3     = str2double( fieldsv3{13} );
            line     = fgetl( fidv3 );
            fieldsv3 = split( strtrim( line ) );
            [ x0_v3, fieldsv3 ] = getuvect( fieldsv3, 'x0', fidv3, active, strcmp( v3, 'fortran' ), unfixed );
            L0_v3    = str2double( fieldsv3{3} );
            line     = fgetl( fidv3 );
            fieldsv3 = split( strtrim( line ) );
            [ g0_v3, fieldsv3 ] = getuvect( fieldsv3, 'g0', fidv3, active, strcmp( v3, 'fortran' ), unfixed );
            [ Hv_v3, fieldsv3 ] = getuvect( fieldsv3, 'Hv', fidv3, active, strcmp( v3, 'fortran' ), unfixed );
         else
            disp( [ 'Out of sync: can''t find problem ', pname, ' in ', v3file,' !' ] )
            return
         end
      end
   end

   %   Compare the different versions.

   dx0v0v1a = norm(x0_v0-x0_v1);
   dx0v0v1r = dx0v0v1a /( 1.e-15+norm(x0_v0) );
   dL0v0v1a = norm(L0_v0-L0_v1);
   dL0v0v1r = dL0v0v1a /( 1.e-15+norm(L0_v0) );
   dg0v0v1a = norm(g0_v0-g0_v1);
   dg0v0v1r = dg0v0v1a /( 1.e-15+norm(g0_v0) );
   dHvv0v1a = norm(Hv_v0-Hv_v1);
   dHvv0v1r = dHvv0v1a /( 1.e-15+norm(Hv_v0) );

   if ( nv== 1 )
      fprintf( [ '%4d %11s ', fmt, '\n' ], ...
                  pbind, pname, dx0v0v1a, dx0v0v1r, dL0v0v1a, dL0v0v1r, dg0v0v1a, dg0v0v1r, dHvv0v1a, dHvv0v1r );
   else
      if ( nv > 1 )
         dx0v0v2a = norm(x0_v0-x0_v2);
         dx0v0v2r = dx0v0v2a /( 1.e-15+norm(x0_v0) );
         dL0v0v2a = norm(L0_v0-L0_v2);
         dL0v0v2r = dL0v0v2a /( 1.e-15+norm(L0_v0) );
         dg0v0v2a = norm(g0_v0-g0_v2);
         dg0v0v2r = dg0v0v2a /( 1.e-15+norm(g0_v0) );
         dHvv0v2a = norm(Hv_v0-Hv_v2);
         dHvv0v2r = dHvv0v2a /( 1.e-15+norm(Hv_v0) );
         if ( nv == 2 )
            fprintf( [ '%4d %11s ', fmt, ' | ', fmt, '\n' ], ...
                     pbind, pname, dx0v0v1a, dx0v0v1r, dL0v0v1a, dL0v0v1r, dg0v0v1a, dg0v0v1r, dHvv0v1a, dHvv0v1r, ...
                                   dx0v0v2a, dx0v0v2r, dL0v0v2a, dL0v0v2r, dg0v0v2a, dg0v0v2r, dHvv0v2a, dHvv0v2r );
         else
            dx0v0v3a = norm(x0_v0-x0_v3);
            dx0v0v3r = dx0v0v3a /( 1.e-15+norm(x0_v0) );
            dL0v0v3a = norm(L0_v0-L0_v3);
            dL0v0v3r = dL0v0v3a /( 1.e-15+norm(L0_v0) );
            dg0v0v3a = norm(g0_v0-g0_v3);
            dg0v0v3r = dg0v0v3a /( 1.e-15+norm(g0_v0) );
            dHvv0v3a = norm(Hv_v0-Hv_v3);
            dHvv0v3r = dHvv0v3a /( 1.e-15+norm(Hv_v0) );
            fprintf( [ '%4d %11s ', fmt, ' | ', fmt, ' | ', fmt, '\n' ], ...
                     pbind, pname, dx0v0v1a, dx0v0v1r, dL0v0v1a, dL0v0v1r, dg0v0v1a, dg0v0v1r, dHvv0v1a, dHvv0v1r, ...
                                   dx0v0v2a, dx0v0v2r, dL0v0v2a, dL0v0v2r, dg0v0v2a, dg0v0v2r, dHvv0v2a, dHvv0v2r, ...
                                   dx0v0v3a, dx0v0v3r, dL0v0v3a, dL0v0v3r, dg0v0v3a, dg0v0v3r, dHvv0v3a, dHvv0v3r );
         end
      end
   end

   
end

fclose( fidv0)
switch ( nv )
case 1
   fclose( fidv1);
case 2
   fclose( fidv1);
case 2
   fclose( fidv1);
   fclose( fidv2);
case 2
   fclose( fidv1);
   fclose( fidv2);
   fclose( fidv3);
end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ vect, fields ] = getvect( fields, label, fid )

k = 0;  % the number of components read
vect = [];
for nline = 1:1000
   nf = length( fields ) - 2;
   for j = 1:nf
      if ( contains( fields{ 2+j }, 'NaN' ) )
          vect( k+j ) = Inf;
      elseif ( strcmp( label, 'fix' ) )
          vect( k+j ) = round(str2double( fields{ 2+j } ));
      else
          vect( k+j ) = str2double( fields{ 2+j } );
      end
   end
   k = k + nf;
   if ( ~feof( fid ) )
      line   = fgetl( fid );
      fields = split( strtrim( line ) );
      if ( ~strcmp( fields{1}, label ) )
         break
      end
   else
      fields = {};
      break
   end
end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ vect, fields ] = getuvect( fields, label, fid, active, is_fortran, unfixed )

[ vect, fields ] = getvect( fields, label, fid );
if ( unfixed && length( vect ) > length( active ) )
   vect = vect( active );
end

return

end


