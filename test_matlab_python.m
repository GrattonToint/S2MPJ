function test_matlab_python(filename, istart)

diary('results_matlab_python');

fid = fopen(filename, 'r' );
i=1;
while( ~feof(fid) )
   line  = fgetl(fid);
   if ( line(1) == '#' )
      continue
   end
   pblist{i} = textscan( line, '%s %s' );
   i = i+1;
end
fclose(fid);

addpath('/Users/sgratton/S2X/S2X/matlab_problems');

% Chemin du répertoire contenant le module Python
modulePath = '/Users/sgratton/S2X/S2X/python_problems';

% Ajouter le répertoire au sys.path
py.sys.path().append(string(modulePath));



% Parcours du dictionnaire

for i = istart:length(pblist)
  name   = pblist{i}{1}{1}; % i-ème clé
  valeur = pblist{i}{2}{1}; % i-ème valeur
  %
  disp(name)
  disp('********   ');
        command = ['wc -l < ', name, '.py 2>/dev/null'] ; ;
        [status, cmdout] = system(command) ;
        numLines = str2double(cmdout);
        if numLines > 40000
            error('File  "%s.py" has more than 40 000 lines.', name);
        end

        % Handle Python files
        
        options.language = 'python';

        eval(['ModLoaded = py.importlib.import_module(''', name, ''');']) % Load the module
        tic
        eval(['Prob_py = ModLoaded.',name, '(',valeur,');'])              % initialize problem
        setuptime = toc;
        tic
        y = ones(double(Prob_py.pb.m),1);
        C     = Prob_py.Lgxy(Prob_py.pb.x0,py.numpy.array(y));                               % compute 
        fxpy = double(C{1}); gxpy=double(C{2});                           % unpack
        evaltime1 = toc;
        tic
        C     = Prob_py.Lgxy(Prob_py.pb.x0,py.numpy.array(y));                               % compute 
        fxpy = double(C{1}); gxpy=double(C{2});                           % unpack
        evaltime2 = toc;
        fprintf('Python:  %.4f %.4f %.4f  s\n', setuptime, evaltime1, evaltime2 );

        % Handle Matlab files

        tic
        eval(['Prob_ma =',name,'(''setup'',', valeur ');'])        % Initialize
        setuptime = toc;
        y = ones(Prob_ma.m,1);
        tic
        eval([' [ fxma, gxma ] =',name,'(''Lgxy'',Prob_ma.x0,y);'])       % compute'
        evaltime1 = toc;
        tic
        eval([' [ fxma, gxma ] =',name,'(''Lgxy'',Prob_ma.x0,y);'])       % compute'
        evaltime2 = toc;
        fprintf('Matlab:  %.4f %.4f %.4f  s\n', setuptime, evaltime1, evaltime2 );
        
        disp(['   function difference ', num2str(abs(fxma-fxpy)/(1e-12+abs(fxma)))])
        disp(['   gradient difference ', num2str(norm(gxma-gxpy)/(1e-12+norm(gxma)))])
        disp('')
end 

diary off;
