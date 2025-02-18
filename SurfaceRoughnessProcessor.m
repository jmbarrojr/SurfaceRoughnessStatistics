classdef SurfaceRoughnessProcessor
    methods (Static)
        function Surface = getSurfStatistics(filename)
            Surface = SurfaceRoughnessProcessor.loadSurface(filename);
            Surface = SurfaceRoughnessProcessor.determineSurfaceType(Surface);
            Surface = SurfaceRoughnessProcessor.determineXandYdir(Surface);
            Surface = SurfaceRoughnessProcessor.roughPhysicalProp(Surface);
            Surface = SurfaceRoughnessProcessor.roughnessStats(Surface);
        end
        
        function SurfStruct = loadSurface(filename)
            ext = filename(end-3:end);
            switch ext
                case '.mat'
                    matObj = matfile(filename, 'Writable', false);
                    SurfStruct.obj = matObj;
                    SurfStruct.varProps = whos(matObj);
                case {'.asc','.ASC','.dat','.DAT','.csv'}
                    SurfStruct = SurfaceRoughnessProcessor.importSurfaceTxtData(filename);
                case {'.xls','xlsx'}
                    SurfStruct = SurfaceRoughnessProcessor.importSurfaceExcel(filename);
                otherwise
                    error('File format might not be yet implemented')
            end
            N = length(SurfStruct.varProps);
            c = 1;
            for n = 1:N
                cls = SurfStruct.varProps(n).class;
                name = SurfStruct.varProps(n).name;
                if strcmp(cls, 'double')
                    varNames{c} = name; %#ok<AGROW>
                    c = c + 1;
                end
            end
            SurfStruct.varNames = varNames;
        end
        
        function SurfStruct = importSurfaceTxtData(filename)
            A = importdata(filename);
            L = size(A, 2);
            if L == 2
                X = A(:, 1);
                Z = A(:, 2);
                save('temp.mat', 'X', 'Z')
            elseif L == 3
                A = sortrows(A, [1, 2]);
                J = find(A(1, 1) - A(:, 1) == 0, 1, 'last');
                A = sortrows(A, [2, 1]);
                I = find(A(1, 2) - A(:, 2) == 0, 1, 'last');
                X = reshape(A(:, 1), I, J)';
                Y = reshape(A(:, 2), I, J)';
                Z = reshape(A(:, 3), I, J)';
                save('temp.mat', 'X', 'Y', 'Z')
            end
            SurfStruct = SurfaceRoughnessProcessor.loadSurface('temp.mat');
        end
        
        function SurfStruct = importSurfaceExcel(filename)
            A = readmatrix(filename);
            L = size(A, 2);
            if L == 2
                X = A(:, 1);
                Z = A(:, 2);
                save('temp.mat', 'X', 'Z')
            elseif L == 3
                A = sortrows(A, [1, 2]);
                J = find(A(1, 1) - A(:, 1) == 0, 1, 'last');
                A = sortrows(A, [2, 1]);
                I = find(A(1, 2) - A(:, 2) == 0, 1, 'last');
                X = reshape(A(:, 1), I, J)';
                Y = reshape(A(:, 2), I, J)';
                Z = reshape(A(:, 3), I, J)';
                save('temp.mat', 'X', 'Y', 'Z')
            end
            SurfStruct = SurfaceRoughnessProcessor.loadSurface('temp.mat');
        end
        
        function SurfStruct = determineSurfaceType(SurfStruct)
            L = length(SurfStruct.varNames);
            if L == 2
                type = '1D-profile';
            elseif L == 3
                type = '2D-surface';
            else
                error('Could not determine the surface type')
            end
            SurfStruct.type = type;
        end
        
        function SurfStruct = determineXandYdir(SurfStruct)
            zname = SurfStruct.varNames{end};
            z = SurfStruct.obj.(zname);
            [sr, ~] = size(z);
            switch SurfStruct.type
                case '1D-profile'
                    if sr == 1
                        Xdir = 2;
                    else
                        Xdir = 1;
                    end
                    Ydir = [];
                case '2D-surface'
                    xname = SurfStruct.varNames{1};
                    x = SurfStruct.obj.(xname);
                    if abs(x(1, 2) - x(1, 1)) == 0
                        Xdir = 1;
                        Ydir = 2;
                    else
                        Xdir = 2;
                        Ydir = 1;
                    end
            end
            SurfStruct.Xdir = Xdir;
            SurfStruct.Ydir = Ydir;
        end
        
        function SurfStruct = roughPhysicalProp(SurfStruct)
            xname = SurfStruct.varNames{1};
            x = SurfStruct.obj.(xname);
            switch SurfStruct.type
                case '1D-profile'
                    Lx = max(x) - min(x);
                    Ly = [];
                case '2D-surface'
                    Lx = max(x(:)) - min(x(:));
                    yname = SurfStruct.varNames{2};
                    y = SurfStruct.obj.(yname);
                    Ly = max(y(:)) - min(y(:));
            end
            SurfStruct.Lx = Lx;
            SurfStruct.Ly = Ly;
        end
        
        function SurfStruct = roughnessStats(SurfStruct)
            zname = SurfStruct.varNames{end};
            z = SurfStruct.obj.(zname);
            kmin = min(z(:));
            SurfStruct.kmin = kmin;
            kmax = max(z(:));
            SurfStruct.kmax = kmax;
            SurfStruct.kt = kmax - kmin;
            MaxZ = sortrows(z(:), 'descend');
            MinZ = sortrows(z(:), 'ascend');
            SurfStruct.kz5 = mean(MaxZ(1:5)) - mean(MinZ(1:5));
            SurfStruct.kbar = mean(z(:));
            SurfStruct.ka = mean(abs(z(:)));
            SurfStruct.krms_z = rms(z(:));
            z = z - mean(z(:));
            SurfStruct.krms = rms(z(:));
            SurfStruct.Sk = skewness(z(:));
            SurfStruct.Fl = kurtosis(z(:));
            xname = SurfStruct.varNames{1};
            x = SurfStruct.obj.(xname);
            SurfStruct.ESx = SurfaceRoughnessProcessor.EffectiveSlope(x, z, SurfStruct.Xdir, SurfStruct.Lx);
            SurfStruct.Rlx = SurfaceRoughnessProcessor.CorrelationLenght(x, z, SurfStruct.Xdir);
            switch SurfStruct.type
                case '2D-surface'
                    yname = SurfStruct.varNames{2};
                    y = SurfStruct.obj.(yname);
                    SurfStruct.ESy = SurfaceRoughnessProcessor.EffectiveSlope(y, z, SurfStruct.Ydir, SurfStruct.Ly);
                    SurfStruct.Rly = SurfaceRoughnessProcessor.CorrelationLenght(y, z, SurfStruct.Ydir);
                otherwise
                    SurfStruct.ESy = [];
                    SurfStruct.Rly = [];
            end
        end
        
        function Es = EffectiveSlope(X, Z, dir, L)
            if dir == 1
                dx = X(2, 1) - X(1, 1);
                mean_dir = 2;
            else
                dx = X(1, 2) - X(1, 1);
                mean_dir = 1;
            end
            dzdx = diff(Z, 1, dir) ./ diff(X, 1, dir);
            Es = 1/L .* mean(trapz(abs(dzdx), dir) .* dx, mean_dir);
        end
        
        function Rlx = CorrelationLenght(X, Z, dir)
            if dir == 1
                dx = X(2, 1) - X(1, 1);
            else
                dx = X(1, 2) - X(1, 1);
            end
            s = size(Z, dir);
            [lags, Zcorr] = SurfaceRoughnessProcessor.MeanAutoCorr_FFT(Z, dir);
            lags = lags .* dx;
            ind = SurfaceRoughnessProcessor.findSlopeCorr(Zcorr);
            Rlx = interp1(Zcorr(s:s+ind-1), lags(s:s+ind-1), 1./exp(1), 'linear');
        end
        
        function [lags, C] = MeanAutoCorr_FFT(A, dir)
            S = size(A);
            s = S(dir);
            nfft = 2*s-1;
            A = A - mean(A(:));
            Afft = fft(A, nfft, dir);
            Afft = Afft ./ s;
            C = Afft .* conj(Afft);
            C = fftshift(ifft(C, [], dir), dir);
            if S(1) > 1 || S(2) > 1
                if dir == 1
                    C = mean(C, dir+1);
                elseif dir == 2
                    C = mean(C, dir-1);
                end
            end
            C = C ./ max(C(:));
            lags = [-linspace(s-1, 1, s-1) 0 linspace(1, s-1, s-1)];
        end
        
        function ind = findSlopeCorr(C)
            l = length(C);
            C = C((l+1)/2:end);
            ind = find(C-min(C) == 0, 1, 'first');
        end
        
        function S = cleanUpStruct(SurfStruct)
            S = struct();
            fields = fieldnames(SurfStruct);
            for n = 4:length(fields)
                f = fields{n};
                S.(f) = SurfStruct.(f);
            end
        end
        
        function SurfAnswers = RoughnessQuestionnaire(batch)
            switch batch
                case true
                    SurfAnswers = SurfaceRoughnessProcessor.loadQuestionnaire('Questionnaire_Batch.txt');
                otherwise
                    prompt = 'Is the surface Homogeneous or Heterogeneous? ';
                    q1 = 'Hom'; q2 = 'Het';
                    SurfAnswers.Q1 = questdlg(prompt, 'Roughness Information', q1, q2, q1);
                    SurfaceRoughnessProcessor.checkAnswer(SurfAnswers.Q1);
                    
                    prompt = 'Is the roughness Regular or Irregular?';
                    q1 = 'Reg'; q2 = 'Irreg';
                    SurfAnswers.Q2 = questdlg(prompt, 'Roughness Information', q1, q2, q1);
                    SurfaceRoughnessProcessor.checkAnswer(SurfAnswers.Q2);
                    
                    prompt = 'Are results for a ...';
                    q1 = 'TBL'; q2 = 'Pipe'; q3 = 'Channel';
                    SurfAnswers.Q3 = questdlg(prompt, 'Roughness Information', q1, q2, q3, q1);
                    SurfaceRoughnessProcessor.checkAnswer(SurfAnswers.Q3);
                    
                    prompt = 'Are results from Experiments or Simulations?';
                    q1 = 'Exp'; q2 = 'Sim';
                    SurfAnswers.Q4 = questdlg(prompt, 'Roughness Information', q1, q2, q1);
                    SurfaceRoughnessProcessor.checkAnswer(SurfAnswers.Q4);
                    
                    prompt = 'What is the general descriptor of this surface, .i.e. "Sandgrain"?';
                    SurfAnswers.Q5 = inputdlg(prompt, 'Roughness Information', [1 50]);
                    SurfAnswers.Q5 = SurfAnswers.Q5{1};
                    SurfaceRoughnessProcessor.checkAnswer(SurfAnswers.Q5);
                    
                    prompt1 = 'What is the last name of the lead author of the study? ';
                    prompt2 = 'What year were the results published? ';
                    prompt3 = 'What is the identifying name of this surface, i.e. "220Grit"? ';
                    prompt4 = 'Include the doi of the publication';
                    temp = inputdlg({prompt1, prompt2, prompt3, prompt4}, 'Roughness Information', [1 50; 1 50; 1 50; 1 50]);
                    SurfAnswers.Q6 = temp{1};
                    SurfaceRoughnessProcessor.checkAnswer(SurfAnswers.Q6)
                    SurfAnswers.Q7 = temp{2};
                    SurfaceRoughnessProcessor.checkAnswer(SurfAnswers.Q7)
                    SurfAnswers.Q8 = temp{3};
                    SurfaceRoughnessProcessor.checkAnswer(SurfAnswers.Q8)
                    SurfAnswers.Q9 = temp{4};
                    SurfaceRoughnessProcessor.checkAnswer(SurfAnswers.Q9)
            end
        end
        
        function ScannerAnswers = ScannerQuestionnaire(SurfAnswers, batch)
            if strcmp(SurfAnswers.Q4, 'Exp')
                switch batch
                    case true
                        ScannerAnswers = SurfaceRoughnessProcessor.loadQuestionnaire('Profiler_Batch.txt');
                    otherwise
                        prompt1 = 'What is the name and model of the profiler/scanner? ';
                        prompt2 = 'What is the uncertainty in the measurement of surface heights in microns? ';
                        prompt3 = ['What is the unit of measurements of the scan(file),'...
                                  'e.g. mm (milimeters), um (microns), in(inches)? '];
                        temp = inputdlg({prompt1, prompt2, prompt3}, 'Scanner Information', [1 50; 1 50; 1 50]);
                        ScannerAnswers.Q1 = temp{1};
                        ScannerAnswers.Q2 = temp{2};
                        ScannerAnswers.Q3 = temp{3};
                end
            else
                ScannerAnswers.Q1 = 'N/A';
                ScannerAnswers.Q2 = 'N/A';
                ScannerAnswers.Q3 = 'N/A';
            end
        end
        
        function SurfAnswers = loadQuestionnaire(file)
            fid = fopen(file, 'r');
            n = 1;
            while ~feof(fid)
                line = strtrim(fgetl(fid));
                if isempty(line) || all(isspace(line)) || strncmp(line, '#', 1)
                else
                    Q = ['Q' num2str(n)];
                    SurfAnswers.(Q) = line;
                    n = n + 1;
                end
            end
            fclose(fid);
        end
        
        function checkAnswer(S)
            if iscell(S) && ~isempty(S)
                for n = 1:length(S)
                    if isempty(S{n})
                        error('Please, select/type a proper answer and/or don''t close the dialog')
                    end
                end
            else
                if isempty(S)
                    error('Please, select/type a proper answer and/or don''t close the dialog')
                end
            end
        end
        
        function exportSurfaceStatistics(SurfStruct, SurfAnswers, ScannerAnswers)
            S = SurfaceRoughnessProcessor.cleanUpStruct(SurfStruct);
            Results = struct2cell(S);
            Fields = fieldnames(S);
            Results = Results(4:end);
            Fields = Fields(4:end);
            [dirName, fileName] = SurfaceRoughnessProcessor.SurfAnswers2filename(SurfAnswers);
            ScannerName = {'Scanner name and model'; 'Scanner uncertainty (microns)';
                           'Unit of scan file/statistics'};
            ScannerInfo = {ScannerAnswers.Q1; ScannerAnswers.Q2; ScannerAnswers.Q3};
            SurfaceName = {'Kind'; 'Type'; 'Flow Type'; 'Results'; 'Descriptor'; 'Author';...
                           'Year'; 'Identifier'; 'doi URL'};
            SurfaceInfo = {SurfAnswers.Q1; SurfAnswers.Q2; SurfAnswers.Q3; SurfAnswers.Q4;...
                           SurfAnswers.Q5; SurfAnswers.Q6; SurfAnswers.Q7; SurfAnswers.Q8; SurfAnswers.Q9};
            varNames = {'Streamwise length of scan';
                        'Spanwise length of scan';
                        'Minimum roughness height';
                        'Maximum roughness height';
                        'Peak-to-trough roughness height';
                        'Average 5 peak-to-trough roughness height';
                        'Average roughness height';
                        'Average of absolute value of the height fluctuations';
                        'Root-mean-square of the total height';
                        'Root-mean-square of the height fluctuations';
                        'Skewness of the height fluctuations';
                        'Flatness of the height fluctuations';
                        'Effective Slope in the streamwise direction';
                        'Correlation length in the streamwise direction';
                        'Effective Slope in the spanwise direction';
                        'Correlation length in the spanwise direction';};
            C = [varNames Fields Results];
            C2 = [ScannerName ScannerInfo];
            C3 = [SurfaceName SurfaceInfo];
            SurfaceRoughnessProcessor.writeExcelResults(C, C2, C3, dirName, fileName);
            SurfaceRoughnessProcessor.exportSurfStats2mat(fullfile(dirName, fileName), S, SurfaceName, SurfAnswers, ScannerName, ScannerAnswers);
            SurfaceRoughnessProcessor.exportSurfaceData2mat(dirName, fileName, SurfStruct);
        end
        
        function writeExcelResults(C, C2, C3, pathname, filename)
            pathfile = fullfile(pathname, filename);
            if isfile(pathfile)
                Nf = dir([pathfile(1:end-5) '*.' pathfile(end-3:end)]);
                N = length(Nf);
                N = N + 1;
                filename = [filename(1:end-5) '(' num2str(N) ').' filename(end-3:end)];
            end
            SurfaceRoughnessProcessor.writeResults(C, pathname, filename);
            SurfaceRoughnessProcessor.writeResults(C2, pathname, filename, 'E1:F3')
            SurfaceRoughnessProcessor.writeResults(C3, pathname, filename, 'E5:F13', true)
        end
        
        function writeResults(C, pathname, filename, Range, flag_rename)
            if ~exist('flag_rename', 'var')
                flag_rename = false;
            end
            if ~ispc
                filename = fullfile(pathname, filename);
                if exist('Range', 'var')
                    writecell(C, filename, 'Range', Range);
                else
                    writecell(C, filename);
                end
            else
                cd(pathname)
                if exist('Range', 'var')
                    xlswrite('temp.xlsx', C, 'Sheet1', Range);
                else
                    xlswrite('temp.xlsx', C, 'Sheet1', 'A1');
                end
                if flag_rename == true
                    movefile('temp.xlsx', filename)
                end
            end
        end
        
        function [dirName, fileName] = SurfAnswers2filename(SurfAnswers)
            pattern = [SurfAnswers.Q1 '_'...
                       SurfAnswers.Q2 '_'...
                       SurfAnswers.Q3 '_'...
                       SurfAnswers.Q4 '_'...
                       SurfAnswers.Q5 '_'...
                       SurfAnswers.Q6 '_'...
                       SurfAnswers.Q7];
            folders = {'Surfaces', 'Flow documentation', 'Paper'};
            dirName = fullfile(pwd, pattern);
            if ~isfolder(dirName)
                mkdir(dirName)
                for n = 1:length(folders)
                    dir = fullfile(dirName, folders{n});
                    if ~isfolder(dir)
                        mkdir(dir)
                    end
                end
            end
            dirName = fullfile(dirName, folders{1});
            fileName = ['SurfaceStatistics_'...
                        SurfAnswers.Q6 '_'...
                        SurfAnswers.Q7 '_'...
                        SurfAnswers.Q8...
                        '.xlsx'];
        end
        
        function exportSurfStats2mat(filename, SurfStruct, SurfaceName, SurfAnswers, ScannerName, ScannerAnswers)
            fields = fieldnames(SurfAnswers);
            for n = 1:length(SurfaceName)
                f = fields{n};
                var = genvarname(SurfaceName{n});
                Surface.(var) = SurfAnswers.(f);
            end
            if strcmp(SurfAnswers.Q4, 'Experiments')
                fields = fieldnames(ScannerAnswers);
                for n = 1:length(ScannerName)
                    f = fields{n};
                    var = genvarname(ScannerName{n});
                    Surface.(var) = ScannerAnswers.(f);
                end
            end
            fields = fieldnames(SurfStruct);
            for n = 1:length(fields)
                f = fields{n};
                Surface.(f) = SurfStruct.(f);
            end
            ext = '.mat';
            filename = [filename(1:end-5) ext];
            if isfile(filename)
                Nf = dir([filename(1:end-4) '*' ext]);
                N = length(Nf);
                N = N + 1;
                filename = [filename(1:end-4) '(' num2str(N) ')' ext];
            end
            save(filename, 'Surface')
        end
        
        function exportSurfaceData2mat(pathname, filename, SurfStruct)
            varNames = SurfStruct.varNames;
            N = length(varNames);
            for n = 1:N
                var = varNames{n};
                SurfaceData.(var) = SurfStruct.obj.(var);
            end
            prt = 'Statistics'; l = length(prt);
            ind = strfind(filename, prt);
            ext = '.mat';
            filename = [filename(1:ind-1) 'Data' filename(ind+l:end-5) ext];
            filename = fullfile(pathname, filename);
            if isfile(filename)
                Nf = dir([filename(1:end-4) '*' ext]);
                N = length(Nf);
                N = N + 1;
                filename = [filename(1:end-4) '(' num2str(N) ')' ext];
            end
            save(filename, 'SurfaceData')
        end
        
        function visualizeResults(SurfStruct)
            type = SurfStruct.type;
            figure(1);
            vars = SurfStruct.varNames;
            switch type
                case '1D-profile'
                    vars = SurfStruct.varNames;
                    X = SurfStruct.obj.(vars{1});
                    Z = SurfStruct.obj.(vars{2});
                    p = plot(X, Z);
                    p.LineWidth = 1.5;
                    xlabel('x [mm]'), ylabel('z [mm]')
                    set(gca, 'FontName', 'Times', 'FontSize', 12)
                case '2D-surface'
                    X = SurfStruct.obj.(vars{1});
                    Y = SurfStruct.obj.(vars{2});
                    Z = SurfStruct.obj.(vars{3});
                    contourf(X, Y, Z)
                    axis equal tight
                    xlabel('x [mm]'), ylabel('y [mm]')
                    set(gca, 'FontName', 'Times', 'FontSize', 12)
                    c = colorbar;
                    c.Label.String = 'z [mm]';
            end
        end
        
        function displayResults(SurfStruct)
            S = SurfaceRoughnessProcessor.cleanUpStruct(SurfStruct);
            DataSet = struct2table(S, 'AsArray', true);
            disp('Roughness Statistics')
            disp(DataSet)
        end
    end
end