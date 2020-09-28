% This script processes surface roughness statistics to support the
% roughness database hosted by the University of Southampton.
%
% The development idea of this script is to require minimum users
% interface when running the code. In other words, there is no need for
% code manipulation or changes to correctly export the roughness statistics.
%
% Two type of surfaces are supported: 1D-line profiles, or 2D-surface
% scans.
%
% The input format supported are MATLAB (*.mat), Excel (*.xls) and 
% ASCII (*.csv or tab-delimited *.txt, or *.dat).
%
% We however ask the users to comply with the standards set for the input
% files containing the roughness height information.
%
% For MATLAB files, put the roughness information in data into variables
% X,Y,Z and use the save(...) function to save the roughness information.
%
% For ASCII/CSV files, the X,Y,Z information should be format in each column.
%
% The supporting functions are nested in the bottom of this script.
% Changing them may likely break the code and/or calculations.
%
% If a bug is detected, either open an issue on github at 
% https://github.com/jmbarrojr/SurfaceRoughnessStatistics/issues
% or send an email to julio.barros@gmail.com using "[BUG]" prefix in the
% email's subject.
%
%
% Authors: Julio Barros (OIST) and Karen Flack (USNA) - 2020

clc, clear, close all

%% PARAMETER INPUTS
% This is useful if you are processing multiple scans/surface in which most 
% of the answers are the same. View 'Questionnaire_Batch.txt' and
% 'Profiler_Batch.txt' for more information.
batch = true;

%% MAIN SECTION
% Choose file to analyze
[filename,pathname] = uigetfile({'*.mat;*.csv;*.asc;*.dat;*.xls;*.xlsx',...
                       'Surface Files (*.mat,*.csv,*.asc,*.dat,*.xlx,*.xlsx)';...
                                 '*.*','All Files(*.*)'},...
    'Choose a Matlab data file with x,y,z coordinates.');
if isempty(filename) == 1 || ~ischar(filename)
    error('No file was selected')
end
% Or  paste in the file
% filename = ['Processed_Surface02_8_12_25grit_CURVTILT.mat'];

% Roughness and Scanner Questionare
SurfAnswers = RoghnessQuestionnaire(batch);
ScannerAnswers = ScannerQuestionnaire(SurfAnswers,batch);
        
% Run function to calculate statistics
Surface = getSurfStatistics(fullfile(pathname,filename));

% Export Surface Statistics
exportSurfaceStatistics(Surface,SurfAnswers,ScannerAnswers)

% Visualize Surface
visualizeResults(Surface)

% Display results on command prompt
displayResults(Surface)

%% ---------------- ROUGHNESS NESTED FUNCTIONS ---------------------------
% MAIN FUNCTION -----------------------------------------------------------
function Surface = getSurfStatistics(filename)
Surface = loadSurface(filename);
Surface = determineSurfaceType(Surface);
Surface = determineXandYdir(Surface);
Surface = roughPhysicalProp(Surface);
Surface = roughnessStats(Surface);
%Surface = cleanUpStruct(Surface);
end
% LOADING MATLAB FUNCTIONS ------------------------------------------------
function SurfStruct = loadSurface(filename)
ext = filename(end-3:end);
switch ext
    % MATLAB FORMAT FILE
    case '.mat'
        matObj = matfile(filename, 'Writable', false);
        SurfStruct.obj = matObj;
        %SurfStruct.varNames = who(matObj);
        SurfStruct.varProps = whos(matObj);
    % TEXT BASED FILES    
    case {'.asc','.ASC','.dat','.DAT','.csv'}
        SurfStruct = importSurfaceTxtData(filename);
    case {'.xls','xlsx'}
        SurfStruct = importSurfaceExcel(filename);
    otherwise
        error('File format might not be yet implemented')
end
% TESTING
% Remove non-numerical variables
N = length(SurfStruct.varProps);
c = 1;
for n=1:N
    cls = SurfStruct.varProps(n).class;
    name = SurfStruct.varProps(n).name;
    if strcmp(cls,'double')
        varNames{c} = name; %#ok<AGROW>
        c = c + 1;
    end
end
SurfStruct.varNames = varNames;
end
% IMPORT TEXT FILES -------------------------------------------------------
function SurfStruct = importSurfaceTxtData(filename)
A = importdata(filename);
L = size(A,2);
if L == 2
    X = A(:,1);
    Z = A(:,2);
    save('temp.mat','X','Z')
elseif L == 3
    % Find # of columns and rows
    A = sortrows(A,[1,2]);
    J = find(A(1,1)-A(:,1) == 0, 1 , 'last');
    A = sortrows(A,[2,1]);
    I = find(A(1,2)-A(:,2) == 0, 1 , 'last');
    X = reshape(A(:,1),I,J)';
    Y = reshape(A(:,2),I,J)';
    Z = reshape(A(:,3),I,J)';
    save('temp.mat','X','Y','Z')
end
SurfStruct = loadSurface('temp.mat');
end
% IMPORT EXCEL FILES ------------------------------------------------------
function SurfStruct = importSurfaceExcel(filename)
A = readmatrix(filename);
L = size(A,2);
if L == 2
    X = A(:,1);
    Z = A(:,2);
    save('temp.mat','X','Z')
elseif L == 3
    % Find # of columns and rows
    A = sortrows(A,[1,2]);
    J = find(A(1,1)-A(:,1) == 0, 1 , 'last');
    A = sortrows(A,[2,1]);
    I = find(A(1,2)-A(:,2) == 0, 1 , 'last');
    X = reshape(A(:,1),I,J)';
    Y = reshape(A(:,2),I,J)';
    Z = reshape(A(:,3),I,J)';
    save('temp.mat','X','Y','Z')
end
SurfStruct = loadSurface('temp.mat');
end
% TYPE OF SCAN FUNCTION ---------------------------------------------------
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
% SURFACE RELATED FUNCTIONS -----------------------------------------------
function SurfStruct = determineXandYdir(SurfStruct)
zname = SurfStruct.varNames{end}; % A bit of a strech assuming the
% height information is the last variable.
z = SurfStruct.obj.(zname);
[sr,~] = size(z); % Size of rows and columns
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
        %yname = SurfStruct.varNames{2};
        x = SurfStruct.obj.(xname);
        if abs(x(1,2) - x(1,1)) == 0
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
% -------------------------------------------------------------------------
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
% SURFACE ROUGHNESS STATISTICS --------------------------------------------
function SurfStruct = roughnessStats(SurfStruct)
zname = SurfStruct.varNames{end};

% Height information is the last variable.
z = SurfStruct.obj.(zname);

% Minimum roughness height
kmin = min(z(:));
SurfStruct.kmin = kmin;

% Maximum roughness height
kmax = max(z(:));
SurfStruct.kmax = kmax;

% Peak-to-trough roughness height (kt)
SurfStruct.kt = kmax - kmin;

% Peak-to-trough roughness height average 5 (kz5)
MaxZ = sortrows(z(:),'descend');
MinZ = sortrows(z(:),'ascend');
SurfStruct.kz5 = mean(MaxZ(1:5)) - mean(MinZ(1:5));

% Mean roughness elevation
SurfStruct.kbar = mean(z(:));

% Average roughness heigth
SurfStruct.ka = mean(abs(z(:)));

% Root-mean-square of the height
SurfStruct.krms_z = rms(z(:));

% Average roughness height (h')
z = z - mean(z(:));

% Root-mean-square of the height with mean height removed (krms)
SurfStruct.krms = rms(z(:));

% Skewness of the roughness height fluctuation
SurfStruct.Sk = skewness(z(:));

% Flatness of the roughness height fluctuation
SurfStruct.Fl = kurtosis(z(:));

% Effective slope and correlation length
xname = SurfStruct.varNames{1};
x = SurfStruct.obj.(xname);
SurfStruct.ESx = EffectiveSlope(x,z,SurfStruct.Xdir,SurfStruct.Lx);
SurfStruct.Rlx = CorrelationLenght(x,z,SurfStruct.Xdir);
switch SurfStruct.type
    case '2D-surface'
        yname = SurfStruct.varNames{2};
        y = SurfStruct.obj.(yname);
        SurfStruct.ESy = EffectiveSlope(y,z,SurfStruct.Ydir,SurfStruct.Ly);
        SurfStruct.Rly = CorrelationLenght(y,z,SurfStruct.Ydir);
    otherwise
        SurfStruct.ESy = [];
        SurfStruct.Rly = [];
end
end
% EFFECTIVE SLOPE ---------------------------------------------------------
function Es = EffectiveSlope(X,Z,dir,L)
if dir == 1
    dx = X(2,1) - X(1,1);
    mean_dir = 2;
else
    dx = X(1,2) - X(1,1);
    mean_dir = 1;
end
dzdx = diff(Z,1,dir) ./ diff(X,1,dir);
Es = 1/L .* mean( trapz(abs(dzdx),dir).*dx, mean_dir); % Efective Slope;
end
% CORRELATION LENGHTSCALE -------------------------------------------------
function Rlx = CorrelationLenght(X,Z,dir)
if dir == 1
    dx = X(2,1) - X(1,1);
else
    dx = X(1,2) - X(1,1);
end
s = size(Z,dir);
% Compute the correlation length in the given direction
[lags,Zcorr] = MeanAutoCorr_FFT(Z,dir);
lags = lags.*dx;
% DEGUB
% figure,plot(lags,Zcorr)
% ind = findSlopeCorr(Zcorr);
% Rlx = interp1(Zcorr(s:s+ind),lags(s:s+ind),1./exp(1),'linear');
Rlx = interp1(Zcorr(s:end),lags(s:end),1./exp(1),'linear');
end
% -------------------------------------------------------------------------
function [lags,C] = MeanAutoCorr_FFT(A,dir)
S = size(A);
s = S(dir);                      % Get the size of A in the desired direction
nfft = 2*s-1;                    % Make the FFT size to be 2*s-1 to have lag=0 at the center
A = A - mean(A(:));              % Remove any mean in the surface
Afft = fft(A,nfft,dir);          % Compute the FFT of A
Afft = Afft ./s;                 % To unbiased the correlation magnitude
C = Afft .* conj(Afft);          % Compute the correlation via FFT
C = fftshift(ifft(C,[],dir),dir);% Make the zero lag at the center

% Ensemble average only if it's a 2d surface
if S(1) > 1 || S(2) > 1
    if dir == 1
        C = mean(C,dir+1);
    elseif dir == 2
        C = mean(C,dir-1);
    end
end
C = C ./ max(C(:)); % Normalize the correlation
lags = [-linspace(s-1,1,s-1) 0 linspace(1,s-1,s-1)];
% DEGUB
% figure,plot(C),title('AutoCorr')
end
% -------------------------------------------------------------------------
function ind = findSlopeCorr(C)
l = length(C);
C = C((l+1)/2:end); % Get only half of C
ind = find(diff(C) > 0,1);
end
% -------------------------------------------------------------------------
function S = cleanUpStruct(SurfStruct)
S = struct();
fields = fieldnames(SurfStruct);
for n=4:length(fields)
    f = fields{n};
    S.(f) = SurfStruct.(f);
end
end

%% ----------------- QUESTIONNAIRE FUNCTIONS ------------------------------
function SurfAnswers = RoghnessQuestionnaire(batch)
switch batch
    case true
        SurfAnswers = loadQuestionnaire('Questionnaire_Batch.txt');
    otherwise
        % Get information to name the Excel output file and directory
        % Surface information
        prompt = 'Is the surface Homogeneous or Heterogeneous? ';
        q1 = 'Hom'; q2 = 'Het';
        SurfAnswers.Q1 = questdlg(prompt,'Roughness Information',...
            q1,q2,q1);
        checkAnswer(SurfAnswers.Q1);
        
        prompt = 'Is the roughness Regular or Irregular?';
        q1 = 'Reg'; q2 = 'Irreg';
        SurfAnswers.Q2 = questdlg(prompt,'Roughness Information',...
            q1, q2, q1);
        checkAnswer(SurfAnswers.Q2);
        
        prompt = 'Are results for a ...';
        q1 = 'TBL'; q2 = 'Pipe'; q3 = 'Channel';
        SurfAnswers.Q3 = questdlg(prompt,'Roughness Information',...
            q1, q2, q3, q1);
        checkAnswer(SurfAnswers.Q3);
        
        prompt = 'Are results from Experiments or Simulations?';
        q1 = 'Exp'; q2 = 'Sim';
        SurfAnswers.Q4 = questdlg(prompt,'Roughness Information',...
            q1,q2,q1);
        checkAnswer(SurfAnswers.Q4);
        
        prompt = 'What is the general descriptor of this surface, .i.e. "Sandgrain"?';
        SurfAnswers.Q5 = inputdlg(prompt,'Roughness Information',[1 50]);
        SurfAnswers.Q5 = SurfAnswers.Q5{1};
        checkAnswer(SurfAnswers.Q5);
        
        prompt1 = 'What is the last name of the lead author of the study? ';
        prompt2 = 'What year were the results published? ';
        prompt3 = 'What is the identifying name of this surface, i.e. "220Grit"? ';
        prompt4 = 'Include the doi of the publication';
        temp = inputdlg({prompt1,prompt2,prompt3,prompt4},...
            'Roughness Information',...
            [1 50; 1 50; 1 50; 1 50]);
        SurfAnswers.Q6 = temp{1};
        checkAnswer(SurfAnswers.Q6)
        SurfAnswers.Q7 = temp{2};
        checkAnswer(SurfAnswers.Q7)
        SurfAnswers.Q8 = temp{3};
        checkAnswer(SurfAnswers.Q8)
        SurfAnswers.Q9 = temp{4};
        checkAnswer(SurfAnswers.Q9)
end
end
% -------------------------------------------------------------------------
function ScannerAnswers = ScannerQuestionnaire(SurfAnswers,batch)
if strcmp(SurfAnswers.Q4,'Exp')
    switch batch
        case true
            ScannerAnswers = loadQuestionnaire('Profiler_Batch.txt');
        otherwise
            prompt1 = 'What is the name and model of the profiler/scanner? ';
            prompt2 = 'What is the uncertainty in the measurement of surface heights in microns? ';
            prompt3 = ['What is the unit of measurements of the scan(file),'...
                      'e.g. mm (milimeters), um (microns), in(inches)? '];
            temp = inputdlg({prompt1,prompt2,prompt3},...
                'Scanner Information',[1 50;1 50;1 50]);
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
% -------------------------------------------------------------------------
function SurfAnswers = loadQuestionnaire(file)
% Open the questionnaire
fid = fopen(file,'r');
n=1;
% Look for the answers in the questionnaire
while ~feof(fid)
    line = strtrim(fgetl(fid));
    
    if isempty(line) || all(isspace(line)) || strncmp(line, '#', 1)
    else
        Q = ['Q' num2str(n)];
        SurfAnswers.(Q) = line;
        n = n+1;
    end
end
% Close the file
fclose(fid);
end
% -------------------------------------------------------------------------
function checkAnswer(S)
if iscell(S) && ~isempty(S)
    for n=1:length(S)
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

%% -------------------- EXPORT FUNCTIONS ----------------------------------
function exportSurfaceStatistics(SurfStruct,SurfAnswers,ScannerAnswers)
% Cleanup Struct for Excel file
S = cleanUpStruct(SurfStruct);

% Put surface statistics in an Excel file structure
Results = struct2cell(S);
Fields = fieldnames(S);
% Chop unnecessary information for Excel file
Results = Results(4:end);
Fields = Fields(4:end);

% Prepare base directory and filename
[dirName,fileName] = SurfAnswers2filename(SurfAnswers);

% Prepare Scanner information
ScannerName = {'Scanner name and model';'Scanner uncertainty (microns)';
               'Unit of scan file/statistics'};
ScannerInfo = {ScannerAnswers.Q1;ScannerAnswers.Q2;ScannerAnswers.Q3};

% Prepare Surface information
SurfaceName = {'Kind';'Type';'Flow Type';'Results';'Descriptor';'Author';...
               'Year';'Identifier';'doi URL'};
SurfaceInfo = {SurfAnswers.Q1;SurfAnswers.Q2;SurfAnswers.Q3;SurfAnswers.Q4;...
SurfAnswers.Q5;SurfAnswers.Q6;SurfAnswers.Q7;SurfAnswers.Q8;SurfAnswers.Q9};

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

% Concatenate cells for Excel file
C = [varNames Fields Results];
C2 = [ScannerName ScannerInfo];
C3 = [SurfaceName SurfaceInfo];

% Write data to Excel spreadsheet
writeExcelResults(C,C2,C3,dirName,fileName)

% Write Surface Statistics to MATLAB file
exportSurfStats2mat(fullfile(dirName,fileName),S,...
                        SurfaceName,SurfAnswers,ScannerName,ScannerAnswers)

% Write Suface Data to MATLAB file
exportSurfaceData2mat(dirName,fileName,SurfStruct)
end
% Wrapper for Saving Statistics In Excel ----------------------------------
function writeExcelResults(C,C2,C3,pathname,filename)
pathfile = fullfile(pathname,filename);
if isfile(pathfile)
    Nf = dir([pathfile(1:end-5) '*.' pathfile(end-3:end)]);
    N = length(Nf);
    N = N+1;
    filename = [filename(1:end-5) '(' num2str(N) ').' filename(end-3:end)];
end
writeResults(C,pathname,filename);
writeResults(C2,pathname,filename,'E1:F3')
writeResults(C3,pathname,filename,'E5:F13',true)
end
% Excel export function ---------------------------------------------------
function writeResults(C,pathname,filename,Range,flag_rename)
if ~exist('flag_rename','var')
    flag_rename = false;
end
if ~ispc
    filename = fullfile(pathname,filename);
    if exist('Range','var')
        writecell(C,filename,'Range',Range);
    else
        writecell(C,filename);
    end
else
    % This might fix a Windows issue
    cd(pathname)
    if exist('Range','var')
        xlswrite('temp.xlsx',C,'Sheet1',Range);
        %writecell(C,'temp.xls','Range',Range,'UseExcel',false);
    else
        xlswrite('temp.xlsx',C,'Sheet1','A1');
        %writecell(C,'temp.xls','UseExcel',false);
    end
    % This might fix a windows issue
    if flag_rename == true
        movefile('temp.xlsx',filename)
    end
end
end
% -------------------------------------------------------------------------
function [dirName,fileName] = SurfAnswers2filename(SurfAnswers)
pattern = [SurfAnswers.Q1 '_'...
    SurfAnswers.Q2 '_'...
    SurfAnswers.Q3 '_'...
    SurfAnswers.Q4 '_'...
    SurfAnswers.Q5 '_'...
    SurfAnswers.Q6 '_'...
    SurfAnswers.Q7 '_'...
    SurfAnswers.Q8];

% Make the directory structure
folders = {'Surfaces','Flow documentation','Paper'};
dirName = fullfile(pwd,pattern);
if ~isfolder(dirName)
    mkdir(dirName)
    for n=1:length(folders)
        dir = fullfile(dirName,folders{n});
        if ~isfolder(dir)
            mkdir(dir)
        end
    end
end
% This is the directory in which the surface related files will be saved
dirName = fullfile(dirName,folders{1});

fileName = ['SurfaceStatistics_'...
...%    SurfAnswers.Q1 '_'...
...%    SurfAnswers.Q2 '_'...
...%    SurfAnswers.Q3 '_'...
...%    SurfAnswers.Q4 '_'...
...%    SurfAnswers.Q5 '_'...
    SurfAnswers.Q6 '_'...
    SurfAnswers.Q7 '_'...
    SurfAnswers.Q8...
    '.xlsx'];
%fileName = 'SurfaceStatistics.xls';
end
% EXPORT SURFACE STATS TO MATLAB ------------------------------------------
function exportSurfStats2mat(filename,SurfStruct,SurfaceName,SurfAnswers,...
                                                ScannerName,ScannerAnswers)
% Pass Surface Questionnaire
fields = fieldnames(SurfAnswers);
for n=1:length(SurfaceName)
    f = fields{n};
    var = genvarname(SurfaceName{n});
    Surface.(var) = SurfAnswers.(f);
end
%Surface.Author = SurfAnswers.Q6;
%Surface.year = SurfAnswers.Q7;

% Pass Scanner Questionnaire
if strcmp(SurfAnswers.Q4,'Experiments')
    fields = fieldnames(ScannerAnswers);
    for n=1:length(ScannerName)
        f = fields{n};
        var = genvarname(ScannerName{n});
        Surface.(var) = ScannerAnswers.(f);
    end
end

% Pass Surface Statistics
fields = fieldnames(SurfStruct);
for n=1:length(fields)
    f = fields{n};
    Surface.(f) = SurfStruct.(f);
end

% Make MATLAB filename
ext = '.mat';
filename = [filename(1:end-5) ext];

% Check if file alredy exist
if isfile(filename)
    Nf = dir([filename(1:end-4) '*' ext]);
    N = length(Nf);
    N = N + 1;
    filename = [filename(1:end-4) '(' num2str(N) ')' ext];
end
save(filename,'Surface')
end

% EXPORT SURFACE DATA TO MATLAB -------------------------------------------
function exportSurfaceData2mat(pathname,filename,SurfStruct)
varNames = SurfStruct.varNames;
N = length(varNames);
for n=1:N
    var = varNames{n};
    SurfaceData.(var) = SurfStruct.obj.(var);
end
prt = 'Statistics'; l = length(prt);
ind = strfind(filename,prt);
% Make MATLAB file name
ext = '.mat';
filename = [filename(1:ind-1) 'Data' filename(ind+l:end-5) ext];
filename = fullfile(pathname,filename);
% Check if file alredy exist
if isfile(filename)
    Nf = dir([filename(1:end-4) '*' ext]);
    N = length(Nf);
    N = N + 1;
    filename = [filename(1:end-4) '(' num2str(N) ')' ext];
end
save(filename,'SurfaceData')
end

%% ------------------- VISUALIZE RESULTS ----------------------------------
function visualizeResults(SurfStruct)
type = SurfStruct.type;
figure(1);
vars = SurfStruct.varNames;
switch type
    case '1D-profile'
        vars = SurfStruct.varNames;
        X = SurfStruct.obj.(vars{1});
        Z = SurfStruct.obj.(vars{2});
        p = plot(X,Z);
        p.LineWidth = 1.5;
        xlabel('x [mm]'),ylabel('z [mm]')
        set(gca,'FontName','Times','FontSize',12)
    case '2D-surface'
        X = SurfStruct.obj.(vars{1});
        Y = SurfStruct.obj.(vars{2});
        Z = SurfStruct.obj.(vars{3});
        contourf(X,Y,Z)
        axis equal tight
        xlabel('x [mm]'),ylabel('y [mm]')
        set(gca,'FontName','Times','FontSize',12)
        c = colorbar;
        c.Label.String = 'z [mm]';
end
end
% -------------------------------------------------------------------------
function displayResults(SurfStruct)
S = cleanUpStruct(SurfStruct);
DataSet = struct2dataset(S); %#ok<STRUCTDTSET>
disp('Roughness Statistics')
disp(DataSet)
end