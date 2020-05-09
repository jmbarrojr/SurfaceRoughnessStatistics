% This script processes surface roughness statistics to support the
% roughness database hosted by the University of Southampton.
%
% The development idea of this script is to require minimum users
% interface when running the code. In other words, there is no need for
% code manipulation or changes to fully export the roughness statistics.
%
% Two type of surfaces are supported: 1D line profiles, or 3D surface
% scans.
%
% However, we ask the user to comply with the standard set of the input
% file containing the roughness height information. The input formart
% supported are MATLAB (*.mat), Excel (*.xls) and ASCII (*.csv or
% tab-delimeted *.txt, or *.dat).
%
% For MATLAB files, put the roughness infomation in data into variables
% X,Y,Z and use the save(...) function to save the roughness information.
%
% For ASCII files, the X,Y,Z information should be format in each column.
%
% The supporting functions are nested in the bottom of this script.
% Changing them may like break the code. If a bug was detected, either
% open an issue github at https://github.com/jmbarrojr/SurfaceRoughnessStatistics/issues
% or send an email to julio.barros@gmail.com using "[BUG]" prefix in the
% email's subject.
%
%
% Authors: Julio Barros (OIST) and Karen Flack (USNA)

clc, clear

%% MAIN SECTION
% Choose file to analyze
[filename,pathname] = uigetfile({'*.mat';'*.txt';'*.csv'},...
                      'Choose a Matlab data file with x,y,z coordinates.');
if isempty(filename) == 1 || ~ischar(filename)
    error('No file was selected')
end
% Or  paste in the file
% filename = ['Processed_Surface02_8_12_25grit_CURVTILT.mat'];

% Roughness and Scanner Questionare
SurfAnswers = RoghnessQuestionnaire();
ScannerAnswers = ScannerQuestionnaire(SurfAnswers);

% Run function to calculate statistics
Surface = getSurfProperties(fullfile(pathname,filename));

% Export Surface Statistics
exportSurfaceStatistics(Surface,SurfAnswers,ScannerAnswers)

%% SUPPORTING FUNCTIONS
% MAIN FUNCTION -----------------------------------------------------------
function Surface = getSurfProperties(filename)
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
    case '.mat'
        matObj = matfile(filename, 'Writable', false);
        SurfStruct.obj = matObj;
        %SurfStruct.varNames = who(matObj);
        SurfStruct.varProps = whos(matObj);
    otherwise
        error('Just matlab file loader implemented')
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

% Peak-to-trough roughness height
SurfStruct.kp = kmax - kmin;

% Peak-to-trough roughness height average 5
MaxZ = sortrows(z(:),'descend');
MinZ = sortrows(z(:),'ascend');
SurfStruct.kp5 = mean(MaxZ(1:5)) - mean(MinZ(1:5));

% Average roughness height over entire surface
SurfStruct.kbar = mean(z(:));

% Root-mean-square of the height
SurfStruct.ka = mean(abs(z(:)));

% Remove mean height
SurfStruct.krms = rms(z(:));

% Average roughness height
z = z - mean(z(:));

% Root-mean-square of the height with mean height removed
SurfStruct.krmsp = rms(z(:));

% Skewness of the roughness height fluctuation
SurfStruct.Sk = skewness(z(:));

% Flatness of the roughness height fluctuation
SurfStruct.Fl = kurtosis(z(:));

% Effective slope and correlation length
xname = SurfStruct.varNames{1};
x = SurfStruct.obj.(xname);
SurfStruct.Esx = EffectiveSlope(x,z,SurfStruct.Xdir);
SurfStruct.Rlx = CorrelationLenght(x,z,SurfStruct.Xdir);
switch SurfStruct.type
    case '2D-surface'
        yname = SurfStruct.varNames{2};
        y = SurfStruct.obj.(yname);
        SurfStruct.Esy = EffectiveSlope(y,z,SurfStruct.Ydir);
        SurfStruct.Rlz = CorrelationLenght(y,z,SurfStruct.Ydir);
end

end
% -------------------------------------------------------------------------
function Es = EffectiveSlope(X,Z,dir)
if dir == 1
    dx = X(2,1) - X(1,1);
    mean_dir = 2;
else
    dx = X(1,2) - X(1,1);
    mean_dir = 1;
end 
L = max(X(:)) - min(X(:));
[dzdx,~] = gradient(Z,dx);
Es = 1/L .* mean( trapz(abs(dzdx),dir).*dx, mean_dir); % Efective Slope;
end
% -------------------------------------------------------------------------
function Rlx = CorrelationLenght(X,Z,dir)
if dir == 1
    dx = X(2,1) - X(1,1);
else
    dx = X(1,2) - X(1,1);
end 
s = size(Z,dir);
% Compute the correlation length in the Y-dir
[lags,Zcorr] = MeanAutoCorr_FFT(Z,dir);
lags = lags.*dx;
ind = findSlopeCorr(Zcorr);
Rlx = interp1(Zcorr(s:s+ind),lags(s:s+ind),1./exp(1),'linear');
end
% -------------------------------------------------------------------------
function [lags,C] = MeanAutoCorr_FFT(A,dir)
S = size(A);                 
s = S(dir);                      % Get the size of A in the desired direction
nfft = 2*s-1;                    % Make the FFT size to be 2*s-1 to make the center as 0
A = A - mean(A(:));              % Remove any mean in the surface
Afft = fft(A,nfft,dir);          % Compute the FFT of A
Afft = Afft ./s;                 % To unbiase the correlation
C = Afft .* conj(Afft);          % Compute the correlation via FFT
C = fftshift(ifft(C,[],dir),dir);% Make the zero lag in the center

% Ensemble average only if it's a 2d surface
if S(1) > 1 || S(2) > 1
    if dir == 1
        C = mean(C,dir+1);
    elseif dir == 2
        C = mean(C,dir-1);
    end
end
C = C ./ max(C(:)); % Normalize the correlation
lags = [linspace(-(s-1),1,s-1) 0 linspace(1,s-1,s-1)];
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

% QUESTIONNAIRE FUNCTIONS ---------------------------------------------------
function SurfAnswers = RoghnessQuestionnaire()
% Get information to name the Excel output file
% Surface information
prompt = 'What kind of suface is it? ';
q1 = 'Homogeneous'; q2 = 'Heterogeneous';
SurfAnswers.Q1 = questdlg(prompt,'Roughness Information',...
                     q1,q2,q1);
checkAnswer(SurfAnswers.Q1);

prompt = 'Is the roughness ...';
q1 = 'Regular'; q2 = 'Irregular';
SurfAnswers.Q2 = questdlg(prompt,'Roughness Information',...
                     q1, q2, q1);
checkAnswer(SurfAnswers.Q2);

prompt = 'Are results for a ...';
q1 = 'TBL'; q2 = 'Pipe'; q3 = 'Channel';
SurfAnswers.Q3 = questdlg(prompt,'Roughness Information',...
                     q1, q2, q3, q1);
checkAnswer(SurfAnswers.Q3);

prompt = 'Are results from ...';
q1 = 'Experiments'; q2 = 'Simulations';
SurfAnswers.Q4 = questdlg(prompt,'Roughness Information',...
                     q1,q2,q1);
checkAnswer(SurfAnswers.Q4);

prompt = 'What is the general descriptor of this surface, .i.e. "Sandgrain"?';
SurfAnswers.Q5 = inputdlg(prompt,'Roughness Information',[1 50]);
checkAnswer(SurfAnswers.Q5);

prompt1 = 'What is the last name of the lead author of the study? ';
prompt2 = 'What year were the results published? ';
prompt3 = 'What is the identifying name of this surface, i.e. "220Grit"? ';
SurfAnswers.Q678 = inputdlg({prompt1,prompt2,prompt3},...
                              'Roughness Information', [1 50; 1 50; 1 50]);
checkAnswer(SurfAnswers.Q678)
end
% -------------------------------------------------------------------------
function ScannerAnswers = ScannerQuestionnaire(SurfAnswers)
if strcmp(SurfAnswers.Q4,'Experiments')
prompt1 = 'What is the name and model of the profiler/scanner? ';
prompt2 = 'What is the uncertainty in the measurement of surface heights in microns? ';
ScannerAnswers.Q1 = inputdlg({prompt1,prompt2},...
                             'Scanner Information',[1 50;1 50]);
else
    ScannerAnswers.Q1 = {'N/A';'N/A'};
end
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

% EXPORT FUNCTIONS --------------------------------------------------------
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
ScannerName = {'Scanner name and model';'Scanner uncertainty (microns)'};
ScannerInfo = ScannerAnswers.Q1;

varNames = {'Streamwise length of scan (mm)';
    'Spanwise length of scan (mm)';
    'Minimum roughness height (mm)';
    'Maximum roughness height (mm)';
    'Peak-to-trough roughness height (mm)';
    'Average 5 peak-to-trough roughness height (mm)';
    'Average roughness height (mm)';
    'Average of absolute value of the height fluctuations (mm)';
    'Root-mean-square of the total height (mm)';
    'Root-mean-square of the height fluctuations (mm)';
    'Skewness of the height fluctuations';
    'Flatness of the height fluctuations';
    'Effective Slope in the steamwise direction';
    'Correlation length in the steamwise direction';
    'Effective Slope in the spanwise direction';
    'Correlation length in the spanwise direction';};

% Concatenate cells for Excel file
C = [varNames Fields Results];
C2 = [ScannerName ScannerInfo];
% Write data to Excel spreadsheet
writecell(C,fullfile(dirName,fileName));
writecell(C2,fullfile(dirName,fileName),'Range','E1:F2')

% Write Surface Statistics to MATLAB file
exportSurfStats2mat(dirName,fileName,SurfStruct,SurfAnswers,ScannerAnswers)

% Write Suface Data to MATLAB file
exportSurfaceData2mat(dirName,fileName,SurfStruct)
end
% -------------------------------------------------------------------------
function [dirName,fileName] = SurfAnswers2filename(SurfAnswers)
pattern = [SurfAnswers.Q1 '_'...
           SurfAnswers.Q2 '_'...
           SurfAnswers.Q3 '_'...
           SurfAnswers.Q4 '_'...
           SurfAnswers.Q5{1} '_'...
           SurfAnswers.Q678{1} '_'...
           SurfAnswers.Q678{2} '_'...
           SurfAnswers.Q678{3}];
       
dirName = fullfile(pwd,pattern);
if ~isfolder(dirName)
    mkdir(dirName)
end
fileName = ['SurfaceStatistics_'...
           SurfAnswers.Q1 '_'...
           SurfAnswers.Q2 '_'...
           SurfAnswers.Q3 '_'...
           SurfAnswers.Q4 '_'...
           SurfAnswers.Q5{1} '_'...
           SurfAnswers.Q678{1} '_'...
           SurfAnswers.Q678{2} '_'...
           SurfAnswers.Q678{3} '.xls'];
end
% -------------------------------------------------------------------------
function exportSurfStats2mat(pathname,filename,SurfStruct,SurfAnswers,ScannerAnswers)
Surface.Author = SurfAnswers.Q678{1};
Surface.year = SurfAnswers.Q678{2};
fields = fieldnames(SurfStruct);
for n=1:length(fields)
    f = fields{n};
    Surface.(f) = SurfStruct.(f);
end
filename = [filename(1:end-4) '.mat'];
save(fullfile(pathname,filename),'Surface')
end
% -------------------------------------------------------------------------
function exportSurfaceData2mat(dirName,fileName,SurfStruct)
type = SurfStruct.type;
switch type
    case '1D-profile'
        SurfaceData.X = SurfStruct.obj.X;
        SurfaceData.Z = SurfStruct.obj.Z;
    case '2D-surface'
        SurfaceData.X = SurfStruct.obj.X;
        SurfaceData.Y = SurfStruct.obj.Y;
        SurfaceData.Z = SurfStruct.obj.Z;
end
prt = 'Statistics'; l = length(prt);
ind = strfind(fileName,prt);
fileName = [fileName(1:ind-1) 'Data' fileName(ind+l:end-4) '.mat'];
save(fullfile(dirName,fileName),'SurfaceData')
end