%This program was written by Julio Barros (OIST) and Karen Flack (USNA)
%Last updated 8 APR 2020
clc
clear

%Choose file to analyze
disp('Choose a Matlab data file with x,y,z coordinates.')
[file,path] = uigetfile('*.mat');
filename=file;

%Or  paste in the file
%filename = ['Filename.mat'];

%Get information to name the Excel output file
disp('--------------------------------------------------------------------------')
disp('Please enter the following information for a consistent naming convention.')
disp('Your response should be the word in "quotes".')
disp(' ')
%Surface information
prompt = 'Is the roughness "homogeneous" or "heterogeneous"? ';
SurfTypeOne = input(prompt,'s');
disp(' ')
prompt = 'Is the roughness "regular" or "irregular"? ';
SurfTypeTwo = input(prompt,'s');
disp(' ')
prompt = 'Are results for a "TBL", "pipe" or "channel"? ';
SurfTypeThree = input(prompt,'s');
disp(' ')
prompt = 'Are results from "experiments" or "simulations"? ';
SurfTypeFour = input(prompt,'s');
disp(' ')
prompt = 'What is the last name of the lead author of the study? ';
SurfTypeFive = input(prompt,'s');
disp(' ')
prompt = 'What year were the results published? ';
SurfTypeSix = input(prompt,'s');
disp(' ')
prompt = 'What is the identifying name of this surface, i.e. "Sandpaper_2"? ';
SurfTypeSeven = input(prompt,'s');
disp('--------------------------------------------------------------------------')

SurfName=append(SurfTypeOne,"_",SurfTypeTwo,"_",SurfTypeThree,"_",...
    SurfTypeFour,"_",SurfTypeFive,"_",SurfTypeSix,"_",SurfTypeSeven,".xls");

%Scanner information for Excel file
disp('Please enter the following information about the profiler/scanner. ')
disp(' ')
prompt = 'What is the name and model of the profiler/scanner? ';
ScannerName = input(prompt,'s');
disp(' ')
prompt = 'What is the uncertainty in the measurement of surface heights in microns? ';
Uncertainty = input(prompt,'s');
Scanner={'Scanner name and model';'Scanner uncertainty (microns)'};
ScannerInfo={ScannerName;Uncertainty};

% Run function to calculate statistics
Surface = getSurfProperties(filename);

% Put information in an excel file
Results = struct2cell(Surface);
Fields=fieldnames(Surface);
Names = {'Streamwise length of scan (mm)';
    'Spanwise length of scan (mm)';
    'Minimum roughness height (mm)';
    'Minimum roughness height (mm)';
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

%Write data to Excel spreadsheet
xlswrite(SurfName,Scanner,'A1:A2');
xlswrite(SurfName,ScannerInfo,'B1:B2');
xlswrite(SurfName,Names,'D1:D16');
xlswrite(SurfName,Fields(7:22),'E1:E16');
xlswrite(SurfName,Results(7:22),'F1:F16');


%% NESTED FUNCTIONS
function Surface = getSurfProperties(filename)
Surface = loadSurface(filename);
Surface = determineSurfaceType(Surface);
Surface = determineXandYdir(Surface);
Surface = roughPhysicalProp(Surface);
Surface = roughnessStats(Surface);
end

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
        varNames{c} = name;
        c = c + 1;
    end
end
SurfStruct.varNames = varNames;
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

function SurfStruct = roughPhysicalProp(SurfStruct)
xname = SurfStruct.varNames{1};
x = SurfStruct.obj.(xname);
switch SurfStruct.type
    case '1D-profile'
        Lx = max(x) - min(x);
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
% height information is the last variable.
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

function ind = findSlopeCorr(C)
l = length(C);
C = C((l+1)/2:end); % Get only half of C
ind = find(diff(C) > 0,1);
end