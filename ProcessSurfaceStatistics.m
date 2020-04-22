clc
clear

filename = 'Processed_GaussianPSD_L=11xkrms_krms=100um_Sk=1_CURVTILT.mat';
Surface = getSurfProperties(filename);

filename = 'Processed_Krms350skew+1tile1tiltcurve.mat';
Surface2 = getSurfProperties(filename);

filename = ['Generated_Surface_gaussian_p1=3.85_p2=0.35_rms=350um_'...
        'Sk=0.97498_ES=0.40465_279.4mm x 190.5mm_100umRes_26-Feb-2018.mat'];
Surface3 = getSurfProperties(filename);

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
zname = SurfStruct.varNames{end}; % A bit of a strech assuming the
% height information is the last variable.
z = SurfStruct.obj.(zname);
kmin = min(z(:));
kmax = max(z(:));
SurfStruct.kp = kmax - kmin;
%kp5 = 
SurfStruct.km = mean(z(:));
SurfStruct.ka = mean(abs(z(:)));
SurfStruct.krms = rms(z(:));
SurfStruct.kstd = std(z(:));
SurfStruct.kSk = skewness(z(:));
SurfStruct.kKu = kurtosis(z(:));

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