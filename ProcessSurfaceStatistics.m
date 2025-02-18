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
[filename, pathname] = uigetfile({'*.mat;*.csv;*.asc;*.dat;*.xls;*.xlsx',...
                       'Surface Files (*.mat,*.csv,*.asc,*.dat,*.xls,*.xlsx)';...
                                 '*.*','All Files(*.*)'},...
    'Choose a Matlab data file with x, y, z coordinates.');
if isempty(filename) == 1 || ~ischar(filename)
    error('No file was selected')
end
% Or  paste in the file
% filename = ['Processed_Surface02_8_12_25grit_CURVTILT.mat'];

% Roughness and Scanner Questionare
SurfAnswers = SurfaceRoughnessProcessor.RoughnessQuestionnaire(batch);
ScannerAnswers = SurfaceRoughnessProcessor.ScannerQuestionnaire(SurfAnswers, batch);
        
% Run function to calculate statistics
Surface = SurfaceRoughnessProcessor.getSurfStatistics(fullfile(pathname, filename));

% Export Surface Statistics
SurfaceRoughnessProcessor.exportSurfaceStatistics(Surface, SurfAnswers, ScannerAnswers)

% Visualize Surface
SurfaceRoughnessProcessor.visualizeResults(Surface)

% Display results on command prompt
SurfaceRoughnessProcessor.displayResults(Surface)