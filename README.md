# Surface Roughness Statistics Code for Roughness Database
The MATLAB script processes surface roughness statistics to support the roughness database hosted by the University of Southampton.

The development idea of this script is to require minimum user interface when running the code. In other words, there is no need for code manipulation or changes to correctly export the roughness statistics.

The script will generate the correct folder structure required by the database.
The structures should look like this:
> 
**Het_Irreg_TBL_turbine-blade_Barros_2014**
> 

## Requirements
A working version of MATLAB. Should work on both Macs and PCs.

Additionally, some basic information about the roughness is advisable to have at hand before running the script. 
These are:
- If roughness is _Homogeneous_ or _Heterogeneous_;
- If roughness is _Regular_ or _Irregular_;
- Are results from a _TBL_, _Pipe_, _Channel_;
- Are results from _Experiments_, _Simulations_;
- General descriptor of this surface, i.e. 'Sandgrain';
- Last name of the lead author of the study;
- Year when the results were published;
- Identifying name of this surface, i.e. '220Grit';
- DOI of the publication.

Additionally, if the results are from _Experiments_, the profiler information is required:
- Name and model of the profiler/scanner;
- Uncertainty in the measurement of surface heights in microns.

## Getting Started
First, visit [Southampton's Roughness Database](http://roughnessdatabase.org) and request user access to upload and contribute your roughness to the community.
Once access is granted, download the MATLAB script via GitHub (cloning or download ZIP).

__You can also download the ZIP container from the [Roughness Database](https://sotonac.sharepoint.com/teams/roughnessdatabase), once user access is granted (Possibly, preferred and simplest way).__

## Installation
Simply clone the repository or put the downloaded files into a preferred folder.

## Usage
Simply click the _Run_ button on MATLAB, and all calculations and database requirements will be done automatically.
A basic roughness and surface profiler questionnaire will need to be answered. 

There are two ways to run the scripts: __User-Prompt__ or __Batch__ questionnaire.
The simplest way is the __User-Prompt__. In this mode, dialog boxes will pop up in which the roughness and work information can be filled.
__To run in this mode, simply change the variable in the script to:__
```
batch = false;
```

### Multiple scans for the same surface/work
If this is the case, then using __Batch__ questionnaire mode is recommended.
Make sure that the below variable is set to:
```
batch = true;
```
Two additional files can be found:
- __Profiler_Batch.txt__
- __Questionnaire_Batch.txt__

Follow the instructions found in these files and run the script.

## File Structure Requirements
Two types of surfaces are supported: 1D-line profiles, or 2D-surface scans.

The input formats supported are MATLAB (*.mat), Excel (*.xls), and ASCII (*.csv or tab-delimited *.txt, or *.dat).

_We however ask the users to comply with the standards set for the input files containing the roughness height information._

__All roughness scan information (height and coordinates) must be in millimeters [mm].__

For MATLAB files, put the roughness information into variables X, Y, Z and use the save(...) function to save the roughness information.

For ASCII/CSV files, the X, Y, Z information should be formatted in each column.

The supporting functions are packed as a MATLAB class in a separate file. Changing them may likely break the code and/or calculations.

If a bug is detected, either open an issue on GitHub at 
https://github.com/jmbarrojr/SurfaceRoughnessStatistics/issues or send an email to julio.barros@gmail.com using "[BUG]" prefix in the email's subject.

## Calculated roughness statistics
- Streamwise length of scan (mm)
- Spanwise length of scan (mm)
- Minimum roughness height (mm)
- Maximum roughness height (mm)
- Peak-to-trough roughness height (mm)
- Average 5 peak-to-trough roughness height (mm)
- Average roughness height (mm)
- Average of absolute value of the height fluctuations (mm)
- Root-mean-square of the total height (mm)
- Root-mean-square of the height fluctuations (mm)
- Skewness of the height fluctuations
- Flatness of the height fluctuations
- Effective Slope in the streamwise direction (mm)
- Correlation length in the streamwise direction
- Effective Slope in the spanwise direction
- Correlation length in the spanwise direction (mm)

## Upload Folder Structure to Roughness Database
Once logged in, go to __Useful files__ and click on the __Create entry__ folder.

Simply drag-and-drop the parent folder created by the MATLAB script.
