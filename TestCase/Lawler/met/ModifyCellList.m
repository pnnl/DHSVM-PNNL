%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program calculate the longitude and latitude of all grid cells 
% (of 1/16th degree) in Livneh forcing set within the defined boundary. 

clear all;clc;

cd('WRF_2008_Morr_NARR');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% min(lower left) lat 
minLat = 37.75000;
% min lon
minLon = -119.75000; 

% row and col number of a bounding box that contains all grids
row = 9;
col = 9;
int = 1/16;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% max lat 
maxLat = minLat + (row-1)*int;
% max lon
maxLon = minLon + (col-1)*int;
% bounding box
BoundLat = minLat : 1/16 : maxLat;
BoundLon = minLon: 1/16 : maxLon;

files = dir('dhsvm_all_var_input_NARR_Morr_*');
for inFile = files'
    [year, lat, lon] = strread(inFile.name, '%*s %*s %*s %*s %*s %*s %d %f %f','delimiter','_');
    tmp = abs(BoundLat-lat);
    [idx idx] = min(tmp); %index of closest value
    temp1 = BoundLat(1, idx); %closest value
    
    tmp = abs(BoundLon-lon);
    [idx idx] = min(tmp); %index of closest value
    temp2 = BoundLon(1, idx); %closest value
    
    fprintf('%.5f %.5f\n', temp1, temp2);
    
    %make a copy of file
    newfilename = strcat('dhsvm_all_var_input_NARR_Morr_',num2str(year),'_',num2str(temp1,'%.5f'),'_',num2str(temp2,'%.5f'));
    movefile(inFile.name, newfilename);
end



