%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program calculate the longitude and latitude of all grid cells 
% (of 1/16th degree) in Livneh forcing set within the defined boundary. 

clear all;clc;

cd('WRF_2014_Morr_NARR');


files = dir('dhsvm_all_var_input_NARR_Morr_*');
for inFile = files'
    [year, lat, lon] = strread(inFile.name, '%*s %*s %*s %*s %*s %*s %d %f %f','delimiter','_');
    %make a copy of file
    newfilename = strcat('dhsvm_all_var_input_NARR_Morr_',num2str(year),'_',num2str(lat,'%.5f'),'_',num2str(lon,'%.5f'));
    movefile(inFile.name, newfilename);
end



