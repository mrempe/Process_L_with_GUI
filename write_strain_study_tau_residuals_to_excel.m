% script to write tau values and residuals from the strain data to an excel spreadsheet  

addpath 'C:\Users\wisorlab\Documents\MATLAB\Brennecke\matlab-pipeline\Matlab\etc\matlab-utils\'

AKR_files = dir(strcat(AKR_directory,'*.txt'));
BA_files  = dir(strcat(BA_directory,'*.txt'));
BL_files  = dir(strcat(BL_directory,'*.txt'));
CDJ_files = dir(strcat(CDJ_directory,'*.txt'));
DBA_files = dir(strcat(DBA_directory,'*.txt'));
OLA_files = dir(strcat(OLA_directory,'*.txt'));


Filenames = [ {AKR_files.name}'; {BA_files.name}'; {BL_files.name}'; {CDJ_files.name}'; {DBA_files.name}'; {OLA_files.name}'];


three_state_Taui = [Taui_AKR3; Taui_BA3; Taui_BL3; Taui_CDJ3; Taui_DBA3; Taui_OLA3];
three_state_Taud = [TauD_AKR3; TauD_BA3; TauD_BL3; TauD_CDJ3; TauD_DBA3; TauD_OLA3];
three_state_residuals = [residual_AKR3'; residual_BA3'; residual_BL3'; residual_CDJ3'; residual_DBA3'; residual_OLA3'];

five_state_Taui = [Taui_AKR5; Taui_BA5; Taui_BL5; Taui_CDJ5; Taui_DBA5; Taui_OLA5];
five_state_Taud = [TauD_AKR5; TauD_BA5; TauD_BL5; TauD_CDJ5; TauD_DBA5; TauD_OLA5];
five_state_residuals = [residual_AKR5'; residual_BA5'; residual_BL5'; residual_CDJ5'; residual_DBA5'; residual_OLA5'];



xl=XL;
%sheet = xl.Sheets.Item(1);
%[numcols,numrows] = xl.sheetSize(sheet);


sheet_output = xl.addSheets({'Output'});



xl.setCells(sheet_output{1},[1,1],{'Files'},'669999');
xl.setCells(sheet_output{1},[2,1],{'3 state Ti'},'669999','true');
xl.setCells(sheet_output{1},[3,1],{'3 state Td'},'669999','true');
xl.setCells(sheet_output{1},[4,1],{'3 state residuals'},'669999','true');
xl.setCells(sheet_output{1},[5,1],{'5 state Ti'},'669999','true');
xl.setCells(sheet_output{1},[6,1],{'5 state Td'},'669999','true');
xl.setCells(sheet_output{1},[7,1],{'5 state residuals'},'669999','true');
xl.setCells(sheet_output{1},[1,2],Filenames,'false','true');
xl.setCells(sheet_output{1},[2,2],three_state_Taui);
xl.setCells(sheet_output{1},[3,2],three_state_Taud);
xl.setCells(sheet_output{1},[4,2],three_state_residuals);
xl.setCells(sheet_output{1},[5,2],five_state_Taui);
xl.setCells(sheet_output{1},[6,2],five_state_Taud);
xl.setCells(sheet_output{1},[7,2],five_state_residuals);



directory = '\\FS1\WisorData\Rempe\Data\strain_study_data\';
xl.saveAs('strain_study_simulation_output_1700starttime.xlsx',directory);
fclose('all');