% Script to run PROCESSLBATCHMODE.m on all the strain data,
% collecting residuals and tau values for both the 5-state model and
% the 3-state model


AKR_directory = '\\FS1\WisorData\Rempe\Data\strain_study_data\AKR\long_files\';
BA_directory  = '\\FS1\WisorData\Rempe\Data\strain_study_data\BA\Long_Files\Autoscore_output_09.08.2015.13.51\';
BL_directory  = '\\FS1\WisorData\Rempe\Data\strain_study_data\BL\Long_Files\Autoscore_output_09.08.2015.14.07\';
CDJ_directory = '\\FS1\WisorData\Rempe\Data\strain_study_data\CDJ\';
DBA_directory = '\\FS1\WisorData\Rempe\Data\strain_study_data\DBA\long_files\';
OLA_directory = '\\FS1\WisorData\Rempe\Data\strain_study_data\OLA\OLA.txt\';


%AKR (3-state model and 5-state model)
[signal_data_AKR3,state_data_AKR3,timestampvec_AKR3,residual_AKR3,best_S_AKR3,UppA_AKR3,LowA_AKR3,dynamic_range_AKR3,Timer_AKR3,Taui_AKR3,TauD_AKR3]=PROCESSLBATCHMODE_NOGUI(AKR_directory,'lactate','3state','17:00:00',1,'END',1);
%);
[signal_data_AKR5,state_data_AKR5,timestampvec_AKR5,residual_AKR5,best_S_AKR5,UppA_AKR5,LowA_AKR5,dynamic_range_AKR5,Timer_AKR5,Taui_AKR5,TauD_AKR5]=PROCESSLBATCHMODE_NOGUI(AKR_directory,'lactate','5state','17:00:00',1,'END',1);

save strain_study_lactate_tau_residuals.mat

%BA (3-state model and 5-state model)
[signal_data_BA3,state_data_BA3,timestampvec_BA3,residual_BA3,best_S_BA3,UppA_BA3,LowA_BA3,dynamic_range_BA3,Timer_BA3,Taui_BA3,TauD_BA3]=PROCESSLBATCHMODE_NOGUI(BA_directory,'lactate','3state','17:00:00',1,'END',1);
[signal_data_BA5,state_data_BA5,timestampvec_BA5,residual_BA5,best_S_BA5,UppA_BA5,LowA_BA5,dynamic_range_BA5,Timer_BA5,Taui_BA5,TauD_BA5]=PROCESSLBATCHMODE_NOGUI(BA_directory,'lactate','5state','17:00:00',1,'END',1);

save strain_study_lactate_tau_residuals.mat

%BL (3-state model and 5-state model)
[signal_data_BL3,state_data_BL3,timestampvec_BL3,residual_BL3,best_S_BL3,UppA_BL3,LowA_BL3,dynamic_range_BL3,Timer_BL3,Taui_BL3,TauD_BL3]=PROCESSLBATCHMODE_NOGUI(BL_directory,'lactate','3state','17:00:00',1,'END',1);
[signal_data_BL5,state_data_BL5,timestampvec_BL5,residual_BL5,best_S_BL5,UppA_BL5,LowA_BL5,dynamic_range_BL5,Timer_BL5,Taui_BL5,TauD_BL5]=PROCESSLBATCHMODE_NOGUI(BL_directory,'lactate','5state','17:00:00',1,'END',1);

save strain_study_lactate_tau_residuals.mat

%CDJ (3-state model and 5-state model)
[signal_data_CDJ3,state_data_CDJ3,timestampvec_CDJ3,residual_CDJ3,best_S_CDJ3,UppA_CDJ3,LowA_CDJ3,dynamic_range_CDJ3,Timer_CDJ3,Taui_CDJ3,TauD_CDJ3]=PROCESSLBATCHMODE_NOGUI(CDJ_directory,'lactate','3state','17:00:00',1,'END',1);
[signal_data_CDJ5,state_data_CDJ5,timestampvec_CDJ5,residual_CDJ5,best_S_CDJ5,UppA_CDJ5,LowA_CDJ5,dynamic_range_CDJ5,Timer_CDJ5,Taui_CDJ5,TauD_CDJ5]=PROCESSLBATCHMODE_NOGUI(CDJ_directory,'lactate','5state','17:00:00',1,'END',1);

save strain_study_lactate_tau_residuals.mat

%DBA (3-state model and 5-state model)
[signal_data_DBA3,state_data_DBA3,timestampvec_DBA3,residual_DBA3,best_S_DBA3,UppA_DBA3,LowA_DBA3,dynamic_range_DBA3,Timer_DBA3,Taui_DBA3,TauD_DBA3]=PROCESSLBATCHMODE_NOGUI(DBA_directory,'lactate','3state','17:00:00',1,'END',1);
[signal_data_DBA5,state_data_DBA5,timestampvec_DBA5,residual_DBA5,best_S_DBA5,UppA_DBA5,LowA_DBA5,dynamic_range_DBA5,Timer_DBA5,Taui_DBA5,TauD_DBA5]=PROCESSLBATCHMODE_NOGUI(DBA_directory,'lactate','5state','17:00:00',1,'END',1);

save strain_study_lactate_tau_residuals.mat

%OLA (3-state model and 5-state model)
[signal_data_OLA3,state_data_OLA3,timestampvec_OLA3,residual_OLA3,best_S_OLA3,UppA_OLA3,LowA_OLA3,dynamic_range_OLA3,Timer_OLA3,Taui_OLA3,TauD_OLA3]=PROCESSLBATCHMODE_NOGUI(OLA_directory,'lactate','3state','17:00:00',1,'END',1);
[signal_data_OLA5,state_data_OLA5,timestampvec_OLA5,residual_OLA5,best_S_OLA5,UppA_OLA5,LowA_OLA5,dynamic_range_OLA5,Timer_OLA5,Taui_OLA5,TauD_OLA5]=PROCESSLBATCHMODE_NOGUI(OLA_directory,'lactate','5state','17:00:00',1,'END',1);

save strain_study_lactate_tau_residuals_starting_at17_00.mat


% you may want to do all the excel-file writing right here in the script (or in a separate post-processing sc)
% to get a list of all the .txt files you can do dir

