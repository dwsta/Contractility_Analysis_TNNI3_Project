function callCellpose(imdir,cfg_data)
imdir_pyformat = strrep(imdir,'\','/'); % Need to replace Windows \ by / for python.
imdir_pyformat = ['"',imdir_pyformat,'"']; % and also add " at the beginning and end so that spaces are not an issue

condacall = cfg_data.SingleCell.CellposeEnvActivate; % e.g.'conda activate cellpose'
pythoncall = cfg_data.SingleCell.CellposePythonCall; % e.g 'python -m cellpose'
model = cfg_data.SingleCell.CellposeModel;           % e.g. cyto
cellposeOtherArgs = cfg_data.SingleCell.CellposeOtherArgs; % e.g. --save_png --save_flows --flow_threshold 1.0

cellposeMandatoryArgs = [' --dir ',imdir_pyformat, ' --pretrained_model ',model];
% cellposeargs = ['--dir ',imdir_pyformat, cellposeMandatoryArgs,' ', cellposeOtherArgs];
diary(fullfile(imdir,'cellpose_run.log'))
command = [condacall,' & ',pythoncall,' ',cellposeMandatoryArgs,' ',cellposeOtherArgs, ' & conda deactivate']
system(command)
diary off
