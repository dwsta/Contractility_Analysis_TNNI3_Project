function singleCellSegmentation(imdir,cfg_data)
% If an image needs some sort of processing before cellpose, it will go
% here

% Call cellpose
if cfg_data.SingleCell.UseSingleCell % in case we add some other way to segment
% str = 'G:\Shared drives\Stanford UCSD_sharedDrive\contractility_analysis_rel_4_0\2021_12_13_Pretreatment_Middle_Well__C_008_Tritc\test';
    if cfg_data.SingleCell.UseCellpose
        callCellpose(imdir,cfg_data)
    end
end