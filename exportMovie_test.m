
%% Load data
MainDir = 'G:\Shared drives\Stanford UCSD_sharedDrive\contractility_analysis_rel_4_0\TaxolFullMovies\contractility_run_20220712_164803_multipass\output\Pretreatment_Middle_Well__E_007_Tritc';
RawVideoLocation = 'G:\Shared drives\Stanford UCSD_sharedDrive\contractility_analysis_rel_4_0\TaxolFullMovies\Pretreatment_Middle\Well__E_007\Tritc';
PivDir = [MainDir,'\piv'];
PivFile = 'deformations_pass4.bin';
fname = fullfile(PivDir,PivFile);
[xdata,ydata,tdata,qdatax,qdatay] = readPIV_bin(fname);
quiverdata = struct;
quiverdata.xvec = xdata;
quiverdata.yvec = ydata;
quiverdata.tvec = tdata;
quiverdata.U = qdatax;
quiverdata.V = qdatay;
quiverdata.scale = 1;
quiverdata.color = [1 0 1];

cdata = struct;
cdata.plottype = 'video';
Movie = imageLoader(RawVideoLocation,'tif','noparallel');
cdata.XData = [1 size(Movie,2)];
cdata.YData = [1 size(Movie,1)];
cdata.tvec  = [1:size(Movie,3)];
cdata.Data = Movie;
cdata.CLim = [100 1000]; 
cdata.colormap = gray;
                
im = imread(fullfile(MainDir,'mask','cp_input_Tritc__E_007_r_0004_c_0004_t_00000000_z_0000_cp_masks.png'));

mask = boundarymask(im);
hboundaries = imagesc(mask);
hboundaries.AlphaData = (mask==1);
cellboundaries = struct();
cellboundaries.CData = hboundaries.CData;
cellboundaries.AlphaData = hboundaries.AlphaData;
cellboundaries.colormap = [0 0 0; 1 0 0];
cellboundaries.CLim = [0 1];

outputdir = MainDir;
outputname = 'video';
outputres = [256 256];
framerate = 30;
xl = [1 size(Movie,2)]; % In case it is zoomed in
yl = [1 size(Movie,1)];

exportMovie(quiverdata,cdata,cellboundaries,outputres,framerate,outputdir,outputname,xl,yl)
