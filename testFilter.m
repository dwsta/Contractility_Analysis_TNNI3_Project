% maindir = 'G:\Shared drives\Stanford UCSD_sharedDrive\contractility_analysis_rel_4_0\TaxolFullMovies\contractility_run_20220711_121055\output\Pretreatment_Middle_Well__E_007_Tritc';
% binfile = 'deformations_pass1.bin';

maindir = 'G:\Shared drives\Stanford UCSD_sharedDrive\contractility_analysis_rel_4_0\TaxolFullMovies\contractility_run_20220712_164803_multipass\output\Pretreatment_Middle_Well__E_007_Tritc';
binfile = 'deformations_pass4.bin';

[xvec,yvec,tvec,Uraw,Vraw] = readPIV_bin(fullfile(maindir,'piv',binfile));
[X,Y,T] = meshgrid(xvec,yvec,tvec);

% Fill up nans
U = fill_nans(Uraw);
V = fill_nans(Vraw);

h = ones(5,1);
H = h*h';
H(3,3) = 0;
H = H/sum(H,'all');

[Umed,Vmed,Diff,NeighborDiff]=deal(nan(size(U)));
bad_vecs = false(size(U));
tic
for iii = 1 : size(U,3)
    %     Umed (:,:,iii) = medfilt2(U(:,:,iii),[5 5]);
    %     Vmed (:,:,iii) = medfilt2(V(:,:,iii),[5 5]);
    %     Diff(:,:,iii)  = sqrt( (U(:,:,iii)-Umed(:,:,iii)).^2 + (V(:,:,iii)-Vmed(:,:,iii)).^2);
    %     NeighborDiff(:,:,iii) = medfilt2(Diff(:,:,iii),[5 5]);
    Umed = medfilt2(U(:,:,iii),[5 5]);
    Vmed = medfilt2(V(:,:,iii),[5 5]);
    Diff = sqrt( (U(:,:,iii)-Umed).^2 + (V(:,:,iii)-Vmed).^2);
    NeighborDiff = medfilt2(Diff,[5 5]);
    bad_vecs(:,:,iii) = Diff./(NeighborDiff+0.1)>1.5;
end
toc
% tic
% Urec=griddata(X(~bad_vecs),Y(~bad_vecs),T(~bad_vecs),Uraw(~bad_vecs),X,Y,T);
% Vrec=griddata(X(~bad_vecs),Y(~bad_vecs),T(~bad_vecs),Vraw(~bad_vecs),X,Y,T);
% toc

tic
Uhole = Uraw;
Uhole(bad_vecs) = nan;
Urec2 = fill_nans(Uhole);
Vhole = Vraw;
Vhole(bad_vecs) = nan;
Vrec2 = fill_nans(Vhole);
toc

% figure(1);
% clf
% quiver(X(:,:,22),Y(:,:,22),Uraw(:,:,22)*scl,Vraw(:,:,22)*scl,0,'b');
% hold on
% quiver(X(:,:,22),Y(:,:,22),Uraw(:,:,22)*scl.*bad_vecs(:,:,22),Vraw(:,:,22).*bad_vecs(:,:,22)*scl,0,'r');
% 
% figure(2);
% clf
% quiver(X(:,:,22),Y(:,:,22),Urec(:,:,22)*scl,Vrec(:,:,22)*scl,0,'k');


% Ustd = stdfilt(U,ones(5));
% Vstd = stdfilt(V,ones(5));

% Umed = nanmedian_3Dneighbors(Uraw);
% Vmed = nanmedian_3Dneighbors(Vraw);
% mag = sqrt((U-Umed).^2 + (V-Vmed).^2);
% rmed = nanmedian_3Dneighbors(mag,[3 3 3]);
% eps = 0.1;
% % thresh =
% bad_vec = (mag./(rmed+eps))