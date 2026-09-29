reference = imread('E:\20240514_Dongju_Pattern\ToAnalyze\dws_4channel_6well_pattern_20240514140824\Rearranged_Files\Well__B_003\r_0022_c_0029\RED_BEADS_ref\RED_BEADS__B_003_r_0022_c_0029_t_00000000_z_0000.tif');
movie = imread('E:\20240514_Dongju_Pattern\ToAnalyze\dws_4channel_6well_pattern_20240514140824\Rearranged_Files\Well__B_003\r_0022_c_0029\RED_BEADS\RED_BEADS__B_003_r_0022_c_0029_t_00000001_z_0000.tif');

imgpath = 'E:\20240514_Dongju_Pattern\ToAnalyze\dws_4channel_6well_pattern_20240514140824\Rearranged_Files\Well__B_003\r_0022_c_0029\RED_BEADS\';
ext = 'tif'
tic
[IMAGES, IMAGE_NAMES] = imageLoader(imgpath,ext);
toc
reference = double(reference);

reference = (reference-min(reference(:)))/(max(reference(:)-min(reference(:)))) * 65535;
reference = uint16(reference);
%consider replacing with ORB features for open source
ptsRreference = detectSURFFeatures(reference);
[featuresReference,validPtsReference] = extractFeatures(reference,ptsRreference,"Method","SURF");

tic
for img = 1:size(IMAGES,3)
    %disp(['number ',num2str(img)])
    movie = double(IMAGES(:,:,img));
    movie = (movie-min(movie(:)))/(max(movie(:)-min(movie(:)))) * 65535;
    

    movie = uint16(movie);
    
    % f1 = figure;
    % %imshow(reference);
    % imshowpair(reference,movie,"falsecolor","Scaling","joint")
    
    %ptsRreference = detectSURFFeatures(reference);
    %ptsMovie = detectSURFFeatures(movie);
    
    ptsMovie = detectSURFFeatures(movie);
    
    [featuresMovie,validPtsMovie] = extractFeatures(movie,ptsMovie,"Method","SURF");

    indexPairs = matchFeatures(featuresReference,featuresMovie);
    
    matchedReference = validPtsReference(indexPairs(:,1));
    matchedMovie = validPtsMovie(indexPairs(:,2));
    
    % f2 = figure
    % showMatchedFeatures(reference,movie,matchedReference,matchedMovie);
    % title('Putatively matched points (including outliers)');
    
    [tform, inlierIdx] = estgeotform2d(matchedMovie,matchedReference,'rigid');
    inlierMovie = matchedMovie(inlierIdx,:);
    inlierReference = matchedReference(inlierIdx,:);
    % 
    % f3 = figure;
    % showMatchedFeatures(reference,movie,inlierReference,inlierMovie);
    % title('Matching points (inliers only)');
    % legend('ptsReference','ptsMovie');
    
    %convert rotation
    invTform = invert(tform);
    Ainv = invTform.A;
    
    ss = Ainv(1,2);
    sc = Ainv(1,1);
    scaleRecovered = hypot(ss,sc);
    %disp(['Recovered translation: ',num2str(tform.Translation)])
    %disp(['Recovered scale: ', num2str(scaleRecovered)])
    
    % Recover the rotation in which a positive value represents a rotation in
    % the clockwise direction.
    thetaRecovered = atan2d(-ss,sc);
    %disp(['Recovered theta: ', num2str(thetaRecovered)])
    
    %tform.RotationAngle=0;
    
    outputView = imref2d(size(reference));
    recovered = imwarp(movie,tform,OutputView=outputView);
    
    % f4 = figure
    % imshowpair(reference,recovered,"falsecolor","Scaling","joint")
end

toc
