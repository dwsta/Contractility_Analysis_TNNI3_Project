function [U,V] = filter_PIV(Uraw,Vraw)

% Fill up nans
U = fill_nans(Uraw);
V = fill_nans(Vraw);
bad_vecs = false(size(U));

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
U(bad_vecs) = nan;
U = fill_nans(U);
V(bad_vecs) = nan;
V = fill_nans(V);

end
