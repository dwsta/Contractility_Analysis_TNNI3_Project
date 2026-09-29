function DIV = divergence_rik(X,Y,U,V)
numwindowsx = size(U,2);
numwindowsy = size(U,1);
N_FRAMES = size(U,3);
XU = unique(X(:));
YU = unique(Y(:));
dx = XU(2)-XU(1);
dy = YU(2)-YU(1);
interior_pt = reshape([1:numwindowsx*numwindowsy],numwindowsy,numwindowsx);
interior_pt = interior_pt(2:end-1,2:end-1);
interior_pt = interior_pt(:);
U = reshape(U,numwindowsx*numwindowsy,N_FRAMES);
V = reshape(V,numwindowsx*numwindowsy,N_FRAMES);
DIV = nan(size(U));
DIV(interior_pt,:) = (U(interior_pt+numwindowsy,:) - U(interior_pt-numwindowsy,:))/(2*dx) + ...
    (V(interior_pt+1,:) - V(interior_pt-1,:))/(2*dy);
DIV = reshape(DIV,numwindowsy,numwindowsx,N_FRAMES);
end
