function [X,Y,T,SXX,SYY,SXY,Xdrift,Ydrift] = readMSM_bin(filename)
fid = fopen(filename);
ver = fread(fid,2,'uchar')';
str = char(ver);
switch str
    case 'v3'
        A = fread(fid,4,'single');
        nx = A(1);
        ny = A(2);
        nf = A(3);
        nr = A(4);
        N = nx*ny*nf*nr;
        X = fread(fid,nx,'single');
        Y = fread(fid,ny,'single');
        T = fread(fid,nf,'single');
        repvec = fread(fid,nr,'single');
        Xdrift = fread(fid,nf,'single');
        Ydrift = fread(fid,nf,'single');
        SXX = reshape(fread(fid,N,'single'),ny,nx,nf,nr);
        SYY = reshape(fread(fid,N,'single'),ny,nx,nf,nr);
        SXY = reshape(fread(fid,N,'single'),ny,nx,nf,nr);        
    case 'v2'
        A = fread(fid,4,'single');
        nx = A(1);
        ny = A(2);
        nf = A(3);
        nr = A(4);
        N = nx*ny*nf*nr;
        X = fread(fid,nx,'single');
        Y = fread(fid,ny,'single');
        T = fread(fid,nf,'single');
        repvec = fread(fid,nr,'single');
        Xdrift = 0*T;
        Ydrift = 0*T;
        SXX = reshape(fread(fid,N,'single'),ny,nx,nf,nr);
        SYY = reshape(fread(fid,N,'single'),ny,nx,nf,nr);
        SXY = reshape(fread(fid,N,'single'),ny,nx,nf,nr);        
end
fclose(fid);