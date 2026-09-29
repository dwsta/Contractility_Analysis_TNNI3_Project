function [X,Y,T,U,V,Xdrift,Ydrift] = readPIV_bin(filename)
fid = fopen(filename);
ver = fread(fid,2,'uchar')';
str = char(ver);
switch str
    case 'v3' % Now including x-y drift (will be 0 if drift correction was off)
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
        U = reshape(fread(fid,N,'single'),ny,nx,nf,nr);
        V = reshape(fread(fid,N,'single'),ny,nx,nf,nr);
    
    case 'v2' % X,Y,T are now vectors instead of dense matrices to save memory
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
        U = reshape(fread(fid,N,'single'),ny,nx,nf,nr);
        V = reshape(fread(fid,N,'single'),ny,nx,nf,nr);
        Xdrift = 0*T;
        Ydrift = 0*T;
%         [X,Y,T,~] = meshgrid(xvec,yvec,tvec,repvec);
    otherwise
        fclose(fid);
        fid = fopen(filename);
        A = fread(fid,4,'single');
        nx = A(1);
        ny = A(2);
        nf = A(3);
        nr = A(4);
        N = nx*ny*nf*nr;
        X = reshape(fread(fid,N,'single'),ny,nx,nf,nr);
        Y = reshape(fread(fid,N,'single'),ny,nx,nf,nr);
        T = reshape(fread(fid,N,'single'),ny,nx,nf,nr);
        U = reshape(fread(fid,N,'single'),ny,nx,nf,nr);
        V = reshape(fread(fid,N,'single'),ny,nx,nf,nr);
        Xdrift = 0*T;
        Ydrift = 0*T;
end
fclose(fid);