function [Xdrift,Ydrift] = readDrift_bin(drift_file)
fid = fopen(drift_file);
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
end
fclose(fid);
end