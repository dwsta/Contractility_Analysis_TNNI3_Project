cycles = 100;
xvec = 1:63; yvec = 1:63; tvec = 1:300;
[X,Y,T] = meshgrid(xvec,yvec,tvec);
p = 100;
U = sin(2*pi*T/p);
V = cos(2*pi*T/p);

scl = 1;
fig = figure;
hquiver = quiver(X(:,:,1),Y(:,:,1),scl*U(:,:,1),scl*V(:,:,1),0);

iframe = 10;

available_frames = [1:3:300];
nFrames = length(available_frames);
U = U(:,:,available_frames);
V = V(:,:,available_frames);


t0 = tic;
frameRate = 30; % target frame rate
for icycles = 1 : cycles*nFrames
    telapsed = toc(t0);
    icycles
    theo_frame = round(1 + mod(telapsed*frameRate,nFrames))
    [~,idx] = min(abs(theo_frame-available_frames));
    iframe = available_frames(idx)
    u = U(:,:,iframe);
    v = V(:,:,iframe);

    set(hquiver,'UData',scl*u);
    set(hquiver,'VData',scl*v);
    drawnow
end
