function hh = cquiver_rik(u,v)
ncols = 8;
cols = jet(8);
cax = newplot();

h = matlab.graphics.chart.primitive.Quiver;
mag = sqrt(u.^2+v.^2);
levels = [0,linspace(min(mag(:)),max(mag(:)),ncols)];
[x,y]=meshgrid(1:size(u,1));
for icol = 1:ncols
    sel = mag>levels(icol) & levels(icol+1) > mag;
    set(h,'Parent',cax,'Color_I',cols(icol,:),'XData',x(sel),'YData',y(sel),'Udata',u(sel),'VData',v(sel));
    hold(cax,'on')
end
h.assignSeriesIndex();

% call CreateFcn explicitly immediately following object
% creation from the M point of view.
% this was the last line of @Quiver/Quiver.m
% There is no obvious place to have the CreateFcn
% automatically executed on the MCOS side, so we call it here
h.CreateFcn;



if nargout>0, hh = h; end
