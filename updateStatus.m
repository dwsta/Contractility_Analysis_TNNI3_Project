function updateStatus(app,Row,newStatus)
app.UITableStatus.Data(Row,:).Status = newStatus;
drawnow
figure(app.UIFigure)
end