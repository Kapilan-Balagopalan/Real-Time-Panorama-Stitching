% url = 'http://192.168.8.106:8080/photo.jpg';
% ss  = imread(url);
% 
% fh = image(ss);
% while(1)
%     ss  = imread(url);
%     set(fh,'CData',ss);
%     drawnow;
% end
     url = 'http://192.168.8.106:8080/video';
%     ss  = implay(url);

vid = videoinput(url);
src = getselectedsource(vid);


 
 vid.FramesPerTrigger = 100;
preview(vid);