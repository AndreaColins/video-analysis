function data = dat2mp4(loadfile,savefile)
% dat2mp4.m
% 2015 version of batch_dat2avi.m
% M.Evans 03.2015. Original file by rsp 080713
% Function to convert mikotron .dat files into .mp4 format using
% videowriter
% A new script was needed (over dat2avi) as some old functions (avifile)
% have been superseded by VideoWriter.m in newer Matlab versions.
% videoWriter.m also handles files larger than 2GB (which aviFile.m
% couldn't)
% This function is called by batch_dat2mp4.m

% load mikrotron .dat files directly into matlab
% rsp 020113

fid = fopen([loadfile '.dat'],'r');

% Load in header information
data = read_mikrotron_datfile_header(fid);

data.loadfile = loadfile;
data.savefile = savefile;

% data.profile = 'MPEG-4'; % This is tiny. Good for talks. Cannot be read by Janelia Whisker Tracker
% data.profile = 'Uncompressed AVI'; % This is huge. DO NOT USE
data.profile = 'Motion JPEG AVI'; % This is the default


%%
% Create mov file to write individual images (img)
fseek(fid,data.offset,-1);



mp4obj = VideoWriter(savefile,data.profile);
% aviobj = avifile(savefile,'compression','none','fps',round(data.framerate/20));

mp4obj.FrameRate = data.framerate/20;

open(mp4obj)

count = 0;
stop = 0;
fprintf('Processing %d frames:\n',data.nframes)
while(~stop)
    if rem(count,100)==0
        fprintf('%d ', count)
    end
    count = count+1;
    imagecounter = fread(fid,2,'uint8');
    junk = fread(fid,2,'uint8');
    tickcounter = fread(fid,2,'uint32');
    digitalinputs = fread(fid,1,'uint8');
    junk = fread(fid,11,'uint8');
    frame = fread(fid,data.imagesize-24,'uint8=>uint8');
    stop = feof(fid);
    if ~stop
        digitalbit0(count) = bitand(digitalinputs,1);
                % disp(fprintf('Doing frame %d:(%d,%d)\n',count,tickcounter))
        frame_res = reshape([zeros(24,1); frame],data.width,data.height)'; % AE reshapes frame in the same way that imagesc is displayed
                 %imagesc(frame_res);
        img = repmat(frame_res, [1,1,3]); % AE duplicates the 2-D matrix into 3-D, colour values stay the same, just a formating thing
        %         writeVideo(mov,img); % AE
        %         aviobj = addframe(aviobj,img);
         
         if count>=data.startframe
             %imagesc(img);
        %pause(0.001)
        writeVideo(mp4obj,img);
         end
    end
    
    
    %title(ftell(fid))
    %      pause
    
end

fprintf('done\n')
close(mp4obj)
% close(mov); % AE necessary??
% aviobj = close(aviobj);
% mov file is now in format that can be passed to 'traceWhiskers.m'

% data.digitalbit0 = digitalbit0; % trigger input

return