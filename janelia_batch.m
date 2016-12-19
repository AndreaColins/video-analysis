% janelia_batch.m
%
% Minimal script for converting files to .avi then running the janelia whisker tracker
%
% Use chimera.m to clean up the data, then extract kappa/angle variables using a mask, quadratic fit and whikerman functions.
%
% Also has instructions for running the tracker at the command line
%
% First videos need to be converted into .avi format with dat2mp4.m
% [~] = dat2mp4(FILENAME,FILENAME);
% Or cd into the video directory and run batch_dat2mp4.m (editted to
% include the appropriate session names.
%
% Second, run (in turn):
% trace FILENAME.avi FILENAME.whiskers
% measure --face bottom FILENAME.whiskers FILENAME.measurements
% classify FILENAME.measurements FILENAME.measurements bottom --px2mm 0.057 -n 1
% reclassify FILENAME.measurements FILENAME.measurements -n 0

%% List of directories to convert/track
clear
close all
main_dir = 'B:\Dario\Behavioral_movies';
% movie_dirs = {};
% Some movies are in Dario's folder, so will need to be tracked from his PC
% or moved.

cd(main_dir)

% Animal 32
dirs_32 = dir(['32\']);
dirs_32 = dirs_32(3:end);
for i = 1:numel(dirs_32);
    subdir = dir(['32\',dirs_32(i).name,'\']);
    if numel(subdir) == 3;
        subdir = subdir(3).name;
        fulldir{i} = ['32\',dirs_32(i).name,'\',subdir,'\'];
    end
end
num_sofar = numel(fulldir);

% Animal 33
dirs_33 = dir(['33\']);
dirs_33 = dirs_33(3:end);
for i = 1:numel(dirs_33);
    subdir = dir(['33\',dirs_33(i).name,'\']);
    if numel(subdir) == 3;
        subdir = subdir(3).name;
        fulldir{num_sofar+i} = ['33\',dirs_33(i).name,'\',subdir,'\'];
    end
end
num_sofar = numel(fulldir);

% Animal 34
dirs_34 = dir(['34\']);
dirs_34 = dirs_34(3:end);
for i = 1:numel(dirs_34);
    subdir = dir(['34\',dirs_34(i).name,'\']);
    if numel(subdir) == 3;
        subdir = subdir(3).name;
        fulldir{num_sofar+i} = ['34\',dirs_34(i).name,'\',subdir,'\'];
    end
end

% Added manually due to folder structure irregularity
fulldir{184} = '33\33\030215_33a\2015_02_03\';

num_sofar = numel(fulldir);

main_dir = 'I:\Dario Campagner\BEHAVIORAL MOVIES';
cd(main_dir)
% Animal 34
dirs_34 = dir(['34\']);
dirs_34 = dirs_34(3:end);
for i = 1:numel(dirs_34);
    subdir = dir(['34\',dirs_34(i).name,'\']);
    if numel(subdir) == 3;
        subdir = subdir(3).name;
        fulldir{num_sofar+i} = ['34\',dirs_34(i).name,'\',subdir,'\'];
    end
end

num_sofar = numel(fulldir);
% Animal 32
dirs_32 = dir(['32\']);
dirs_32 = dirs_32(3:end);
for i = 1:numel(dirs_32);
    subdir = dir(['32\',dirs_32(i).name,'\']);
    if numel(subdir) == 3;
        subdir = subdir(3).name;
        fulldir{num_sofar+i} = ['32\',dirs_32(i).name,'\',subdir,'\'];
    end
end

num_sofar = numel(fulldir);
% Animal 33
dirs_33 = dir(['33\']);
dirs_33 = dirs_33(3:end);
for i = 1:numel(dirs_33);
    subdir = dir(['33\',dirs_33(i).name,'\']);
    if numel(subdir) == 3;
        subdir = subdir(3).name;
        fulldir{num_sofar+i} = ['33\',dirs_33(i).name,'\',subdir,'\'];
    end
end

num_sofar = numel(fulldir);
% Animal 38
dirs_38 = dir(['38\']);
dirs_38 = dirs_38(3:end);
for i = 1:numel(dirs_38);
    subdir = dir(['38\',dirs_38(i).name,'\']);
    if numel(subdir) == 3;
        subdir = subdir(3).name;
        fulldir{num_sofar+i} = ['38\',dirs_38(i).name,'\',subdir,'\'];
    end
end

% main_dir = 'I:\Dario Campagner\BEHAVIORAL MOVIES';

% main_dir = 'B:\Dario\Behavioral_movies';

% Single whisker sessions only
keepers = [1,5,7,16,27,30,32,35,37,40,42,44,45,50,51,52,54,56,57,58, ...
    60,61,62,64,67,72,75,79,82,86,89,69,73,76,77,...
    78,80,81,83,84,85,87,90,93,107,108,109,110,111,...
    114,117,184,141,143,147,150,154,122,135,138,145,...
    148,151,155,160,174,176,178,125,128,129,131,132,136,140,142,146,149,156,152,163,166,169,171,172,180,...
    186 187 188 189 190 191 192 193 194 195 196 197 198 215 219 220 221 223 224 229 231 232 234,...
   238 240 199 201 203 205 207 210];
dirs_1_whisker = {};
for i = 1:numel(keepers);
    dirs_1_whisker{i} = fulldir{keepers(i)};
end

%% Main conversion loop
% parfor
for i = 25%:numel(dirs_1_whisker); % For all chosen sessions (above)
    
    %% CD to a video directory, then convert all .dat files to .avi
    
    cd(main_dir)
    dir_name = dirs_1_whisker{i};
    cd(dir_name)
    
    files = dir('*.dat');
    names = files(1);
    names = names.name;
    session = names(1:10);
    vid = names(12:19);
    
    for fidx = 315:length(files)
    loadfile = files(fidx).name(1:end-4);
    savefile = ['C:\Users\Dario\Desktop\BEHAVIORAL_DATA_ANALYSIS\janelia tracker\janelia_tracker_scratch',files(fidx).name(1:end-4) '.avi'];
    savefile_nas = [files(fidx).name(1:end-4) '.avi'];
    
    disp(sprintf('Loading %s...',loadfile))
    [~] = dat2mp4(loadfile,savefile);
    disp(sprintf('Saved %s',loadfile))
    
    copyfile(savefile,savefile_nas);
    delete(savefile);
    
    end
    
    
    %% Now run Janelia tracker on the videos
    
    files = dir([session '_' vid '_*.avi']);
    
    disp(['Tracking session ',session])
    
    for fidx = 315:length(files)
        loadfile = files(fidx).name(1:end-4);
        savefile = [loadfile,'.avi'];
        
        disp(['Tracing vid ',loadfile])
        str = ['trace ',savefile,' ',loadfile,'.whiskers'];
        system(str);
        
        disp(['Measuring vid ',loadfile])
        str = ['measure --face bottom ',loadfile,'.whiskers ',loadfile,'.measurements'];
        system(str);
        
        disp(['Classifying vid ',loadfile])
        str = ['classify ',loadfile,'.measurements ',loadfile,'.measurements bottom --px2mm 0.057 -n 1'];
        system(str);
        
        disp(['Re-classifying vid ',loadfile])
        str = ['reclassify ',loadfile,'.measurements ',loadfile,'.measurements -n 1'];
        system(str);
        
    end
    
%     
%     for fidx =1:length(files)% 1:length(files)
%         loadfile = files(fidx).name(1:end-4);
%         [best_m,best_w,kappa_w,theta_w,r_base,dropped_frames] = chimera(loadfile);
        chimera
%         save([loadfile,'_clean.mat'],'kappa_w','theta_w','r_base','dropped_frames','-append');
%     end
    % Then run chimera.m (currently doing this separately, eventually it
    % would happen here, such that the pipeline is run end-to-end on a
    % given session using only one script
    
end