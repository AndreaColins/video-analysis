% batch_dat2mp4 
% 2015 version of batch_dat2avi.m
% M.Evans 03.2015. Original file by rsp 080713
% 
% Process mikrotron dat files in current directory using dat2mp4.m
% Creates an .mp4 file for tracking with the Janelia whisker tracker
%
% Needs to be provided with a 'session' and 'vid' string (usually set by
% janelia_batch.m)
% 
% 
% close all
% clear all

% Process all dat files in the current directly that match [session '_' vid '*.dat']

% session = '050514f';
% vid = '20140505';

% Process all *.dat files in current directory that
files = dir(['Z:\Pole\Video data\2016_08_24\TT3\TT3_20160824_121539','*.dat']);

%disp(sprintf('Processing session %s',session))

%% Main loop
for fidx = 1:length(files)
    loadfile = files(fidx).name(1:end-4);
    savefile = [files(fidx).name(1:end-4) '.avi'];
    
    disp(sprintf('Loading %s...',loadfile))
    [~] = dat2mp4(loadfile,savefile);
    disp(sprintf('Saved %s',loadfile))
    
end
clear fidx files loadfile savefile
