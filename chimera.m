% function [best_m,best_w,kappa_w,theta_w,r_base,dropped_frames] = chimera(loadfile);
% chimera.m
%
% Minimal script for taking janelia whisker tracker data and cleaning up the data,
% then extracting kappa/angle variables using a mask, quadratic fit and whikerman functions.
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
%
% Third, load the data into matlab:
% w = LoadWhiskers('FILENAME.whiskers');
% m = LoadMeasurements('FILENAME.measurements');
%
% Fourth, clean up data using code found below (should be a function call)
%
% Fifth, fit whikerman bezier curve to data (should be a function call)


%% CD to a video directory. Could be done already in an outer janelia_batch.m loop
%cd Z:\Dario\Behavioral_movies\32\010415_32a\2015_04_01\
%clear all;


% Work out file names to track
% files = dir('*.avi');
% names = files(1);
% names = names.name;
% session = names(1:10);
% vid = names(12:19);

%plotting = 1; % set to 1 if you want to plot as you go, good for debugging

%% Load measurements and whiskers files for a given loadfile
clear all
files = dir('*.avi');
for fidx = 315:numel(files); %
    plotting = 1;
    tic    
    loadfile = files(fidx).name(1:end-4); % was -4 for .dat or .avi, -3 for .tr
    % loadfile = [session,'_',vid,'_',file];
    
    disp(['Processing trial ',loadfile])
    
    
    
    whiskerfile = [loadfile,'.whiskers'];
    measurementfile = [loadfile,'.measurements'];
    
    w = LoadWhiskers(whiskerfile);
    m = LoadMeasurements(measurementfile);
    
    
    %% Loop over frames and work out best_w, whisker file equvalent to best_m, the best whisker tracked in that frame
    % Checks whiski's own solution first, making sure the answer is of only one
    % whisker per frame, then if that whisker is the right length.
    % If no tracked whiskers are long enough, or frame was dropped all
    % together, set best_w to 0.
    
    disp(['Looping over frames to determine best_m'])
    
    best_w = zeros(1,numel(unique([m(:).fid]))); % Not tracked frames are assigned zeros
    best_m = zeros(1,numel(unique([m(:).fid])));
    
    dropped = 0;
    for i = unique([m(:).fid]);
        candidate_Ids = find([m(:).fid] == i);
        candidates_m = m(candidate_Ids);
        candidates_w = w(candidate_Ids);
        
        lengths = [candidates_m(:).length];
        [ml,mx] = max(lengths);
        
        if ml<=100 % discard short whiskers
            best_m(i+1) = 0;
            best_w(i+1) = 0;
        else
            
            best_m(i+1) = candidate_Ids(mx);
            
            % now work out which best_w corresponds to best_m
            candidate_wid = candidates_m(mx).wid;
            wx = find([candidates_w.id] == candidate_wid);
            best_w(i+1) = candidate_Ids(wx);
        end
        %     end
        clear real_whisker candidates_m candidates_w candidate_Ids candidate_wid mx wx
    end
    
    clear dropped i lengths ml
    
    %% Check whisker end point is within a reasonable distance from face
    % TO DO use hypotenuse distance from mean base position, not just y axis
    
    x = find(best_m);
    h = [[m(best_m(x)).follicle_x];[m(best_m(x)).follicle_y]];
    % y axis
    base_thresh = mean(h(2,:)) - std(h(2,:))*5; % 5 std as acceptable range
    discard = find(h(2,:)<base_thresh);
    
    % x axis
   	width_thresh = [mean(h(1,:)) - std(h(1,:))*3, mean(h(1,:)) + std(h(1,:))*3]; % 3 std as acceptable range
    discard_2 = [find(h(1,:) < width_thresh(1)) , find(h(1,:) > width_thresh(2))];
    
    discard = [discard, discard_2];
    best_m(x(discard)) = 0;
    best_w(x(discard)) = 0;
    
    clear discard discard_2 base_thresh h x
    
    %% Save whisker backbone pixels for best candidates
    disp(['Save whisker backbone pixels for best candidates'])
    wx = best_w(find(best_w));
    
    % determine distribution of lengths for determining size of
    % whisker_backbone matrix
    l = zeros(1,numel(wx));
    for i = 1:numel(wx);
        l(i) = numel(w(wx(i)).x);
    end
    
    whisker_backbone = zeros(2,max(l),numel(wx));
    for i = 1:numel(wx);
        n = numel(w(wx(i)).x);
        whisker_backbone(1,1:n,i) = w(wx(i)).x;
        whisker_backbone(2,1:n,i) = w(wx(i)).y;
    end
    
    clear n l i
    
    save([loadfile,'_clean.mat'],'whisker_backbone','best_m','best_w','width_thresh');
    
    %% Plot whisker backbones
    if plotting
        figure(1);clf;
        for i = 1: size(whisker_backbone,3); %4464;
            I = find(whisker_backbone(1,:,i));
            plot(whisker_backbone(1,I,i),whisker_backbone(2,I,i),'k');
            %             % in 3D
            %             plot3(i*ones(1,numel(I)),squeeze(whisker_backbone(1,I,i)),squeeze(whisker_backbone(2,I,i)),'color',[.5 .5 .5]);
            %             drawnow;
            hold all;
%             drawnow;plotting = 1
        end;
    end
    
    clear I i
    
    %% CHIMERA_2 method
    % Compute snout mask for each frame
    % Arclength (s) of each whisker backbone point
    % Compute x(s) and y(s) for whole whisker (5th order)
    % Compute x(s) and y(s) for mask segment (2nd order)
    % Use polyfit values to determine curvature at s(0)
    % Extrapolate to snout contour and find contour_intercept
    % Plug s(contour_intercept) into x(s) and y(s) to compute kappa at base
    % Plot comparison (circles) of curvature at base
    
    %%    If running this code after having saved 'best' whisker backbones, uncomment below to re-load
    
    %     loadfile = files(fidx).name(1:end-3);
    %     whiskerfile = [loadfile,'.whiskers'];
    %     measurementfile = [loadfile,'.measurements'];
    %
    %     w = LoadWhiskers(whiskerfile);
    %     m = LoadMeasurements(measurementfile);
    %
    %     load([loadfile,'_clean.mat'],'best_m','whisker_backbone','-mat');
    %     load([loadfile,'.tr'],'rall','kappa_all','theta_all','-mat');
    
    %%
    x = find(best_m);             % non-dropped frames
    %     time_per_loop = zeros(1,numel(x));
    %kappa_base = zeros(1,numel(x));
    kappa_w = zeros(1,numel(x));
    theta_w = zeros(1,numel(x));
    r_base = zeros(numel(x),2,3);
    
    % Set some parameters
%     x = find(best_m);             % non-dropped frames
    % roi for snout countour
    roi.x = 1:400; roi.y = 100:289; % Needs adjusting from session to session
    
    % Mask fit range, must be long enough to include interpolation of all
    % whiskers to face, but no longer as the fit to the face can degrade
%    mask_xrange = 190:290; % Likely to need changing from session to session
    mask_xrange = floor(width_thresh(1)):ceil(width_thresh(2));
    quad_fit_length = 100; % Amount of whisker to fit, in pixels
    z = 30; % Distance of mask from face contour
    
    % Plot mask_xrange boundaries on top of whisker backbones to see if the
    % range will do the job
    if plotting
        plot([mask_xrange(1),mask_xrange(1)],[0,300],'r');
        plot([mask_xrange(end),mask_xrange(end)],[0,300],'r');
        title('If any whisker base looks very close to one of the red lines, consider widening mask_ xrange')
        xlim([0,400]);
        drawnow;
    end
    
    % Snout mask fit params, probably don't need changing
    sigma = 6;
    hsize = 3*[sigma sigma];
    filter = fspecial('gaussian', hsize, sigma);
    face_z = 5;
    
    %% execution options
    plotting = 1; % plot as you go
    pausing = 0;  % pause after each frame is tracked
    zooming = 0;  % zoom in around the follicle
    debugging =1;
    %
    tic
   % clf
    for i = 1: numel(x); % 1570; %1500
        i
        
%            tic
        
        %% Find face gradient with whikerman code
        if debugging
        disp('Fitting snout mask');
        end
        [~,image] = poletracker([loadfile,'.dat'],15,x(i),0); % frame(i)
        
        ims.raw = image(roi.y,roi.x,1);
        ims.w = size(ims.raw,2);
        ims.h = size(ims.raw,1);
        ims.mfilt = medfilt2(ims.raw,[5 5]);
        
        ims.gfilt = imfilter(ims.mfilt,filter,'replicate');
        [ims.Del1.x,ims.Del1.y] = gradient(ims.gfilt);
        nvec = [1; 1]/sqrt(2);
        ims.grad = nvec(1)*ims.Del1.x+nvec(2)*ims.Del1.y;
        snout.x = 1:length(roi.x);
        [~,snout.y] = min(ims.grad);
        
        % shift contour in direction of gradient slightly
        snout.x = snout.x + face_z*nvec(1,:);
        snout.y = snout.y + face_z*nvec(2,:);
        
        % shift to account for roi mask
        s_contour(1,:) = snout.x + roi.x(1);
        s_contour(2,:) = snout.y + roi.y(1);
        
        if plotting
            clf;
            imagesc(image); colormap gray
            hold all;
            plot(s_contour(1,:),s_contour(2,:),'b');
        end
        
        clear snout ims nvec
        
        %% Compute arclength for each point of whisker_backbone
        if debugging
        disp('Compute arc length for each point of whisker backbone')
        end
        %         for i = 1000:2000;
        I = find(whisker_backbone(1,:,i));
        if I
            whisker_x = whisker_backbone(1,I,i);
            whisker_y = whisker_backbone(2,I,i);
            
            % Need to find a good way of ensuring whisker_y(1) is near the
            % face
            [~,endpoint] = max([whisker_y(1) , whisker_y(end)]);
            if endpoint == 2;
                whisker_x = fliplr(whisker_x);
                whisker_y = fliplr(whisker_y);
            end
            clear endpoint
            %             [whisker_y,sort_idx] = sort(whisker_y,2,'descend');
            %             whisker_x = whisker_x(sort_idx);
            
            len_s = sqrt(diff(whisker_x).^2 + diff(whisker_y).^2);
            arclen_s = [0 , cumsum(len_s)]; % zero for first entry
        end
        
        if plotting
            plot(whisker_x+1,whisker_y+1,'k');
            plot(whisker_x(1)+1,whisker_y(1)+1,'.w');
        end
        
        %         % Dubug plot
        %         plot(whisker_x,whisker_y,'color',[.5,.5,.5])
        %         hold all
        %         plot(whisker_x(1),whisker_y(1),'c.')
        %         plot(whisker_x(end),whisker_y(end),'y.')
        %         end
        
        clear len_s I
        
        %% Compute x(s) and y(s)
        if debugging
        disp('Compute x(s) and y(s)');
        end
        % Consider smoothing x and y first
        norm_s = arclen_s./max(arclen_s);
        
        [x_s] = polyfit(norm_s,whisker_x,5); % ,[],mu_x_s % Originally scaled but need correct coefficients for kappa computation later
        [y_s] = polyfit(norm_s,whisker_y,5);
        
        s_range = linspace(0,1,100); % linspace(min(arclen_s),max(arclen_s),100);
        
        whisker_poly_s(1,:) = polyval(x_s,s_range); % ,[],mu_x_s
        whisker_poly_s(2,:) = polyval(y_s,s_range); % ,[],mu_y_s
        
        if plotting
            plot(whisker_poly_s(1,:)+1,whisker_poly_s(2,:)+1,'m');
        end
        
        %% Extrapolate 5th order poly to face
%         if debugging
%         disp('Extrapolate to face and find intersection (pixels)');
%         end
%         tneg = linspace(-0.4,0.1,200);
%         w_ex(1,:) = polyval(x_s,tneg); %,[],mu_x_s);
%         w_ex(2,:) = polyval(y_s,tneg); %,[],mu_y_s);
%         
%         % Find intersection point between snout contour and whisker
%         % extrapolation. Whikerman style.
%         
%         % restrict s_contour to mask_xrange
%         face_points = s_contour(:,mask_xrange);
%         dist = zeros(size(w_ex,2),size(face_points,2));
%         
%         for a = 1:size(w_ex,2);
%             for b = 1:size(face_points,2);
%                 dist(a,b) = sum((w_ex(:,a) - face_points(:,b)).^2);
%             end
%         end
%         
%         % Find fp as min dist between contour and whisker
%         [min_c,c_idx] = min(dist');
%         [~,s_idx] = min(min_c);
%         fp_px = [w_ex(1,s_idx), w_ex(2,s_idx)];
%         
%         if plotting
%             plot(w_ex(1,:)+1,w_ex(2,:)+1,'color',[.5 .5 .5]);
%             plot(face_points(1,:),face_points(2,:),'color',[.5 .5 .5])
%             plot(fp_px(1)+1,fp_px(2)+1,'wo');
%         end
        
        %         %% Fit straight line to points of snout/whisker around fp and find
        %         % intersection of those points with fzero
        %         disp('Find whisker/face intersection with fzero')
        %         s_idx = min([size(w_ex,2)-1, s_idx]);
        %         c_idx = c_idx(s_idx);
        %
        %         p1 = polyfit(w_ex(1,s_idx-1:s_idx+1),w_ex(2,s_idx-1:s_idx+1),1);
        %
        %         p2 = polyfit(face_points(1,c_idx-3:c_idx+3),face_points(2,c_idx-3:c_idx+3),1); % More than 3 points to deal with pixelisation
        %
        %         fp(1) = fzero(@(h) polyval(p1-p2,h),0); % intersection of two lines in x
        %         fp(2) = polyval(p2,fp(1));
        %
        %         if plotting
        %             plot(fp(1)+1,fp(2)+1,'w+');
        %             %    drawnow;
        %         end
        %
        %         %% Work out arclength from fp to s(0)
        %         % Need to check sign of angle between fp estimates
        %         disp('Compute arc length from fp to s(0)')
        %         diff_fp = fp_px - fp;
        %         s_fp = tneg(s_idx) - (sqrt(diff_fp(1).^2 + diff_fp(2).^2))./max(arclen_s).*sign(atan2(diff_fp(1),diff_fp(2)));
        %
        clear dist w_ex min_c c_idx s_idx fp a b p1 p2 diff_fp
        
        %% Fit quadratic curve to whisker with snout mask + ROI chosen by distance from this mask
        if debugging
        disp(['Fit quadratic curve to snout contour'])
        end
        
        % Linear fit to snout to determine gradient to shift mask
        % Maybe use fp(1)-10:fp(1)+10 as mask_xrange?
        [g,~,mug] = polyfit(s_contour(1,mask_xrange),s_contour(2,mask_xrange),1);
        gv = polyval(g,1:2,[],mug);
        grad = gv(1) - gv(2);
        
        % Quadratic fit to face
        [p,~,mup] = polyfit(s_contour(1,mask_xrange),s_contour(2,mask_xrange),1);
        pv = polyval(p,s_contour(1,mask_xrange),[],mup);
        
        % Mask - Should be arclength, using distance from fp
        % z is distance out from the face in pixels, set at the start of the loop
        mask_1 = [(s_contour(1,mask_xrange))-z*grad; pv-z];
        [pm1,~,mum1] = polyfit(mask_1(1,:),mask_1(2,:),1);
        
        
        %% Set this_poly to x(s) and y(s)
        if debugging
        disp(['Work out mask and quadratic fit ROI'])
        end
        % remember     s_range = linspace(0,1,100);
        %              whisker_poly_s(1,:) = polyval(x_s,s_range);
        
        this_poly = whisker_poly_s;
        xrange = this_poly(1,:);
        
        % Quadratic fit to snout contour
        % mask 1
        pv_m1 = polyval(pm1,xrange,[],mum1);
        
        if plotting
            plot(xrange+1,pv_m1+1,'c')
        end
        
        % Find length of this_poly
        len = sqrt((diff(this_poly(1,:))).^2 + (diff(this_poly(2,:))).^2);
        % then compute cumulative sum of len
        cumlen = cumsum(len);
        
        % Find closest point on this_poly to the mask
        [~,idx_m1] = min((this_poly(2,:) - pv_m1).^2);
        
        % reset length to zero at mask intersection point
        length_at_mask = cumlen(idx_m1);
        
        cumlen_reset = cumlen-length_at_mask;
        
        % find element of whisker_poly arclen > quad_fit_length (nominally 100)
        % This returns a range of numbers, so we pick the smallest
        mask_int = find(abs(cumlen_reset)>=quad_fit_length); % 50

        
        if max(floor(cumlen_reset)) <= quad_fit_length
            disp('Changing mask roi length as whisker is short')
            quad_fit_length = floor(max(cumlen_reset));
            mask_int = find(abs(cumlen_reset)>=quad_fit_length)
            % reset quad_fit_length to 100;
            quad_fit_length = 100;
        end
        
        [~,idx_end] = min(abs(cumlen_reset(mask_int)));
        idx_end = mask_int(idx_end);
        
        % Make sure order is preserved
        idx_s = sort([idx_m1 idx_end]);
        
        % Region of whisker to fit is therefore defined as:
        fit_points = [this_poly(1,idx_s(1):idx_s(2));this_poly(2,idx_s(1):idx_s(2))];
        
        %% Compute x(s) and y(s) at mask with quadratic fit
        if debugging
        disp(['Compute x(s) and y(s) for quadratic fit w/ mask'])
        end
        
        s_range = linspace(0,1,100);
        fit_points_s = s_range(idx_s(1):idx_s(2));
        
        x_sq = polyfit(fit_points_s, fit_points(1,:),2);
        y_sq = polyfit(fit_points_s, fit_points(2,:),2);
        %         [wp,~,wmu] = polyfit(fit_points(2,:),fit_points(1,:),2);
        
        sq_range = linspace(min(fit_points_s),max(fit_points_s),21); % need to choose on the fly
        
        whisker_poly_sq(1,:) = polyval(x_sq,sq_range); % ,[],mu_x_s
        whisker_poly_sq(2,:) = polyval(y_sq,sq_range); % ,[],mu_y_s
        
        
        if plotting
            plot(whisker_poly_sq(1,:)+1,whisker_poly_sq(2,:)+1,'g');
        end
        %     drawnow
        
        %% Compute kappa at base by evaluating poly around base and plugging numbers into whiskerman functions
        if debugging
        disp('Compute kappa with whiskerman functions from 2nd order solution')
        end
        s_fm = sq_range(1); % s at mask
        tbase = [s_fm-0.01, s_fm, s_fm+0.01];
        r_base(i,1,:) = polyval(x_sq,tbase);
        r_base(i,2,:) = polyval(y_sq,tbase);
        
        if plotting
            plot(squeeze(r_base(i,1,:))+1,squeeze(r_base(i,2,:))+1,'y.')
            temp = WhikerMan_functions({squeeze(r_base(i,:,:)),[0:0.01:1]},8); b_tmp = temp{:};
            plot(b_tmp(1,:)+1,b_tmp(2,:)+1,'y');
            drawnow;
        end
        
        p0 = r_base(i,:,1);
        p2 = r_base(i,:,3);
        b = r_base(i,:,2);
        
        % rearrange B(t) = (1-t).^2 *p0 + 2*t*(1-t)*p1 + t.^2 *p2
        p1 = 2.*b - 0.5.*p0 - 0.5.*p2;
        
        r_base(i,:,2) = p1;
        
        temp = WhikerMan_functions({squeeze(r_base(i,:,:)),0},6);
        kappa_w(i) = temp{:} + eps;
        temp = WhikerMan_functions({squeeze(r_base(i,:,:)),0},5);
        theta_w(i) = temp{:} + eps;
        
        %% Save whisker_poly_2
        % This is slow so good to avoid if possible
%         save([loadfile,'_clean.mat'],'whisker_poly_s','whisker_poly_sq','-append');
        
        clear p0 p1 p2 b tbase s_fm
        
        
        %% Cleaning up
        clear x_prime y_prime x_2prime y_2prime s_fp s diff_fp whisker_x whisker_y x_s y_s temp p0 p1 p2 b
        if zooming
            xlim([200,300]); ylim([50,280]);
            drawnow;
        end
        
        if pausing
            pause
        end
        
%        toc
    end
    
    %% Fill in dropped frames by interpollating frames either side
    d = 1:numel(best_m);
    d(x) = 0;
    dropped_frames = find(d);
    
    kappa = zeros(1,numel(best_m));
    theta = zeros(1,numel(best_m));
    kappa(x) = kappa_w;
    theta(x) = theta_w;
    dropped_frames = find(kappa==0);
    
    disp(['Filling in ',num2str(numel(dropped_frames)),' dropped frames']);
    for j = dropped_frames;
        if j == 1;
           kappa(j) = kappa(end);
           theta(j) = theta(end);
        else
        kappa(j) = kappa(j-1);
        theta(j) = theta(j-1);
        end
    end
    
    kappa_w = kappa;
    theta_w = theta;
    
   save([loadfile,'_clean.mat'],'kappa_w','theta_w','r_base','dropped_frames','-append');
    
    %%
    clear a x k whisker_backbone whisker_poly whisker_poly_2 kappa_all theta_all 
    clear arclen_s b_tmp circles cumlen cumlen_reset face_z filter fit_points fit_points_s fp_px g grad gv hsize
    clear idx_end idx_m1 idx_s kappa_base kappa_baseq kappa_bw kappa_bw2 kappa_bw5 len length_at_mask mask_1 mask_int
    clear mask_xrange measurementfile mug mup norm_s p pc pv_m1 r_base roi s_range sigma sort_idx sq_range pv whiskerfile
    clear tbase theta_bw theta_bw2 theta_bw5 this_poly time_per_loop tneg whisker_poly_s whisker_poly_sq x_sq xrange y_sq z loadfile
    
    clear ax i image s_contour mum1 pm1 r rall temp wp face_points
    
    
    %% Plot theta/kappa as sanity check
    % Always plot here
    
    disp(['Sanity check plotting'])
    % Load .tr file data
    
    x = find(best_m);             % non-dropped frames
    a = [m(best_m(x)).angle];     % whiski angle
    k = [m(best_m(x)).curvature]; % whiski curvature
    
    % Plotting
    figure(3);clf;
    ax(1) = subplot(2,1,1);
    plot(x,-a) ;% angle from janelia tracker has different reference frame
    hold all;
    plot(theta_w);

    title('Theta');
    legend('Janelia','Chimera');
    
    ax(2) = subplot(2,1,2);
    plot(x,k)
    hold all
    plot(kappa_w)
    plot(dropped_frames,kappa_w(dropped_frames),'k.');
    title('Kappa');
    legend('Janelia','Chimera','Dropped')%
    
    linkaxes(ax,'x');
    % end
    drawnow;
    % pause;
    
    %%
    clear a x k whisker_backbone whisker_poly whisker_poly_2 kappa_all theta_all theta_w kappa_w w m best_m best_w
    clear ax i image s_contour mum1 pm1 r rall temp wp
    toc
end

