function [output] = WhikerMan_functions(input,f_id)
% WhikerMan_functions
% output = WhikerMan_functions(input,f_id)
%    Standalone function that runs private functions of WhikerMan. 'input' and 'output' 
%    are cell arrays containing the necessary inputs and output for that
%    function. Help for individual functions can be found by digging into
%    the script for this function of WhikerMan itself.
% 
%    Each function can be called independently by choosing f_id, it's function ID, 
%    as follows:
%
%    1. Auto initialise 
%
%    2. Load frame 
%
%    3. Predict contour: r = predict_contour(rold,ddr,handles.miumax)
%    Predict where the whisker will be in the next frame
%
%    4. Post processing
%
%    5. Base angle 
%
%    6. Curvature
%
%    7. Bimfun: 
%    Evaluate error along whisker fit
%    Also callable with standalone Bimfun.m, allowing for optimisation with fminunc
%
%    8. bezierval: bezierval(rnew,0:.01:1)
%    Evaluate position of bezier curve 'bez' at points t . Included in
%    Bimfun.m too.
%
%    9. bezierdtval: tvec = bezierdtval(r,[0 1]) 
%    Evaluate first derivative of bezier curve 'bez' at points t 
%
%    10. bezierdt2val: 
%    Evaluate second derivative of bezier curve 'bez' at points t 
%
%    11. fit_curve 
%
%    12. find follicle: 
%    Locate follicle based on follicle mask and extrapolated whisker curve 
%
%    13. tmpfun 
%
%    14. snout_segment: 
%    Find follicle mask along snout contour 

%% List of functions, in generic input-output format

%%%% 1. Auto initialise %%%%
if f_id == 1;
    trfile = input{1};
    image = input{2};
    handles = input{3};
    [rinit, goodsolution] = auto_initialise(trfile, image, handles);
    ouput{1} = rinit;
    ouput{2} = goodsolution;
end

%%%% 2. Load frame %%%%
if f_id == 2;
    videopmtrs = input{1};
    frameidx = input{2};
    frame = load_frame(videopmtrs,frameidx);
    output{1} = frame;
end

%%%% 3. Predict contour %%%%
if f_id == 3;
    rold = input{1};
    ddr = input{2};
    miumax = input{3};
    r = predict_contour(rold, ddr, miumax);
    output{1} = r;
end

%%%% 4. Post processing %%%%
if f_id == 4;
    r = input{1};
    folliclemask = input{2};
    s0 = input{3};
    dt = input{4};
    sstar = input{5};
    [rnew, fp, fpidx, tfp, sfp, s, tstar, sp2] = postproc(r, folliclemask, s0, dt, sstar);
    output{1} = rnew; 
    output{2} = fp; 
    output{3} = fpidx;
    output{4} = tfp;
    output{5} = sfp;
    output{6} = s;
    output{7} = tstar; 
    output{8} = sp2;
end

%%%% 5. Base angle %%%%
if f_id == 5;
    r = input{1};
    t = input{2};
    theta = base_angle(r,t);
    output{1} = theta;
end

%%%% 6. Curvature %%%%
if f_id == 6;
    r = input{1};
    t = input{2};
    kappa = curvature(r,t);
    output{1} = kappa;
end

%%%% 7. Bimfun - evaluate error along whisker fit %%%%
if f_id == 7;
    z = input{1};
    r = input{2};
    n = input{3};
    im = input{4};
    dt = input{5};
    sigma = input{6};
    [E,gradE] = Bimfun(z, r, n, im, dt, sigma);
    output{1} = E;
    output{2} = gradE;
end

%%%% 8. bezierval - evaluate position of bezier curve 'bez' at points t %%%%
if f_id == 8;
    bez = input{1};
    t = input{2};
    b = bezierval(bez,t);
    output{1} = b;
end

%%%% 9. bezierdtval - evaluate first derivative of bezier curve 'bez' at points t %%%%
if f_id == 9;
    bez = input{1};
    t = input{2};
    dBdt = bezierdtval(bez,t);
    output{1} = dBdt;
end

%%%% 10. bezierdt2val - evaluate second derivative of bezier curve 'bez' at points t %%%%
if f_id == 10;
    bez = input{1};
    t = input{2};
    dB2dt2 = bezierdt2val(bez,t);
    output{1} = dB2dt2;
end

%%%% 11. fit_curve %%%%
if f_id == 11;
    r = input{1};
    order = input{2};
    [p,s] = fit_curve(r,order);
    output{1} = p;
    output{2} = s;
end

%%%% 12. find follicle - locate follicle based on follicle mask and extrapolated whisker curve %%%%
if f_id == 12;
    r = input{1};
    dt = input{2};
    folliclemask = input{3};
    [fp,fidx,tfp] = find_follicle(r,dt,folliclemask);
    output{1} = fp;
    output{2} = fidx;
    output{3} = tfp;
end

%%%% 13. tmpfun %%%%
if f_id == 13;
    t = input{1};
    r = input{2};
    fp = input{3};
    E = tmpfun(t,r,fp);
    output{1} = E;
end

%%%% 14. snout_segment - find follicle mask along snout contour %%%%
if f_id == 14;
    image = input{1};
    roi = input{2};
    folliclemask_old = input{3};
    fpidx_old = input{4};
    theta = input{5};
    [folliclemask] = snout_segment(image, roi, folliclemask_old, fpidx_old, theta);
    output{1} = folliclemask;
end





%% FUNCTIONS TAKEN DIRECTLY FROM WHIKERMAN

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 1. Auto initialise %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rinit, goodsolution] = auto_initialise(trfile, image, handles)

load(trfile,'-mat','theta_all','rall')

% choose a representative subset of frames:
[theta,frame] = sort(theta_all);
theta = theta(2:end-1);
frame = frame(2:end-1);
n = length(theta);
% Crude way to pick the subset, since biased to the most common angle.
% But, so long as the subset is big enough, it should satisfice.
interval = 2;%max([floor(n/100),1]);
frame_subset = frame(1:interval:end);

% for each selected frame, compute line integral (E) of contour on
% current image (high E is bad):
dt = 0.01;
t = 0:dt:1;
E = zeros(1,length(frame_subset));
figure, 
% set(gca,'ydir','reverse')
imagesc(image), colormap gray, hold on
for fidx = 1:length(frame_subset)
    frame = frame_subset(fidx);
    bez = squeeze(rall(frame,:,:));
    r = bezierval(bez,t);
    plot(r(1,:),r(2,:),'r')
    E(fidx) = dt*sum(interp2(image,r(1,:),r(2,:)));
end
clear fidx frame bez r

% find the contour/frame with the lowest E:
[Emin,Eidx] = min(E);
xlabel(sprintf('min(E)=%f',Emin))
if Emin < handles.energy_threshold
    % we've found a decent initial condition.  Use it:
    rinit = squeeze(rall(frame_subset(Eidx),:,:));
    goodsolution = true;
else
    % we've not.  Flag this and force manual approach.
    rinit = zeros(size(rall,2),size(rall,3));
    goodsolution = false;
    title('Press any key to continue with manual initialisation')
    pause
end
clear E Eidx
close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 2. Load frame %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function frame = load_frame(videopmtrs,frameidx)

switch videopmtrs.type
    case 'avi'
        video = read(videopmtrs.vObj, frameidx);
        frame = video(:,:,:,1);
    case 'dat'
        offset = videopmtrs.header.imagesize * (frameidx-1) + videopmtrs.offset;
        fseek(videopmtrs.fid,offset,-1);
        tmp = fread(videopmtrs.fid,videopmtrs.header.imagesize-24,'uint8=>uint8');
        tmp = reshape([tmp; zeros(24,1)],videopmtrs.width,videopmtrs.height)';
        frame = uint8(zeros(videopmtrs.height,videopmtrs.width,3));
        frame(:,:,1) = tmp;
        frame(:,:,2) = tmp;
        frame(:,:,3) = tmp;
        clear tmp
    otherwise
        error('Unhandled video file type')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 3. Predict contour %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function r = predict_contour(rold, ddr, miumax)
% The predictor is a linear combination of two terms.
% Term 1: P0 & P2 come from r = rold + ddr
%         'P1prev' is set so that the contour has the same shape
%         (curvature) as in the previous frame
% Term 2: P0 & P2 as above.
%         'P1linear' is the midpoint of P0 and P2.  The point of this is
%         that, when the contour is linear, the P1 coordinate is
%         ill-defined (any point along the P0-P2 line works).
%         Using the midpoint gives more stable extrapolation
%         beyond the range t=[0,1] and avoids instability.
% P1prev and P1linear are combined according to how linear the contour
% was in the previous frame.  If it was near-linear, P1linear is
% weighted most; if it was curved, P1prev is weighted most.
% Computation of term 1:
P0old = rold(:,1);
P1old = rold(:,2);
P2old = rold(:,3);
% P0 and P2 are easy:
P0 = P0old + ddr(:,1);
P2 = P2old + ddr(:,3);
% To set P1prev, the idea is that the shape of the quadratic Bezier
% curve is defined by the triangle P0-P1-P2.  So, once we've moved
% P0 and P2, we need to move P1 such that the triangle
% P0-P1-P2 is the same size and shape (but not orientation or
% position) as the triangle P0old-P1old-P2old.
% So, compute angle of line P0-P1 wrt the line P0-P2 so that we can
% find P1prev by rotation wrt the new P0-P2 line:
tmp = ((P2old-P0old)'*(P1old-P0old))/(norm(P2old-P0old)*norm(P1old-P0old));
% dirty hack 300114 copied 04/03/2014
tmp = min([tmp 1]);
theta = acos(tmp);

n = [0 -1; 1 0]*(P2old-P0old);  %normal to P0-P2
if (P1old-P0old)'*n > 0
    % sign of angle is correct
else
    theta = -theta;
end
clear tmp n
% Set P1 by rotation wrt the (new) P0-P2 and rescaling wrt the
% length of the line P0old-P1old:
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
P1prev = P0 + norm(P1old-P0old)*R*(P2-P0)/norm(P2-P0);
clear theta R
% Computation of term 2:
P1linear = (P0+P2)/2;
% Linearly combine P1prev and P1linear.  To determine the balance
% between them, quantify how curved the curve of the previous frame
% was.
% 'miu' is length of the normal from P1 to P0-P2, divided by the
% length of P0-P2.  mui = 0 for a straight line.  miu>0 for a
% curve:
miu = norm(P0old-2*P1old+P2old)/norm(P2old-P0old);
% Introduce a sensitivity parameter.  When miu>=miumax, P1=P1prev.
P1 = min([miu/miumax,1])*P1prev + (1-min([miu/miumax,1]))*P1linear;
clear P1prev P1linear miu miumax
r = [P0 P1 P2];
clear P0old P1old P2old P0 P1 P2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 4. Post processing %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rnew, fp, fpidx, tfp, sfp, s, tstar, sp2] = postproc(r, folliclemask, s0, dt, sstar)

% compute arc length parameter 's':
t = -1:dt:1.5;
idx = find(t==0);
b = bezierval(r,t);
db = diff(b,1,2);
ds = sqrt(sum(db.^2,1));
s = [0 cumsum(ds)]; % s(t=-1) = 0
s = s - s(idx); % now s(t=0) = 0
idx = find(t==1);
sp2 = s(idx);
clear db ds idx

% Find intersection between whisker polynomial and folliclemask
[fp, fpidx, tfp] = find_follicle(r,dt,folliclemask);
% if B(s) is the (x,y) point at 's', then fp = B(s(fpidx)).  this is used in snout detection.
[~,idx] = min((t-tfp).^2);
if length(idx)~=1
    error('')
end
sfp = s(idx);
clear idx

% Renormalise contour using polynomial (non-Bezier) approach:
% Relocate P0 and P2 so that there are at standardised distances from fp
% (s0) and adjust P1 to preserve shape.
if ~isempty(s0)
    % s0(1) is location of P0 wrt fp on the reference contour (first_frame)
    rnew = zeros(size(r));
%             rnew(:,1) = [polyval(p(1,:),sfp+s0(1));polyval(p(2,:),sfp+s0(1))];
%             rnew(:,3) = [polyval(p(1,:),sfp+s0(2));polyval(p(2,:),sfp+s0(2))];
%     idx = find(t==tfp);
    sprime = s - sfp;   % now sprime(t=tfp) = 0
    d = (sprime-s0(1)).^2;
    [~,idx] = min(d);
    rnew(:,1) = b(:,idx);
    clear d idx
    d = (sprime-s0(2)).^2;
    [~,idx] = min(d);
    rnew(:,3) = b(:,idx);
    clear d idx
    % find P1:
    % trying out another way 170414:
    % old way:
%     rnew(:,2) = 2*(bezierval(r,.5) - .25*rnew(:,1) - .25*rnew(:,3));
    % new way
    rnew(:,2) = 0.5*(rnew(:,1)-r(:,1)) + 0.5*(rnew(:,3)-r(:,3)) + r(:,2);
    
else
    % it's the first frame, so using this to define s0
    rnew = r;
end

[~,idx] = min((s-sfp-sstar).^2);
tstar = t(idx)
clear idx
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 5. Base angle %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function theta = base_angle(r,t)

% order = size(p,2)-1;
% pds = zeros(2,order);
% pds(1,:) = p(1,1:order).*(order:-1:1);
% pds(2,:) = p(2,1:order).*(order:-1:1);
% 
% dxds = polyval(pds(1,:),s);
% dyds = polyval(pds(2,:),s);
% theta = atan2(-dyds,dxds)*(180/pi);
% clear order pds dxds dyds

dBdt = bezierdtval(r,t);
theta = atan2(-dBdt(2,:),dBdt(1,:))*(180/pi);
clear dBdt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 6. Curvature %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function kappa = curvature(r,t)

% order = size(p,2)-1;
% pds = zeros(2,order);
% pds(1,:) = p(1,1:order).*(order:-1:1);
% pds(2,:) = p(2,1:order).*(order:-1:1);
% pds2(1,:) = pds(1,1:order-1).*(order-1:-1:1);
% pds2(2,:) = pds(2,1:order-1).*(order-1:-1:1);
% 
% dxds = polyval(pds(1,:),s);
% dyds = polyval(pds(2,:),s);
% d2xds2 = polyval(pds2(1,:),s);
% d2yds2 = polyval(pds2(2,:),s);
% 
% kappa = dxds.*d2yds2 - dyds.*d2xds2;
% kappa = kappa ./ (dxds.^2+dyds.^2).^(3/2);
% 
% clear order pds pds2 dxds dyds d2xds2 d2yds2

dBdt = bezierdtval(r,t);
d2Bdt2 = bezierdt2val(r,t);
kappa = (dBdt(1,:).*d2Bdt2(2,:)-dBdt(2,:).*d2Bdt2(1,:)) ./ (dBdt(1,:).^2+dBdt(2,:).^2).^(3/2);
clear dBdt d2Bt2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function bez = bezierfit(r,order)
% % Fit bezier curve of order 'order' to control points r
% % Initial conditions:
% switch order
%     case 1
%         bez0 = zeros(2,2);
%         bez0(:,1) = r(:,1);
%         bez0(:,2) = r(:,end);
%     case 2
%         bez0 = zeros(2,3);
%         bez0(:,1) = r(:,1);
%         bez0(:,3) = r(:,end);
%         m = ceil(size(r,2)/2);
%         bez0(:,2) = r(:,m);
%     otherwise
%         error('write more code!')
% end
% clear m
% % Minimise regression error:
% [bez,Emin,exitflag,output] = fminunc(@(bez) Bfun(bez,r),bez0,....
%                 optimset('Display','off'));%'TolX',.1,'TolFun',.01));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 7. Bimfun - evaluate error along whisker fit %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [E,gradE] = Bimfun(z, r, n, im, dt, sigma)

% Assume order 2 bezier curve
bez = zeros(2,3);
bez(:,1) = r(:,1) + z(1)*n(:,1);
bez(:,2) = z(2:3)';
bez(:,3) = r(:,3) + z(4)*n(:,2);
Ereg = 0.5 * sigma * sum(sum((bez-r).^2));
gradEreg = zeros(4,1);
gradEreg(1) = sigma * (bez(:,1)-r(:,1))'*n(:,1);
gradEreg(2:3) = sigma * (bez(:,2)-r(:,2));
gradEreg(4) = sigma * (bez(:,3)-r(:,3))'*n(:,2);
clear r
t = 0:dt:1;
r = bezierval(bez,t);
% dr = diff(r,1,2);
% ds = sqrt(sum(dr.^2,1));
% s = [0 cumsum(ds)];
% snew = linspace(s(2),s(end-1),length(s));
% clear dr ds
% % figure, plot(s,t,'.',snew,0,'.'),pause
% tnew = interp1(s,t,snew);
% r = bezierval(bez,tnew);
% t = tnew;
% clear s snew tnew

% find the points on the bezier curve that are internal to the image:
w = size(im.s,2);
h = size(im.s,1);
gdpt = find((r(1,:)>=1)&(r(1,:)<=w)&(r(2,:)>=1)&(r(2,:)<=h));
clear w h

%im.s
E = dt*sum(interp2(im.s,r(1,gdpt),r(2,gdpt)));

% fprintf('Eim=%.2f Ereg=%.2f Etot=%.2f\n',E,Ereg,E+Ereg)
E = E + Ereg;
gradE = zeros(4,1);
idx = interp2(im.dx,r(1,gdpt),r(2,gdpt));
idy = interp2(im.dy,r(1,gdpt),r(2,gdpt));
gradE(1) = dt*sum(idx.*((1-t(gdpt)).^2)*n(1,1)+idy.*((1-t(gdpt)).^2)*n(2,1));
gradE(2) = dt*sum(idx.*(2*(1-t(gdpt)).*t(gdpt)));
gradE(3) = dt*sum(idy.*(2*(1-t(gdpt)).*t(gdpt)));
gradE(4) = dt*sum(idx.*(t(gdpt).^2)*n(1,2)+idy.*(t(gdpt).^2)*n(2,2));
clear gdpt idx idy
gradE = gradE + gradEreg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 8. bezierval - evaluate position of bezier curve 'bez' at points t %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function b = bezierval(bez,t)
% Evaluate bezier curve with parameters bez(:,1),bez(:,2),... at points t
% size(t) = [1,N], where N is number of points at which to evaluate
% function.  Typically, t = [0,1]
% size(bez) = [2,order+1]
order = size(bez,2)-1;
switch order
    case 1
        p0 = bez(:,1);
        p1 = bez(:,2);
        b = p0 + (p1-p0)*t;
    case 2
        p0 = bez(:,1);
        p1 = bez(:,2);
        p2 = bez(:,3);
        b = p0*(1-t).^2 + 2*p1*((1-t).*t) + p2*t.^2;
    otherwise
        error('Write more code!')
end
clear p0 p1 p2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 9. bezierdtval - evaluate first derivative of bezier curve 'bez' at points t %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dBdt = bezierdtval(bez,t)
% Evaluate 1st deriv of bezier curve wrt t, with parameters bez(:,1),bez(:,2),... at points t
% size(t) = [1,N], where N is number of points at which to evaluate function.  Typically, t = [0,1]
% size(bez) = [2,order+1]
order = size(bez,2)-1;
switch order
    case 2
        p0 = bez(:,1);
        p1 = bez(:,2);
        p2 = bez(:,3);
%         dBdt = -p0*2*(1-t) + 2*p1*(1-2*t) + 2*p2*t;
        dBdt = 2*(p0-2*p1+p2)*t + 2*(-p0+p1)*ones(size(t));
    otherwise
        error('Write more code!')
end
clear p0 p1 p2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 10. bezierdt2val - evaluate second derivative of bezier curve 'bez' at points t %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dB2dt2 = bezierdt2val(bez,t)
% Evaluate 2nd deriv of bezier curve wrt t, with parameters bez(:,1),bez(:,2),... at points t
% size(t) = [1,N], where N is number of points at which to evaluate function.  Typically, t = [0,1]
% size(bez) = [2,order+1]
order = size(bez,2)-1;
switch order
    case 2
        p0 = bez(:,1);
        p1 = bez(:,2);
        p2 = bez(:,3);
        dB2dt2 = 2*(p0-2*p1+p2) * ones(size(t));
    otherwise
        error('Write more code!')
end
clear p0 p1 p2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 11. fit_curve %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p,s] = fit_curve(r,order)
% Fit polynomial curve (parameterised by 's') to the control points:
dr = sqrt(sum(diff(r,1,2).^2,1));
s = [0 cumsum(dr)];
p = zeros(2,order+1);
p(1,:) = polyfit(s,r(1,:),order);
p(2,:) = polyfit(s,r(2,:),order);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 12. find follicle - locate follicle based on follicle mask and extrapolated whisker curve %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fp,fidx,tfp] = find_follicle(r,dt,folliclemask)

% compute arc length parameter 's':
% Extrapolate the whisker back towards the follicle. 
% Old version (-1:dt:0) breaks if segment is small and far from the follicle, 
% such that the fit-length-extrapolation doesn't reach the follicle mask
tneg = -2:dt:0;      

w = bezierval(r,tneg);
dw = diff(w,1,2);      % MHE these don't do anything
ds = sqrt(sum(dw.^2,1));
% sneg = [0 cumsum(ds)] - sum(ds); % s(1)=0 ~ P0
clear dw ds

% w = zeros(2,length(tneg));
% w(1,:) = polyval(p(1,:),sneg);
% w(2,:) = polyval(p(2,:),sneg);

%keyboard
dist = zeros(length(tneg),size(folliclemask,2));
for i = 1:length(tneg)
    for j = 1:size(folliclemask,2)
        dist(i,j) = sum((w(:,i)-folliclemask(:,j)).^2,1);
    end
end
clear i j
[m,fidx] = min(dist');  % min over folliclemask for fixed s
% locate first minimum (in order of decreasing t)
% sidx = find(diff(m)>0,1);
[~,sidx] = min(m);      % min of the min (over s)
%     fidx = fidx(sidx);      % pt on f closest to w
% now we have Order(0) follicle position estimate:
fp = [w(1,sidx); w(2,sidx)];
clear dist m tmp
 
% keyboard
% use linear interpolation to increase accuracy of fp estimate.
% find line that locally fits folliclemask (p2) and line that locally fits
% whisker contour (p1):
sidx = min([size(w,2)-1 sidx]);
% keyboard
p1 = polyfit(w(1,sidx-1:sidx+1),w(2,sidx-1:sidx+1),1);
fidx = fidx(sidx);
p2 = polyfit(folliclemask(1,fidx-1:fidx+1),folliclemask(2,fidx-1:fidx+1),1);
% now solve for the point where p1 and p2 intersect:
M = [p1(1) -1; p2(1) -1];
fp = -inv(M)*[p1(2);p2(2)];
clear p1 p2 M

% find the 't' value of fp.  Do this by minimising ||f(t)-fp||^2 wrt t,
% starting from initial condition 't0':
t0 = tneg(sidx);
tfp = fminunc(@(t) tmpfun(t, r, fp), t0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 13. tmpfun %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function E = tmpfun(t,r,fp)
fphat = zeros(2,1);
% fphat(1) = polyval(p(1,:),s);
% fphat(2) = polyval(p(2,:),s);
fphat = bezierval(r,t);
E = sum((fp-fphat).^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 14. snout_segment - find follicle mask along snout contour %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [folliclemask] = snout_segment(image, roi, folliclemask_old, fpidx_old, theta)

% Filter and differentiate:
ims.raw = image(roi.y,roi.x,1);
ims.w = size(ims.raw,2);
ims.h = size(ims.raw,1);
ims.mfilt = medfilt2(ims.raw,[5 5]);
sigma = 6;
hsize = 3*[sigma sigma];
filter = fspecial('gaussian', hsize, sigma);
ims.gfilt = imfilter(ims.mfilt,filter,'replicate');
clear sigma hsize filter
% ims.Del1.x = [zeros(size(ims.gfilt,1),1) diff(ims.gfilt,1,2)];
% ims.Del1.y = [zeros(1,size(ims.gfilt,2)); diff(ims.gfilt,1,1)];
% ims.Del2.x = [zeros(size(ims.gfilt,1),1) diff(ims.gfilt,2,2) zeros(size(ims.gfilt,1),1)];
% ims.Del2.x = [zeros(1,size(ims.gfilt,2)); diff(ims.gfilt,2,1); zeros(1,size(ims.gfilt,2))];
[ims.Del1.x,ims.Del1.y] = gradient(ims.gfilt);

% Specify the desired gradient measure:
if ~isempty(folliclemask_old)
% Revised code (290414)
% Go back to using normal to snout mask instead of tangent to whisker
% solution:
%     % use normal to snout mask in previous frame to define gradient
%     % direction:
    tvec = folliclemask_old(:,fpidx_old) - folliclemask_old(:,fpidx_old+1);
    tvec = tvec / norm(tvec);
    R = [0 1;-1 0];
    nvec = R*tvec;
    clear tvec R
% So following lines (290414) commented out:
%     % use tangent to whisker solution in previous frame to define gradient
%     % direction
%     nvec = [-cosd(theta); sind(theta)];
%     nvec = nvec/norm(nvec);
    
else
    % no previous frame, so use default:
    nvec = [1; 1]/sqrt(2);
end 
ims.grad = nvec(1)*ims.Del1.x+nvec(2)*ims.Del1.y;

snout.x = 1:length(roi.x);
[~,snout.y] = min(ims.grad);
% Take every Nth pixel:
N = 5;
snout.x = snout.x(1:N:end);
snout.y = snout.y(1:N:end);
Npts = length(snout.x);
clear N

% Quality control on snout detection, debugging:
%     figure
%     subplot(2,2,1), imagesc(ims.raw), colormap gray
%     hold on, plot(snout.x,snout.y,'r.-')
%     subplot(2,2,2), imagesc(ims.gfilt)
%     subplot(2,2,3), imagesc(ims.mfilt)
%     subplot(2,2,4), imagesc(ims.grad)
%     pause

% The above produces a snout estimate that is outside the face
% slightly.  Warp the contour in a bit in direction of local normals

z = 5;
folliclemask(1,:) = snout.x + z*nvec(1,:);
folliclemask(2,:) = snout.y + z*nvec(2,:);
clear nvec

%  Debug:
%     figure
%     imagesc(ims.raw), colormap gray, colorbar
%     hold on,
%     plot(folliclemask(1,:),folliclemask(2,:),'m-',trackmask(1,:),trackmask(2,:),'g-')

% Align masks to complete image:
folliclemask(1,:) = folliclemask(1,:) + roi.x(1);
folliclemask(2,:) = folliclemask(2,:) + roi.y(1);
clear roi

% The extreme points are probably junk.  Cut them:
L = 1;
folliclemask(:,1:L) = [];
folliclemask(:,end-L+1:end) = [];
clear L

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%