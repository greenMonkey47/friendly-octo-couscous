clc
clear all
close all

%% Cut Detection and Path Optimization

disp('Loading Data....')

load AllData
outwidth_all = [400 350 350 570 496 572 496 534 350 490 490 490 572];

idx_video = 9
fps = AllData{idx_video}.fps;
Data = AllData{idx_video}.Data;
resolution = AllData{idx_video}.resolution;

AR = 4/3;  % Required Aspect Ratio

l = min([length(Data{1}(:,1)) length(Data{2}(:,1)) length(Data{3}(:,1)) length(Data{4}(:,1)) length(Data{5}(:,1)) length(Data{6}(:,1)) ]);
L = AllData{idx_video}.NoofFrames;
delta = floor(abs(l-L))/2;
st = 1+7;
ed = l;
% st=2001;ed=3000;
AData = [Data{1}(st:ed,1) Data{2}(st:ed,1) Data{3}(st:ed,1) Data{4}(st:ed,1) Data{5}(st:ed,1) Data{6}(st:ed,1)];
tAData = AData;
% AData(1260:1366,:) = 890;

l = length(AData);
AData = AData(:,[1 3 4]);
data = median(AData');
data = data';

figure;
clr = 'rgbkyc'
for i=1:size(AData,2)
    plot(AData(:,i),clr(i));
    hold on
end
legend('m','n','o','c','d','u' );
axis([0 l 0 1366])
% out_width = AllData{idx_video}.resolution(2) * AR;
out_width = outwidth_all(idx_video) * AR;
size_data = size(data,1);
cut_dist = round(0.8*(out_width)); %was 150

%% Finding cuts in gaze data
% disp('Cut Detection using CVX...')
% 
N = size(data,1);
s_skip = 3;      % s-skip distance
fixtime = 24;    % Fixation time
k=24 ;      % no more than 1 cut in k frames
% tic
cuts_cvx = cut_detect_cvx(data,cut_dist,s_skip,fixtime,k);
% toc
figure,
%subplot(211)
plot(data,'.b')
hold on;
scatter(cuts_cvx,data(cuts_cvx),20,'*r')
% plot(abs(D1*data))
% plot(1:3600,ones(3600,1)*l2,'-k');
% plot(n3+1:N-n3,abs(D4*data.*(x(n3+1:N-n3))),'-m');
% plot(n3+1:N-n3,abs(D4*data),'-g');
% plot(n3+1:N-n3,abs(D4*data).*(x(n3+1:N-n3)) + (1-x(n3+1:N-n3))*(l2)-l2,'-g')
axis([0 l 0 1366])


%% DP
disp('Cut Detection using DP ....! ')
tic
% flag:{1-DP on all data, 0-DP on median}
[cuts_dp,dp_output,img,c1,cuts21] = cut_detect_DP(data,out_width,k*1,15,200,cut_dist,tAData,1);
%[cuts_dp_2,dp_output_2,img_2,c2,cuts22] = cut_detect_DP(data,out_width,k,15,100,cut_dist,AData,0);

toc



scatter(cuts_dp,dp_output(cuts_dp),20,'ok');
% scatter(cuts21,dp_output(cuts21),20,'+r');
% scatter(c1,dp_output(c1),20,'*r');

% scatter(cuts_dp_2,dp_output_2(cuts_dp_2),20,'or');
% scatter(cuts22,dp_output_2(cuts22),20,'+k');
% scatter(c2,dp_output_2(c2),20,'*b');

plot(dp_output,'-k');
legend('with all users')
% plot(dp_output_2,'-g');
% legend('with Median')

% plot(temp(:,1),'-g')
% plot(temp(:,2),'-m');
% plot(temp(:,3),'-k');

img2 = img(size(img,1):-1:1,:);
img2 = imresize(img2,[1366 l]);
figure,imshow(img2,[])

per_frameVar = [];
for i=1:N
    bg = max(i-s_skip,1);
    ed = min(i+s_skip,N);
    t = AData(bg:ed,:);
    per_frameVar(i,1) = var(t(:));
end
%per_frameVar = var(AData');

mini = min(per_frameVar);
per_frameVar = per_frameVar-mini;
maxi = max(per_frameVar);

per_frameVar = 1 - ((1-(per_frameVar/maxi))*0.3)' ;

% std
per_frameStd = [];
for i=1:N
    bg = max(i-5,1);
    ed = min(i+5,N);
    t = AData(bg:ed,:);
    per_frameStd(i,1) = std(t(:));
end
%per_frameVar = var(AData');

mini = min(per_frameStd);
per_frameStd = per_frameStd-mini;
maxi = max(per_frameStd);
W = outwidth_all(idx_video)*AR;

per_frameStd = W*(1 - ((1-(per_frameStd/maxi))*0.3)') ;


%% Path Optimization

% load original cuts
A = importdata(['./Videos/Original_Cuts/' AllData{idx_video}.filename(1:end-4) '_shots.txt'], ' ');
cuts_org = A(:,1);
% cuts_org(1:4) = [];
cuts = cuts_dp;
for i=length(cuts):-1:1
    if min(abs(cuts(i)-cuts_org)) < k
        cuts(i) = [];
    end
end



bool = ones(N,4);
tbool = ones(N,1);
%% generate bool variables
% for original cuts

% [dp_output_smooth,L2_cuts]= smooth_DP(cuts21,c1,dp_output,out_width,bool) ; 


bool(cuts_org(1):cuts_org(1)+s_skip,1) = 0;
tbool(cuts_org(1)) = 1;
for i=2:length(cuts_org)
        bool(cuts_org(i):cuts_org(i)+s_skip,1) = 0;
        bool(cuts_org(i):cuts_org(i),2) = 0;
        bool(cuts_org(i)-1:cuts_org(i),3) = 0;
        bool(cuts_org(i)-2:cuts_org(i),4) = 0;
        tbool(cuts_org(i)) = 1;
end
% for detected cuts
for i=1:length(cuts)
    if cuts(i)>3
        bool(cuts(i):cuts(i)+s_skip,1) = 0;
        bool(cuts(i):cuts(i),2) = 0;
        bool(cuts(i)-1:cuts(i),3) = 0;
        bool(cuts(i)-2:cuts(i),4) = 0;
        tbool(cuts(i)) = 1;
    end
end
bool(1:10,1) = 0;
bool(l-3:l,[0 2 3]+1) = 0;
bool(l-3:l,1) = 1;
bool = bool(1:l,:);

% bool variables for zoom : tbool
tbool(cuts_org(1):cuts_org(1)+s_skip,1) = 0;
for i=2:length(cuts_org)
        tbool(cuts_org(i):cuts_org(i)+s_skip,1) = 0;
        tbool(cuts_org(i):cuts_org(i),2) = 0;
        tbool(cuts_org(i)-1:cuts_org(i),3) = 0;
        tbool(cuts_org(i)-2:cuts_org(i),4) = 0;
end
% for detected cuts
for i=1:length(cuts)
    if cuts(i)>15
        tbool(cuts(i)-15:cuts(i)+15,1) = 0;
        tbool(cuts(i):cuts(i),2) = 0;
        tbool(cuts(i)-1:cuts(i),3) = 0;
        tbool(cuts(i)-2:cuts(i),4) = 0;
    end
end
tbool(1:10,1) = 0;
tbool(l-3:l,[0 2 3]+1) = 0;
tbool(l-3:l,1) = 1;
tbool = tbool(1:l,:);




lambda0 = 0.005;
lambda1 = 50;
lambda2 = 0;
lambda3 = 30;
vc1 = 3;
vc2 = 3;

thresh = out_width*0.05;
% tic
% [opt_data,temp1, temp2]=path_optimization_cvx(data,bool,tbool,lambda0,lambda1,lambda2,lambda3,vc1,vc2,thresh,out_width);
% toc
tic
[opt_data_dp,temp1_dp, temp2_dp,zoom_dp]=path_optimization_cvx(dp_output,bool,tbool,lambda0,lambda1,lambda2,lambda3,vc1,vc2,thresh,out_width,per_frameStd',1);
toc
z = zoom_dp/W;
% z = z + (1-max(z));

% generate w,h,x,y
Y = 271; W = outwidth_all(idx_video)*AR; H = outwidth_all(idx_video);
h = []; x = []; y = [];
g = opt_data_dp;

for i=1:N
   h(i,1) = round(H*z(i));
   w(i,1) = round(W*z(i));
   x(i,1) = round(g(i) -w(i)/2);
   y(i,1) = round(Y + ((H-h(i))/2));
end

figure,
%subplot(211)
plot(data,'.b')
hold on;
plot(dp_output,'-k')
clr = 'rgbkyc'
for i=1:3
    plot(AData(:,i),['.' clr(i)]);
    hold on
end
legend('m','n','o','c','d','u' );
scatter(cuts_dp,dp_output(cuts_dp),20,'*r');
scatter(cuts_org,dp_output(cuts_org),20,'*m');
% plot(opt_data,'-g')
plot(opt_data_dp,'-r')
plot(x,'-g')
plot(x+w,'-g')
plot(opt_data_dp+W/2,'-g')
plot(opt_data_dp-W/2,'-g')
legend('Gaze Data', 'Track')

figure
subplot(211)
for i=1:3
    plot(AData(:,i),['.' clr(i)]);
    hold on
end
plot(opt_data_dp,'-k')

legend('m','n','o','c','d','u' );
subplot(212)
plot(per_frameStd,'-b')
hold on
plot(z*W,'-r')

fileID = fopen(['./generateCroppedVideo/opt_path/9_zoom_x.txt'],'w');
fprintf(fileID,'%f\n',[x]');
fclose(fileID);
fileID = fopen(['./generateCroppedVideo/opt_path/9_zoom_y.txt'],'w');
fprintf(fileID,'%f\n',[y]');
fclose(fileID);
fileID = fopen(['./generateCroppedVideo/opt_path/9_zoom_w.txt'],'w');
fprintf(fileID,'%f\n',[w]');
fclose(fileID);
fileID = fopen(['./generateCroppedVideo/opt_path/9_zoom_h.txt'],'w');
fprintf(fileID,'%f\n',[h]');
fclose(fileID);
% fileID = fopen(['./generateCroppedVideo/opt_path/' AllData{idx_video}.filename(1:end-4) '_optpath.txt'],'w');
% fprintf(fileID,'%f \n',opt_data_dp');
% fclose(fileID);