% DP method to find best x - position from the gaze
%  gaze_data is x-coordinate od the gazedata, original_width is width of
%  the input video , out_width is the width of the output video
% sigma is the variance of the gaussian
% duration = 30
% sigma with value 30 then e(-1/2) occurs at 30 frame on either side

function [backtrack,scenes,cuts,gaze_data,weight_map,cuts2]=DP_findcut2(gaze_data,original_width,out_width,cut_dist,sigma,duration,AllData,flag_data)
    
    dis =10;
    edges = 1:dis:original_width; % downsample the data
    
%     [a,b]=find(AllData>original_width);
%     
%     for i=1:length(a)
%         AllData(a(i),b(i)) = original_width;
%     end
%     [a,b]=find(AllData<0);
%     
%     for i=1:length(a)
%         AllData(a(i),b(i)) = 1;
%     end

    
    [gaze_data] = discretize(gaze_data,edges);
   % AllData(3918,:)
    AllData = discretize(AllData,edges);
    %AllData(3918,:)
    
    % gaze data is changed discrete data
    
    n_frames = size(gaze_data,1);
    
    original_width = length(edges); %change original width
    out_width = out_width/dis;
    cut_dist = cut_dist/dis;
    sigma = round(sigma/dis);
    
    weight_map = zeros(original_width,n_frames);    
    % size(weight_map)
    temp = weight_map;
    frames = 1:n_frames;
    
    scenes  = cell(original_width,1);
    scenes_temp  = cell(original_width,1);
    
    scenes2  = cell(original_width,1);
    scenes_temp2  = cell(original_width,1);
    
    for i =1:original_width
        scenes{i,1}=[1];
        scenes2{i,1}=[1];
    end
%     size(AllData)
    % max(AllData(:))
    % pause
    list = [1 4 3]
    if flag_data
        for i=1:n_frames
            for j=1:3 %this need to be "n" users 
                x = round(AllData(i,list(j))); % gaze at first frame.
%                  [i j x AllData(i,j)]
                if ~isnan(x) 
                    weight_map(x,:) = weight_map(x,:)-exp(-((frames-i).^2)/(2*sigma^2));
                    temp(x,i) = 1;
                end
            end
        end
        weight_map = weight_map/3;
    else
        for i=1:n_frames
            x = round(gaze_data(i)); % gaze at first frame.
            %[i x gaze_data(i)]
            if ~isnan(x)
            weight_map(x,:) = weight_map(x,:)-exp(-((frames-i).^2)/(2*sigma^2));
            temp(x,i) = 1;
            end
        end
    end
    
    filt = fspecial('gaussian',sigma+1,sigma);
    weight_map = conv2(temp,filt,'same');
    [a,b]=find(weight_map>0);
    for i=1:length(a)
        weight_map(a(i),b(i)) = -weight_map(a(i),b(i));
    end
    %weight_map = -weight_map;
    
    filt = fspecial('gaussian',5,10);
    temp = conv2(temp,filt,'same');
%     figure,imshow(weight_map,[]);
%     pause
%     
    
    DP_mat = zeros(original_width,n_frames);
    indi = zeros(original_width,n_frames);
    backtrack = zeros(n_frames,1);
    
    DP_mat(:,1) = weight_map(:,1);
    su = sum(weight_map,2);
    zindi = find(su==0);
    
    for i= 2:n_frames %frame_number
%         [i n_frames]
        for j=1:original_width            
             DP_mat(j,i) = inf;
             if(temp(j,i)==0)
                 continue;
             end
            for k =1:original_width
                if(temp(k,i-1)==0)
                 continue;
                end
                %if i>=388
            %    [i j k]
%                [gaze_data(j),gaze_data(k),i,duration,dis]
%                temp(j,i),temp(k,i)
                %end
                [cost,flag] = cost_travel(j,k,out_width,cut_dist,scenes{k},scenes2{k},i,duration,dis); %scenes
                if(cost+DP_mat(k,i-1)+weight_map(j,i) < DP_mat(j,i))
                    
                    DP_mat(j,i) = cost+weight_map(j,i)+DP_mat(k,i-1);
                    indi(j,i) = k;
                    
                    if(flag==1)
                        scenes_temp{j} = [scenes{k},i];
                    else
                        scenes_temp{j} = scenes{k};
                    end
                    x1 = round(j*dis+dis/2);
                    x2 = round(k*dis+dis/2);
                    
                    if( (abs(x1-x2)> 1 && abs(x1-x2) < cut_dist) && flag~=1)                       
                        scenes_temp2{j} = [scenes2{k},i];
                    else
                        scenes_temp2{j} = scenes2{k};
                    end
                    
                else
                    DP_mat(j,i) = DP_mat(j,i);
                end
            end
        end
        scenes = scenes_temp;
        scenes2 = scenes_temp2;
    end
    
    [~,ind] = min(DP_mat(:,n_frames));
    backtrack(n_frames)=ind ;
    cuts = scenes{ind} ;
    cuts2 = scenes2{ind};
    
    for i=1:n_frames-1
        ind = indi(ind,n_frames+1-i);
        backtrack(n_frames-i) = ind;
    end
    
    backtrack = backtrack * dis;
end

function [cost,flag] = cost_travel(x1,x2,out_width,cut_dist,scenes,scenes2,i,duration,dis)
    flag=0;
     
    x1 = round(x1*dis+dis/2);
    x2 = round(x2*dis+dis/2);
    
    % if the distance is more than w , then there is no transition cost but
    % there is a scene cost.
    if(abs(x1-x2)>= cut_dist)
        flag=1;
%         cost = 10*exp(-(i-scenes(end))/(duration/2)) ;  % cost e_{-time/(duration/2)} with value reaching e_{-1} at duration/2   
        cost = 2 + 2*exp(-(i-scenes(end))/(duration)) ;  % cost e_{-time/(duration/2)} with value reaching e_{-1} at duration/2   

    else
%         cost = 20*(1 - exp(-abs(x1-x2)/(2*((out_width/8))))); % cost is 1-e_{-distance/(W/8)} with value reaching 0.9 at w/2
%         cost = cost + 2*exp(-(i-scenes2(end))/(duration/2));
        cost = 2*(1 - exp(-abs(x1-x2)/(2*((out_width/8))))); % cost is 1-e_{-distance/(W/8)} with value reaching 0.9 at w/2
    end
    
end