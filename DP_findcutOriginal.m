% DP method to find best x - position from the gaze
%  gaze_data is x-coordinate od the gazedata, original_width is width of
%  the input video , out_width is the width of the output video
% sigma is the variance of the gaussian
% duration = 30
% sigma with value 30 then e(-1/2) occurs at 30 frame on either side

function [backtrack,scenes,cuts,gaze_data]=DP_findcut2(gaze_data,original_width,out_width,sigma,duration)
    tic
    
    dis =10;
    edges = 1:dis:original_width; % downsample the data
    [Y] = discretize(gaze_data,edges);
   
    gaze_data = Y;
    
    n_frames = size(gaze_data,1);
    
    original_width = length(edges) %change original width
    out_width = out_width/dis;
    
    weight_map = zeros(original_width,n_frames);    
    frames = 1:n_frames;
    
    scenes  = cell(original_width,1);
    scenes_temp  = cell(original_width,1);
    for i =1:original_width
        scenes{i,1}=[1];
        scenes_temp  = cell(original_width,1);
    end
    
    for i=1:n_frames
        x = round(gaze_data(i)); % gaze at first frame.  
        %[i x gaze_data(i)]
        weight_map(x,:) = weight_map(x,:)-exp(-((frames-i).^2)/(2*sigma^2));
        %weight_map(x,i) = 1;
    end
    
    %figure,imshow(weight_map,[]);
    %pause
    
    toc
    DP_mat = zeros(original_width,n_frames);
    indi = zeros(original_width,n_frames);
    backtrack = zeros(n_frames,1);
    
    DP_mat(:,1) = weight_map(:,1);
    su = sum(weight_map,2);
    zindi = find(su==0);
    
    for i= 2:n_frames %frame_number
        i
        for j=1:original_width            
             DP_mat(j,i) = inf;
             if(find(zindi==j))
                 continue;
             end
            for k =1:original_width

                if(find(zindi==k))
                 continue;
                end
                [j k]
                [cost,flag] = cost_travel(gaze_data(j),gaze_data(k),out_width,scenes{k},i,duration,dis); %scenes
                if(cost+DP_mat(k,i-1)+weight_map(j,i) < DP_mat(j,i))
                    DP_mat(j,i) = cost+weight_map(j,i)+DP_mat(k,i-1);
                    indi(j,i) = k;
                    if(flag==1)
                        scenes_temp{j} = [scenes{k},i];
                    else
                        scenes_temp{j} = scenes{k};
                    end
                else
                    DP_mat(j,i) = DP_mat(j,i);
                end
            end
                    pause

        end
        scenes = scenes_temp;
    end
    
    [~,ind] = min(DP_mat(:,n_frames));
    backtrack(n_frames)=ind ;
    cuts = scenes{ind} ;
    
    for i=1:n_frames-1
        ind = indi(ind,n_frames+1-i);
        backtrack(n_frames-i) = ind;
    end
    toc
end

function [cost,flag] = cost_travel(x1,x2,out_width,scenes,i,duration,dis)
    flag=0;
     
    x1 = round(x1*dis+dis/2);
    x2 = round(x2*dis+dis/2);
    
    % if the distance is more than w , then there is no transition cost but
    % there is a scene cost.
    if(abs(x1-x2)>= out_width)
        flag=1;
        cost = exp(-(i-scenes(end))/(duration/2)) ;  % cost e_{-time/(duration/2)} with value reaching e_{-1} at duration/2   
    else
        cost = 2*(1 - exp(-abs(x1-x2)/(2*((out_width/8))))); % cost is 1-e_{-distance/(W/8)} with value reaching 0.9 at w/2
    end
    
end