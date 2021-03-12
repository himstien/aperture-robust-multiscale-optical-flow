%

clear
close all

% dirName = '~/POST_DOC/DATA/ATIS_PROPHESEE_DATA/';
dirName = '~/POST_DOC/DATA/atisData/bar_square/';

% dataSet = 'slide_depth_events_';
% dataSet = 'dual_pendulum_6_5e6';
% dataSet = 'shapes_004_fixed_cleaned_cross_'
dataSet = 'multiPattern1_fixed_'
% 'cityday_fixed_20mill_bestTheta_';
% 'peripheric1_17-12-04_15-53-19_cut_td_10mill_events_fixed_test';

disp('Loading raw data');

% raw_data = load([dirName 'peripheric1_17-12-04_15-53-19_cut_td_60mill_events_fixed_test' '.txt']); %load([dirName '/' dataSet '.txt']);
raw_data = load([dirName '/' dataSet '.txt']);
raw_data(:,1) = raw_data(:,1) + 1;
raw_data(:,2) = raw_data(:,2) + 1;
raw_data(:,3) = raw_data(:,3) - raw_data(1,3);

% disp('Loading grayscale data');
% rawGrayEvents = load_atis_aps([dirName '/shapes_004_aps.dat']);


%%

disp('Loading optical flow data');

LMFlow_matlab = false;

if(LMFlow_matlab)
    OpticalFlow_data = load([dirName '/' dataSet '_spatialFiltered_flowSurface_spatialSize_5_5_20_v1.mat']);
    OpticalFlow_data = OpticalFlow_data.OpticalFlow_data;
else
        OpticalFlow_data = load([dirName '/' dataSet '_FARMSOut_batch.txt']);
%     OpticalFlow_data = load([dirName '/' dataSet '_FARMSOut_bench_500us.txt']);
    OpticalFlow_data = OpticalFlow_data(OpticalFlow_data(:,5)~=0 & OpticalFlow_data(:,5)~=0, :);
end


OpticalFlow_data(:,1) = OpticalFlow_data(:,1) + 1;
OpticalFlow_data(:,2) = OpticalFlow_data(:,2) + 1;

[tempTheta, tempRad] = cart2pol(OpticalFlow_data(:,7), OpticalFlow_data(:,8));
tempTheta = rad2deg(wrapTo2Pi(tempTheta));

% allAngles_ = atan2(OpticalFlow_data(:,4), OpticalFlow_data(:,3));
% allAngles_ = wrapTo2Pi(allAngles_-pi/2);
%%

addGrayLevel = false;

close all

colMaps = hsv(25);
colorIndices = ceil( (tempTheta+0.001)/180);
% OpticalFlow_data(:,4) = OpticalFlow_data(:,4) - 1;

% OpticalFlow_data(:,2) = 127 - OpticalFlow_data(:,2) +1;

t_last = OpticalFlow_data(1, 3);

size_1 = 320;
size_2 = 240;

imageToShow = 0.5*ones(size_1, size_2);

quiverToShowCorr = [];
quiverToShowOrig = [];

gcf = figure(1);
% gcf.Position = [0, 0, 2000, 1000];
gcf.Visible = 'on';

window = 10000;
eventJump = 5;
count = 0;
startTime = 0; % 48000 for cross data
startGrayTime = 16.43e6; %4800000;

tic

if(addGrayLevel)
    eventsGrayInit = find(rawGrayEvents.ts > 0 & rawGrayEvents.ts < startGrayTime );
    eventsToShowGray = [rawGrayEvents.x(eventsGrayInit)' rawGrayEvents.y(eventsGrayInit)' rawGrayEvents.ts(eventsGrayInit)' rawGrayEvents.gray(eventsGrayInit)'];
    
    [~, indxS] = sort(eventsToShowGray(:,3), 'descend');
    eventsToShowGray = eventsToShowGray(indxS, :);
    
    eventsToShowGray(:,5) = eventsToShowGray(:,1)*320+eventsToShowGray(:,2);
    [C, uniq_i] = unique(eventsToShowGray(:,5), 'stable');
    eventsToShowGray = eventsToShowGray(uniq_i,:);
    
    imageGray = 0.5*ones(size_1, size_2);
    
    c=0;
    
    for in = 1:size(eventsToShowGray, 1)
        i_ = eventsToShowGray(in, 2);
        j_ = eventsToShowGray(in, 1);
        
        imageGray(i_+1, j_+1) = eventsToShowGray(in, 4)*4048;
    end
end

% imshow(imageGray);
toc

saveImages = true;

if(saveImages)
    mkdir([dirName '_rect_quiver/' dataSet '_rect_FARMS_histogram_images_' num2str(window) 'events_both_cpp_']);
end

makeImagesCorrected = true;


for i=startTime:window:OpticalFlow_data(end,3)-window
    
    imageExtra = 0.5*ones(size_1, size_2);
    count=count+1;
    
    clc
    disp(num2str((100*i/size(OpticalFlow_data,1))));
    %     if( abs(OpticalFlow_data(i, 3) - t_last) < 2000)
    
    indxEvents = find(raw_data(:,3)> i & raw_data(:,3) < i + window);
    rawEvents = raw_data(indxEvents, :);
    
    if(addGrayLevel)
        eventsGrayInWindow = find(rawGrayEvents.ts > startGrayTime+ i & rawGrayEvents.ts <= (i+window + startGrayTime) );
        eventsToShowGray = [rawGrayEvents.x(eventsGrayInWindow)' rawGrayEvents.y(eventsGrayInWindow)' rawGrayEvents.ts(eventsGrayInWindow)' rawGrayEvents.gray(eventsGrayInWindow)'];
    end
    
    indxEventsFlow = find(OpticalFlow_data(:,3)> i & OpticalFlow_data(:,3) < i + window);
    flowEvents = OpticalFlow_data(indxEventsFlow,:);
    flowEvents(:,5) = 1;
    flowEvents(:,9) = 1;
    
    for k = 1:numel(indxEvents)
        imageToShow(rawEvents(k,2), rawEvents(k,1)) = (128 + 64*rawEvents(k,4))/255;
    end
    
    % % %         for j=1:eventJump:numel(eventsGrayInWindow)
    % % %             imageGray(eventsToShowGray(j,2)+1, eventsToShowGray(j,1)+1) = eventsToShowGray(j, 4)*4048;
    % % %         end
    
    %         imageToShow(OpticalFlow_data(j,2), OpticalFlow_data(j,1)) = (128 + 64*OpticalFlow_data(j,4))/255;
    
    if(LMFlow_matlab)
        quiverToShowOrig = [flowEvents(1:eventJump:end, 2), flowEvents(1:eventJump:end,1), flowEvents(1:eventJump:end,6), flowEvents(1:eventJump:end,5)];
        quiverToShowCorr = [flowEvents(1:eventJump:end, 2), flowEvents(1:eventJump:end,1), flowEvents(1:eventJump:end,8), flowEvents(1:eventJump:end,7)];
    else
        quiverToShowOrig = [flowEvents(1:eventJump:end, 2), flowEvents(1:eventJump:end,1), flowEvents(1:eventJump:end,9).*cos(flowEvents(1:eventJump:end,10)), flowEvents(1:eventJump:end,9).*sin(flowEvents(1:eventJump:end,10))];
        quiverToShowCorr = [flowEvents(1:eventJump:end, 2), flowEvents(1:eventJump:end,1), flowEvents(1:eventJump:end,5).*cos(flowEvents(1:eventJump:end,6)), flowEvents(1:eventJump:end,5).*sin(flowEvents(1:eventJump:end,6))];
        %             spatialScales = flowEvents(1:eventJump:end, 11);
    end
    
    %     else
    
    % %         subplot(4, 4, [1 2 5 6]);
    % %         imshow(imageToShow);
    % %         hold on;
    % %         if(addGrayLevel)
    % %             imshow(imageGray);
    % %         end
    
    %         drawnow;
    if(addGrayLevel)
        imshow(imageGray);
    end
    %         drawnow;
    
    if(~isempty(flowEvents))
        
        
        
        t_last = flowEvents(1,3);
        
        % %             subplot(4, 4, [1 2 5 6]);
        
        % %             hold on;
        
        if(LMFlow_matlab)
            angle_ = atan2(quiverToShowOrig(:,4), quiverToShowOrig(:,3));
            angle = wrapTo2Pi(angle_-pi/2);
            angle_ = wrapTo2Pi(angle_-pi/2);
            
            for j = 1:size(quiverToShowOrig,1)
                quiver(quiverToShowOrig(j,2), quiverToShowOrig(j,1), quiverToShowOrig(j,4), quiverToShowOrig(j,3), 15, 'Color', colMaps(ceil(angle_(j)/(pi/12)),:));
            end
            
        else
            angle_ = flowEvents(1:eventJump:end,10);  %  atan2(quiverToShowOrig(:,3), quiverToShowOrig(:,4));
            angle = wrapTo2Pi(flowEvents(1:end,10));
            angle_ = wrapTo2Pi(angle_);
            
            colMaps_ = [];
            %                 pointsToScatter = zeros(size(quiverToShowOrig,1), 3);
            for j = 1:size(quiverToShowOrig,1)
                %                     pointsToScatter(j, :) = [quiverToShowOrig(j,1), quiverToShowOrig(j,2), ceil(angle_(j)/(pi/12))];
                colMaps_(j,:) = colMaps(ceil(angle_(j)/(pi/12))+1,:);
                
                imageExtra(quiverToShowOrig(j,1), quiverToShowOrig(j,2)) = ceil(angle_(j)/(pi/12));
                
                
                %                     hs = scatter(quiverToShowOrig(j,2), quiverToShowOrig(j,1), 20, 'MarkerFaceColor', colMaps(ceil(angle_(j)/(pi/12)),:), 'MarkerEdgeColor', colMaps(ceil(angle_(j)/(pi/12)),:));
                %                     alpha(hs, 0.5);
                %                     plot(quiverToShowOrig(j,2), quiverToShowOrig(j,1), '.', 'MarkerSize', 5, 'Color', colMaps(ceil(angle_(j)/(pi/12))+1,:));
                %                     drawnow;
            end
            
            
            
            %                 drawnow;
            %                 [sortedD, sortIndex] = sort(pointsToScatter(:,3));
            %                 pointsToScatter = pointsToScatter(sortIndex,:);
            %                 colMaps_ = hsv(size(pointsToScatter,1));
            subplot(4, 4, [1 2 5 6]);
            imshow(imageToShow);
            hold on;
            
%             t = annotation('textbox', [0.4, 0.885, 0.1, 0.1], 'String', [num2str(round ((t_last ) /1000)) ' msec'], 'EdgeColor', 'none');
%             t.FontSize = 12;
%             t.FontWeight = 'bold';
            
            hs = scatter(quiverToShowOrig(:,2), quiverToShowOrig(:,1), 3, colMaps_, 'Filled');
            
            for j = 1:size(quiverToShowCorr,1)
                quiver(quiverToShowCorr(j,2), quiverToShowCorr(j,1), 10*sin(angle_(j)), 10*cos(angle_(j)), 1, 'Color', colMaps_(j,:));
                hold on;
            end
            
            title('LOCAL FLOW');
%             ax = gca;
%             ax.Position = [0 0.5 0.5 0.5];
            %                 alpha(hs, 0.5);
            
            
            
            %         subplot(3, 4, [3 4 7 8]);
            %         imshow(imageToShow);
            %         pcolor(imageExtra);
            
            % % %                 title('Local Flow ');
        end
        
        %             angles = atan2(quiverToShowOrig(j,4), quiverToShowOrig(j,3));
        %             angles = wrapTo2Pi(angles-pi/2) ;
        %             angles = rad2deg(angles)  ;
        
        [hisogramDataOrig(count,:), xaxis_] = hist(angle, 0:pi/50:2*pi); %0:15:360);
        
        
        subplot(4, 4, [9 10 13 14]);        
        polarplot(xaxis_-pi/2, hisogramDataOrig(count,:)/sum(hisogramDataOrig(count,:)), 'LineWidth', 3);
%     rlim([0 0.5]);

%         ax = gca;
%             ax.Position = [0.5 0.5 0.5 0.5];
        %         title(['Time: ' num2str(round ((t_last ) /1000)) ' msec']);
        %             plot(xaxis_, hisogramDataOrig(count,:), 'LineWidth', 3);
        %             ylim([0 185]);
        %             xlim([0 2*pi]);
        %             xticks([0:pi/4:2*pi]);
        %             xtickformat('%.2f');
        %             title(['Distribution flow directions']);
        %             xlabel('angle [radian]');
        %             ylabel('frequency');
        
        
        %             subplot(4, 4, [1 2 3 4 5 6 7 8]); %correctyed
        %             hold on;
        %             drawnow;
        if(LMFlow_matlab)
            angle_ = atan2(quiverToShowCorr(:,3), quiverToShowCorr(:,4));
            angle = wrapTo2Pi(angle_-pi/2);
            angle_ = wrapTo2Pi(angle_-pi/2);
            
            for j = 1:size(quiverToShowCorr,1)
                quiver(quiverToShowCorr(j,2), quiverToShowCorr(j,1), quiverToShowCorr(j,3), -quiverToShowCorr(j,4), 15, 'Color', colMaps(ceil(angle_(j)/(pi/12))+1,:));
            end
            
        else
            angle_ = flowEvents(1:eventJump:end,6);% atan2(quiverToShowCorr(:,3), quiverToShowCorr(:,4));
            angle = wrapTo2Pi(flowEvents(1:end,6));
            angle_ = wrapTo2Pi(angle_);
            
            %                 for j = 1:size(quiverToShowOrig,1)
            %                     quiver(quiverToShowCorr(j,2), quiverToShowCorr(j,1), quiverToShowCorr(j,3), quiverToShowCorr(j,4), 15,'Color', colMaps(ceil(angle_(j)/(pi/12)),:));
            %                 end
            
            pointsToScatter = zeros(size(quiverToShowCorr,1), 3);
            
            for j = 1:size(quiverToShowCorr,1)
                pointsToScatter(j, :) = [quiverToShowCorr(j,1), quiverToShowCorr(j,2), ceil(angle_(j)/(pi/12))];
                
                colMaps_(j,:) = colMaps(ceil(angle_(j)/(pi/12))+1,:);
                
                imageExtra(quiverToShowCorr(j,1), quiverToShowCorr(j,2)) = ceil(angle_(j)/(pi/12));
                
                
                %                     hs = scatter(quiverToShowOrig(j,2), quiverToShowOrig(j,1), 20, 'MarkerFaceColor', colMaps(ceil(angle_(j)/(pi/12)),:), 'MarkerEdgeColor', colMaps(ceil(angle_(j)/(pi/12)),:));
                %                     alpha(hs, 0.5);
                %                     plot(quiverToShowOrig(j,2), quiverToShowOrig(j,1), '.', 'MarkerSize', 5, 'Color', colMaps(ceil(angle_(j)/(pi/12))+1,:));
                %                     drawnow;
            end
            %                 [sortedD, sortIndex] = sort(pointsToScatter(:,3));
            %                 pointsToScatter = pointsToScatter(sortIndex,:);
            %                 colMaps_ = hsv(size(pointsToScatter,1));
            
            subplot(4, 4, [3 4 7 8]);
            imshow(imageToShow);
            hold on;
            
            

            hs = scatter(quiverToShowCorr(:,2), quiverToShowCorr(:,1), 3, colMaps_, 'Filled');
            
            indxRect = [];
            
            indxRect = find(quiverToShowOrig(:,2)>=380 & quiverToShowOrig(:,2)<=590 & quiverToShowOrig(:,1)>=140 & quiverToShowOrig(:,1)<=245);
            
            for j = 1:size(quiverToShowCorr,1)
                quiver(quiverToShowCorr(j,2), quiverToShowCorr(j,1), 10*sin(angle_(j)), 10*cos(angle_(j)), 1, 'Color', colMaps_(j,:));
                hold on;
            end
            
            
            
            %                 for j = 1:numel(indxRect)
            %                     hold on;
            %                     r = rectangle('Position', [quiverToShowOrig(indxRect(j),2)-spatialScales(indxRect(j))-1, quiverToShowOrig(indxRect(j),1)-spatialScales(indxRect(j))-1, spatialScales(j)+1, spatialScales(j)+1], 'Curvature', 0.3, 'LineWidth', 2);
            %                     r.EdgeColor = [0 1 1];
            %                 end
            
            title(['ARMS FLOW (Corrected)']);
            % %             title(['Time: ' num2str(round ((t_last ) /1000)) ' msec']);
            
        end
        
        %             angles = atan2(quiverToShowCorr(j,4), quiverToShowCorr(j,3));
        %             angles = wrapTo2Pi(angles-pi/2) ;
        %             angles = rad2deg(angles)  ;
        [hisogramDataCorr(count,:), xaxis_] = hist(angle, 0:pi/50:2*pi); %0:15:360);
        
        [maxxx, indxMax(count)] = max(hisogramDataCorr(count,:));
        
        subplot(4, 4, [11 12 15 16]);
        %             hold on;
        %             plot(xaxis_, hisogramDataCorr(count,:), 'LineWidth', 3, 'Color', 'r');
        polarplot(xaxis_-pi/2, hisogramDataCorr(count,:)/sum(hisogramDataOrig(count,:)), 'LineWidth', 3, 'Color', 'r');
%         rlim([0 0.35]);
        
        %             ylim([0 185]);
        %             xlim([0 2*pi]);
        %             xticks([0:pi/4:2*pi]);
        %             xtickformat('%.2f');
        %             xlabel('angle [radian]');
        %             title(['Distribution flow directions']);
        %             ylabel('frequency');
        
        % % %         title(['Time: ' num2str(round ((t_last ) /1000)) ' msec']);
        
        
        
        drawnow;
        
        hold off;
        
        imageToShow = 0.5*ones(size_1, size_2);
        quiverToShowCorr = [];
        quiverToShowOrig = [];
        
        
        if(saveImages)
            %             saveas(gcf, [dirName '_1/' dataSet '_1_FARMS_histogram_images_' num2str(window) 'events_both_cpp_' '/opticalFlow_all_scatter_image_colored_' num2str(count) '.png']);
            
            saveas(gcf, [dirName '_rect_quiver/' dataSet '_rect_FARMS_histogram_images_' num2str(window) 'events_both_cpp_' '/opticalFlow_all_scatter_image_colored_' num2str(count) ], 'png');
            
            
            %             saveas(gcf, [dirName '_rect/' dataSet '_rect_FARMS_histogram_images_' num2str(window) 'events_both_cpp_' '/opticalFlow_all_scatter_image_colored_' num2str(count) ], 'epsc');
            %             saveas(gcf, [dirName '/' dataSet '_histogram_images_' num2str(window) 'events_both_cpp_' '/opticalFlow_all_scatter_image_colored_' num2str(count) '.pdf']);
            %             saveas(gcf, [dirName '/' dataSet '_histogram_images_' num2str(window) 'events_both_cpp_' '/opticalFlow_all_scatter_image_colored_' num2str(count) '.eps']);
            %             saveas(gcf, [dirName '/' dataSet '_histogram_images_' num2str(window) 'events_both_cpp_' '/opticalFlow_all_scatter_image_colored_' num2str(count) '.fig']);
        end
        %         pause(0.1);
        clf(gcf);
    end
    
end

close all;

