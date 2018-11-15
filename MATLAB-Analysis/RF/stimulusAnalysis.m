clear xyname xypath xyTrajTrue xytraj VAR
[xyname, xypath] = uigetfile('*.*','Select stimGL framevar (xytraj) file:');
VAR=importdata([xypath,xyname]);
xytraj = VAR.data(:,5:6); 


xyTrajTrue(:,1) = xytraj(:,2); % x is second row
xyTrajTrue(:,2) = xytraj(:,1); % y is first row

disp(['Stimulus Duration is :', num2str(size(xyTrajTrue,1)/(120*3)), 'seconds'])
        % /(120*3) - 3 frames of colour and 120 fps projector

response_offset = 0;
degreesTraj = atan2( diff(xyTrajTrue(response_offset+1:end,2),1), ...
                     diff(xyTrajTrue(response_offset+1:end,1),1) ); %dy/dx
degreesTraj = (degreesTraj)*180/pi;

degreesTraj( find(diff(xyTrajTrue(response_offset+1:end,2),1)== 0 & ...
                  diff(xyTrajTrue(response_offset+1:end,1),1)== 0) ,1 ) = NaN;
    % when there is no change in x or y, we set it to not a number

degreeBin = 45;
    
% ------>   0 degrees
leftToRight = (degreesTraj>0-degreeBin & degreesTraj<0+degreeBin);
    xyLtoR = xyTrajTrue(leftToRight,:);    % Find xy traj of this motion
% <------   -180 OR +180 degrees
rightToLeft = (degreesTraj>+180-degreeBin & degreesTraj<=180 | degreesTraj<-180+degreeBin & degreesTraj>=-180);
    xyRtoL = xyTrajTrue(rightToLeft,:);    % Find xy traj of this motion
% /\        +90 degrees
downToUp = (degreesTraj>90-degreeBin & degreesTraj<90+degreeBin);
    xyDtoU = xyTrajTrue(downToUp,:);    % Find xy traj of this motion
% \/        -90 degrees
upToDown = (degreesTraj>-90-degreeBin & degreesTraj<-90+degreeBin);
    xyUtoD = xyTrajTrue(upToDown,:);    % Find xy traj of this motion


figure('Name','Stimulus Trajectories')
subplot(2,2,1)
    scatter(xyTrajTrue(leftToRight,1), xyTrajTrue(leftToRight,2), 10,[0.8 0.8 0.8],'filled')
    xlim([ min(xyTrajTrue(:,1)), max(xyTrajTrue(:,1)) ]); ylim([ min(xyTrajTrue(:,2)), max(xyTrajTrue(:,2)) ]);
    title('Left to Right')
subplot(2,2,3)
    scatter(xyTrajTrue(rightToLeft,1), xyTrajTrue(rightToLeft,2), 10,[0.8 0.8 0.8],'filled')
    xlim([ min(xyTrajTrue(:,1)), max(xyTrajTrue(:,1)) ]); ylim([ min(xyTrajTrue(:,2)), max(xyTrajTrue(:,2)) ]);
    title('Right to Left')
subplot(2,2,2)
    scatter(xyTrajTrue(downToUp,1), xyTrajTrue(downToUp,2), 10,[0.8 0.8 0.8],'filled')
    xlim([ min(xyTrajTrue(:,1)), max(xyTrajTrue(:,1)) ]); ylim([ min(xyTrajTrue(:,2)), max(xyTrajTrue(:,2)) ]);
    title('Down to Up')
subplot(2,2,4)
    scatter(xyTrajTrue(upToDown,1), xyTrajTrue(upToDown,2), 10,[0.8 0.8 0.8],'filled')
    xlim([ min(xyTrajTrue(:,1)), max(xyTrajTrue(:,1)) ]); ylim([ min(xyTrajTrue(:,2)), max(xyTrajTrue(:,2)) ]);
    title('Up to Down')
    

%% Heat Maps 
figure('Name','Stimulus Heat Map')


colormap(hot) % number inside indicates number of colours to use
numBins = round([ max(xyTrajTrue(:,1))/50, max(xyTrajTrue(:,2))/60]); % Number of bins for each dimensions
                % essentially [480/N, 640/N]        for a 480x640 projection
gaussianSigma = 0; % standard deviation for gaussian for smoothing
                      % set to 0 for no change
                
% imagesc() has all inputs transposed: columns go horizontally - hence are x co-ordinates
               
                % manually added elements are to force histogram to fill plot area
subplot(2,2,1)
[N2,c2] = hist3([[xyTrajTrue(leftToRight,1);0;max(xyTrajTrue(:,1))], ...   % Trajectoriess Histogram - i.e. how often is the trajectory in the same bin as we used for spikes
                 [xyTrajTrue(leftToRight,2);0;max(xyTrajTrue(:,2))]], ...
                 numBins); N2(1,1)=N2(1,1)-1; N2(end,end) = N2(end,end)-1; % subtract elements we added
              if gaussianSigma<=0
                imagesc(c2{1}([1 end]),c2{2}([1 end]),N2');     % c2 - pixel centres, N - pixel values (NOTE: TRANSPOSE)
              else
                imagesc(c2{1}([1 end]),c2{2}([1 end]),imgaussfilt(N2',gaussianSigma));
              end             
              axis xy % ensure y axis points up
              colorbar
              title('Left to Right')
              xlim([ min(xyTrajTrue(:,1)), max(xyTrajTrue(:,1)) ]); ylim([ min(xyTrajTrue(:,2)), max(xyTrajTrue(:,2)) ]);
subplot(2,2,3)
[N2,c2] = hist3([[xyTrajTrue(rightToLeft,1);0;max(xyTrajTrue(:,1))], ...   % Trajectoriess Histogram - i.e. how often is the trajectory in the same bin as we used for spikes
                 [xyTrajTrue(rightToLeft,2);0;max(xyTrajTrue(:,2))]], ...
                 numBins); N2(1,1)=N2(1,1)-1; N2(end,end) = N2(end,end)-1; % subtract elements we added
              if gaussianSigma<=0
                imagesc(c2{1}([1 end]),c2{2}([1 end]),N2');     % c2 - pixel centres, N - pixel values (NOTE: TRANSPOSE)
              else
                imagesc(c2{1}([1 end]),c2{2}([1 end]),imgaussfilt(N2',gaussianSigma));
              end
              axis xy % ensure y axis points up
              colorbar
              title('Right to Left')
              xlim([ min(xyTrajTrue(:,1)), max(xyTrajTrue(:,1)) ]); ylim([ min(xyTrajTrue(:,2)), max(xyTrajTrue(:,2)) ]);

subplot(2,2,2)
[N2,c2] = hist3([[xyTrajTrue(downToUp,1);0;max(xyTrajTrue(:,1))], ...   % Trajectoriess Histogram - i.e. how often is the trajectory in the same bin as we used for spikes
                 [xyTrajTrue(downToUp,2);0;max(xyTrajTrue(:,2))]], ...
                 numBins); N2(1,1)=N2(1,1)-1; N2(end,end) = N2(end,end)-1; % subtract elements we added
              if gaussianSigma<=0
                imagesc(c2{1}([1 end]),c2{2}([1 end]),N2');     % c2 - pixel centres, N - pixel values (NOTE: TRANSPOSE)
              else
                imagesc(c2{1}([1 end]),c2{2}([1 end]),imgaussfilt(N2',gaussianSigma));
              end
              axis xy % ensure y axis points up
              colorbar
              title('Down to Up')
              xlim([ min(xyTrajTrue(:,1)), max(xyTrajTrue(:,1)) ]); ylim([ min(xyTrajTrue(:,2)), max(xyTrajTrue(:,2)) ]);
subplot(2,2,4)
[N2,c2] = hist3([[xyTrajTrue(upToDown,1);0;max(xyTrajTrue(:,1))], ...   % Trajectoriess Histogram - i.e. how often is the trajectory in the same bin as we used for spikes
                 [xyTrajTrue(upToDown,2);0;max(xyTrajTrue(:,2))]], ...
                 numBins); N2(1,1)=N2(1,1)-1; N2(end,end) = N2(end,end)-1; % subtract elements we added
              if gaussianSigma<=0
                imagesc(c2{1}([1 end]),c2{2}([1 end]),N2');     % c2 - pixel centres, N - pixel values (NOTE: TRANSPOSE)
              else
                imagesc(c2{1}([1 end]),c2{2}([1 end]),imgaussfilt(N2',gaussianSigma));
              end
              axis xy % ensure y axis points up
              colorbar
              title('Up to Down')
              xlim([ min(xyTrajTrue(:,1)), max(xyTrajTrue(:,1)) ]); ylim([ min(xyTrajTrue(:,2)), max(xyTrajTrue(:,2)) ]);