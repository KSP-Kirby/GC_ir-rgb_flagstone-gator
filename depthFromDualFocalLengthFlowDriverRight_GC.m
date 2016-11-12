% version 24 August 2016
% for dual flea images
% used calibration from 24 Aug 2016
% uses right image as reference


addpath('C:\Users\Richard\Documents\MATLAB\flow-code-matlab')     % this is flow to color
addpath('C:\Users\Richard\Documents\MATLAB\gco-v3.0\matlab')      % Graph Cuts
%load('flow_flagstone.mat')
load('flow_rect.mat')
imageSet = 1;

usePreviousDataCost = 0;

scaleImage = 1;
vl = -uv_vl{imageSet}(:,:,1);
ir = -uv_ir{imageSet}(:,:,1);


params.f_l = 3.9676;        % left camera is ir
params.f_r = 7.8545;         % right camera is RGB
params.b =  3*25.4;         % stereo baseline
params.d =  78.4954;         % dual focal length baseline
params.minZ = 600;
params.maxZ = 3000;
params.pixelDimL = .0048;
params.pixelDimR = .006;
startSeq = 1;
[zEst0, deltaXest] = centerEstimate(ir, vl, params);

[m,n] = size(vl);
for i = 1:m
    zEst(i,:) = (deltaXest*params.f_r)./vl(i,:);
end
                
initCostFactor = .005;

for smoothingFactor = [21]                                        % this smooths between depth labels 
    for dataFactor = [90]                                        % for nice smooth along a single line, 15, 400, 40, 0 works well
        for horzNeighborMaskWeight = [5]
            for verticalNeighborMaskWeight = [2]
                %smoothingFactor = 10;    % higher is smoother, zero = no smoothing
                %dataFactor = 90;
                %horzNeighborMaskWeight = 5;
                %verticalNeighborMaskWeight = 2;


                minLabel = 70;
                maxLabel = 150;
                f_l = 3.9676;        % left camera is ir
                f_r = 7.8545;         % right camera is RGB
                b =  25.4*3;         % stereo baseline
                d =  78.4954;
                
                pixelDimR = .006;
                pixelDimL = .0048;  

                [numRowsR,numColsR] = size(vl);
                [numRowsL,numColsL] = size(ir);
                optimalLabellingOut = [];

                numSites = numRowsR*numColsR;

                % labels are in mm
                % min label for the stereo camera setup is about 850
                % max label is about 2500
                labels = (minLabel:1:maxLabel)*10;  %*10 converts from cm to mm
                numLabels = length(labels);
                h = GCO_Create(numSites,numLabels);

                if usePreviousDataCost == 1
                    load('dataCost.mat')
                else
                    dataCost = zeros(numLabels, numSites);

                    site = 1;

                    h1 = waitbar(0, 'Constructing data cost matrix');
                    for k = 1:numRowsR
                        for i = 1:numColsR
                            for j = 1: numLabels
                                Z = labels(j);
                                %if i + disparity <= numCols
                                    m = (pixelDimL/pixelDimR)*(f_r/f_l)*(Z+d)/(Z);      % Z is for the right camera
                                    rightPixel_1 = i;
                                    x_r = (rightPixel_1-numColsR/2)*pixelDimR;
                                    
                                    X_r = Z*x_r/f_r;
                                    X_l = X_r+b;
                                    x_l1 = f_l*X_l/(Z+d);
        
                                    x_l = (x_r*Z*f_l+b*f_r*f_l)/(f_r*Z+f_r*d);
                                    leftPixel_1 = numColsL/2+x_l/pixelDimL;        % this is a fractional pixel
                                    %rp(j) = rightPixel_1;
                                    if (leftPixel_1 < 1)
                                        dataCost(j, site) = 10;
                                    elseif (leftPixel_1 > numColsL)
                                        dataCost(j, site) = 10;
                                    else
                                        %wb_intrp = wb(k, floor(scaledPixelLocation)) + ((wb(k,ceil(scaledPixelLocation))-wb(k,floor(scaledPixelLocation))) * (scaledPixelLocation-floor(scaledPixelLocation)));                      % this interpolates between the two interger wb values
                                        %wb_intrp = wb(k,round(scaledPixelLocation));
                                        k_l = round((k - numRowsR/2)/m + numRowsL/2 - 16);
                                        dataCostInit = initCostFactor*abs(Z - zEst(k, rightPixel_1));
                                        dataCostMatch = abs(m*ir(k_l, round(leftPixel_1)) - vl(k,rightPixel_1));
                                        dataCost(j, site) = abs(m*ir(k_l, round(leftPixel_1)) - vl(k,rightPixel_1)) + dataCostInit;

                                    end
                                
                                
                                
%                                     % This is if left image is reference
%                                     m = (f_r/f_l)*(Z+d)/(Z);
%                                     leftPixel_1 = i;
%                                     x_l = (leftPixel_1-numCols/2)*pixelDim;
%                                     x_r = (-b*f_r*f_l+f_r*x_l*Z+f_r*x_l*d)/(Z*f_l);
%                                     rightPixel_1 = numCols/2+x_r/pixelDim;        % this is a fractional pixel
%                                     if (rightPixel_1 < 1)
%                                         dataCost(j, site) = 10;
%                                     elseif (rightPixel_1 > numCols)
%                                         dataCost(j, site) = 10;
%                                     else
%                                         dataCost(j, site) = abs(m*ir(k,leftPixel_1) - vl(k,round(rightPixel_1)));

                                    
                                %else
                                %    dataCost(j,site) = NaN;
                                %end   
                            end
                            site = site + 1;
                        end
                        waitbar(k/numRowsR)
                    end
                    close(h1)
                    
%                     figure
%                     plot(labels, energy)
%                     figure         
%                     plot(vl)
%                     hold all
%                     plot(ir)
%                     [val,index] = min(energy);
%                     plot(leftPixel_1,ir(leftPixel_1), '*r')
%                     plot(rp(index),vl(round(rp(index))), '*r')
%                     hold off

                    minDataCost = min(min(dataCost));
                    maxDataCost = max(max(dataCost));

                    %scale data cost between 1 and 1000

                    scaleFactor =50/maxDataCost;
                    dataCost = dataFactor*scaleFactor*dataCost+1;
                end


                GCO_SetDataCost(h,cast(dataCost,'int32'));

                % Smooth cost is number of labels by number of labels
                smoothCost = zeros(numLabels, numLabels);
                for i = 1:numLabels
                    smoothCost(i,:) = abs((1:1:numLabels) - i);
                end

                smoothCost = smoothCost*smoothingFactor;

                GCO_SetSmoothCost(h,cast(smoothCost, 'int32'));

                neighbors = sparse(numSites,numSites);

                colCount = 1;
                for i = 1:numSites-1
                    if i ~=numColsR
                        neighbors(i,i+1) = horzNeighborMaskWeight;      % must be upper triangle
                        colCount = colCount + 1;
                    else
                        colCount = 1;       % this prevents wrapping from the end of one roll to the start of the next
                    end
                end

                % Epipolar line to radial line neighbors
                for i = 1:numSites-numColsR - 1
                    neighbors(i,i+numColsR) = verticalNeighborMaskWeight;
                end

                GCO_SetNeighbors(h,neighbors);

                GCO_Expansion(h);

                [E, D, S] = GCO_ComputeEnergy(h); 

                optimalLabelling = cast(GCO_GetLabeling(h),'double');

                %figure
                %plot(optimalLabelling);

                GCO_Delete(h);
                
%                 figure
%                 plot(optimalLabelling)
%                 depthMap = optimalLabelling*10 + minLabel*10;
%                 m = (f_r/f_l)*(depthMap+d)./(depthMap);
%                 leftPixel_1 = 1:1:640;
%                 x_l = (leftPixel_1-320)*pixelDim;
%                 x_r = (-b*f_r*f_l+f_r*x_l.*depthMap'+f_r*x_l*d)./(depthMap'*f_l);
%                 rightPixel_1 = 320+x_r/pixelDim;
%                 
%                 for i = 1:640
%                     if (rightPixel_1(i) >= 1) && (rightPixel_1(i) <=640)
%                         vl_adj(i) = (1/m(i))*vl(round(rightPixel_1(i)));
%                     end
%                 end
%                 
%                 plot(vl)
%                 hold all
%                 plot(ir)
%                 plot(vl_adj)
%                 hold off
%                 
%                 title('Flow along epipolar line adjusted by estimated depth')
%                 legend('Flow1','Flow2', '(1/m)*Flow1(x+h)')
%                 xlabel('Pixel')
%                 ylabel('Flow (pixels/frame)')
%                 
%                 figure
%                 plot(vl)
%                 hold all
%                 plot(ir)    
% 
%                 leftPixel_1 = 340;
%                 plot(leftPixel_1,ir(leftPixel_1),'*')
%                 x_l = (leftPixel_1-320)*pixelDim;
%                 for i = 1:length(labels)
%                     z = labels(i);
%                     x_r = (-b*f_r*f_l+f_r*x_l*z+f_r*x_l*d)/(z*f_l);
%                     rightPixel_1 = 320+x_r/pixelDim;
%                     m = (f_r/f_l)*(z+d)/(z);
%                     v1_flow(i) = (m)*ir(leftPixel_1);
%                     pix(i) = rightPixel_1;
%                     %plot(rightPixel_1,v1_flow, '*r')
%                 end
%                 plot(pix,v1_flow)
% 
%                 axis([180,405,5,35])
%                 title('Energy minimization search path')
%                 xlabel('Pixels')
%                 ylabel('Flow (pixels/frame)')
%                 legend('Flow 1', 'Flow 2','Eval Point', 'Eval Energy')

                depthMap = reshape(optimalLabelling,numColsR,length(optimalLabelling)/numColsR);

                depthMap = depthMap';
                %imtool(depthMap/241)
                figure
                imshow(depthMap/length(labels))
                %imageSizer
                filename=strcat('depthMap_',num2str(smoothingFactor),'_',num2str(dataFactor),'_',num2str(horzNeighborMaskWeight),'_',num2str(verticalNeighborMaskWeight));
                save(filename,'depthMap');  
            end
        end
    end
end
