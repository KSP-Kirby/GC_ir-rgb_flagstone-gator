function [zEst0, deltaXest] = centerEstimate(w0_left, w0_right, params )
% version 15 October 2015
% this only deals with the first estimate where there is no delta Z
    
    tic
    [rowsL, colsL] = size(w0_left);
    [rowsR, colsR] = size(w0_right);
    kernelWidth = 6;            % for front camera, TODO: put this in params

    kernel=ones(1,kernelWidth)/kernelWidth;
    
    % Add padding of one smoothing kernel width
    w0_left = [ones(1,kernelWidth)*w0_left(rowsL/2,1),w0_left(rowsL/2,:),ones(1,kernelWidth)*w0_left(rowsL/2,end)];
    w0_right = [ones(1,kernelWidth)*w0_right(rowsR/2,1),w0_right(rowsR/2,:),ones(1,kernelWidth)*w0_right(rowsR/2,end)];
    
    %smooth
    w0_left=conv(w0_left,kernel,'same'); 
    w0_right=conv(w0_right,kernel,'same'); 
    
    %remove padding
    w0_left = w0_left(kernelWidth+1:end-kernelWidth);
    w0_right = w0_right(kernelWidth+1:end-kernelWidth);
    
    [numRows, numCols] = size(w0_left);
    centerCol = round(numCols/2);
    centerRow = round(numRows/2);
    if w0_left(centerCol) < 0
        w0_left = -w0_left;
        w0_right = -w0_right;
    end
    
    % find initial offset
    plot(w0_left)
    hold all
    plot(w0_right)
    
    for Z = params.minZ:10:params.maxZ
        m = (params.pixelDimL/params.pixelDimR)*(params.f_r/params.f_l)*(Z+params.d)/(Z);
        x_r = 0;
        
        X_r = Z*x_r/params.f_r;
        X_l = X_r+params.b;
        x_l = params.f_l*X_l/(Z+params.d);
        
        leftPixel = x_l/params.pixelDimL + colsL/2; 
        plot(round(leftPixel), w0_right(320)/m, '*r')
        

        flowLeft = w0_right(320)/m;
        plot(round(leftPixel),flowLeft,'*')
        
        crossingSign = sign(flowLeft - w0_left(centerRow,round(leftPixel)));
        if Z ~= params.minZ
            if lastCrossingSign ~= crossingSign;
                break;
            end  
        end
        lastCrossingSign = crossingSign;
            
%        title(num2str(crossingSign))
  
    end

    deltaXest = Z*w0_right(320)/params.f_r;
    
    zEst0 = (deltaXest*params.f_r)./w0_right;


end

