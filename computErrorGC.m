function [ rmsErr ] = computErrorGC( frontFlow, depthMap, minLabel, labels )
% 24 August 2016 Dual Flea Images Sandstone
% Run depthFromDualFocalLengthFlowDriver_GC
% Call computErrorGC(vl, depthMap, minLabel, labels )

left = 50;
right = 640-50;
top = 50;
bottom = 480-50;

fl = 7.8545/.006;      %fl scaled by same factor image was scaled at

depthMapXY = depthMap;

imshow(depthMapXY/length(labels));


hold on
plot(left:right,top:top,'r')
plot(right:right,top:bottom,'r')
plot(left:right,bottom:bottom,'r')
plot(left:left,top:bottom,'r')
hold off;

stepSize = 1;

depth = [];
flow = [];
Vel = [];
if stepSize == 1
    index = 1;
    for i = left:right
        for j = top:bottom
            depth(index) = ((depthMapXY(j,i) - 1) + minLabel)*10;
            flow(index) = frontFlow(j,i);
            Vel(index) = flow(index).*depth(index)/fl ;
            index = index + 1;
        end
    end
else
    index = 1;
    for i = left:right
        for j = top:bottom
            depth(index) = ((depthMapXY(j,i) - 1)*10 + minLabel)*10;
            flow(index) = frontFlow(j,i);
            Vel(index) = flow(index).*depth(index)/fl;
            index = index + 1;
        end
    end
end
    

vMean = 20;
rmsErr = sqrt(mean((Vel-vMean).*conj(Vel-vMean)));
disp(strcat('RMS Error:',num2str((rmsErr/20)*100),'%'))

vMean = mean(Vel);
rmsErr = sqrt(mean((Vel-vMean).*conj(Vel-vMean)));
disp(strcat('RMS Error:',num2str((rmsErr/20)*100),'%'))

end

