% computes flow using brox
% 2016_10_29 data, IR-RGB-Flagstone images
% requires that RGB images be converted to GS first


addpath 'C:\Users\Richard\Documents\MATLAB\pami2010Matlab';
addpath 'C:\Users\Richard\Documents\MATLAB\flow-code-matlab';

h = waitbar(0,'Processing Data')
images = [1, 2, 3, 4, 5];

vl = 'vl';
ir = 'ir';

%extn = '.tif'; 
extn = '_rect.tif';
framesBetweenFlowComp = 4;

    for j = 1:length(images)-framesBetweenFlowComp
        i = images(j);
        i1 = images(j+framesBetweenFlowComp);
        if i < 10
            vl_1=imread(strcat(vl,'0',num2str(i),extn));
            ir_1=(imread(strcat(ir,'0',num2str(i),extn)));
            titlef1 = strcat(vl,'0',num2str(i));
            titleb1 = strcat(ir,'0',num2str(i));
        else
            vl=imread(strcat(vl,num2str(i),extn));
            ir=imread(strcat(ir,num2str(i),extn));
            titlef1 = strcat(vl,num2str(i));
            titleb1 = strcat(ir,num2str(i));
        end

        if i1 < 10
            vl_2=imread(strcat(vl,'0',num2str(i1),extn));
            ir_2=imread(strcat(ir,'0',num2str(i1),extn));
            titlef2 = strcat(vl,'0',num2str(i1));
            titleb2 = strcat(ir,'0',num2str(i1));
        else
            vl_2=imread(strcat(vl,num2str(i1),extn));
            ir_2=imread(strcat(ir,num2str(i1),extn));
            titlef2 = strcat(vl,num2str(i1));
            titleb2 = strcat(ir,num2str(i1));
        end
        
        if (strcmp(extn,'_rect.tif'))
            vl_1 = convert1Dto3D(vl_1);
            vl_2 = convert1Dto3D(vl_2);
            ir_1 = convert1Dto3D(ir_1);
            ir_2 = convert1Dto3D(ir_2);
        end

        plotTitle = {strcat('From 9 Aug 2016 Images:', titlef1,'-',titlef2,'<<->>', titleb1, '-', titleb2);...
            strcat('method: Brox')};
        
        uv_vl{j} = mex_LDOF(double(vl_1), double(vl_2));
        waitbar((j-.5)/length(images))

        uvi_vl{j} = uint8(flowToColor(uv_vl{j}));

        uv_ir{j} = mex_LDOF(double(ir_1), double(ir_2));
        uvi_ir{j} = uint8(flowToColor(uv_ir{j}));

        horzLine_l = 240;
        horzLine_r = 512;
        curFigure = figure;
        plot(uv_vl{j}(horzLine_l,:,1),'LineWidth',2)
        hold all
        plot(uv_ir{j}(horzLine_r,:,1),'LineWidth',2)
        title(plotTitle);
        legend('RGB', 'IR')
        ylabel('pixels')
        xlabel('Flow in pixels along center horizontal line')
        %axis([0,450,0,20])
        hold off
        filename = strcat('brox',num2str(images(j)),'.jpg');
        saveas(curFigure,filename)
        filename = strcat('brox',num2str(images(j)),'.fig');
        saveas(curFigure,filename)

        waitbar((j)/length(images))
    end
    
%end

close(h)

save('flow','uv_ir','uv_vl','uvi_ir','uvi_vl')
clear 'uv_ir'
clear 'uv_vl'
clear 'uvi_ir'
clear 'uvi_vl'
close all
