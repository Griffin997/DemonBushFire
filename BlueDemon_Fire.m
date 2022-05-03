%%% Lil' Blue Demons Tryin to Put Out a Grassland Fire
clc;
clear;
close all;

%Side length of grassland
sL = 20;
mF = 5000;
fireProb = 0.666;

%The list of all bushes on fire
fireList = [];
acctFire = []; %Fire's that have been accounted for

bushMoments = 0;
demonSteps = 0;

GL_pos_color = zeros(sL^2,5);

%Set number of demons
dN = 6;
demonStat = zeros(dN,6); %x, y, moving, Dx, Dy
xRand = sL*rand(dN,1);%randi([1 sL],dN,1);
yRand = sL*rand(dN,1);%randi([1 sL],dN,1);
demonStat(:,1) = xRand;
demonStat(:,2) = yRand;

v = VideoWriter('DemonCollective.avi');
    
v.FrameRate = 3;

open(v);


for i = 1:sL
    
    for j = 1:sL
        
        row = (i-1)*sL + j;
        
        GL_pos_color(row,1) = i; %x position
        GL_pos_color(row,2) = j; %y position
        GL_pos_color(row,3) = 0; %on fire status
        
    end
end


%Initial Figure
figure('Renderer', 'painters', 'Position', [50 50 1400 900]);
scatter(GL_pos_color(:,1),GL_pos_color(:,2),floor(mF/sL),'o','markerfacecolor',[0,0.5,0]...
    ,'markeredgecolor',[0.1,0.1,0.1]);
hold on
scatter(demonStat(:,1),demonStat(:,2),floor(mF/sL)*2,'o','markerfacecolor',[0.2,0,0.8]...
    ,'markeredgecolor',[0.1,0.1,0.1]);
set(gca,'color',[0.6 1 0.6])
xlim([0,sL+1])
ylim([0,sL+1])
set(gca,'XColor', 'none','YColor','none')
title('Initial Bush Field')
set(gca,'FontSiz',18)
hold off

frame = getframe(gcf);
writeVideo(v,frame);
writeVideo(v,frame);

% pause(1);


for ii = 1:666
    
    bushMoments = bushMoments + length(fireList);
    
    oldFire = [];
    newFire = [];
    oldFire = fireList;
    fireList = setOnFire(sL,fireList,floor(sL*1.33),fireProb);
    newFire = setdiff(fireList,oldFire);
    
%     if ~isempty(fireList)
%            dViable = find(demonStat(:,3)==0);
%            
%            fireAvail = setdiff(fireList,acctFire);
%            
%            fPosit = GL_pos_color(fireAvail,1:2);
%            
%            for k = 1:min(length(dViable),length(fireAvail))
%                
%                dOcc = find(demonStat(:,3)==1);
%                OccFire = demonStat(dOcc,6);
%                OccFire = OccFire(OccFire~=0);
%                fireOpt = setdiff(fireAvail,OccFire);
%                fPosit = GL_pos_color(fireOpt,1:2);
%                
%                dV = dViable(k);
%                
%                x1 = demonStat(dV,1);
%                y1 = demonStat(dV,2);
%                
%                
%                x2 = fPosit(:,1);
%                y2 = fPosit(:,2);
%                
%                [distMat,dX,dY] = calcD(x1,y1,x2,y2);
%                
%                [~,closeF] = min(distMat);
%                
%                fA = fireOpt(closeF);
%                
%                acctFire = [acctFire fA];
%                
%                demonStat(dV,3) = 1;
%                demonStat(dV,4) = dX(closeF)/distMat(closeF);
%                demonStat(dV,5) = dY(closeF)/distMat(closeF);
%                demonStat(dV,6) = fireOpt(closeF);
%               
%            end
%            
%            
%            
%     end
    
    
    
    %Replotting Figure
    scatter(GL_pos_color(:,1),GL_pos_color(:,2),floor(mF/sL),'o','markerfacecolor',[0,0.5,0]...
        ,'markeredgecolor',[0.1,0.1,0.1]);
    hold on
    scatter(GL_pos_color(fireList,1),GL_pos_color(fireList,2),floor(mF/sL),'o','markerfacecolor',[0.7,0,0]...
        ,'markeredgecolor',[0.1,0.1,0.1]);
    scatter(demonStat(:,1),demonStat(:,2),floor(mF/sL)*2,'o','markerfacecolor',[0.2,0,0.8]...
        ,'markeredgecolor',[0.1,0.1,0.1]);
    scatter(GL_pos_color(newFire,1),GL_pos_color(newFire,2),floor(mF/sL*3),'s','markerfacecolor',[0.91, 0.41, 0.17]...
                ,'markeredgecolor',[0.1,0.1,0.1]);
    xlim([0,sL+1])
    ylim([0,sL+1])
    set(gca,'color',[0.6 1 0.6])
    title({strcat('Iteration:',string(ii),' - Bushes on fire:', ...
        string(length(fireList)));strcat('Bush Fire Moments:',string(bushMoments) ...
        ,' - Demon Steps:',string(demonSteps)  )});
    set(gca,'XColor', 'none','YColor','none')
    set(gca,'FontSiz',18)
    hold off;
    
    frame = getframe(gcf);
    writeVideo(v,frame);
    
%     pause(0.5);
    
    for jj = 1:6
        
        if ~isempty(fireList)
            [demonStat,acctFire] = rerouteDemon(demonStat,fireList,GL_pos_color,acctFire);
        end
        
        demonStat = checkStats(demonStat);        
        
        demonStat(:,1:2) = demonStat(:,1:2) + demonStat(:,4:5); %Have to figure out a better way to do this so that they go effectively...
        
        moveMat = (demonStat(:,4)~=0) + (demonStat(:,5)~=0);
        moveMat = moveMat>0;
        
        demonSteps = demonSteps + sum(moveMat);
        
        scatter(GL_pos_color(:,1),GL_pos_color(:,2),floor(mF/sL),'o','markerfacecolor',[0,0.5,0]...
            ,'markeredgecolor',[0.1,0.1,0.1]);
        hold on
        scatter(GL_pos_color(fireList,1),GL_pos_color(fireList,2),floor(mF/sL),'o','markerfacecolor',[0.7,0,0]...
            ,'markeredgecolor',[0.1,0.1,0.1]);
        scatter(demonStat(:,1),demonStat(:,2),floor(mF/sL)*2,'o','markerfacecolor',[0.2,0,0.8]...
            ,'markeredgecolor',[0.1,0.1,0.1]);
        xlim([0,sL+1])
        ylim([0,sL+1])
        set(gca,'color',[0.6 1 0.6])
        title({strcat('Iteration:',string(ii),' - Bushes on fire:', ...
            string(length(fireList)));strcat('Bush Fire Moments:',string(bushMoments) ...
            ,' - Demon Steps:',string(demonSteps)  )});
        set(gca,'FontSiz',18)
        
        [fireList,acctFire,demonStat] = putOutFire(fireList,acctFire,demonStat,GL_pos_color,mF,sL);
        
        set(gca,'XColor', 'none','YColor','none')
        hold off;
        frame = getframe(gcf);
        writeVideo(v,frame);
        
%         pause(0.5);
        
        
    end

%     pause(2);
    
    
end

scatter(GL_pos_color(:,1),GL_pos_color(:,2),floor(mF/sL),'o','markerfacecolor',[0,0.5,0]...
    ,'markeredgecolor',[0.1,0.1,0.1]);
hold on
scatter(GL_pos_color(fireList,1),GL_pos_color(fireList,2),floor(mF/sL),'o','markerfacecolor',[0.7,0,0]...
    ,'markeredgecolor',[0.1,0.1,0.1]);
scatter(demonStat(:,1),demonStat(:,2),floor(mF/sL)*2,'o','markerfacecolor',[0.2,0,0.8]...
    ,'markeredgecolor',[0.1,0.1,0.1]);
xlim([0,sL+1])
ylim([0,sL+1])
set(gca,'color',[0.6 1 0.6])
title({strcat('Iteration:',string(ii),' - Bushes on fire:', ...
    string(length(fireList)));strcat('Bush Fire Moments:',string(bushMoments) ...
    ,' - Demon Steps:',string(demonSteps)  )});
set(gca,'XColor', 'none','YColor','none')
set(gca,'FontSiz',18)
hold off

frame = getframe(gcf);
writeVideo(v,frame);

close(v);


function fireList = setOnFire(sL,fireList,pos,p)
    randP = rand(1,pos);
    succ = (randP<p);
    count_succ = sum(succ);
    newFire = randi([1 sL^2],1,count_succ);
    fireList = union(fireList, newFire);
end

function [dist,Dx,Dy] = calcD(x1,y1,x2,y2)

    Dx = x2 - x1;
    
    Dy = y2 - y1;
    
    dist = (Dx.^2 + Dy.^2).^(1/2);

end

function [fHL,fHA,dS] = putOutFire(fL,fA,dS,gL,mF,sL)

    dViable = find(dS(:,3)==1);
           
%     fPosit = gL(fA,1:2);
    
    
    fHL = fL;
        
    fHA = fA;
           

    for i = 1:length(dViable)
        
        dV = dViable(i);
        
        x1 = dS(dV,1);
        y1 = dS(dV,2);
        
        fF = dS(dV,6);

        x2 = gL(dS(dV,6),1);
        y2 = gL(dS(dV,6),2);
        
        [dist,~,~] = calcD(x1,y1,x2,y2);
        
%         fFind = find(dist<0.5);
        
        
        
        if dist<1
            
%             fF = fA(fFind);
            
            fHL = setdiff(fHL,fF);
            fHA = setdiff(fHA,fF);
            
            dS(dV,3) = 0;
            dS(dV,4) = 0;
            dS(dV,5) = 0;
            
            scatter(gL(fF,1),gL(fF,2),floor(mF/sL*3),'s','markerfacecolor',[0.4,0.9,0.9]...
                ,'markeredgecolor',[0.1,0.1,0.1]);
            
        end
        
    end

end

function dS = checkStats(dS)

    dOcc = find(dS(:,3)==1);
    fireFocus = dS(dOcc,6);
    if numel(fireFocus) ~= numel(unique(fireFocus))
        error('Overlap Detected');
    end

end

function [dS,aF] = rerouteDemon(dS,fL,gL,aF);

    %dS = demonStats
    %fL = fireList
    %gL = GL_pos_color
    %aF = acctFire

    dViable = find(dS(:,3)==0);
    
    fA = setdiff(fL,aF);

    for k = 1:min(length(dViable),length(fA))

       dOcc = find(dS(:,3)==1);
       OccFire = dS(dOcc,6);
       OccFire = OccFire(OccFire~=0);
       fireOpt = setdiff(fA,OccFire);
       fPosit = gL(fireOpt,1:2);

       dV = dViable(k);

       x1 = dS(dV,1);
       y1 = dS(dV,2);

       x2 = fPosit(:,1);
       y2 = fPosit(:,2);

       [distMat,dX,dY] = calcD(x1,y1,x2,y2);

       [~,closeF] = min(distMat);

       fAcc = fireOpt(closeF);

       aF = [aF fAcc];

       dS(dV,3) = 1;
       dS(dV,4) = dX(closeF)/distMat(closeF);
       dS(dV,5) = dY(closeF)/distMat(closeF);
       dS(dV,6) = fireOpt(closeF);

    end 
    
end