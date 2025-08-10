function [] = polarwitherrorbarRan(theta,rho,thetaError,rhoError,color)
%% modified from polarwitherrorbar 
% Ran 20210223
switch color
    case 1
        colormap = [1 0 0];
    case 2
        colormap = [0 1 0];    
    case 3
        colormap = [0 0 1]; 
end
fake = polar(theta,max(rho+rhoError)); set(fake,'Visible','off'); hold on; 
h=polar(theta,rho,'-sb');
set(h,'color',colormap);
% polar(theta*ones(1,3),[0, rho(ni), rho(ni)+error(ni)],'-r'); 
h=polar(theta*ones(1,3),[rho-rhoError, rho, rho+rhoError],'-g'); % rho error
set(h,'color',colormap);
temp = theta-thetaError:0.01:theta+thetaError;
h=polar(temp,rho*ones(1,length(temp)),'-g');
set(h,'color',colormap);
end