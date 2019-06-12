%% Plots planets
function [h] = plotPlanet(planetName,solarSystem,gridCheck)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% - planetName   name of desired planet [string]
%
% OUTPUTS
% -              plot of desired planet [figure]
%
% Author - Vivek Desai
%
% Last Update - 03/29/2019
%
% References
% - original author: Ryan Gray
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if gridCheck == 1
    h = figure('Color','k','units','normalized','outerposition',[.5 .2 .5 .8]);
    hold on;
    rotate3d on;
    set(gca, 'NextPlot','add', 'Visible','off');
    grid on
else
    h = figure('Color', 'k');
    hold on;
    rotate3d on;
    set(gca, 'NextPlot','add', 'Visible','off');
    axis equal;
    axis auto;
end

axis vis3d;
npanels = 180;   % Number of globe panels around the equator deg/panel = 360/npanels
alpha   = 0.8; % globe transparency level, 1 = opaque, through 0 = invisible

if planetName(1:3) == 'Mer'
    erad    = 2439.7e3; % equatorial radius (meters)
    prad    = 2439.7e3; % polar radius (meters)
    R = 0.387;
elseif planetName(1:3) == 'Ven'
    erad    = 6051.8e3; % equatorial radius (meters)
    prad    = 6051.8e3; % polar radius (meters)
    R = 0.722;
elseif planetName(1:3) == 'Ear'
    erad    = 6378.137; % equatorial radius (kmeters)
    prad    = 6356.752; % polar radius (kmeters)
    R = 1;
elseif planetName(1:3) == 'Mar'
    erad    = 3396.2e3; % equatorial radius (meters)
    prad    = 3376.2e3; % polar radius (meters)
    R = 1.52;
elseif planetName(1:3) == 'Jup'
    erad    = 71492e3; % equatorial radius (meters)
    prad    = 66854e3; % polar radius (meters)
    R = 5.2;
elseif planetName(1:3) == 'Sat'
    erad    = 60268e3; % equatorial radius (meters)
    prad    = 54364e3; % polar radius (meters)
    R = 9.58;
elseif planetName(1:3) == 'Ura'
    erad    = 25559e3; % equatorial radius (meters)
    prad    = 24973e3; % polar radius (meters)
    R = 19.2;
elseif planetName(1:3) == 'Nep'
    erad    = 24764e3; % equatorial radius (meters)
    prad    = 24341e3; % polar radius (meters)
    R = 30.1;
elseif planetName(1:3) == 'Sun'
    erad    = 695700e3;% equatorial radius (meters)
    prad    = 695700e3;% polar radius (meters)
    R = 0;
elseif planetName(1:3) == 'Moo'
    erad    = 1738.1e3;% equatorial radius (meters)
    prad    = 1736.0e3;% polar radius (meters)
    R = 1.002570038465126;
end

if solarSystem == 1
    AU = 149597870700;
    R = R*AU/10;
else
    R = 0;
end

[x, y, z] = ellipsoid(R, 0, 0, erad, erad, prad, npanels);
globe = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
fileName = strcat(planetName,'Map.jpg');
cdata = imread(fileName);
set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');