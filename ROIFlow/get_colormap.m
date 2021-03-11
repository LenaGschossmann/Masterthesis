function [cmap] = get_colormap(r, y, g, nelem)
% Create colormap
cmap=[linspace(g(1),y(1),nelem);linspace(g(2),y(2),nelem);linspace(g(3),y(3),nelem)];
cmap(:,end+1:end+nelem)=[linspace(y(1),r(1),nelem);linspace(y(2),r(2),nelem);linspace(y(3),r(3),nelem)];
end