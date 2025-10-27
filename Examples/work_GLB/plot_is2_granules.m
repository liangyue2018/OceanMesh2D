function plot_is2_granules(datestr, save_flag)
% Visualize the locations of ICESat-2 ATL03 granules over coastal buffer zones.
%

% Check input:
if nargin < 2
    save_flag = false;
end
assert(ischar(datestr), 'Error: Invalid date string.');

% load data
load('./coast_polygon_c.mat', 'shape', 'shapeout');
load("./S_" + datestr + ".mat",'S');

% set figure
fig = figure(Units="inches");
fig.Position(3) = 10;
fig.Position(4) = 5;
fsz = 12;
latlim = [-60 66];
lonlim = [-180 180];

% plot
gx = geoaxes;
hold on;
geoplot(gx, shape, FaceColor='none', ...
                   FaceAlpha=.25, ...
                   EdgeColor='k', ...
                   EdgeAlpha=.25);
geoplot(gx, shapeout, FaceColor='b', ...
                      FaceAlpha=.5, ...
                      EdgeColor='k', ...
                      EdgeAlpha=.5);
for n = 1:numel(S)
    orb = S(n).orbit_info;
    geoplot(gx, orb.bnd_poly, FaceColor='none', ...
                              FaceAlpha=.5, ...
                              EdgeColor='g', ...
                              EdgeAlpha=.5, ...
                              LineWidth=2);
    geoplot(gx, S(n).shape, 'r.', MarkerSize=5);
end
geobasemap(gx, 'grayland');

% set ticks
%gx.LatitudeAxis.TickValues = [];
%gx.LongitudeAxis.TickValues = [];
geotickformat(gx, 'dd');
geolimits(gx, latlim, lonlim);
gx.LatitudeAxis.FontSize = fsz;
gx.LongitudeAxis.FontSize = fsz;
gx.LatitudeLabel.String = '';
gx.LongitudeLabel.String = '';
gx.Scalebar.Visible = 'off';

if save_flag
    sfnm = ['1_is2_granules_' datestr '.png'];
    exportgraphics(fig, sfnm, Resolution=900);
end

end