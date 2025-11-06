function plot_is2_swh(datestr, version, save_flag)
% Visualize the coastal SWH records derived from ICESat-2 altimetry.
%
% Check input:
if nargin < 2
    version = '006';
    save_flag = false;
end
assert(ischar(datestr) && length(datestr) == 6, 'Error: Invalid date string.');

% load timetable
fnm = sprintf('icesat2_atl03_%s_v%s.csv', datestr, version);
opts = detectImportOptions(fnm);
opts = setvartype(opts, {'beam' 'rcs'}, {'string' 'string'});
opts = setvaropts(opts, 'Time', TimeZone='UTC');
TT = readtimetable(fnm, opts);
loc = TT.missFrac == 0;
TT = TT(loc, :);
[~, idx] = max(TT.swh);

% load coastal polygons
load('coast_polygon_c.mat', 'shape', 'shapeout');

% set figure
fig = figure(Units="inches");
fig.Position(3) = 10;
fig.Position(4) = 5;
fsz = 12;
ms = 10;
latlim = [-60 66];
lonlim = [-180 180];

% plot
gx = geoaxes;
hold on;
geoplot(gx, shape, FaceColor='none', ...
                   FaceAlpha=.25, ...
                   EdgeColor='k', ...
                   EdgeAlpha=.25);
geoplot(gx, shapeout, FaceColor='none', ...
                      FaceAlpha=.5, ...
                      EdgeColor='k', ...
                      EdgeAlpha=.5);
geoscatter(gx, TT.Lat, TT.Lon, ms, TT.swh, "filled", MarkerFaceAlpha=.5);
geoscatter(gx, TT.Lat(idx), TT.Lon(idx), 5*ms, TT.swh(idx), 's', LineWidth=1.5, MarkerEdgeAlpha=.5);
colormap(gx, 'jet');
c = colorbar(gx);
ylabel(c, '$H_{s}\ \mathrm{[m]}$', Interpreter='latex');
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
    sfnm = sprintf('1_is2_swh_%s.png', datestr);
    exportgraphics(fig, sfnm, Resolution=900);
end

end