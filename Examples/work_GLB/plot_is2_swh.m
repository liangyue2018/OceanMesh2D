function TT = plot_is2_swh(datestr, version, save_flag)
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
ms = 5;
limit = [0 round(prctile(TT.swh, 99))];
cmap = 'jet';
clbl = '$H_{s}\ \mathrm{[m]}$';
xtks = -180:60:180;
ytks = -60:30:66;
lonlim = [min(xtks) max(xtks)];
latlim = [min(ytks) max(ytks)];

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
geoscatter(gx, TT.Lat, TT.Lon, ms, TT.swh, "filled", MarkerFaceAlpha=.25);
geoscatter(gx, TT.Lat(idx), TT.Lon(idx), 15*ms, TT.swh(idx), 'p', LineWidth=1.5, MarkerEdgeAlpha=1);
geobasemap(gx, 'grayland');
clim(gx, limit);
colormap(gx, cmap);
c = colorbar(gx);
ylabel(c, clbl, Interpreter='latex');

% set ticks
gx.LongitudeAxis.TickValues = xtks;
gx.LatitudeAxis.TickValues = ytks;
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