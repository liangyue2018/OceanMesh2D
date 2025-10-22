%[text] # 1 ICESat-2 ATL03 Data Access
%[text] ## 1.1 Get coastal polygons
% load shapefile
reso = 'c';
fnm = "D:\soft\OceanMesh2D\datasets\GSHHS_shp\" + reso + '\GSHHS_' + reso + '_L1.shp';
GT = readgeotable(fnm);

% subset
latlim = [-60 67];
lonlim = [-180 180];
shape = geoclip(GT.Shape, latlim, lonlim);
idx = shape.NumRegions ~= 0;
shape = shape(idx);

% filter with area
a = area(shape) / 1e6;% km2
shape = shape(a > 100);

% buffer
bufwidth = {3 .25}; % degree
direction = {'out' 'in'};
shp = cell(1,2);
for i = 1:2
    shp_b = buffer(shape, bufwidth{i}, Direction=direction{i});
    shp_b = geoclip(shp_b, latlim, lonlim);
    idx = shp_b.NumRegions ~= 0;
    shp_b = shp_b(idx);
    shp{i} = union(shp_b);
end
shapeout = subtract(shp{1}, shp{2});

% get lonb & latb
GT = table(shapeout, VariableNames={'Shape'});
T = geotable2table(GT, ["Lat","Lon"]);
[latb, lonb] = polyjoin(T.Lat, T.Lon);

% simplify
[latb, lonb, err, tol] = reducem(latb, lonb, .1);
fprintf('Simplify coastal polygon (Npts=%d) with tol=%.2f deg, err=%.1f%%.\n', ... %[output:group:3878cf5f] %[output:3ddc1177]
    numel(latb), tol, err*100); %[output:group:3878cf5f] %[output:3ddc1177]
shapeout = geopolyshape(latb, lonb);
shapeout.GeographicCRS = geocrs(4326);

bbox = {[-180 180;27 66]
        [-180 180;-27 27]
        [-180 180;-90 -27]};
shapeoutC = cell(3,1);
for i = 1:numel(bbox)
    bnd = bbox2poly(bbox{i});
    reg = geopolyshape(bnd(:,2), bnd(:,1));
    reg.GeographicCRS = geocrs(4326);
    shapeoutC{i} = intersect(shapeout, reg);
    %shapeoutC{i} = geoclip(shapeout, bbox{i}(2,:), bbox{i}(1,:));
end

% save and write shapefile
fnm = "./coast_polygon_" + reso;
delete(fnm + '*');
save(fnm,'shape','latb','lonb','shapeout','shapeoutC');
for i = 1:numel(shapeoutC)
    GT = table(shapeoutC{i}, VariableNames={'Shape'});
    sfnm = fnm + '_' + num2str(i);
    shapewrite(GT, sfnm);
    zip(sfnm, sfnm + '*');
end
%%
%[text] ## 1.2 Download ICESat-2 ATL03 data
% using python script `download_icesat2.py`
%%
%[text] ## 1.3 Read ICESat-2 ATL03 data
% define filepath (Gyyyy/MM/dd/)
datestr = '20181014';
filepath = "./G" + datestr(1:4) + '/' ...
           + datestr(5:6) + '/' ...
           + datestr(7:8) + '/';

% get filename
fList = dir(filepath + '*.h5');
filename = {fList.name}';
fnm = fullfile(filepath, filename);

% read h5file
gtx = {'gt1r' 'gt2r' 'gt3r'};
kw = {1 Inf 100};
S = struct;
for n = 1:numel(fnm) %[output:group:8cc1cbd2]
    fprintf('Processing file %d of %d: %s\n', n, numel(fnm), fnm{n}); %[output:61549d4c]
    for rn = gtx
        try
            lon_ph = h5read(fnm{n}, "/" + rn{1} + '/heights/lon_ph', kw{:});
            lat_ph = h5read(fnm{n}, "/" + rn{1} + '/heights/lat_ph', kw{:});
            S(n,1).(rn{1}).shape = geopointshape(lat_ph, lon_ph);
            S(n,1).(rn{1}).GeographicCRS = geocrs(4326);
        catch
            fprintf('    Error reading %s: %s\n', fnm{n}, rn{1});
        end
    end
end %[output:group:8cc1cbd2]

% save result
sfnm = ['is2_granules_' datestr '.mat'];
save(sfnm, 'S');
%%
%[text] ## 1.4 Visualization
% load data
datestr = '20181014';
load('coast_polygon_c.mat');
load(['is2_granules_' datestr '.mat'],'S');

fig = figure(Units="inches");
fig.Position(3) = 10;
fig.Position(4) = 8;

gx = geoaxes;
hold on;
geoplot(gx, shape, FaceColor='k', ...
                   FaceAlpha=.5, ...
                   EdgeColor='k', ...
                   EdgeAlpha=.5);
geoplot(gx, shapeout, FaceColor='b', ...
                      FaceAlpha=.5, ...
                      EdgeColor='k', ...
                      EdgeAlpha=.5);
for n = 1:numel(S)
    for rn = {'gt1r' 'gt2r' 'gt3r'}
        try
            geoplot(gx, S(n,1).(rn{1}).shape, 'r.');
        catch
        end
    end
end

%sfnm = ['is2_granules_' datestr '.png'];
%exportgraphics(fig, sfnm, Resolution=900);

%[appendix]{"version":"1.0"}
%---
%[metadata:styles]
%   data: {"code":{"fontSize":"12"},"heading1":{"fontSize":"18"},"heading2":{"fontSize":"14"},"heading3":{"fontSize":"12"},"normal":{"fontFamily":"Arial","fontSize":"12"},"title":{"fontSize":"24"}}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
%[output:3ddc1177]
%   data: {"dataType":"text","outputData":{"text":"Simplify coastal polygon (Npts=4067) with tol=0.10 deg, err=1.1%.\n","truncated":false}}
%---
%[output:61549d4c]
%   data: {"dataType":"text","outputData":{"text":"Processing file 1 of 3: .\\G2018\\10\\14\\ATL03_20181014001049_02350102_006_02_subsetted.h5\nProcessing file 2 of 3: .\\G2018\\10\\14\\ATL03_20181014013805_02360101_006_02_subsetted.h5\nProcessing file 3 of 3: .\\G2018\\10\\14\\ATL03_20181014054647_02380110_006_02_subsetted.h5\n","truncated":false}}
%---
