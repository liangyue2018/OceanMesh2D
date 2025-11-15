function shapeoutC = make_coast_polygon(reso, bbox, bufwidth, tol, varargin)
% Create coastal buffer polygons from GSHHS L1 data.
% Input:
%   reso: 'c' (crude), 'l' (low), 'i' (intermediate), 'h' (high) or 'f' (full)
%   bbox: [lonmin lonmax; latmin latmax]
%   bufwidth: buffer width for {outer inner} in degrees
%   tol: tolerance in degrees
%
% Optional inputs:
%  'min_area': minimum area (km^2) to keep land polygons (default: 100)
%  'latq', 'lonq': query points to remove land polygons (default: [])
%  'clip_lat': latitudes to clip geopolyshape into several parts (default: [])
%  'save_flag': whether to save output files (default: false)
%
% Output:
%   shape: coastline
%   latb, lonb: vector data of coastal buffer polygons
%   shapeout: coastal buffer polygons (geopolyshape)
%   shapeoutC: cell array of geopolyshapes if clipped by latitudes
%
% check and parse inputs
assert(ischar(reso) && ismember(reso, {'c','l','i','h','f'}), 'Error: Invalid reso.');
assert(isnumeric(bbox) && size(bbox,1) == 2 && size(bbox,2) == 2, 'Error: Invalid bbox.');
assert(iscell(bufwidth) && numel(bufwidth) == 2, 'Error: bufwidth should be a cell array with 2 elements.');
assert(isnumeric(tol) && tol > 0, 'Error: Invalid tol.');

p = inputParser;
addParameter(p, 'min_area', 100);
addParameter(p, 'latq', []);
addParameter(p, 'lonq', []);
addParameter(p, 'clip_lat', []);
addParameter(p, 'save_flag', false);
parse(p, varargin{:});
inp = p.Results;

% load shapefile
fnm = "D:\soft\OceanMesh2D\datasets\GSHHS_shp\" + reso + '\GSHHS_' + reso + '_L1.shp';
GT = readgeotable(fnm);

% subset
latlim = bbox(2,:);
lonlim = bbox(1,:);
shape = geoclip(GT.Shape, [latlim(1)-.5 latlim(2)+.5], lonlim);
idx = shape.NumRegions ~= 0;
shape = shape(idx);

% filter with area
a = area(shape) / 1e6;% km2
shape = shape(a > inp.min_area);

% buffer
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

% remove some land polygons by query points
if ~isempty(inp.latq) && ~isempty(inp.lonq)
    latq = inp.latq(:);
    lonq = inp.lonq(:);
    [latcells, loncells] = polysplit(latb, lonb);
    ploc = true(size(latcells));
    for i = 1:length(latcells)
        [in, on] = inpolygon(lonq, latq, loncells{i}, latcells{i});
        if any(in)
            ploc(i) = false;
        end
        assert(nnz(on) == 0, 'Error: Query points are on the edge');
    end
    assert(nnz(~ploc) == length(latq), 'Error: Check query points');
    latcells = latcells(ploc);
    loncells = loncells(ploc);
    [latb, lonb] = polyjoin(latcells, loncells);
end

% simplify
[latb, lonb, err, tol] = reducem(latb, lonb, tol);
fprintf('Simplify coastal polygon (Npts=%d) with tol=%.2f deg, err=%.1f%%.\n', ... %[output:group:3232e88e] %[output:193b5451]
    numel(latb), tol, err*100); %[output:group:3232e88e] %[output:193b5451]
shapeout = geopolyshape(latb, lonb);
shapeout.GeographicCRS = geocrs(4326);

% clip by latitudes
if ~isempty(inp.clip_lat)
    latbnd = [inp.clip_lat(:)' latlim];
    latbnd = sort(latbnd, "descend");
    n = numel(latbnd)-1;
    bboxc = cell(n,1);
    for i = 1:n
        bboxc{i} = [lonlim;latbnd(i+1) latbnd(i)];
    end
    shapeoutC = cell(n,1);
    for i = 1:n
        bnd = bbox2poly(bboxc{i});
        reg = geopolyshape(bnd(:,2), bnd(:,1));
        reg.GeographicCRS = geocrs(4326);
        shapeoutC{i} = intersect(shapeout, reg);
    end
else
    shapeoutC = {shapeout};
end

% save and write shapefile
if inp.save_flag
    fnm = "./coast_polygon_" + reso;
    delete(fnm + '*');
    save(fnm,'shape','latb','lonb','shapeout','shapeoutC');
    for i = 1:numel(shapeoutC)
        GT = table(shapeoutC{i}, VariableNames={'Shape'});
        sfnm = fnm + '_' + num2str(i);
        shapewrite(GT, sfnm);
        zip(sfnm, sfnm + '*');
    end
end