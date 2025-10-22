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
for n = 1:numel(fnm)
    fprintf('Processing file %d of %d: %s\n', n, numel(fnm), fnm{n});
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
end

% save result
sfnm = ['is2_granules_' datestr '.mat'];
save(sfnm, 'S');
