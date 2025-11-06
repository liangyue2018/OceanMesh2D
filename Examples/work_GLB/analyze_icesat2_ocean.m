function TT = analyze_icesat2_ocean(datestr, varargin)
% Read ICESat-2 ATL03 geolocated photon data, calculate ocean surface
% wave parameters using strong beam profiles and write output.
% Usage:
%   TT = analyze_icesat2_ocean('20181014');
%   TT = analyze_icesat2_ocean('201810**');
%
% Inputs:
%   datestr - 'yyyyMMdd' or 'yyyyMM**'
%
% Optional inputs:
%   'oceanSegLen' - ocean segment length to analyze (default: 10.24 km)
%   'depth_threshold' - depth threshold used to filter geosegs (default: -10 m)
%   'distance_threshold' - distance threshold from the coastline (km) used to filter geosegs (default: 3 km)
%   'version' - ATL03 version, '006' or '007' (default: '006')
%   'weight_threshold' - Threshold of photon weight parameter (default: 0.5)
%   'bufWin' - half the height of the ocean surface buffer window (default: 20 m)
%   'interval' - along-track spacing for resampling (default: 5 m)
%   'method' - interpolation method for resampling (default: 'linear')
%   'lackwidth' - lackwidth for checking data gaps (default: 60 m)
%   'missing_threshold' - maximum fraction of missing data allowed (default: 0.2)
%   'psd_flag' - true or false (default: false)
%   'save_flag' - 0 or 1 (default: 1)
%
% Outputs:
%   TT - timetable containing wave parameters (SWH) retrieved from ICESat-2 ATL03 geolocated photon data
%
% check and parse inputs
assert(ischar(datestr) && length(datestr) == 8, 'ERROR: Invalid date string.');

p = inputParser;
addParameter(p, 'oceanSegLen', 10.24, @(x) x > 0);
addParameter(p, 'depth_threshold', -10, @(x) x <= 0);
addParameter(p, 'distance_threshold', 3);
addParameter(p, 'version', '006', @(x) ismember(x, {'006', '007'}));
addParameter(p, 'weight_threshold', 0.2, @(x) x > 0 && x < 1);
addParameter(p, 'bufWin', 20, @(x) x > 0);
addParameter(p, 'interval', 5, @(x) x > 0);
addParameter(p, 'method', 'linear', @(x) ismember(x, {'linear', 'pchip'}));
addParameter(p, 'lackwidth', 60, @(x) x > 0);
addParameter(p, 'missing_threshold', 0.2, @(x) x > 0 && x < 1);
addParameter(p, 'psd_flag', false, @(x) islogical(x));
addParameter(p, 'save_flag', 1, @(x) ismember(x, [0 1]));

parse(p, varargin{:});
inp = p.Results;

% define filepath
if all(isstrprop(datestr, 'digit')) % yyyyMMdd
    filepath = fullfile('data', "ATL03_v" + inp.version, "G" + datestr(1:4), datestr(5:6), datestr(7:8));
    output = sprintf('icesat2_atl03_%s_v%s.csv', datestr, inp.version);
elseif strcmp(datestr(7:8), '**') % yyyyMM**
    filepath = fullfile('data', "ATL03_v" + inp.version, "G" + datestr(1:4), datestr(5:6), '*');
    output = sprintf('icesat2_atl03_%s_v%s.csv', datestr(1:6), inp.version);
else
    error('ERROR: Invalid datestr format.');
end

% get filenames
fstruct = dir(fullfile(filepath, '*.h5'));
h5file = fullfile({fstruct.folder}', {fstruct.name}');
assert(~isempty(h5file), 'ERROR: No h5 files in %s', filepath);

% read h5files
TCell = cell(3*length(h5file),1);
count = 0;
for i = 1:length(h5file)
    fprintf('Granule %d of %d: %s ', i, length(h5file), h5file{i}(end-60:end));

    % Get orbit information
    orb = read_granule_info(h5file{i});
    if isempty(orb.beam_list)
        fprintf('* Skip due to [1] sparse data.\n');
        continue
    end
    if orb.qa_flag == 1
        fprintf('* Skip due to [2] quality assessment failure.\n');
        continue
    end

    % Read strong beams
    numRows = 0;
    for j = 1:numel(orb.beam_list)
        beam = orb.beam_list{j};

        % Read geolocation and geophys_corr data
        TG = read_gtx_geodata(h5file{i}, beam);
        if isempty(TG)
            continue
        end

        % Read heights data
        TH = read_gtx_heights(h5file{i}, beam);
        if isempty(TH)
            continue
        end

        % calculate
        T = calc_wave_params(TH);
        if isempty(T)
            continue
        end
        T.beam = repmat(string(beam), height(T), 1);
        T.rcs = repmat(string(orb.rcs), height(T), 1);
        numRows = numRows + height(T);
        count = count + 1;
        TCell{count} = T;
    end
    if numRows == 0
        fprintf('* Skip due to [3] no enough ocean segments.\n');
        continue
    end
    fprintf('- Processed %d records.\n', numRows);
end
assert(count > 0, 'ERROR: No valid records found.');
TCell = TCell(1:count);
T = vertcat(TCell{:});
TT = table2timetable(T);

% save results
if inp.save_flag
    fprintf('Saving results to %s\n', output);
    writetimetable(TT, output);
end
if nargout == 0
    clear TT;
end

function orb = read_granule_info(h5file)
% Read granule information from ATL03 HDF5 file
%
% get rgt cycle segment [ttttccss]
[~, filename, ~] = fileparts(h5file);
tok = regexp(filename, '_(\d{8})_', 'tokens');
orb.rcs = tok{1}{1};

% get strong beams and check validity
beam_list_l = {'gt1l' 'gt2l' 'gt3l'};
beam_list_r = {'gt1r' 'gt2r' 'gt3r'};
orb.sc_orient = h5read(h5file, '/orbit_info/sc_orient'); % [0 1 2] [backward forward transition]
switch orb.sc_orient
    case 0
        orb.beam_list = beam_list_l; % strong lead weak
    case 1
        orb.beam_list = beam_list_r; % weak lead strong
    case 2
        orb.beam_list = {};
        return
end
for gtx = orb.beam_list
    try
        % check if height dataset exists
        info = h5info(h5file, "/" + gtx{1} + '/heights/h_ph');
        assert(~isempty(info));

        % check if ocean segment length is desired
        segment_length = h5read(h5file, "/" + gtx{1} + '/geolocation/segment_length'); % [meters] DOUBLE
        segment_ph_cnt = h5read(h5file, "/" + gtx{1} + '/geolocation/segment_ph_cnt'); % _FillValue=0
        surf_type = h5read(h5file, "/" + gtx + '/geolocation/surf_type', [2 1], [1 Inf])'; % Nx1 [0 1] [not-type is-type]
        sloc = segment_ph_cnt > 0 & surf_type == 1;
        assert(check_length(segment_length(sloc), inp.oceanSegLen));
    catch
        orb.beam_list = setdiff(orb.beam_list, gtx);
    end
end
if isempty(orb.beam_list)
    return
end

% get ascending/descending flag
start_region = h5read(h5file, '/ancillary_data/start_region'); % INTEGER_4
end_region = h5read(h5file, '/ancillary_data/end_region'); % INTEGER_4
assert(start_region == end_region, 'Start region (%d) and end region (%d) are different.', start_region, end_region);
orb.region = sprintf('%02d', start_region);
if ismember(orb.region, {'01','02','03','12','13','14'})
    orb.gt_orient = 1; % ascending
elseif ismember(orb.region, {'05','06','07','08','09','10'})
    orb.gt_orient = 0; % descending
else
    error('Invalid segment number: %s', orb.region);
end

% get bounding polygon
orb.bnd_poly_lat = h5read(h5file, '/orbit_info/bounding_polygon_lat1'); % [-90 90] FLOAT
orb.bnd_poly_lon = h5read(h5file, '/orbit_info/bounding_polygon_lon1'); % [-180 180] FLOAT

% get granule time info
t1 = h5read(h5file, '/ancillary_data/granule_start_utc'); % STRING
t2 = h5read(h5file, '/ancillary_data/granule_end_utc'); % STRING
orb.gran_start_utc = datetime(t1, InputFormat='yyyy-MM-dd''T''HH:mm:ss.SSSSSSZ', TimeZone='UTC');
orb.gran_end_utc = datetime(t2, InputFormat='yyyy-MM-dd''T''HH:mm:ss.SSSSSSZ', TimeZone='UTC');

% get quality assessment flag
orb.qa_flag = h5read(h5file, '/quality_assessment/qa_granule_pass_fail'); % [0 1] [pass fail]
end

function TG = read_gtx_geodata(h5file, gtx)
% Read geolocation and geophys_corr group for a given beam (gtx)
%
% Create timetable with UTC time
delta_time = h5read(h5file, "/" + gtx + '/geolocation/delta_time'); % seconds since 2018-01-01 DOUBLE
utc_time = gps2utc(delta_time);
TG = timetable(utc_time);
TG.Properties.Description = "Geolocation parameters and geophysical corrections posted at " + newline + ...
                           "~20m along-track segment interval for " + upper(gtx);

% Read datasets in /gtx/geolocation group
TG.ph_idx_beg = h5read(h5file, "/" + gtx + '/geolocation/ph_index_beg'); % _FillValue=0
TG.podppd_flag = h5read(h5file, "/" + gtx + '/geolocation/podppd_flag'); % [0 1 2 3 4 5 6 7] [NOMINAL DEGRAGE ...]
TG.ref_ph_idx = h5read(h5file, "/" + gtx + '/geolocation/reference_photon_index'); % _FillValue=0
TG.ref_ph_lat = h5read(h5file, "/" + gtx + '/geolocation/reference_photon_lat'); % [-90 90] DOUBLE
TG.ref_ph_lon = h5read(h5file, "/" + gtx + '/geolocation/reference_photon_lon'); % [-180 180] DOUBLE
TG.segment_dist_x = h5read(h5file, "/" + gtx + '/geolocation/segment_dist_x'); % [meters] DOUBLE
TG.segment_id = h5read(h5file, "/" + gtx + '/geolocation/segment_id');
TG.segment_length = h5read(h5file, "/" + gtx + '/geolocation/segment_length'); % [meters] DOUBLE
TG.segment_ph_cnt = h5read(h5file, "/" + gtx + '/geolocation/segment_ph_cnt'); % _FillValue=0
TG.surf_type = h5read(h5file, "/" + gtx + '/geolocation/surf_type', [2 1], [1 Inf])'; % Nx1 [0 1] [not-type is-type] (five type: land ocean sea-ice land-ice inland-water)

% Read datasets in /gtx/geophys_corr group
TG.dac = h5read(h5file, "/" + gtx + '/geophys_corr/dac'); % [meters] FLOAT _FillValue=3.4028235E38
TG.tide_equilibrium = h5read(h5file, "/" + gtx + '/geophys_corr/tide_equilibrium'); % [meters] FLOAT _FillValue=3.4028235E38
TG.tide_ocean = h5read(h5file, "/" + gtx + '/geophys_corr/tide_ocean'); % [meters] FLOAT _FillValue=3.4028235E38
TG.geoid = h5read(h5file, "/" + gtx + '/geophys_corr/geoid'); % [meters] FLOAT _FillValue=3.4028235E38
TG.geoid_free2mean = h5read(h5file, "/" + gtx + '/geophys_corr/geoid_free2mean'); % [meters] FLOAT _FillValue=3.4028235E38

% Add dataset descriptions
TG.Properties.VariableDescriptions = {'Index of the first photon in a given segment', ...
                                      'Flag indicates the quality of ATL03 geo-segments', ...
                                      'Index of the reference photon within a segment', ...
                                      'Latitude of the reference photon', ...
                                      'Longitude of the reference photon', ...
                                      'Along-track distance (meters) from the equator crossing to the segment start', ...
                                      'Unique number for the segment', ...
                                      'Along-track length (meters) of the segment', ...
                                      'Number of photons in a given along-track segment', ...
                                      'Surf type (5xN, land, ocean, sea-ice, land-ice and inland-water) for each segment', ...
                                      'Dynamic atmospheric correction (±0.2m) from MOG2D (6-hourly, 0.25°)', ...
                                      'Long period equilibrium tide (meters)', ...
                                      'Short period ocean tides (diurnal and semi-diurnal, ±4m)', ...
                                      'Geoid height (meters) above WGS84 reference ellipsoid in the tide-free system', ...
                                      'Geoid free-to-mean conversion (meters)'};

% Coarse selection for ocean segments
sloc = TG.segment_ph_cnt > 0 & ...
       (TG.podppd_flag == 0 | TG.podppd_flag == 4) & ...
       (TG.surf_type == 1);

% Exclude Fillvalues (no tides) for ocean segments
chkVars = {'dac', 'tide_equilibrium', 'tide_ocean', 'geoid', 'geoid_free2mean'};
fillVal = realmax('single'); % 3.4028235e+38
tf = ismissing(TG{:,chkVars}, fillVal);
sloc = sloc & ~any(tf, 2);
TG = TG(sloc, :);
if ~check_length(TG.segment_length, inp.oceanSegLen)
    TG = [];
    return
end
assert(nnz(TG.ph_idx_beg == 0) == 0 && nnz(TG.ref_ph_idx == 0) == 0, 'Error: Found fill values (0)');
tf = TG{:,chkVars} > 1e38;
if any(tf, 'all')
    error('Error: Found fill values (%e) in variable(s): %s', fillVal, strjoin(chkVars(any(tf, 1)), ', '));
end

% Fine selection based on water depth and distance from coastline (if specified)
TG.depth_ocn_seg = grdtrack(TG.ref_ph_lon, TG.ref_ph_lat, 'GEBCO_2025'); % 15 arc-second
sloc = TG.depth_ocn_seg <= inp.depth_threshold;
if ~isempty(inp.distance_threshold)
    TG.dist_ocn_seg = grdtrack(TG.ref_ph_lon, TG.ref_ph_lat, 'earth_dist_01m'); % 1 arc-minute
    TG.dist_ocn_seg = TG.dist_ocn_seg * (-1); % >0: ocean to coastline
    sloc = sloc & TG.dist_ocn_seg >= inp.distance_threshold;
end
TG = TG(sloc, :);
if ~check_length(TG.segment_length, inp.oceanSegLen)
    TG = [];
    return
end

% Get MSS for each ocean segment
TG.mss = grdtrack(TG.ref_ph_lon, TG.ref_ph_lat, 'DTU21'); % 1 arc-minute
end

function TH = read_gtx_heights(h5file, gtx)
% Read heights group for a given beam (gtx)
%
% Create timetable with UTC time
delta_time = h5read(h5file, "/" + gtx + '/heights/delta_time'); % seconds since 2018-01-01 DOUBLE
utc_time = gps2utc(delta_time);
TH = timetable(utc_time);
TH.Properties.Description = "Contains array of the parameters stored at" + newline + ...
                            "the photon detection rate for " + upper(gtx);

% Filter extra geosegs in some cases
begIdx = double(TG.ph_idx_beg); % unique to a specific granule
endIdx = begIdx + double(TG.segment_ph_cnt) - 1;
%refIdx = begIdx + double(TG.ref_ph_idx) - 1;
sloc = begIdx > 0 & endIdx <= height(TH);
TG = TG(sloc, :);
if ~check_length(TG.segment_length, inp.oceanSegLen)
    return
end
begIdx = max(1, begIdx(sloc));
endIdx = min(endIdx(sloc), height(TH));

% Get photon mask from the indicies of selected geosegs
loc = false(height(TH), 1);
for k = 1:height(TG)
    b = begIdx(k);
    e = endIdx(k);
    loc(b:e) = true;
end
TH = TH(loc, :);

% Read datasets in /gtx/heights group
TH.dist_ph_along = h5read_block(h5file, "/" + gtx + '/heights/dist_ph_along', loc); % [meters] FLOAT
TH.h_ph = h5read_block(h5file, "/" + gtx + '/heights/h_ph', loc); % [meters] FLOAT
TH.lat_ph = h5read_block(h5file, "/" + gtx + '/heights/lat_ph', loc); % [-90 90] FLOAT
TH.lon_ph = h5read_block(h5file, "/" + gtx + '/heights/lon_ph', loc); % [-180 180] FLOAT
TH.quality_ph = h5read_block(h5file, "/" + gtx + '/heights/quality_ph', loc); % '006': [0 1 2 3] [nominal ...], '007': [0 3 4 5 ...] [nominal tep noise-burst noise-streak]
TH.signal_conf_ph = h5read_block(h5file, "/" + gtx + '/heights/signal_conf_ph', {2 loc}); % Nx1 [-2 -1 0 1 2 3 4] [TEP not-type noise background low med high]
TH.weight_ph = h5read_block(h5file, "/" + gtx + '/heights/weight_ph', loc); % '006': [0 255], '007': [0 65535]
if strcmp(inp.version, '007')
    TH.weight_ph = double(TH.weight_ph) / 65535;
else
    TH.weight_ph = double(TH.weight_ph) / 255;
end

% Add dataset descriptions
TH.Properties.VariableDescriptions = {'Along-track distance (meters) from the segment start to each photon', ...
                                      'Height of each received photon above WGS84 reference ellipsoid (meters)', ...
                                      'Latitude of each received photon', ...
                                      'Longitude of each received photon', ...
                                      '* Photon quality as to saturation and after-pulse (updated in 007)', ...
                                      'Confidence level associated with each photon event selected as signal', ...
                                      '* Calculated weight of each received photon (updated in 007)'};

% Update signal_class_ph for version '007'
if strcmp(inp.version, '007')
    TH.signal_class_ph = h5read_block(h5file, "/" + gtx + '/heights/signal_class_ph', loc); % [-1 0 1 2 3 4 5] [ignored likely-noise likely-signal ...]
    TH.Properties.VariableDescriptions{'signal_class_ph'} = '* Photon-rate flag based on evaluation of weight_ph';
end

% Asign geoseg-rate parameters to each photon
TH.dist_ph_along = double(TH.dist_ph_along) + repelem(TG.segment_length .* double(TG.segment_id - TG.segment_id(1)), TG.segment_ph_cnt);
TH = renamevars(TH, 'dist_ph_along', 'cum_dist_along');
TH.h_ph = TH.h_ph - repelem(TG.dac + TG.tide_equilibrium + TG.tide_ocean, TG.segment_ph_cnt);
TH.geoid_mean = repelem(TG.geoid + TG.geoid_free2mean, TG.segment_ph_cnt);
TH.depth_ocn_seg = repelem(TG.depth_ocn_seg, TG.segment_ph_cnt);
TH.dist_ocn_seg = repelem(TG.dist_ocn_seg, TG.segment_ph_cnt);
TH.mss = repelem(TG.mss, TG.segment_ph_cnt);

% Select photon heights for med and high confidence, no saturation or at least signal-below (for version '007' only)
loc = TH.signal_conf_ph >= 3 & TH.quality_ph == 0 & TH.weight_ph > inp.weight_threshold;
if strcmp(inp.version, '007')
    loc = loc & (TH.signal_class_ph >= 2);
end
TH = TH(loc, :);
if isempty(TH)
    return
end
chkVars = {'geoid_mean', 'depth_ocn_seg', 'dist_ocn_seg', 'mss'};
assert(nnz(isnan(TH{:, chkVars})) == 0, 'Some variables contain NaN values');

% Filter photon heights with a priori estimate of the surface elevation (MSS)
loc = abs(TH.h_ph - TH.mss) <= inp.bufWin;
TH = TH(loc, :);
if isempty(TH)
    return
end

% Get unique utc_time and cum_dist_along
TH.h_ph = TH.h_ph - TH.geoid_mean; % convert to orthometric height
TH = renamevars(TH, 'h_ph', 'h_ortho');
TH = TH(:, {'cum_dist_along', 'h_ortho', 'lat_ph', 'lon_ph', 'depth_ocn_seg', 'dist_ocn_seg'});
[~, ia] = unique(TH.utc_time);
TH = TH(ia, :);
[~, ia] = unique(TH.cum_dist_along);
TH = TH(ia, :);
if height(TH) < inp.oceanSegLen * 1000 / 0.7
    TH = [];
    return
end
end

function T = calc_wave_params(TH)
% Calculate ocean surface wave parameters from photon height profile
%
% Initialize a table for output
names = {'Time' 'Lat' 'Lon' 'missFrac' 'depth' 'distance' 'swh'};
types = {'datetime' 'double' 'double' 'double' 'double' 'double' 'double'};
if inp.psd_flag
    names = [names 'pwl'];
    types = [types 'double'];
end
T = table('Size', [0 numel(names)], 'VariableTypes', types, 'VariableNames', names);
T.Time.TimeZone = 'UTC';

% Define resampling grid
xq = TH.cum_dist_along(1):inp.interval:TH.cum_dist_along(end);
oceanSegPts = floor(inp.oceanSegLen * 1000 / inp.interval) + 1;
stepPts = floor((oceanSegPts - 1) / 2);
starts = 1:stepPts:(numel(xq) - oceanSegPts + 1);
ends = starts + oceanSegPts - 1;
nSeg = numel(starts);
if nSeg < 1
    return
end
qloc = check_gaps(TH.cum_dist_along, xq, inp.lackwidth);

% calculate wave parameters for each overlapping ocean segment
for k = 1:nSeg
    s = starts(k);
    e = ends(k);
    xq_seg = xq(s:e); % 5m

    % set missing threshold
    missFrac = nnz(qloc(s:e)) / numel(xq_seg);
    if missFrac > inp.missing_threshold
        continue
    end

    % Remove outliers and use S-G filter
    loc = TH.cum_dist_along >= xq_seg(1) & TH.cum_dist_along < xq_seg(end);
    x_seg = TH.cum_dist_along(loc); % 0.7m
    v_seg = double(TH.h_ortho(loc));
    tf = isoutlier(v_seg, "movmedian", 51, "SamplePoints", x_seg);
    x_seg = x_seg(~tf);
    v_seg = v_seg(~tf);
    v_seg = sgolayfilt(v_seg, 5, 51);

    % Resample the along-track photon height profile
    vq_seg = interp1(x_seg, v_seg, xq_seg, inp.method, 0);
    assert(nnz(isnan(vq_seg)) == 0, 'Error: Interpolation resulted in NaN values');

    % get the mean time and location
    mtime = mean(posixtime(TH.utc_time(loc)));
    mtime = datetime(mtime, ConvertFrom='posixtime', TimeZone='UTC');
    mlon = mean(TH.lon_ph(loc));
    mlat = mean(TH.lat_ph(loc));
    mdepth = min(TH.depth_ocn_seg(loc));
    mdist = min(TH.dist_ocn_seg(loc));

    % S-G filter and detrend
    vq_seg = detrend(vq_seg - mean(vq_seg));
    swh = 4 * std(vq_seg); % significant wave height

    % calculate wavenumber spectrum if desired
    if inp.psd_flag
        N = numel(vq_seg);
        nsc = floor(N / 8);
        nov = floor(nsc / 2);
        nff = max(256, 2^nextpow2(nsc));
        fs = 1 / inp.interval; % [1/m]
        [pxx, f] = pwelch(vq_seg, hamming(nsc), nov, nff, fs);
        [fp, IKM] = calc_fp(f, pxx); % 0:fs/N:fs/2 [nff/2+1 1]
        if IKM < 3
            pwl = nan;
        else
            pwl = 1 / fp; % peak wavelength
        end
    end

    % Store results
    if inp.psd_flag
        T(end+1, :) = {mtime, mlat, mlon, missFrac, mdepth, mdist, swh, pwl}; %#ok<AGROW>
    else
        T(end+1, :) = {mtime, mlat, mlon, missFrac, mdepth, mdist, swh}; %#ok<AGROW>
    end
end

end

end

function utc_time = gps2utc(delta_time)
% Convert ICESat-2 GPS time (seconds since 1980-01-06 00:00:00 UTC) to UTC datetime
% UTC time = GPS time - leap seconds
% Input:
%   delta_time - seconds since 2018-01-01 00:00:00 UTC 
% Output:
%   utc_time - datetime array in UTC
%
atlas_sdp_gps_epoch = 1198800000;
gps_time = delta_time + atlas_sdp_gps_epoch;
leap_seconds = 18; % as of 2017-01-01
t0 = datetime(1980, 1, 6, 0, 0, 0, TimeZone='UTC');
utc_time = t0 + seconds(gps_time - leap_seconds);
end

function ok = check_length(segment_length, threshold)
% Check if the segment length is above the threshold
assert(isvector(segment_length) && isscalar(threshold), 'Error: Invalid inputs.');
totalSegLen = sum(double(segment_length)) / 1000; % [km]
ok = totalSegLen >= double(threshold);
end

function out = h5read_block(h5file, dataset, loc)
% Read HDF5 dataset using logical indexing array
%
% check dimensions
info = h5info(h5file, dataset);
dims = info.Dataspace.Size;
if isscalar(dims)
    assert(islogical(loc) && numel(loc) == dims, 'Error: Dimensions do not match.');
    idx1 = find(loc);
elseif numel(dims) == 2
    assert(iscell(loc) && numel(loc) == 2, 'Error: For 2D dataset, loc must be a cell array with 2 elements.');
    assert(loc{1} <= dims(1) && islogical(loc{2}) && numel(loc{2}) == dims(2), 'Error: Dimensions do not match.');
    idx1 = find(loc{2});
    idx2 = loc{1};
else
    error('Error: Only 1D or 2D datasets are supported.');
end
assert(~isempty(idx1), 'Error: No valid indices found.');

% Get continuous block
d = diff(idx1) ~= 1;
starts = idx1([true; d]);
ends = idx1([d; true]);

blocks = cell(numel(starts), 1);
for k = 1:numel(starts)
    s = starts(k);
    e = ends(k);
    cnt = e - s + 1;
    if isscalar(dims)
        % 1D dataset
        start = s;
        count = cnt;
    else
        % 2D dataset
        start = [idx2 s];
        count = [1 cnt];
    end
    blk = h5read(h5file, dataset, start, count);
    blocks{k} = blk(:);
end

out = vertcat(blocks{:});
assert(iscolumn(out) && length(out) == length(idx1), 'Error: Output vector does not match expected length.');
end

function qloc = check_gaps(x, xq, lackwidth)
% Check gaps for input data when resampling
%
assert(isvector(x) && isvector(xq) && isscalar(lackwidth), 'Error: Invalid inputs.');
dfx = diff(x);
gap_idx = find(dfx > lackwidth);
qloc = false(size(xq));
for g = gap_idx(:)'
    gstart = x(g);
    gend = x(g+1);
    qloc = qloc | (xq > gstart & xq < gend);
end
end