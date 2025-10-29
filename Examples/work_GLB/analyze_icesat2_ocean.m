function S = analyze_icesat2_ocean(datestr, varargin)
% Read ICESat-2 ATL03 geolocated photon data, calculate ocean surface
% wave parameters using strong beam profiles and write output.
% Usage:
%   S = analyze_icesat2_ocean('20181014', 1);
%   S = analyze_icesat2_ocean('201810**', 1);
%
% Inputs:
%   datestr - 'yyyyMMdd' or 'yyyyMM**'
%
% Optional inputs:
%   'oceanSegLen' - ocean segment length to analyze (default: 10.24 km)
%   'depth' - depth limit to filter geosegs (default: -10 m)
%   'distance' - distance from the coastline (km) used to filter geosegs (default: [])
%   'save_flag' - 0 or 1 (default: 1)
%
% Outputs:
%   S - nx1 struct, n is the number of granules
%   S.bnd_poly - geopolyshape of granule bounding polygon
%   S.shape - geopointshape of ref_ph_lat/ref_ph_lon
%   To be added...
%
% check and parse inputs
assert(ischar(datestr) && length(datestr) == 8, 'ERROR: Invalid date string.');

p = inputParser;
addParameter(p, 'oceanSegLen', 10.24, @(x) x > 0);
addParameter(p, 'depth', -10, @(x) x <= 0);
addParameter(p, 'distance', []);
addParameter(p, 'save_flag', 1, @(x) ismember(x, [0 1]));

parse(p, varargin{:});
inp = p.Results;

% define filepath
if all(isstrprop(datestr, 'digit')) % yyyyMMdd
    filepath = "./data/G" + datestr(1:4) + '/' + datestr(5:6) + '/' + datestr(7:8) + '/';
    output = "S_" + datestr + '.mat';
elseif strcmp(datestr(7:8), '**') % yyyyMM**
    filepath = "./data/G" + datestr(1:4) + '/' + datestr(5:6) + '/*/';
    output = "S_" + datestr(1:6) + '.mat';
else
    error('ERROR: Invalid datestr format.');
end

% get filename
fstruct = dir(filepath + '*.h5');
h5file = fullfile({fstruct.folder}', {fstruct.name}');
assert(~isempty(h5file), 'ERROR: No h5 files in %s', filepath);

% read h5file
S = struct();
tic;
fprintf('+-----------------------------------------------------------------------------+\n');
fprintf('Start processing %d granules for %s...\n', length(h5file), datestr);
for n = 1:length(h5file)
    fprintf('Granule %d of %d: %s\n', n, length(h5file), h5file{n}(end-60:end));

    % Get orbit information
    orb = read_granule_info(h5file{n});
    if isempty(orb.beam_list)
        fprintf('>>>>Skip granule %d due to empty beam list.\n', n);
        continue
    end
    if orb.qa_flag == 1
        fprintf('>>>>Skip granule %d due to quality assessment failure.\n', n);
        continue
    end
    S(n).orbit_info = orb;

    % Read strong beams
    ref_ph_lat = cell(numel(orb.beam_list),1);
    ref_ph_lon = cell(numel(orb.beam_list),1);
    for i = 1:numel(orb.beam_list)
        beam = orb.beam_list{i};
        GT = read_gtx_geodata(h5file{n}, beam);
        if isempty(GT)
            continue
        end

        ref_ph_idx = h5read(h5file{n}, "/" + beam + '/geolocation/reference_photon_index'); % _FillValue=0



        ref_ph_lat{i}(1,:) = h5read(h5file{n}, "/" + beam + '/geolocation/reference_photon_lat'); % [-90 90] DOUBLE
        ref_ph_lon{i}(1,:) = h5read(h5file{n}, "/" + beam + '/geolocation/reference_photon_lon'); % [-180 180] DOUBLE
    end
    S(n).shape = geopointshape(ref_ph_lat, ref_ph_lon);
    S(n).shape.GeographicCRS = geocrs(4326);
end
elapsedTime = toc;
fprintf('Finished processing %d granules in %s.\n', length(h5file), formatElapsedTime(elapsedTime));
if inp.save_flag
    fprintf('Saving output to %s\n', output);
    save(output, 'S');
end
fprintf('+-----------------------------------------------------------------------------+\n');

if nargout == 0
    clear S;
end

function orb = read_granule_info(h5file)
% Read granule information from ATL03 HDF5 file
orb.atl03_file = h5file;

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
        surf_type = h5read(h5file, "/" + gtx + '/geolocation/surf_type', [2 1], [1 Inf]); % 5xN [0 1] [not-type is-type]
        sloc = segment_ph_cnt > 0 & surf_type == 1;
        totalSegLen = sum(segment_length(sloc)) / 1000; % [km]
        assert(totalSegLen >= inp.oceanSegLen);

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
bnd_poly_lat = h5read(h5file, '/orbit_info/bounding_polygon_lat1'); % [-90 90] FLOAT
bnd_poly_lon = h5read(h5file, '/orbit_info/bounding_polygon_lon1'); % [-180 180] FLOAT
orb.bnd_poly = geopolyshape(bnd_poly_lat, bnd_poly_lon);
orb.bnd_poly.GeographicCRS = geocrs(4326);

% get granule time info
t1 = h5read(h5file, '/ancillary_data/granule_start_utc'); % STRING
t2 = h5read(h5file, '/ancillary_data/granule_end_utc'); % STRING
orb.gran_start_utc = datetime(t1, InputFormat='yyyy-MM-dd''T''HH:mm:ss.SSSSSSZ', TimeZone='UTC');
orb.gran_end_utc = datetime(t2, InputFormat='yyyy-MM-dd''T''HH:mm:ss.SSSSSSZ', TimeZone='UTC');

% get quality assessment flag
orb.qa_flag = h5read(h5file, '/quality_assessment/qa_granule_pass_fail'); % [0 1] [pass fail]
end

function GT = read_gtx_geodata(h5file, gtx)
% Read geolocation and geophys_corr group for a given beam (gtx)
delta_time = h5read(h5file, "/" + gtx + '/geolocation/delta_time'); % seconds since 2018-01-01 DOUBLE
utc_time = gps2utc(delta_time);
GT = timetable(utc_time);
GT.Properties.Description = "Geolocation parameters and geophysical corrections posted at " + newline + ...
                           "~20m along-track segment interval for " + upper(gtx);

% Read datasets in /gtx/geolocation group
GT.ph_index_beg = h5read(h5file, "/" + gtx + '/geolocation/ph_index_beg'); % _FillValue=0
GT.podppd_flag = h5read(h5file, "/" + gtx + '/geolocation/podppd_flag'); % [0 1 2 3 4 5 6 7] [NOMINAL DEGRAGE ...]
GT.ref_ph_idx = h5read(h5file, "/" + gtx + '/geolocation/ref_ph_idx'); % _FillValue=0
GT.ref_ph_lat = h5read(h5file, "/" + gtx + '/geolocation/ref_ph_lat'); % [-90 90] DOUBLE
GT.ref_ph_lon = h5read(h5file, "/" + gtx + '/geolocation/ref_ph_lon'); % [-180 180] DOUBLE
GT.segment_dist_x = h5read(h5file, "/" + gtx + '/geolocation/segment_dist_x'); % [meters] DOUBLE
GT.segment_id = h5read(h5file, "/" + gtx + '/geolocation/segment_id');
GT.segment_length = h5read(h5file, "/" + gtx + '/geolocation/segment_length'); % [meters] DOUBLE
GT.segment_ph_cnt = h5read(h5file, "/" + gtx + '/geolocation/segment_ph_cnt'); % _FillValue=0
GT.surf_type = h5read(h5file, "/" + gtx + '/geolocation/surf_type', [2 1], [1 Inf]); % 5xN [0 1] [not-type is-type]

% Read datasets in /gtx/geophys_corr group
GT.dac = h5read(h5file, "/" + gtx + '/geophys_corr/dac'); % [meters] FLOAT _FillValue=3.4028235E38
GT.tide_equilibrium = h5read(h5file, "/" + gtx + '/geophys_corr/tide_equilibrium'); % [meters] FLOAT _FillValue=3.4028235E38
GT.tide_ocean = h5read(h5file, "/" + gtx + '/geophys_corr/tide_ocean'); % [meters] FLOAT _FillValue=3.4028235E38
GT.geoid = h5read(h5file, "/" + gtx + '/geophys_corr/geoid'); % [meters] FLOAT _FillValue=3.4028235E38
GT.geoid_free2mean = h5read(h5file, "/" + gtx + '/geophys_corr/geoid_free2mean'); % [meters] FLOAT _FillValue=3.4028235E38

% Add dataset descriptions
GT.Properties.VariableDescriptions = {'datetime UTC', ...
                                      'Index of the first photon in a given segment', ...
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
sloc = GT.segment_ph_cnt > 0 & ...
       (GT.podppd_flag == 0 | GT.podppd_flag == 4) & ...
       (GT.surf_type == 1);
GT = GT(sloc, :);
totalSegLen = sum(GT.segment_length) / 1000; % [km]
if totalSegLen < inp.oceanSegLen
    GT = [];
    return
end
assert(nnz(GT.ph_index_beg == 0) == 0 && nnz(GT.ref_ph_idx == 0) == 0, 'Error: Found fill values (0)');

% Fine selection based on water depth and distance from the coastline
GT.depth_ocn_seg = get_lonlat_elevation(GT.ref_ph_lon, GT.ref_ph_lat, 'GEBCO_2025'); % 15 arc-second
sloc = GT.depth_ocn_seg <= inp.depth;
if ~isempty(inp.distance)
    GT.dist_coast_seg = get_lonlat_dist2coast(GT.ref_ph_lon, GT.ref_ph_lat);
    sloc = sloc & GT.dist_coast_seg >= inp.distance;
end
GT = GT(sloc, :);
totalSegLen = sum(GT.segment_length) / 1000; % [km]
if totalSegLen < inp.oceanSegLen
    GT = [];
    return
end
fillval = 3.4028235e+38;
chkVars = {'dac', 'tide_equilibrium', 'tide_ocean', 'geoid', 'geoid_free2mean'};
TF = ismissing(GT{:,chkVars}, fillval);
if any(TF, 'all')
    error('Error: Found fill values (%e) in variable(s): %s', fillval, strjoin(chkVars(any(TF, 1)), ', '));
end

% Get MSS for each ocean segment
GT.mss = get_lonlat_elevation(GT.ref_ph_lon, GT.ref_ph_lat, 'DTU21'); % 1 arc-minute
end

%

end

function utc_time = gps2utc(delta_time)
% Convert ICESat-2 GPS time (seconds since 1980-01-06 00:00:00 UTC) to UTC datetime
% UTC time = GPS time - leap seconds
% Input:
%   delta_time - seconds since 2018-01-01 00:00:00 UTC 
% Output:
%   utc_time - datetime array in UTC

atlas_sdp_gps_epoch = 1198800000;
gps_time = delta_time + atlas_sdp_gps_epoch;
leap_seconds = 18; % as of 2017-01-01
t0 = datetime(1980, 1, 6, 0, 0, 0, TimeZone='UTC');
utc_time = t0 + seconds(gps_time - leap_seconds);
end

function elapsedTimeStr = formatElapsedTime(elapsedTime)
% Format elapsed time (seconds) for display
if elapsedTime < 60
    elapsedTimeStr = sprintf('%.2f seconds', elapsedTime);
elseif elapsedTime < 3600
    elapsedTimeStr = sprintf('%.2f minutes', elapsedTime / 60);
else
    elapsedTimeStr = sprintf('%.2f hours', elapsedTime / 3600);
end
end
