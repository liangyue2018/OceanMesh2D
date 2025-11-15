from harmony import Client, Collection, Request, BBox
from concurrent.futures import as_completed
import h5py
import datetime as dt
import time, calendar
import os

def versionToCMR(version):
    '''
    Convert version string to CMR ID
    '''
    version_dict = {'006': 'C2596864127-NSIDC_CPRD',
                    '007': 'C3326974349-NSIDC_CPRD'}
    if version not in version_dict:
        raise ValueError(f'Unknown version: {version}')
    return version_dict[version]

def check_results(final_json):
    '''
    Check results of Harmony job
    '''
    # check job status
    job_status = final_json['status']
    job_message = final_json['message']
    if job_status.lower() != 'successful':
        raise RuntimeError(f'Harmony job failed with status: {job_status}, message: {job_message}')
    
    # print data size reduction
    numGrans = final_json['numInputGranules']
    original_size = final_json['originalDataSize']
    output_size = final_json['outputDataSize']
    pct_reduction = final_json['dataSizePercentChange']
    print(f"Job completed successfully, achieved a {pct_reduction} ({original_size} -> {output_size})")
    return numGrans

def isGranuleDesired(filename, shapefile):
    '''
    Check if granule is desired based on segment number determined by shapefile
    '''
    # check inputs
    if not filename.startswith('ATL03_'):
        raise ValueError(f'Invalid filename: {filename}')
    
    if shapefile.endswith('_1.zip'):
        segnums = ['02', '03', '05', '06']
    elif shapefile.endswith('_2.zip'):
        segnums = ['14', '01', '07', '08']
    elif shapefile.endswith('_3.zip'):
        segnums = ['09', '10', '12', '13']
    else:
        raise ValueError(f'Unknown shapefile: {shapefile}')
    
    # extract segment number
    ss = filename.split('_')[2][6:8]
    return ss in segnums

def isGranuleValid(fname):
    '''Check if granule is valid
    '''
    try:
        with h5py.File(fname, 'r') as f:
            _ = dict(f.attrs)
        return True
    except:
        return False

def formatElapsedTime(elapsedTime):
    '''
    Format elapsed time in seconds for display
    '''
    if elapsedTime < 60:
        return f"{elapsedTime:.2f} seconds"
    elif elapsedTime < 3600:
        return f"{elapsedTime/60:.1f} minutes"
    else:
        return f"{elapsedTime/3600:.1f} hours"

def download_icesat2_atl03(date_range: dict[str, dt.datetime], shapefile: str = None, version: str = '007', save_dir: str = None, 
                           bbox: list[float] = None, max_results: int = None, granule_name: str = None):
    '''
    Download ICESat-2 ATL03 data for given date range and shape file
    '''
    # validate inputs
    if 'start' not in date_range or 'stop' not in date_range:
        raise ValueError('date_range must have "start" and "stop" keys!')
    if shapefile:
        if not os.path.exists(shapefile):
            raise FileNotFoundError(f'Shape file {shapefile} not found!')
    if bbox:
        if len(bbox) != 4:
            raise ValueError('bbox must be a list of four floats: [west, east, south, north]')
        bbox = BBox(w=bbox[0], e=bbox[1], s=bbox[2], n=bbox[3])
    
    # set env
    os.environ['NUM_REQUESTS_WORKERS'] = '4'

    # create harmony request
    cmr_id = versionToCMR(version)
    harmony_client = Client()
    request = Request(collection=Collection(id=cmr_id),
                      temporal=date_range,
                      shape=shapefile,
                      skip_preview=True,
                      spatial=bbox,
                      max_results=max_results,
                      granule_name=granule_name)
    if not request.is_valid():
        raise ValueError(f"Invalid request parameters: {request.error_messages()}")
    
    # submit request
    job_id = harmony_client.submit(request)
    print(f'Submitted job {job_id}')

    # wait for processing to complete
    harmony_client.wait_for_processing(job_id, show_progress=True)

    # retrieve job results
    final_json = harmony_client.result_json(job_id)
    numGrans = check_results(final_json)
    urls = list(harmony_client.result_urls(job_id))
    if len(urls) != numGrans:
        raise RuntimeError(f"Expected {numGrans} granules, but got {len(urls)}")
    
    # filter granules with segment number determined by shapefile
    if shapefile:
        filtered_urls = []
        for url in urls:
            filename = os.path.basename(url)
            if not isGranuleDesired(filename, shapefile):
                continue
            filtered_urls.append(url)

        urls = filtered_urls
        if not urls:
            raise ValueError(f"No valid granules found for shapefile: {shapefile}")

    # setup save directory
    os.makedirs(save_dir, exist_ok=True)
    ofnm = os.path.join(save_dir, "download_urls.txt")
    with open(ofnm, "w") as f:
        f.writelines(f"{url}\n" for url in urls)
    os.system(f"sort -u {ofnm} -o {ofnm}")

    # download data files
    def _submit_download(url):
        return harmony_client.download(url, directory=save_dir, overwrite=True)
    
    pending = {_submit_download(url): (url, 1) for url in urls}
    downloaded = []
    failed = []
    max_attempts = 5

    while pending:
        for f in as_completed(list(pending.keys())):
            url, attempt = pending.pop(f)
            try:
                result = f.result()

                # check validity
                if not isGranuleValid(result):
                    os.remove(result)
                    raise RuntimeError(f"Invalid HDF5 file: {result}")
                
                # rename files
                fpath = os.path.dirname(result)
                fname = os.path.basename(result)
                fnew = os.path.join(fpath, f"ATL03_{fname.split('ATL03_')[1]}")
                os.rename(result, fnew)

                downloaded.append(fnew)
            except Exception as e:
                print(f"WARNING: Download failed (attempt {attempt}/{max_attempts}): {e}.")
                if attempt < max_attempts:
                    time.sleep(30)
                    pending[_submit_download(url)] = (url, attempt + 1)
                else:
                    failed.append(url)

    # final checks
    if failed:
        raise RuntimeError(f"ERROR: Download failed for {len(failed)} files, need to retry manually."
                           f"Run the following command to retry:"
                           f"   bash download_icesat2.sh {date_range['start'].strftime('%Y%m%d')} {date_range['stop'].strftime('%Y%m%d')} {shapefile} {version} {job_id}")
    else:
        print(f"All downloads completed successfully.")
    
    return

if __name__ == "__main__":
    # define harmony request parameters
    yList = [2023]
    mList = [1]
    shape_list = ['coast_polygon_c_1.zip', 
                  'coast_polygon_c_2.zip', 
                  'coast_polygon_c_3.zip']
    version = '007'

    tm_fmt = '%Y-%m-%d %H:%M:%S'
    for year in yList:
        for month in mList:
            # get number of days in month
            nDays = calendar.monthrange(year, month)[1]

            for day in range(1, nDays + 1):
                # No available data before 2018-10-14
                if year == 2018 and month == 10 and day < 14:
                    continue

                # set date range for request
                date_range = {'start': dt.datetime(year, month, day, 0, 0, 0),
                              'stop': dt.datetime(year, month, day, 0, 0, 0) + dt.timedelta(days=1)}
                save_dir = os.path.join('data', f'ATL03_v{version}', f'G{year:04d}', f'{month:02d}', f'{day:02d}')
                tic = time.time()
                print(f"+{'-' * 84}+")
                print(f"INFO: Start downloading ICESat-2 ATL03 v{version} data for {year}-{month}-{day} at "
                      f"{dt.datetime.now().strftime(tm_fmt)}")

                # looping over shape files
                for shapefile in shape_list:
                    print(f"####Options:\n"
                          f"    time_range = {date_range['start'].strftime(tm_fmt)} to {date_range['stop'].strftime(tm_fmt)}\n"
                          f"    shapefile = {shapefile}\n"
                          f"    save_dir = {save_dir}")
                    
                    download_icesat2_atl03(date_range, 
                                           shapefile, 
                                           version=version, 
                                           save_dir=save_dir)

                print(f"INFO: Finished downloading ICESat-2 ATL03 v{version} data for {year}-{month}-{day} at "
                      f"{dt.datetime.now().strftime(tm_fmt)}")
                toc = time.time()
                print(f"Total time taken: {formatElapsedTime(toc - tic)}")
                print(f"+{'-' * 84}+")
            print(f"* Completed all downloads for {year}-{month}!")