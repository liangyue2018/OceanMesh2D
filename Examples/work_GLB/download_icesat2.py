import harmony
import datetime as dt
import time
import os

def print_results(job_summary):
    '''
    Prints results of Harmony job
    '''
    job_status = job_summary['status']
    if job_status.lower() != 'successful':
        raise RuntimeError('Job did not complete successfully!')
    original_size = job_summary['originalDataSize']
    output_size = job_summary['outputDataSize']
    pct_reduction = job_summary['dataSizePercentChange']
    print(f'Achieved a {pct_reduction} ({original_size} -> {output_size}).')

path = os.path.join('.','G2018','10','14')
os.makedirs(path, exist_ok=True)
harmony_client = harmony.Client()
cmr_id = 'C2596864127-NSIDC_CPRD' # '006'
shape_list = ['./coast_polygon_c_1.zip', 
              './coast_polygon_c_2.zip',
              './coast_polygon_c_3.zip']
date_range = {'start': dt.datetime(2018,10,14,0,0,0),
              'stop': dt.datetime(2018,10,15,0,0,0)}

for shape in shape_list:
    start_time = time.time()
    request = harmony.Request(collection=harmony.Collection(id=cmr_id),
                              temporal=date_range,
                              shape=shape,
                              )
    job_id = harmony_client.submit(request)
    print(f'Job ID: {job_id}')
    harmony_client.wait_for_processing(job_id, show_progress=True)
    job_summary = harmony_client.result_json(job_id)
    print_results(job_summary)
    futures = harmony_client.download_all(job_id, directory=path, overwrite=True)
    filelist = [f.result() for f in futures]
    for fout in filelist:
        fpath = os.path.dirname(fout)
        fname = os.path.basename(fout)
        os.rename(fout, os.path.join(fpath, f"ATL03_{fname.split('ATL03_')[1]}"))
    print(f"Finished download for {shape} using {(time.time() - start_time)/3600:.2f} hours")



