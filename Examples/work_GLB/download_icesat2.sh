#!/bin/bash
# Download ICESat-2 ATL03 data using Harmony API, required: 
# NASA Earthdata Login (~/_netrc)
# curl
user=$(basename "$HOME")
if [[ "$user" = "caolang" ]]; then
	abspath="."
elif [[ "$user" = "ac6vfo7a3a" ]]; then
	abspath="/work/share/ac6vfo7a3a/caolang/data/ICESat-2"
else
	abspath="."
fi

# Step 1: Parse start/end date (yyyyMMdd) and shape file
if [ "$#" -eq 3 ]; then
	start_date="$1"
	end_date="$2"
	shape="$3"
	if [[ ! "$start_date" =~ ^[0-9]{8}$ ]] || [[ ! "$end_date" =~ ^[0-9]{8}$ ]]; then
		echo "Usage: $0 [start_date(yyyyMMdd)] [end_date(yyyyMMdd)] [shapefile]" >&2
		exit 1
	fi
	start_iso="${start_date:0:4}-${start_date:4:2}-${start_date:6:2}T00:00:00.000Z"
	end_iso="${end_date:0:4}-${end_date:4:2}-${end_date:6:2}T00:00:00.000Z"

	# set variables
	time_range="\"$start_iso\":\"$end_iso\""
	shapefile="$shape"
	save_dir="$abspath/G${start_date:0:4}/${start_date:4:2}/${start_date:6:2}"
else
	# default variables
	time_range="\"2018-10-14T00:00:00.000Z\":\"2018-10-15T00:00:00.000Z\""
	shapefile="./coast_polygon_c_1.zip"
	save_dir="$abspath/G2018/10/14"
fi
version="006"

echo "Downloading ICESat-2 ATL03 data using:"
echo "    time_range=$time_range"
echo "    shapefile=$shapefile"
echo "    save_dir=$save_dir"

# Step 2: Build the URL for request
if [ "$version" == "006" ]; then
	cmr_id="C2596864127-NSIDC_CPRD"
elif [ "$version" == "007" ]; then
	cmr_id="C3326974349-NSIDC_CPRD"
else
	echo "Unsupported version: $version"
	exit 1
fi
variable="all"
url="https://harmony.earthdata.nasa.gov/$cmr_id/\
ogc-api-coverages/1.0.0/collections/$variable/coverage/\
rangeset?forceAsync=true&subset=time($time_range)&maxResults=1" #&maxResults=1

# Step 3: Submit the request and get the JSON response
response=$(curl -Lnbj -sS "$url" -F "shapefile=@$shapefile;type=application/shapefile+zip")
#response=$(curl -Lnbj -sS https://harmony.earthdata.nasa.gov/jobs/04078a79-dc2b-4186-9608-34b5f488d91e)
jobID=$(echo "$response" | jq -r '.jobID // empty')
if [ -z "$jobID" ]; then
	echo "ERROR: Job submission failed" >&2
	echo "$response" | jq
	exit 1
fi
echo "Job submitted successfully, jobID=$jobID"
echo "$response" | jq

# construct job status URL
numGrans=$(echo "$response" | jq -r '.numInputGranules')
job_url="https://harmony.earthdata.nasa.gov/jobs/$jobID"

# wait until status is successful
while :; do
	status_json=$(curl -Lnbj -sS "$job_url")
	status=$(echo "$status_json" | jq -r '.status')
	message=$(echo "$status_json" | jq -r '.message')
	progress=$(echo "$status_json" | jq -r '.progress')
	if [ "$status" = "successful" ]; then
		original_size=$(echo "$status_json" | jq -r '.originalDataSize')
		output_size=$(echo "$status_json" | jq -r '.outputDataSize')
		pct_reduction=$(echo "$status_json" | jq -r '.dataSizePercentChange')
		echo "Job completed successfully, achieved a $pct_reduction ($original_size -> $output_size)"
		break
	fi
	if [ "$status" = "failed" ]; then
		echo "$message" >&2
		exit 1
	fi
	echo "Job status: $status with progress: $progress%"
	sleep 60
done

# extract download links whose title or href matches ATL03_*.h5
mapfile -t hrefs < <(echo "$status_json" | \
jq -r '.links | .[] | select(.title | test("ATL03_.*\\.h5")) | .href' | \
tr -d '\r')

if [ "${#hrefs[@]}" -ne "$numGrans" ]; then
	echo "ERROR: Expected $numGrans links, but found ${#hrefs[@]}" >&2
	exit 1
fi

# Step 4: Download the results
mkdir -p "$save_dir"
url_list="$save_dir/download_urls.txt"
printf '%s\n' "${hrefs[@]}" >> "$url_list"
sort -u "$url_list" -o "$url_list"

format="%{http_code} %{time_total} %{size_download}\n"
h5_list=()
for href in "${hrefs[@]}"; do
	h5file="$save_dir/$(basename "$href")"
	h5_list+=("$h5file")
	metrics=$(curl -Lnbj -sS --fail --retry 3 --retry-delay 5 --write-out "$format" "$href" -o "$h5file")
	read -r http_code time_total size_download <<< "$metrics"
	if [ "$http_code" -ne 200 ]; then
		echo "    Failed to download: $h5file" >&2
		rm -f "$h5file"
	fi
	time_min=$(awk -v t="$time_total" 'BEGIN{printf "%.1f", t/60}')
	size_mib=$(awk -v b="$size_download" 'BEGIN{printf "%.1f", b/1024/1024}')
	printf "Downloaded %s (%.1f MiB) in %.1f minutes\n" "$h5file" "$size_mib" "$time_min"
done

# Step 5: Check the missing files
missing=()
for fnm in "${h5_list[@]}"; do
	if ! [ -f "$fnm" ]; then
		missing+=("$fnm")
	fi
done

if [ ${#missing[@]} -ne 0 ]; then
	echo "ERROR: Expected $numGrans files, but missing ${#missing[@]}:" >&2
	printf '    %s\n' "${missing[@]}" >&2
fi