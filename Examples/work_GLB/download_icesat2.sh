#!/usr/bin/env bash
# Download ICESat-2 ATL03 data using Harmony API, required: 
# NASA Earthdata Login (~/_netrc)
# curl
#
# Step 1: Parse start/end date (yyyyMMdd), shape file and version
start_date="$1"
end_date="$2"
shapefile="$3"
version="$4"
if ! [[ "$start_date" =~ ^[0-9]{8}$ ]] || ! [[ "$end_date" =~ ^[0-9]{8}$ ]]; then
	echo "ERROR: Date format: [start_date(yyyyMMdd)] [end_date(yyyyMMdd)]" >&2
	exit 1
fi
if ! [[ -f "$shapefile" ]]; then
	echo "ERROR: Shapefile not found: $shapefile" >&2
	exit 1
fi
if [ "$version" = "006" ]; then
	cmr_id="C2596864127-NSIDC_CPRD"
elif [ "$version" = "007" ]; then
	cmr_id="C3326974349-NSIDC_CPRD"
else
	echo "ERROR: Unsupported version: $version" >&2
	exit 1
fi
start_iso="${start_date:0:4}-${start_date:4:2}-${start_date:6:2}T00:00:00.000Z"
end_iso="${end_date:0:4}-${end_date:4:2}-${end_date:6:2}T00:00:00.000Z"

# set variables
time_range="\"$start_iso\":\"$end_iso\""
save_dir="./data/ATL03_v$version/G${start_date:0:4}/${start_date:4:2}/${start_date:6:2}"

echo "####Options:"
echo "    time_range=$time_range"
echo "    shapefile=$shapefile"
echo "    save_dir=$save_dir"

# Step 2: Build the URL for request
variable="all"
url="https://harmony.earthdata.nasa.gov/$cmr_id/\
ogc-api-coverages/1.0.0/collections/$variable/coverage/\
rangeset?forceAsync=true&subset=time($time_range)&skipPreview=true" #&maxResults=1

if [ "$#" -eq 4 ]; then
	# Step 3: Submit the request and get the JSON response
	response=$(curl -Lnbj -sS "$url" -F "shapefile=@$shapefile;type=application/shapefile+zip")
	sleep 3
	jobID=$(echo "$response" | jq -r '.jobID // empty')
elif [ "$#" -eq 5 ]; then
	# Get jobID from command line argument
	jobID="$5"
	response=$(curl -Lnbj -sS "https://harmony.earthdata.nasa.gov/jobs/$jobID")
	sleep 3
else
	echo "ERROR: Usage: $0 <start_date> <end_date> <shapefile> <version> [<jobID>]" >&2
	exit 1
fi

if [ -z "$jobID" ]; then
	echo "ERROR: Job submission failed" >&2
	echo "$response" | jq
	exit 1
fi
echo "Job submitted successfully, jobID=$jobID"
#echo "$response" | jq

# construct job status URL
numGrans=$(echo "$response" | jq -r '.numInputGranules')
job_url="https://harmony.earthdata.nasa.gov/jobs/$jobID"
set -x

# wait until status is successful
while :; do
	tmp=$(mktemp)
	http_code=$(curl -Lnbj -sS -w "%{http_code}" -o "$tmp" "$job_url")
	curl_rc=$?
	if [ "$http_code" -ne 200 ] || [ $curl_rc -ne 0 ]; then
		sleep 60
		continue
	fi
	status_json=$(cat "$tmp"); rm -f "$tmp"
	status=$(echo "$status_json" | jq -r '.status // empty')
	message=$(echo "$status_json" | jq -r '.message // empty')
	progress=$(echo "$status_json" | jq -r '.progress // empty')
	if [ "$status" = "successful" ]; then
		original_size=$(echo "$status_json" | jq -r '.originalDataSize // empty')
		output_size=$(echo "$status_json" | jq -r '.outputDataSize // empty')
		pct_reduction=$(echo "$status_json" | jq -r '.dataSizePercentChange // empty')
		echo "Job completed successfully, achieved a $pct_reduction ($original_size -> $output_size)"
		break
	fi
	if [ "$status" = "failed" ] || [ "$status" = "complete_with_errors" ]; then
		echo "ERROR: $message" >&2
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

# filter granules with segment number determined by shape file
case "$(basename "$shapefile")" in
	*_1.zip) segnum=("02" "03" "05" "06") ;;
	*_2.zip) segnum=("14" "01" "07" "08") ;;
	*_3.zip) segnum=("09" "10" "12" "13") ;;
	*) echo "ERROR: Unknown shapefile: $shapefile" >&2; exit 1 ;;
esac

filtered_hrefs=()
for href in "${hrefs[@]}"; do
	fnm=$(basename "$href")
	ss=$(awk -F'_' '{print substr($3,7,2)}' <<< "$fnm")
	for sn in "${segnum[@]}"; do
		if [ "$ss" = "$sn" ]; then
			filtered_hrefs+=("$href")
			break
		fi
	done
done
hrefs=("${filtered_hrefs[@]}")
if [ "${#hrefs[@]}" -eq 0 ]; then
	echo "INFO: No valid granules found for $shapefile" >&2
	exit 0
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
	if [ -s "$h5file" ]; then
		continue
	fi
	metrics=$(curl -Lnbj -sS --fail \
		--connect-timeout 30 \
		--retry 5 --retry-delay 10 --retry-max-time 600 \
		--speed-limit 1024 --speed-time 30 \
		--write-out "$format" \
		"$href" -o "$h5file")
	curl_rc=$?
	read -r http_code time_total size_download <<< "$metrics"
	if [ $curl_rc -ne 0 ] || [ "$http_code" -ne 200 ]; then
		echo "ERROR: Failed to download: $h5file" >&2
		rm -f "$h5file"
		continue
	fi
	time_min=$(awk -v t="$time_total" 'BEGIN{printf "%.1f", t/60}')
	size_mib=$(awk -v b="$size_download" 'BEGIN{printf "%.1f", b/1024/1024}')
	printf "Downloaded %s (%.1f MiB) in %.1f minutes\n" "$h5file" "$size_mib" "$time_min"
done

# Step 5: Check the validity of downloaded files
if ! command -v h5ls >/dev/null 2>&1; then
	if [ "$(basename "$HOME")" = "caolang" ]; then
		eval "$(conda shell.bash hook)" && conda activate icepyx_env
	elif [ "$(basename "$HOME")" = "ac6vfo7a3a" ]; then
		module load mathlib/hdf5/1.12.2-intel-21
	fi
fi
invalid=()
for fnm in "${h5_list[@]}"; do
	if ! h5ls "$fnm" >/dev/null 2>&1; then
		invalid+=("$fnm")
		rm -f "$fnm"
	fi
done

if [ ${#invalid[@]} -ne 0 ]; then
	echo "ERROR: Expected ${#hrefs[@]} files, but has ${#invalid[@]} invalid file(s):" >&2
	printf '    %s\n' "${invalid[@]}" >&2
	echo "Run the script to download again:" >&2
	echo "    bash $0 $start_date $end_date $shapefile $version $jobID" >&2
fi
