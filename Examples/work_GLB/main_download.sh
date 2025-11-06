#!/bin/bash
# main script to download ICESat-2 ATL03 data from Oct 2018, required:
# download_icesat2.sh

# for year in {2018..2025}; do
for year in 2018; do
   	# for month in {01..12}; do
   	for month in 10; do
		# for day in {01..31}; do
		for day in {14..31}; do
			# check if the date is valid
			if ! date -d "${year}-${month}-${day}" >/dev/null 2>&1; then
				continue
			fi

			# set start/end date and version
			start_date="${year}${month}${day}"
			end_date=$(date -d "${year}-${month}-${day} + 1 day" +"%Y%m%d")
			version="007"
			tic=$(date +%s)
			echo "------------------------------------------------------------"
			echo "INFO: Start downloading ICESat-2 ATL03 v$version data for $start_date to $end_date at $(date "+%Y-%m-%d %H:%M:%S")"

			# loop over shape files
			for index in {1..3}; do
				shape_file="coast_polygon_c_${index}.zip"
				bash download_icesat2.sh "$start_date" "$end_date" "$shape_file" "$version"
				rc=$?
				if [ $rc -ne 0 ]; then
					exit $rc
				fi
			done

			echo "INFO: Finished downloading ICESat-2 ATL03 v$version data for $start_date to $end_date at $(date "+%Y-%m-%d %H:%M:%S")"
			toc=$(date +%s)
			elapsed=$(awk -v t="$((toc - tic))" 'BEGIN{printf "%.1f", t/3600}')
			echo "Total time taken: $elapsed hours"
			echo "------------------------------------------------------------"
		done
    done
done