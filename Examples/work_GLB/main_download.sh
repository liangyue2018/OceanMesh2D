#!/bin/bash
# main script to download ICESat-2 ATL03 data from Oct 2018, required:
# download_icesat2.sh

# for year in {2018..2025}; do
for year in 2018; do
   	# for month in {01..12}; do
   	for month in 10; do
		# for day in {01..31}; do
		for day in {15..20}; do
			# check if the date is valid
			if ! date -d "${year}-${month}-${day}" >/dev/null 2>&1; then
				continue
			fi

			# limit to Oct 14, 2018 for test
			# if [ "$year" -ne 2018 ] || [ "$month" -ne 10 ] || [ "$day" -ne 14 ]; then
			# 	continue
			# fi

			# specify jobID manually (optional)
			# jobID=("8d5b27bf-7e6d-4313-947f-40e5d32a7ae4" 
			#        "3a6bb08d-4aeb-4aaa-bb3d-184ceb3b89ff" 
			# 	   "b68d23ad-f767-4199-bc0e-26c50edbfb32")

			# set start and end date
			start_date="${year}${month}${day}"
			end_date=$(date -d "${year}-${month}-${day} + 1 day" +"%Y%m%d")
			tic=$(date +%s)
			echo "------------------------------------------------------------"
			echo "INFO: Start downloading ICESat-2 ATL03 data for $start_date to $end_date at $(date "+%Y-%m-%d %H:%M:%S")"

			# loop over shape files
			for index in {1..3}; do
				shape_file="coast_polygon_c_${index}.zip"
				bash download_icesat2.sh "$start_date" "$end_date" "$shape_file"
				#bash download_icesat2.sh "$start_date" "$end_date" "$shape_file" "${jobID[$((index-1))]}"
				rc=$?
				if [ $rc -ne 0 ]; then
					exit $rc
				fi
			done

			echo "INFO: Finished downloading ICESat-2 ATL03 data for $start_date to $end_date at $(date "+%Y-%m-%d %H:%M:%S")"
			toc=$(date +%s)
			elapsed=$(awk -v t="$((toc - tic))" 'BEGIN{printf "%.1f", t/3600}')
			echo "Total time taken: $elapsed hours"
			echo "------------------------------------------------------------"
		done
    done
done