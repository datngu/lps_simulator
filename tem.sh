download_attempts=5
for ((i=1; i<=$download_attempts; i++)); do
    wget ${url_path} -O ${sample_id}.cram

    actual_md5=$(md5sum ${sample_id}.cram | awk '{print $1}')

    if [ "$actual_md5" == "$md5" ]; then
        echo "MD5 sum matches! File downloaded successfully."
        break
    else
        echo "MD5 sum does not match. Retrying..."
        rm -f downloaded_file.zip
    fi
done