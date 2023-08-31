#!/bin/bash

filename=$1
count=0

while IFS= read -r line; do
    if [[ "$line" == *"END"* ]]; then
        ((count++))
    fi
done < "$filename"

echo "Number of lines containing 'END' in $filename: $count"

