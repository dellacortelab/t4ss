#!/bin/bash
#Developed by Ethan Smith + ChatGPT 5/23/23

if [ $# -ne 1 ]; then
  echo "Usage: $0 <input_file>"
  exit 1
fi

input_file="$1"
temp_file="between.pdb"
output_file="output.pdb"

# Remove the fifth letter in the fourth column
sed -E '/^ATOM/s/(^(.{18}).{3}).(.)/\1\3/' "$input_file" > "$temp_file"

# Add a space between the third and fourth letter in the fourth column
sed -E '/^ATOM/s/(^(.{18}).{2})(.)/\1 \3/' "$temp_file" > "$output_file"

# Clean up temporary file
rm "$temp_file"

echo "Output written to $output_file"

