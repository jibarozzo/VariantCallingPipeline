#!/bin/bash

#input_file="fastqc.err"

# Check if the input file is provided as an argument

if [ -z "$1" ]; then
  echo "Usage: $0 <input_file> <output_file>"
  exit 1
fi

## Check if the output file is provided as an argument
if [ -z "$2" ]; then
  echo "Usage: $0 <input_file> <output_file>"
  exit 1
fi

input_file="$1"
output_file="$2"

# Process the input file
while IFS= read -r line; do
  if [[ $line == "Failed to process"* ]]; then
    # Extract the part of the line after the last '/'
    extracted_part="${line##*/}" # Selecting .fq.gz files
   #grep "$extracted_part.*\.fq.gz$" >> "$output_file"     
    echo "$extracted_part" >> "$output_file"
  fi
done < "$input_file"
