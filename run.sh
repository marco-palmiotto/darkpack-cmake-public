#! /bin/bash

# Function to generate a random number in the interval [min, max]
# random_real_in_range() {
#     local min=$1
#     local max=$2
#     awk -v min="$min" -v max="$max" 'BEGIN{srand(); print (min + rand() * (max - min))}'
# }

# # Function to generate a random number in the interval 10^[min, max]
# random_real_in_log_range() {
#     local min=$1
#     local max=$2
#     awk -v min="$min" -v max="$max" 'BEGIN{srand(); print 10^(min + rand() * (max - min))}'
# }

# # Function to generate a random number between $1 and $2 times $3
# random_real_in_range_ref() {
#     local min=$1
#     local max=$2
#     local ref=$3
#     awk -v min="$min" -v max="$max" -v ref="$ref" 'BEGIN{srand(); print ref*(min + rand() * (max - min))}'
# }

# # Intervals
# INTERVAL1_MIN=-4.0
# INTERVAL1_MAX=0.0

# random1=$(random_real_in_log_range $INTERVAL1_MIN $INTERVAL1_MAX)

# INTERVAL2_MIN=0.3
# INTERVAL2_MAX=3.0

# random2=$(random_real_in_range $INTERVAL2_MIN $INTERVAL2_MAX)

# INTERVAL3_FAC_MIN=0.25
# INTERVAL3_FAC_MAX=4.0

# random3=$(random_real_in_range_ref $INTERVAL3_FAC_MIN $INTERVAL3_FAC_MAX $random2)

# echo "./build/src/main/scalar_scan_cma-2dim.x src/models/scalar.lha out/scans/result.out $random1 $random2"

# ./build/src/main/scalar_scan_cma-2dim-loop.x \
#     src/models/scalar.lha \
#     out/scans/result.out \
#     0.3 2.13  \
#     2.3e+3 4.0e+3 10

./build/src/Test/scalar_computation.x \
    src/models/scalar.lha \
    out/scans/input_to_test_all.dat \
    out/scans/output_from_test_all.out