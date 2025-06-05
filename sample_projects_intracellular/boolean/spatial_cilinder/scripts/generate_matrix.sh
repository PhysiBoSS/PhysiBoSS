#!/bin/bash

# filepath: generate_matrix.sh 
# This script generates a grid of points in 3D space with specified intervals.
# The output is formatted as CSV with four columns: x, y, z, and a constant value (0.0).
{
#     for x in $(seq -670.0 20.0 670.0); do
#         for y in $(seq -670.0 20.0 670.0); do
#             # Skip only if both x AND y are between -30 and -10
#             if (( $(echo "$x > -60 && $x < 20 && $y > -60 && $y < 20" | bc -l) )); then
#                 continue
#             fi
#             for z in $(seq -50.0 20.0 30.0); do
#                 printf "%.1f,%.1f,%.1f,0.0\n" $x $y $z
#             done
#         done
#     done
# } > cylinder_unit_cells.csv



    for x in $(seq -330.0 20.0 330.0); do
        for y in $(seq -330.0 20.0 330.0); do
            # Skip only if both x AND y are between -30 and -10
            if (( $(echo "$x > -60 && $x < 20 && $y > -60 && $y < 20" | bc -l) )); then
                continue
            fi
            for z in $(seq -50.0 20.0 30.0); do
                printf "%.1f,%.1f,%.1f,0.0\n" $x $y $z
            done
        done
    done
} > cylinder_unit_cells_360.csv



