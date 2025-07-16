#!/bin/bash

# filepath: generate_matrix.sh 
# This script generates a grid of points in 3D space with specified intervals.
# The output is formatted as CSV with four columns: x, y, z, and a constant value (0.0).
# {
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


# {
#     for x in $(seq -330.0 20.0 330.0); do
#         for y in $(seq -330.0 20.0 330.0); do
#             # Skip only if both x AND y are between -30 and -10
#             if (( $(echo "$x > -60 && $x < 20 && $y > -60 && $y < 20" | bc -l) )); then
#                 continue
#             fi
#             for z in $(seq -10.0 20.0 30.0); do
#                 printf "%.1f,%.1f,%.1f,0\n" $x $y $z
#             done
#         done
#     done
# } > cylinder_unit_cells_360.csv

# {
#     for x in $(seq -330 20 330); do
#         for y in $(seq -330 20 330); do
#             # Skip only if both x AND y are between -30 and -10
#             if (( $(echo "$x > -60 && $x < 20 && $y > -60 && $y < 20" | bc -l) )); then
#                 continue
#             fi
#             printf "%.1f,%.1f,0,0\n" $x $y
#         done
#     done
# } > cylinder_unit_cells_330_2D.csv



size=40  # Matrix size (40x40 for values 0-39)
rim_thickness=2
{
    for x in $(seq 0 1 39); do
        for y in $(seq 0 1 39); do
            # Check if position is within the rim (first/last 2 rows/columns)
            if [ $x -lt $rim_thickness ] || [ $x -ge $(($size - $rim_thickness)) ] || \
               [ $y -lt $rim_thickness ] || [ $y -ge $(($size - $rim_thickness)) ]; then
                printf "%d\t%d\t0\n" $x $y
            fi
        done
    done
} > cylinder_unit_ecm_2D_40_doble.tsv