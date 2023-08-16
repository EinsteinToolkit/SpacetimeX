#! /bin/bash

set -e

echo "Preparing..." >&2

dir="mec-4-t100"
prefix="h.t"
suffix=".ah3.gp"
control="mec-4-t100/BH_diagnostics.ah3.gp"
output="mec-4-t100.horizons"

backslash=\\

rm -f $output.png.*

{
    cat <<EOF
        set terminal png
        set title "Horizons for mec-4"
        set xlabel "x"
        set ylabel "y"
        set zlabel "z"
        set grid
        set size ratio -1
        set view 0, 0
EOF
    shopt -s nullglob
    files="$dir/$prefix?$suffix \
           $dir/$prefix??$suffix \
           $dir/$prefix???$suffix \
           $dir/$prefix????$suffix \
           $dir/$prefix?????$suffix \
           $dir/$prefix??????$suffix \
           $dir/$prefix???????$suffix \
           $dir/$prefix????????$suffix \
           $dir/$prefix?????????$suffix"
    for file in $files; do
        iter=$(echo $file | sed -e "s+$dir/$prefix\(.*\)$suffix+\1+")
        iterf=$(printf "%09d" $iter)
        time=$(echo $iter | awk '{ printf "%.3f", $1 * 0.8 / 524288; }')
        echo "Plotting iteration $iter, time $time..." >&2
        input1=$dir/$prefix$iter.ah1.gp
        input2=$dir/$prefix$iter.ah2.gp
        input3=$dir/$prefix$iter.ah3.gp
        cat <<EOF
            set output "$output.png.$iterf"
            set label 1 "t=$time" at -17,11,0
            splot [-15:15] [-10:10] $backslash
EOF
        sep=""
        if test -f $input1; then
            cat <<EOF
                $sep "$input1" using 4:5:6 title "individual horizon" with lines 1 $backslash
EOF
            sep=","
        fi
        if test -f $input2; then
            cat <<EOF
                $sep "$input2" using 4:5:6 title "" with lines 1 $backslash
EOF
            sep=","
        fi
        if test -f $input3; then
            if (( $iter < 52428800 )); then
                cat <<EOF
                    $sep "$input3" using 4:5:6 title "pretracking horizon" with lines 2 $backslash
EOF
            else
                cat <<EOF
                    $sep "$input3" using 4:5:6 title "common horizon" with lines 3 $backslash
EOF
            fi
            sep=","
        fi
        cat <<EOF
        
EOF
    done
} | gnuplot

echo "Creating movie..." >&2
convert $output.png.* $output.mng

echo "Cleaning up..." >&2
rm $output.png.*

echo "Done." >&2
