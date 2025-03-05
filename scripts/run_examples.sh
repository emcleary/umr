#!/bin/env bash

dir="examples_testing"

if [ -d $dir ]; then
    rm -rfv $dir
fi

mkdir $dir


function run_ex04 {
    local alg="$1"
    local density="$2"
    local output="${dir}/ex04_${alg}_${density}"
    echo $output

    rm -f data*.vtk
    ./bin/examples/ex04_mesh_sources -r $alg -d $density > ${output}.out 2> ${output}.err

    mv data_000_ex04_mesh_sources_none.vtk ${dir}/ex04_none_${alg}_${density}.vtk
    mv data_001_ex04_mesh_sources_points.vtk ${dir}/ex04_points_${alg}_${density}.vtk
    mv data_002_ex04_mesh_sources_lines.vtk ${dir}/ex04_lines_${alg}_${density}.vtk
    mv data_003_ex04_mesh_sources_arc.vtk ${dir}/ex04_arc_${alg}_${density}.vtk
}

function run_ex06 {
    local alg="$1"
    local angle="$2"
    local density="$3"
    local output="${dir}/ex06_${alg}_${angle}_${density}"
    echo $output

    rm -f data*.vtk
    ./bin/examples/ex06_shapes -r $alg -a $angle -d $density > ${output}.out 2> ${output}.err

    mv ex06_shape_circles.vtk ${dir}/ex06_shape_circles_${alg}_${angle}_${density}.vtk
    mv ex06_shape_quads.vtk ${dir}/ex06_shape_quad_${alg}_${angle}_${density}.vtk
    mv ex06_shape_triangles.vtk ${dir}/ex06_shape_triangles_${alg}_${angle}_${density}.vtk
    mv ex06_shape_arbitrary.vtk ${dir}/ex06_shape_arbitrary_${alg}_${angle}_${density}.vtk
}

function run_ex07 {
    local alg="$1"
    local angle="$2"
    local output="${dir}/ex07_${alg}_${angle}"
    echo $output

    rm -f data*.vtk
    ./bin/examples/ex07_corners -r $alg -a $angle > ${output}.out 2> ${output}.err

    mv ex07_corners.vtk ${dir}/ex07_${alg}_${angle}.vtk
}


for density in 0.20 100; do
    run_ex06 chew-uniform 30 $density
    for alg in ruppert chew-nonuniform; do
	for angle in 20 25; do
	    run_ex06 $alg $angle $density
	done
    done
done


run_ex07 chew-uniform 10 # same for all angles
for alg in ruppert chew-nonuniform; do
    for angle in `seq 1 30`; do
	run_ex07 $alg $angle
    done
done


errors=false
for err in "${dir}/ex*.err"; do
    if [ -s "$err" ]; then
	echo "Error! See $err"
	errors=true
    else
	rm $err
    fi
done

if [ $errors == true ]; then
    echo "Errors found"
else
    echo "No errors found!"
fi
