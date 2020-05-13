#!/bin/bash

# PARAMETERS
i_param="false"
e_param=0.5
rmin_param=4
order_param=4
freq_param=25

# FILES
mesh0=sigbee_visu.mesh
mesh1=sigbee_visu.met.mesh

# EXE
mshmet=your/path/to/mshmet
mmg2d=your/path/to/mmg2d

if [[ "$i_param" == "true" ]] ; then
    $mshmet -g -i -e $e_param -rmin $rmin_param -order $order_param -freq $freq_param $mesh0
fi


if [[ "$i_param" == "false" ]] ; then
    $mshmet -g -e $e_param -rmin $rmin_param -order $order_param -freq $freq_param $mesh0
fi

cp $mesh0 $mesh1
$mmg2d $mesh1
