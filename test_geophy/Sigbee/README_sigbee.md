# Explaination to use mshmet for 2D geophysical applications

## Software needed :
* medit (optional for vizualization)
* mmg(2d)

## STEP 0 : Start with visualization mesh and P1 solution associated

Open the sigbee_visu.mesh with medit to vizualize it and (play with  M and L keyboard, Z and shift+Z for Zoom/DeZoom)
```medit sigbee_visu.mesh```

## STEP 1 : Call mshmet with the corresponding option
```mshmet -g (-i) -e XX -rmin XX -order XX -freq XX (-ppw XX) sigbee_visu.mesh```

* -g > use gradient instead of hessien here it is needed for geophysical applications
* -i > optional to keep isotropic mesh or not
* -e > Parameter to highlight the interface component (default = 1)
* -rmin > Reduction factor at interfaces (default = 1 -> no reduction)
* -order > polynomial order desired (default = 1)
* -freq > frequency of the source used (default = 10.0Hz)
* -ppw > Optional parameter to force the number of point per wavelength

/!\ Here the velocity model included in sigbee_visu.sol is used in the code by the fact that the solution and the mesh have the same name.


## STEP 1.5 : Vizualize the metric size map generated

A metric called sigbee_visu.met.sol has been generated.

To visualize it :
```cp sigbee_visu.mesh sigbee_visu.met.mesh```
```medit sigbee_visu.met.mesh```

## STEP 2 : Generate the new mesh
```mmg2d sigbee_visu.met.mesh -hgrad -1```

The final desired mesh has been generated under : sigbee_visu.met.o.mesh. (a .sol file has also been generated just ignore it)

## STEP 2.5 : Visualize the final mesh
```medit sigbee_visu.met.o.mesh```


# Sum up the entire workflow :

```mshmet -g (-i) -e XX -rmin XX -order XX -freq XX (-ppw XX) sigbee_visu.mesh```
```cp sigbee_visu.mesh sigbee_visu.met.mesh```
```mmg2d  sigbee_visu.met.mesh -hgrad -1```

(see sigbee_mesh.bash script)
