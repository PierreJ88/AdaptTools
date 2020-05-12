# Explaination to use mshmet for 2D geophysical applications

## Software needed :
* medit (optional for vizualization)
* mmg(2d)

## STEP 0 : Start with visualization mesh and P1 solution associated

Open the bilayer_visu.mesh with medit to vizualize it and (play with  M and L keyboard, Z and shift+Z for Zoom/DeZoom)
```medit bilayer_visu.mesh```

## STEP 1 : Call mshmet with the corresponding option
```mshmet -g (-i) -e XX -rmin XX -order XX -freq XX (-ppw XX) bilayer_visu.mesh```

* -g > use gradient instead of hessien here it is needed for geophysical applications
* -i > optional to keep isotropic mesh or not
* -e > Parameter to highlight the interface component (default = 1)
* -rmin > Reduction factor at interfaces (default = 1 -> no reduction)
* -order > polynomial order desired (default = 1)
* -freq > frequency of the source used (default = 10.0Hz)
* -ppw > Optional parameter to force the number of point per wavelength

/!\ Here the velocity model included in bilayer_visu.sol is used in the code by the fact that the solution and the mesh have the same name.


## STEP 1.5 : Vizualize the metric size map generated

A metric called bilayer_visu.met.sol has been generated.

To visualize it :
```cp bilayer_visu.mesh bilayer_visu.met.mesh```
```medit bilayer_visu.met.mesh```

## STEP 2 : Generate the new mesh
```mmg2d bilayer_visu.met.mesh -hgrad -1```

The final desired mesh has been generated under : bilayer_visu.met.o.mesh. (a .sol file has also been generated just ignore it)

## STEP 2.5 : Visualize the final mesh
```medit bilayer_visu.met.o.mesh```


# Sum up the entire workflow :

```mshmet -g (-i) -e XX -rmin XX -order XX -freq XX (-ppw XX) bilayer_visu.mesh```
```cp bilayer_visu.mesh bilayer_visu.met.mesh```
```mmg2d  bilayer_visu.met.mesh -hgrad -1```

(see bilayer_mesh.bash script)
