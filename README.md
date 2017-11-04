3D Volumetric Model Analyzer
============================

This program was written to answer the following question:

> Given a 3D mesh model of the blood vessels in a rat's brain, what is the
> distance covered if all the vessels were laid out in a straight line.

Related questions of no particular interest included: how many rats died in the
creation of this model, what crazy person would want to lay out all those blood
vessels in a straight line, and why are your hands covered in blood and rat fur?
At least I hope it's rat fur...

Anyway, the perfect answer would involve calculating and measuring the medial
axis for the entire model (e.g. using the
[shrinking ball](https://vimeo.com/84859998) algorithm). Another approach is to
convert the mesh model into a volume (voxel) representation, and apply a
thinning algorithm until each voxel is left with just two neighbors on average.
Knowing the exact size of each voxel, it is then fairly simple to approximate
the medial axis distance.

`binvoxinfo.go` handles the last part of the problem. It expects a thinned voxel
model as input. Conversion to a volume representation and thinning are handled
by Patrick Min's [binvox](http://www.patrickmin.com/binvox/) and
[thinvox](http://www.patrickmin.com/thinvox/) programs.

Example
-------

```
go build binvoxinfo.go
binvox -d 512 -cb -rotx -aw -c that_poor_rat.stl
thinvox that_poor_rat.binvox
binvoxinfo thinned.binvox
```

Adjust binvox arguments as needed. Higher dimensions will produce more accurate
distance measurements. Ideally, use the maximum of 1024. On Windows, the 32-bit
thinvox.exe fails to allocate memory for dim > 834. Use the 64-bit Linux version
instead or pass `-d 834` to binvox.
