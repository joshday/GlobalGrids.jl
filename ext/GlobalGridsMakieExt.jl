module GlobalGridsMakieExt

using GlobalGrids, Makie, GeoInterface

GeoInterface.@enable_makie Makie GlobalGrids.AbstractGrid
GeoInterface.@enable_makie Makie GlobalGrids.Cell
GeoInterface.@enable_makie Makie GlobalGrids.LonLat

end
