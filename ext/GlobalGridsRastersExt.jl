module GlobalGridsRastersExt

import Rasters: AbstractRaster, X, Y, At, Near, dims
import GlobalGrids: h3cells, h3join, H3Cell, LonLat
import GeoInterface as GI


h3join(r::AbstractRaster{T, 2}, res::Integer = 10) where {T} = h3join((r, r.dims), res)

# function h3cells(r::AbstractRaster{T, 2}, res::Integer = 10) where {T}
#     # Step 1) Initialize output hexagons to cover the extent of the raster
#     # all_cells = h3cells(GI.extent(r), res; containment)
#     out = Dict{H3Cell, Vector{T}}()

#     coords = [LonLat(lon, lat) for lon in r.dims[1].val, lat in r.dims[2].val]

#     # Step 2) Iterate through the raster and populate hexagons
#     for (ll, val) in zip(coords, r)
#         if !ismissing(val)
#             cell = H3Cell(ll, res)
#             v = get!(out, cell, T[])  # get! is required depending on containment mode
#             push!(v, val)
#         end
#     end



#     # for lon in r.dims[1], lat in r.dims[2]
#     #     val = r[X(At(lon)), Y(At(lat))]
#     #     cell = H3Cell(LonLat(lon, lat), res)
#     #     v = get!(out, cell, T[])  # get! is required depending on containment mode
#     #     (!ismissing(val) || !dropmissing) && push!(v, val)
#     # end

#     # # Step 3) For cells without data, find the closest point from the raster
#     # for (k, v) in out
#     #     if isempty(v)
#     #         lon, lat = GI.centroid(k)
#     #         val = r[X(Near(lon)), Y(Near(lat))]
#     #         (!ismissing(val) || !dropmissing) && push!(v, val)
#     #     end
#     # end

#     # Step 4) drop empty vectors (only happens for missing values)
#     # dropmissing && filter!(kv -> !isempty(kv[2]), out)
#     return out
# end

end
