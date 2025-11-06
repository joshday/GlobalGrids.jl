module GlobalGridsRastersExt

import Rasters: AbstractRaster
import GlobalGrids: cells, AbstractCell

cells(::Type{T}, ::Nothing, r::AbstractRaster{<:Any, 2}, res::Integer; kw...) where {T <: AbstractCell} =
    cells(T, nothing, (r, r.dims), res; kw...)

end
