# TODO: remove `earthscatter` and `earthsurface` in favor of a single `earthplot`.
# This uses only the in-place versions `earthscatter!/surface!`.
# TODO: rename `earthplot` to `climplot`

export earthscatter, earthscatter!, earthsurface!, earthsurface!

##########################################################################################
# # Scatter
##########################################################################################
"""
    earthscatter(A::ClimArray; kwargs...) → fig, ax, el, cb
Plot a `ClimArray` with space type `UnstructuredGrid` as a scatter plot with the color
of the points being the value of `A` at these points. This requires that `A` has only
one dimension `Coord`, the coordinates.

Keyword values are propagated to `scatter` and keys of interest are `cmap, vmin, vmax`.
"""
function earthscatter(A::ClimArray; 
        source = "+proj=longlat +datum=WGS84", dest = "+proj=eqearth", 
        colorbar = true, cbkwargs = NamedTuple(), 
        kwargs...
    )
    fig = GeoMakie.Figure()
    ax = GeoMakie.GeoAxis(fig[1,1]; source, dest)
    el = earthscatter!(ax, A; kwargs...)
    if colorbar
        cb = GeoMakie.Colorbar(fig[1, 2], el; cbkwargs...)
    else
        cb = nothing
    end
    return fig, ax, el, cb
end

function earthscatter!(ax, A::ClimArray; title = string(A.name), colormap = :dense, kwargs...)
    ax.title = title
    if hasdim(A, Coord)
        coords = dims(A, Coord)
        lon = [l[1] for l in coords]
        lat = [l[2] for l in coords]
        data = A.data
    elseif dims(A)[1] isa Lon
        londim = dims(A, Lon)
        latdim = dims(A, Lat)
        lon = [l for lat in latdim for l in londim]
        lat = [lat for lat in latdim for l in londim]
        data = vec(A.data)
    else
        error("Unknown spatial dimensions for input.")
    end
    GeoMakie.scatter!(ax, lon, lat; color = vec(A), colormap, kwargs...)
end

##########################################################################################
# # Surface
##########################################################################################
"""
    earthsurface(A::ClimArray; kwargs...) → fig, ax, el, cb
Plot a `ClimArray` with space type `UnstructuredGrid` as a scatter plot with the color
of the points being the value of `A` at these points. This requires that `A` has only
one dimension `Coord`, the coordinates.

Keyword values are propagated to `scatter` and keys of interest are `cmap, vmin, vmax`.
"""
function earthsurface(A::ClimArray; 
        coastlines = true, source = "+proj=longlat", 
        dest = "+proj=eqearth", coastkwargs = NamedTuple(), 
        colorbar = true, cbkwargs = NamedTuple(), 
        kwargs...
    )
    fig = GeoMakie.Figure()
    ax = GeoMakie.GeoAxis(fig[1,1]; coastlines, source, dest, coastkwargs)
    el = earthsurface!(ax, A; kwargs...)
    if colorbar
        cb = Colorbar(fig[1, 2], el; cbkwargs...)
    end
    return fig, ax, el, cb
end

function earthsurface!(ax, A::ClimArray; colormap = :dense, kwargs...)
    if hasdim(A, Coord)
        # TODO: @pkeil this is for you
        error("Surface plots for `Coord` arrays are not supported yet!")
    elseif dims(A)[1] isa Lon
        lon = dims(A, Lon).val
        lat = dims(A, Lat).val
        data = A.data
    else
        error("Unknown spatial dimensions for input.")
    end
    lon, lat, data = _getwrappeddata(lon, lat, data)
    ax.title = string(A.name)
    # TODO: Change this to `contourf`
    GeoMakie.surface!(ax, lon, lat, data; shading = false, colormap, kwargs...)
end

# TODO: extend this to generic central meridian `lon_0`
function _getwrappeddata(lon, lat, s)
    if lon[end] == 359.5 && lon[1] == 0.5
        s = copy(s)
        lon = collect(lon)
        # add one more entry to have perfect overlap in plotting
        x = s[1, :]
        s = vcat(Array(s), x')
        push!(lon, 360.5)
    end
    if any(>(180), lon) # we need to circshift to work well with fucking Proj4.
        wi = findfirst(>(180), lon)
        s = circshift(s, wi)
        lon = circshift(wrap_lon.(lon), wi)
    end
    return lon, lat, s
end
