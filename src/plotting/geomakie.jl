# TODO: remove `climscatter` and `climsurface` in favor of a single `earthplot`.
# This uses only the in-place versions `climscatter!/surface!`.
# TODO: rename `earthplot` to `climplot`

export climscatter, climscatter!, climsurface!, climsurface, climplot

"""
    climplot(A::ClimArray; kwargs...) → fig, ax, el, cb
Main plotting function that dispatches to [`climscatter!`](@ref) if `A` has
an [`UnstructuredGrid`](@ref) dimension, or to [`climsurface!`](@ref) for [`LonLatGrid`](@ref).

Return the figure, axis, plotted element, and colorbar.

Optionally you can provide keyword `scatter = true` to force using the scatterplot.
"""
function climplot(A, args...; scatter = spacestructure(A) == UnstructuredGrid(),
    source = "+proj=longlat +datum=WGS84", dest = "+proj=eqearth", 
    colorbar = true, kwargs...)
    
    # @assert A has only space dim
    
    # vmin = haskey(kwargs, :vmin) ? kwargs[:vmin] : quantile(data, 0.025)
    # vmax = haskey(kwargs, :vmax) ? kwargs[:vmax] : quantile(data, 0.975)
    fig = GeoMakie.Figure()
    ax = GeoMakie.GeoAxis(fig[1,1]; source, dest, title = string(A.name))
    if scatter 
        el = climscatter!(ax, A; kwargs...)
    else
        el = climsurface!(ax, A; kwargs...)
    end

    if colorbar
        cb = GeoMakie.Colorbar(fig[1, 2], el; cbkwargs...)
        if A.attrib isa Dict && haskey(A.attrib, "units")
            cb.label(A.attrib["units"])
        end
    else
        cb = nothing
    end
    return fig, ax, el, cb
end

##########################################################################################
# # Scatter
##########################################################################################
function climscatter!(ax, A::ClimArray; title = string(A.name), colormap = :dense, kwargs...)
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
    climsurface(A::ClimArray; kwargs...) → fig, ax, el, cb
Plot a `ClimArray` with space type `UnstructuredGrid` as a scatter plot with the color
of the points being the value of `A` at these points. This requires that `A` has only
one dimension `Coord`, the coordinates.

Keyword values are propagated to `scatter` and keys of interest are `cmap, vmin, vmax`.
"""
function climsurface(A::ClimArray; 
        coastlines = true, source = "+proj=longlat", 
        dest = "+proj=eqearth", coastkwargs = NamedTuple(), 
        colorbar = true, cbkwargs = NamedTuple(), 
        kwargs...
    )
    fig = GeoMakie.Figure()
    ax = GeoMakie.GeoAxis(fig[1,1]; coastlines, source, dest, coastkwargs)
    el = climsurface!(ax, A; kwargs...)
    if colorbar
        cb = GeoMakie.Colorbar(fig[1, 2], el; cbkwargs...)
    end
    return fig, ax, el, cb
end

function climsurface!(ax, A::ClimArray; colormap = :dense, kwargs...)
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
