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
    colorbar = true, name = string(DimensionalData.name(A)), kwargs...)
    
    # @assert A has only space dim
    # vmin = haskey(kwargs, :vmin) ? kwargs[:vmin] : quantile(data, 0.025)
    # vmax = haskey(kwargs, :vmax) ? kwargs[:vmax] : quantile(data, 0.975)

    fig = GeoMakie.Figure()
    ax = GeoMakie.GeoAxis(fig[1,1]; source, dest)
    if scatter 
        el = climscatter!(ax, A; kwargs...)
    else
        el = climsurface!(ax, A; kwargs...)
    end

    if colorbar
        cb = GeoMakie.Colorbar(fig[1, 2], el; label = name)
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
# Notice that `A` is not declared as `ClimArray`, but assumed to be.
# Duck-typing for Observables.
function climsurface!(ax, A; colormap = :dense, kwargs...)
    if hasdim(A, Coord)
        # TODO: @pkeil this is for you
        error("Surface plots for `Coord` arrays are not supported yet!")
    end
    # TODO: This `180` needs to depend on lon_0
    l = findfirst(>(180), dims(A, Lon).val)
    if !isnothing(l)
        B = GeoMakie.lift(A -> longitude_circshift(A, l), A)
    else
        B = GeoMakie.lift(identity, A)
    end
    lon = dims(B, Lon).val
    lat = dims(B, Lat).val
    data = GeoMakie.lift(B -> B.data, B)
    # TODO: Change this to `contourf`
    GeoMakie.surface!(ax, lon, lat, data; shading = false, colormap, kwargs...)
end


##########################################################################################
# # Observables overloads
##########################################################################################
# These overloads allows me to write generic code that does not care if input
# is observable or not. This leads to simple, clean, small, elegant code for
# animating fields by simply replacing them with their observables.

GeoMakie.lift(f, x::ClimArray) = f(x)
DimensionalData.dims(o::GeoMakie.Observable, args...) = dims(o.val, args...)
DimensionalData.name(o::GeoMakie.Observable, args...) = name(o.val, args...)
spacestructure(A::GeoMakie.Observable) = spacestructure(A.val)