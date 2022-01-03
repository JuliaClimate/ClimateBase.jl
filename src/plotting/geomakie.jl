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
    
    # TODO: Perhaps setting custom colorrange leads to better plots?
    # vmin = haskey(kwargs, :vmin) ? kwargs[:vmin] : quantile(data, 0.025)
    # vmax = haskey(kwargs, :vmax) ? kwargs[:vmax] : quantile(data, 0.975)

    fig = GeoMakie.Figure()
    ax = GeoMakie.GeoAxis(fig[1,1]; source, dest, title = GeoMakie.Observable(""))
    
    coastplot = lines!(ax, GeoMakie.coastlines(); color = :black, overdraw = true, coastkwargs...)
    translate!(coastplot, 0, 0, 99) # ensure they are on top of other plotted elements

    # TODO: @assert A has only space dimension

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
# Notice that `A` is not declared as `ClimArray`, but assumed to be.
# Duck-typing for Observables.
function climscatter!(ax, A; colormap = :dense, kwargs...)
    if hasdim(A, Coord)
        lonlat = dims(A, Coord).val
    elseif dims(A)[1] isa Lon
        londim = dims(A, Lon)
        latdim = dims(A, Lat)
        lonlat = [GeoMakie.Point2f0(l,lat) for lat in latdim for l in londim]
    else
        error("Unknown spatial dimensions for input.")
    end
    data = GeoMakie.lift(A -> vec(A.data), A)
    GeoMakie.scatter!(ax, lonlat; color = data, colormap, kwargs...)
end

##########################################################################################
# # Surface
##########################################################################################
# TODO: Change this to `contourf`

# Notice that `A` is not declared as `ClimArray`, but assumed to be.
# Duck-typing for Observables.
function climsurface!(ax, A; colormap = :dense, kwargs...)
    # TODO: @assert A has only space
    if hasdim(A, Coord)
        # TODO: @pkeil this is for you
        error("Surface plots for `Coord` arrays are not supported yet!")
    end
    # TODO: This needs to depend on lon_0 or use the "meridian cut" function
    B = GeoMakie.lift(A -> longitude_circshift(A), A)
    lon = dims(B, Lon).val
    lat = dims(B, Lat).val
    data = GeoMakie.lift(B -> B.data, B)
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