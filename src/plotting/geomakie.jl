export climscatter, climscatter!, climsurface!, climsurface, climplot

"""
    climplot(A::ClimArray; kwargs...) → fig, ax, el, cb
Main plotting function that dispatches to [`climscatter!`](@ref) if `A` has
an [`CoordinateSpace`](@ref) dimension, or to [`climsurface!`](@ref) for [`OrthogonalSpace`](@ref).

Return the figure, axis, plotted element, and colorbar.

Optionally you can provide keyword `scatter = true` to force using the scatterplot.

Plotting from ClimateBase.jl also works with `Observable`s that enclose a `ClimArray`.
You can update the values of the observable with another `ClimArray` with the same spatial
dimension, and the plot will be updated. See documentation online for examples.
"""
function climplot(A, args...; scatter = spacestructure(A) == CoordinateSpace(),
    source = "+proj=longlat +datum=WGS84", dest = "+proj=eqearth",
    colorbar = true, name = string(DimensionalData.name(A)), kwargs...)

    # TODO: Perhaps setting custom colorrange leads to better plots?
    # vmin = haskey(kwargs, :vmin) ? kwargs[:vmin] : quantile(data, 0.025)
    # vmax = haskey(kwargs, :vmax) ? kwargs[:vmax] : quantile(data, 0.975)

    fig = GeoMakie.Figure()

    ax = GeoMakie.GeoAxis(fig[1,1];
        source, dest, title = GeoMakie.Observable(""), coastlines = true,
        xtickformat = GeoMakie.geoformat_ticklabels,
        ytickformat = GeoMakie.geoformat_ticklabels,
        args...,
    )

    el = scatter ? climscatter!(ax, A; kwargs...) : climsurface!(ax, A; kwargs...)
    if colorbar
        cb = GeoMakie.Colorbar(fig[1, 2], el; label = name)
    else
        cb = nothing
    end
    return fig, ax, el, cb
end

function has_only_space_dim(A)
    d = dims(A)
    if length(d) == 1
        return d[1] isa Coord
    elseif length(d) == 2
        dimtypes = basetypeof.(d)
        return (Lon ∈ dimtypes) && (Lat ∈ dimtypes)
    else
        return false
    end
end

##########################################################################################
# # Scatter
##########################################################################################
# Notice that `A` is not declared as `ClimArray`, but assumed to be.
# Duck-typing for Observables.
function climscatter!(ax, A; colormap = :dense, kwargs...)
    @assert has_only_space_dim(A) "ClimArray must have only spatial dimension!"
    if hasdim(A, Coord)
        lonlat = dims(A, Coord).val
    elseif dims(A)[1] isa Lon
        londim = dims(A, Lon)
        latdim = dims(A, Lat)
        lonlat = [GeoMakie.Point2f0(l,lat) for lat in latdim for l in londim]
    else
        error("Unknown spatial dimensions for input.")
    end
    data = GeoMakie.lift(A -> vec(gnv(A)), A)
    GeoMakie.scatter!(ax, lonlat; color = data, colormap, kwargs...)
end

##########################################################################################
# # Surface
##########################################################################################
# TODO: Change this to `contourf`

# Notice that `A` is not declared as `ClimArray`, but assumed to be.
# Duck-typing for Observables.
function climsurface!(ax, A; colormap = :dense, kwargs...)
    @assert has_only_space_dim(A) "ClimArray must have only spatial dimension!"
    if hasdim(A, Coord)
        # TODO: @pkeil this is for you
        error("Surface plots for `Coord` arrays are not supported yet!")
    end
    # TODO: This needs to depend on lon_0 or use the "meridian cut" function
    B = GeoMakie.lift(A -> longitude_circshift(A), A)
    lon = dims(B, Lon).val
    lat = dims(B, Lat).val
    data = GeoMakie.lift(B -> gnv(B), B)
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