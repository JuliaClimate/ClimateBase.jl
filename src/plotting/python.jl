# Uncomment the following to install the necessary packages for plotting.
# Pkg.add(["Conda", "PyCall"])
# using Conda
# Conda.add("cartopy")

using PyCall, PyPlot, ClimateBase
ccrs = pyimport("cartopy.crs")
LONLAT = ccrs.PlateCarree()
DEFPROJ = ccrs.Mollweide()

"""
    earthscatter(A::ClimArray, projection = ccrs.Mollweide(); kwargs...) → fig, ax, cb
Plot a `ClimArray` with space type `CoordinateSpace` as a scatter plot with the color
of the points being the value of `A` at these points. This requires that `A` has only
one dimension, the coordinates.

Keyword values are propagated to `scatter` and keys of interest are `cmap, vmin, vmax`.
"""
function earthscatter(A::ClimArray, projection = DEFPROJ; kwargs...)
    coords = dims(A, Coord)
    @assert length(dims(A)) == 1 "Input ClimArray must only have one `Coord` dimension."
    eqlon = [l[1] for l in coords]
    eqlat = [l[2] for l in coords]
    fig = figure()
    fig.tight_layout()
    ax = subplot(111, projection=projection)
    ax, cb = earthscatter!(ax, eqlon, eqlat, A.data, projection; kwargs...)
    if A.name ≠ Symbol("")
        ax.set_title(string(A.name))
    end
    return fig, ax, cb
end

function earthscatter!(ax, A, projection = DEFPROJ; kwargs...)
    coords = dims(A, Coord)
    lon = [l[1] for l in coords]
    lat = [l[2] for l in coords]
    earthscatter!(ax, lon, lat, A, projection; kwargs...)
end

function earthscatter!(ax, lon, lat, A, projection = DEFPROJ;
    coast = true, vmin = minimum(A), vmax = maximum(A), levels = 11,
    ticks = nothing, kwargs...)
    lon = Array(lon)
    lat = Array(lat)
    lvls = range(vmin, vmax, length = levels)
    cmap = matplotlib.cm.get_cmap(get(kwargs, :cmap, :viridis), length(lvls)-1)
    norm = matplotlib.colors.Normalize(vmin=lvls[1], vmax=lvls[end])
    cax, kw = matplotlib.colorbar.make_axes(ax,location="right",pad=0.02,shrink=0.8)
    cb = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm)

    ax.scatter(lon, lat; transform = LONLAT, c=A, kwargs..., cmap = cmap, vmin = vmin, vmax = vmax)
    coast && ax.coastlines()
    if projection == LONLAT
        ax.set_xticks(-180:60:180, crs=LONLAT)
        ax.set_yticks(-90:30:90, crs=LONLAT)
    end
    ax.gridlines(alpha = 0.25)
    cb.set_ticks(isnothing(ticks) ? lvls[1:max(1, levels÷5):end] : ticks)
    return ax, cb
end

"""
    earthsurface(A::ClimArray, projection = ccrs.Mollweide(); kwargs...) → fig, ax, cb
Plot a `ClimArray` with space type `OrthogonalSpace` as a surface plot with the color
of the surface being the value of `A`. This requires that `A` has exactly two dimensions,
`Lon, Lat`.

Keyword values are `coast=true` to enable coastlines and also `cmap, vmin, vmax, levels`.
"""
function earthsurface(A, projection = DEFPROJ;
    coast = true, kwargs...)
    fig = figure()
    fig.tight_layout()
    ax = subplot(111, projection=projection)
    ax, cb = earthsurface!(ax, A; kwargs...)
    if projection == LONLAT
        ax.set_xticks(-180:60:180, crs=LONLAT)
        ax.set_yticks(-90:30:90, crs=LONLAT)
    end
    if A.name ≠ Symbol("")
        ax.set_title(string(A.name))
    end
    return fig, ax, cb
end

function earthsurface!(ax, A::ClimArray; kwargs...)
    ax.gridlines(alpha = 0.25)
    lon, lat, s = _getwrappeddata(A)
    earthsurface!(ax, lon, lat, s; kwargs...)
end

function earthsurface!(ax, lon, lat, A;
    vmin = minimum(A), vmax = maximum(A), levels = 11,
    coast = true, specific_contour = nothing, kwargs...)

    lon, lat, A = _getwrappeddata(copy(Array(lon)), copy(Array(lat)), A)
    lvls = range(vmin, vmax, length = levels)

    # TODO: Add colorbar as a separate function
    lvls = range(vmin, vmax, length = levels)
    cmap = matplotlib.cm.get_cmap(get(kwargs, :cmap, :viridis), length(lvls)-1)
    # cmap.set_under("k")
    norm = matplotlib.colors.Normalize(vmin=lvls[1], vmax=lvls[end])
    cax, kw = matplotlib.colorbar.make_axes(ax,location="right",pad=0.02,shrink=0.8)
    cb = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm)
    coast && ax.coastlines(; color = "gray")

    ax.contourf(Array(lon), Array(lat), clamp.(A', vmin, vmax), lvls, cmap = cmap,
                transform=LONLAT, vmin=vmin, vmax = vmax)

    if specific_contour isa Array
        ax.contour(Array(lon), Array(lat), clamp.(A', vmin, vmax);
        levels = specific_contour, colors = "r", transform=LONLAT)
    end
    return ax, cb
end

function _getwrappeddata(A::ClimArray)
    lon = copy(Array(dims(A, Lon).val))
    lat = copy(Array(dims(A, Lat).val))
    s = copy(Array(A.data)) # assumes data have Lon x Lat dimensions
    _getwrappeddata(lon, lat, s)
end

function _getwrappeddata(lon, lat, s)
    if lon[end] == 359.5 && lon[1] == 0.5
        # add one more entry to have perfect overlap in plotting
        x = s[1, :]
        s = vcat(Array(s), x')
        push!(lon, 360.5)
    end
    return lon, lat, s
end
