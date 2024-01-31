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
function climplot end

function climscatter end
function climscatter! end
function climsurface end
function climsurface! end
