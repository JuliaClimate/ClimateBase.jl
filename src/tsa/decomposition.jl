#=
timeseries decompositions, i.e. separating some data into components
=#
# TODO: Finish this
seasonal_decomposition(A::AbDimArray, b) = seasonal_decomposition(dims(A, Time), A, b)

function seasonal_decomposition(t, A::AbDimArray, fs::Vector)
    # make Sinusoidal method
end
