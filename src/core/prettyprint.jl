# printing of a single clim array
function Base.show(io::IO, ::MIME"text/plain", A::ClimArray)
    summary(io, A)
    print(io, "and")
    printstyled(io, " data: "; color=:green)
    dataA = data(A)
    print(io, summary(dataA), "\n")
    x = 2length(dims(A)) + attriblength(A.attrib) + 5
    custom_show(io, data(A), x)
end

attriblength(d::AbstractDict) = length(d)
attriblength(d) = 0

# Define summary
function Base.summary(io::IO, A::ClimArray)
    printstyled(io, "ClimArray"; color=:blue)
    if A.name ≠ Symbol("")
        print(io, " (named ")
        printstyled(io, A.name; color=:blue)
        print(io, ")")
    end

    print(io, " with dimensions:\n")
    for d in dims(A)
        print(io, " ", d, "\n")
    end
    if !isnothing(A.attrib)
        printstyled(io, "attributes: "; color=:magenta)
        show(io, MIME"text/plain"(), A.attrib)
        print(io, '\n')
    end
end

# Thanks to Michael Abbott for the following function
function custom_show(io::IO, A::AbstractArray{T,0}, x) where T
    Base.show(IOContext(io, :compact => true, :limit => true), A)
end
function custom_show(io::IO, A::AbstractArray{T,1}, x) where T
    Base.show(IOContext(io, :compact => true, :limit => true, :displaysize => displaysize(io) .- (x, 0)), A)
end
function custom_show(io::IO, A::AbstractArray{T,2}, x) where T
    Base.print_matrix(IOContext(io, :compact => true, :limit => true, :displaysize => displaysize(io) .- (x, 0)), A)
end
function custom_show(io::IO, A::AbstractArray{T,N}, x) where {T,N}
    o = ones(Int, N-2)
    frame = A[:, :, o...]
    onestring = join(o, ", ")
    println(io, "[:, :, $(onestring)]")
    Base.print_matrix(
        IOContext(io, :compact => true, :limit=>true, :displaysize => displaysize(io) .- (x, 0)),
        frame)
    print(io, "\n[and ", prod(size(A,d) for d=3:N) - 1," more slices...]")
end

# printing of a vector of climarrays
function Base.show(io::IO, A::ClimArray)
    printstyled(io, "ClimArray"; color=:blue)
    if A.name ≠ Symbol("")
        print(io, " (named ")
        printstyled(io, A.name; color=:blue)
        print(io, ")")
    end
    s = join(string.(size(A)), "×")
    d = join([replace(string(basetypeof(d)), "Ti"=>"Time") for d in dims(A)], "×")
    print(io, " with $s $d data")
end
