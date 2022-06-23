#=
Code related with input output (IO) of .nc files directly to/from ClimArrays
utilizing the NCDatasets.jl package and a buttload of convenience code.
An initial version of parts of this code was taken from:
https://github.com/rafaqz/GeoData.jl
=#
using NCDatasets: NCDatasets, NCDataset
export NCDatasets, NCDataset
export nckeys, ncdetails, ncglobalattr, ncsize
export ncread, ncwrite

dim_to_commonname(::Lat) = "lat"
dim_to_commonname(::Lon) = "lon"
dim_to_commonname(::Time) = "time"
dim_to_commonname(::Pre) = "level"
dim_to_commonname(D::Dim) = string(DimensionalData.name(D))
dim_to_commonname(::Coord) = "cell"

const POSSIBLE_CELL_NAMES = ("ncells", "cell", "rgrid", "grid")


"""
    nckeys(file::String)
Return all keys of the `.nc` file in `file`.
"""
function nckeys(path::String)
    NCDataset(path) do ds
        return keys(ds)
    end
end
nckeys(a::NCDataset) = keys(a)

"""
    ncdetails(file, io = stdout)
Print details about the `.nc` file in `file` on `io`.
"""
function ncdetails(file, io = stdout)
    NCDataset(file) do ds
        show(io, MIME"text/plain"(), ds)
    end
end
ncdetails(ds::NCDataset, io = stdout) = show(io, MIME"text/plain"(), ds)

"""
    ncsize(file, var)
Return the size of the variable of the `.nc` file without actually loading any data.
"""
function ncsize(file, var)
    NCDataset(file) do ds
        return size(ds[var])
    end
end


"""
    ncglobalattr(file::String) → Dict
Return the global attributes of the .nc file.
"""
function ncglobalattr(file::String)
    NCDataset(file) do ds
        return Dict(ds.attrib)
    end
end

#########################################################################
# Imports
#########################################################################
include("netcdf_read.jl")
include("netcdf_write.jl")