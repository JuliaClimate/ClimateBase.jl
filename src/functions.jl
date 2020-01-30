"""
    periodmean(C::ClimGrid; startdate::Tuple, enddate::Tuple)

Mean of array data over a given period.
"""
function periodmean(C::ClimGrid; start_date::Tuple=(Inf, ), end_date::Tuple=(Inf,), level=1)

    if start_date != (Inf,) || end_date != (Inf,)
        C = temporalsubset(C, start_date, end_date)
    end

    datain = C.data.data

    # Mean and squeeze
    if ndims(datain) == 2
        dataout = datain
    elseif ndims(datain) == 3
        if size(datain, 3) == 1 # already an average on single value
            dataout = dropdims(datain, dims=3)
        else
            dataout = dropdims(Images.meanfinite(datain, 3), dims=3)
        end
    elseif ndims(datain) == 4
        if size(datain, 4) == 1
            dataout = dropdims(datain[:, :, level, :], dims = 3)
        else
            dataout = dropdims(Images.meanfinite(datain, 3), dims=3)
        end
    end


    # Build output AxisArray
    FD = buildarray_climato(C, dataout)

    # Return climGrid type containing the indice
    return ClimGrid(FD, longrid=C.longrid, latgrid=C.latgrid, msk=C.msk, grid_mapping=C.grid_mapping, dimension_dict=C.dimension_dict, timeattrib=C.timeattrib, model=C.model, frequency=C.frequency, experiment=C.experiment, run=C.run, project=C.project, institute=C.institute, filename=C.filename, dataunits=C.dataunits, latunits=C.latunits, lonunits=C.lonunits, variable="periodmean", typeofvar=C.typeofvar, typeofcal="climatology", varattribs=C.varattribs, globalattribs=C.globalattribs)
end

function buildarray_climato(C::ClimGrid, dataout)
    lonsymbol = Symbol(C.dimension_dict["lon"])
    latsymbol = Symbol(C.dimension_dict["lat"])
    FD = AxisArray(dataout, Axis{lonsymbol}(C[1][Axis{lonsymbol}].val), Axis{latsymbol}(C[1][Axis{latsymbol}].val))
    return FD
end

function buildarrayinterface(axisArraytmp, A)
    latsymbol = Symbol(A.dimension_dict["lat"])
    lonsymbol = Symbol(A.dimension_dict["lon"])
    if ndims(axisArraytmp) == 2
        axisArray = AxisArray(axisArraytmp, Axis{lonsymbol}(A[1][Axis{lonsymbol}].val), Axis{latsymbol}(A[1][Axis{latsymbol}].val))
    elseif ndims(axisArraytmp) == 3
        axisArray = AxisArray(axisArraytmp, Axis{lonsymbol}(A[1][Axis{lonsymbol}].val), Axis{latsymbol}(A[1][Axis{latsymbol}].val), Axis{:time}(A[1][Axis{:time}].val))
    elseif ndims(axisArraytmp) == 4
        axisArray = AxisArray(axisArraytmp, Axis{lonsymbol}(A[1][Axis{lonsymbol}].val), Axis{latsymbol}(A[1][Axis{latsymbol}].val), Axis{:plev}(A[1][Axis{:plev}].val), Axis{:time}(A[1][Axis{:time}].val))
    end
    return axisArray
end

function buildarray_annual(C::ClimGrid, dataout, numYears)
    lonsymbol = Symbol(C.dimension_dict["lon"])
    latsymbol = Symbol(C.dimension_dict["lat"])
    FD = AxisArray(dataout, Axis{lonsymbol}(C[1][Axis{lonsymbol}].val), Axis{latsymbol}(C[1][Axis{latsymbol}].val), Axis{:time}(Dates.year.(DateTime.(numYears))))
    return FD
end



function buildarray_resample(C::ClimGrid, dataout, newtime)
    lonsymbol = Symbol(C.dimension_dict["lon"])
    latsymbol = Symbol(C.dimension_dict["lat"])
    FD = AxisArray(dataout, Axis{lonsymbol}(C[1][Axis{lonsymbol}].val), Axis{latsymbol}(C[1][Axis{latsymbol}].val), Axis{:time}(newtime))
    return FD
end


"""
    function temporalsubset(C::ClimGrid, startdate::Date, enddate::Date)

Returns the temporal subset of ClimGrid C. The temporal subset is defined by a start and end date.

"""
function temporalsubset(C::ClimGrid, datebeg::Tuple, dateend::Tuple)

    T = typeof(get_timevec(C)[1])
    timeV = get_timevec(C)
    idxtimebeg, idxtimeend = timeindex(timeV, datebeg, dateend, T)

    startdate = buildtimetype(datebeg, T)
    enddate = buildtimetype(dateend, T)

    # some checkups
    @argcheck startdate <= enddate

    dataOut = C[1][Axis{:time}(idxtimebeg:idxtimeend)]

    # The following control ensure that a 1-timestep temporal subset returns a 3D Array with time information on the timestep. i.e. startdate == enddate
    if ndims(dataOut) == 2
        timeV = startdate
        latsymbol = Symbol(C.dimension_dict["lat"])
        lonsymbol = Symbol(C.dimension_dict["lon"])
        data2 = fill(NaN, (size(dataOut,1), size(dataOut, 2), 1))
        data2[:,:,1] = dataOut
        dataOut = AxisArray(data2, Axis{lonsymbol}(C[1][Axis{lonsymbol}].val), Axis{latsymbol}(C[1][Axis{latsymbol}].val), Axis{:time}(C[1][Axis{:time}].val))

    end

    return ClimGrid(dataOut, longrid=C.longrid, latgrid=C.latgrid, msk=C.msk, grid_mapping=C.grid_mapping, dimension_dict=C.dimension_dict, timeattrib=C.timeattrib, model=C.model, frequency=C.frequency, experiment=C.experiment, run=C.run, project=C.project, institute=C.institute, filename=C.filename, dataunits=C.dataunits, latunits=C.latunits, lonunits=C.lonunits, variable=C.variable, typeofvar=C.typeofvar, typeofcal=C.typeofcal, varattribs=C.varattribs, globalattribs=C.globalattribs)

end

"""
    get_timevec(C::ClimGrid)

Returns time vector of ClimGrid C.
"""
get_timevec(C::ClimGrid) = C[1][Axis{:time}][:]

"""
    buildtimetype(datetuple, f)

Returns the adequate DateTime for temporal subsetting using DateType *f*
"""
function buildtimetype(date_tuple, f)

    if length(date_tuple) == 1
        dateout = f(date_tuple[1], 01, 01)
    elseif length(date_tuple) == 2
        dateout = f(date_tuple[1], date_tuple[2], 01)
    elseif length(date_tuple) == 3
        dateout = f(date_tuple[1], date_tuple[2], date_tuple[3])
    elseif length(date_tuple) == 4
        dateout = f(date_tuple[1], date_tuple[2], date_tuple[3], date_tuple[4], 00, 00)
    elseif length(date_tuple) == 5
        dateout = f(date_tuple[1], date_tuple[2], date_tuple[3], date_tuple[4], date_tuple[5], 00)
    elseif length(date_tuple) == 6
        dateout = f(date_tuple[1], date_tuple[2], date_tuple[3], date_tuple[4], date_tuple[5], date_tuple[6])
    end

    return dateout
end

"""
    timeindex(timeVec, start_date, end_date, T)

Return the index of time vector specified by start_date and end_date. T is the DateTime type (see NCDatasets.jl documentation).
"""
function timeindex(timeV, datebeg::Tuple, dateend::Tuple, T)

    # Start Date
    if !isinf(datebeg[1])
        # Build DateTime type
        start_date = buildtimetype(datebeg, T)
        # @argcheck start_date >= timeV[1]
        idxtimebeg = findfirst(timeV .>= start_date)[1]
    else
        idxtimebeg = 1
    end
    # End date
    if !isinf(dateend[1])
        # Build DateTime type
        end_date = buildtimetype(dateend, T)
        # @argcheck end_date <= timeV[end]
        idxtimeend = findlast(timeV .<= end_date)[1]
    else
        idxtimeend = length(timeV)
    end

    if !isinf(datebeg[1]) && !isinf(dateend[1])
        @argcheck start_date <= end_date
    end
    return idxtimebeg, idxtimeend
end
