using DrWatson
using Pkg
Pkg.activate(desktop("GeoMakieTests", "Project.toml"))
using ClimateBase
using GLMakie, GeoMakie

Pkg.activate(raw"C:\Users\m300808\ownCloud\Projects\DYAMONDcomparison\Project.toml")
Ru = ncread(datadir("NETCDF_FILES", "CERES_EBAF_gea250.nc"), "toa_sw_all_mon")
Rs = ncread(datadir("NETCDF_FILES", "ngc2001_atm_2d_ml_dailymean.nc"), "rsdt")
A = timemean(Rs)
B = timemean(Ru)


# %%


# test it
fig, ax, el, cb = climplot(B)
display(fig)

# %%
obs =  Observable(Ru[Time(1)])
fig, ax, el, cb = climplot(obs)
display(fig)

for i in 1:12
    obs[] = Ru[Time(i)]
    sleep(0.5)
end

MakieLayout.deactivate_interaction!(ax, :rectanglezoom)
point = select_point(ax.scene)
on(point) do val
    println(val)
end

# %%

fig, ax, el, cb = earthsurface(A)
display(fig)


# %% MWE
# using GLMakie, GeoMakie
