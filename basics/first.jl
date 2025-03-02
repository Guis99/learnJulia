struct Location
    name::String
    lat::Float32
    lon::Float32
end

loc1 = Location("Los Angeles", 34.0522,-118.2437)

loc1.name   # "Los Angeles"
loc1.lat    # 34.0522
loc1.lon    # -118.2437

sites = Location[]
push!(sites, Location("Los Angeles", 34.0522,-118.2437))
push!(sites, Location("Las Vegas", 36.1699,-115.1398))

display(sites)