import netCDF4
import os

def create_nc_file(filename):
    if os.path.exists(filename):
        os.remove(filename)
    
    ds = netCDF4.Dataset(filename, "w", format="NETCDF4")
    
    # Create time dimension
    ds.createDimension("time", 10)
    
    # Create time variable
    time_var = ds.createVariable("time", "f8", ("time",))
    time_var.units = "seconds since 1970-01-01 00:00:00"
    time_var.axis = "T"  # Providing axis='T'
    # time_var.standard_name = "time"  # Intentionally omitted
    
    time_var[:] = range(10)
    
    ds.close()
    print(f"Created {filename}")

if __name__ == "__main__":
    nc_file = "test_issue_1245.nc"
    create_nc_file(nc_file)
