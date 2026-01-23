import netCDF4
import os
from compliance_checker.runner import ComplianceChecker

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

def run_checker(filename):
    # Run CF-1.6 check
    res, errors = ComplianceChecker.run_checker(
        filename,
        ["cf:1.6"],
        verbose=1,
        criteria="normal",
        output_filename="report.txt",
        output_format="text"
    )
    
    print(f"Return code: {res}")
    # Read report to see errors
    with open("report.txt", "r") as f:
        print(f.read())

if __name__ == "__main__":
    nc_file = "test_issue_1245.nc"
    create_nc_file(nc_file)
    run_checker(nc_file)
    # cleanup
    if os.path.exists(nc_file):
        os.remove(nc_file)
    if os.path.exists("report.txt"):
        os.remove("report.txt")
