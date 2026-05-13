# CF Feature Parity Matrix

This document is intended to show where the IOOS Compliance Checker supports CF features as outlined by the conformance lists and conventions documents.
Most CF features are supported, but several are either partially implemented or are not presently implemented, as denoted by the table below.

| Section | Implemented | Notes | GitHub Tracking Issue |
|---------|-------------|-------|-----------------------|
| 2.1. Filename | yes |
| 2.2. Data Types | yes |
| 2.3. Naming Conventions | yes |
| 2.4. Dimensions | yes |
| 2.5. Variables | yes |
| 2.5.1. Missing data, valid and actual range of data |yes |
| 2.6. Attributes | yes |
| 2.6.1. Identification of Conventions | yes |
| 2.6.2. Description of file contents | yes |
| 2.6.3. External Variables | yes |
| 2.7. Groups | no | Handling code exists, not implemented on checker level.  Requires significant code changes to code that references other variables or attributes. | https://github.com/ioos/compliance-checker/issues/945 |
| 2.7.1. Scope | " |
| 2.7.2. Application of attributes | " |
| 3. Description of the Data | yes |
| 3.1. Units | yes |
| 3.2. Long Name  | yes |
| 3.3. Standard Name | yes |
| 3.4. Ancillary Data | yes |
| 3.5. Flags | yes |
| 4. Coordinate Types | yes |
| 4.1. Latitude Coordinate | yes |
| 4.2. Longitude Coordinate | yes |
| 4.3. Vertical (Height or Depth) Coordinate | yes |
| 4.3.1. Dimensional Vertical Coordinate | yes |
| 4.3.2. Dimensionless Vertical Coordinate | yes |
| 4.3.3. Parametric Vertical Coordinate | yes |
| 4.4. Time Coordinate | yes |
| 4.4.1. Calendar | yes |
| 4.5. Discrete Axis | yes |
| 5. Coordinate Systems and Domain | partial | Needs ragged and indexed array support | https://github.com/ioos/compliance-checker/issues/949
| 5.1. Independent Latitude, Longitude, Vertical, and Time Axes | yes
| 5.2. Two-Dimensional Latitude, Longitude, Coordinate Variables | yes
| 5.3. Reduced Horizontal Grid | yes
| 5.4. Timeseries of Station Data | yes
| 5.5. Trajectories | yes
| 5.6. Horizontal Coordinate Reference Systems, Grid Mappings, and Projections | partial | Needs implementation of a few more sections | https://github.com/ioos/compliance-checker/issues/950
| 5.6.1. Use of the CRS Well-known Text Format | yes |
| 5.7. Scalar Coordinate Variables | yes |
| 5.8. Domain Variables | yes |
| 6. Labels and Alternative Coordinates | yes
| 6.1. Labels | yes
| 6.1.1. Geographic Regions | partial | Implicit from handling of standard names - does not check region names
| 6.1.2. Taxon Names and Identifiers | yes
| 6.2. Alternative Coordinates | no? |
| 7. Data Representative of Cells | yes |
| 7.1. Cell Boundaries | yes |
| 7.2. Cell Measures | yes |
| 7.3. Cell Methods | yes |
| 7.3.1. Statistics for more than one axis | yes |
| 7.3.2. Recording the spacing of the original data and other information | yes |
| 7.3.3. Statistics applying to portions of cells | yes |
| 7.3.4. Cell methods when there are no coordinates | ? |
| 7.4. Climatological Statistics | yes |
| 7.5. Geometries | partial | Various portions of the geometry section need to be implemented to bring in line with conformance doc. See issue. | https://github.com/ioos/compliance-checker/issues/955
| 8. Reduction of Dataset Size | yes |
| 8.1. Packed Data | partial |
| 8.2. Lossless Compression by Gathering | yes |
| 8.3. Lossy Compression by Coordinate Subsampling | no | Section requires considerable work to implement and not many files have been found with this interpolation. Presently not implemented even in development branches. | https://github.com/ioos/compliance-checker/issues/957
| 8.3.1. Tie Points and Interpolation Subareas | " |
| 8.3.2. Coordinate Interpolation Attribute | " |
| 8.3.3. Interpolation Variable | " |
| 8.3.4. Subsampled, Interpolated and Non-Interpolated Dimensions | " |
| 8.3.5. Tie Point Mapping Attribute | " |
| 8.3.6. Tie Point Dimension Mapping | " |
| 8.3.7. Tie Point Index Mapping | " |
| 8.3.8. Interpolation Parameters | " |
| 8.3.9. Interpolation of Cell Boundaries | " |
| 8.3.10. Interpolation Method Implementation | " |
| 9. Discrete Sampling Geometries | yes |
| 9.1. Features and feature types | yes |
| 9.2. Collections, instances and elements | yes |
| 9.3. Representations of collections of features in data variables | yes |
| 9.3.1. Orthogonal multidimensional array representation | yes |
| 9.3.2.  Incomplete multidimensional array representation | yes |
| 9.3.3.  Contiguous ragged array representation | yes |
| 9.3.4. Indexed ragged array representation | yes |
| 9.4. The featureType  attribute | yes |
| 9.5. Coordinates and metadata | yes |
| 9.6. Missing Data | ? |
