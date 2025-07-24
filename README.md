# 3dtrees-pc-standard
las/laz file checks and standardization


1. Basic LAS Validation

    When a LAS file is loaded, it shall validate the LAS object using las_check().

    If any errors are returned by las_check(), it shall store them in the error list.

    If the LAS file has fewer than 100 points, it shall raise an error and skip further processing.

2. Bounding Box Consistency

    If a warning indicates that the header's bounding box does not match the points, it shall update the LAS object and store a warning.

    When checking spatial extent, it shall calculate the extend in X, Y, and Z dimensions.

    If any dimension is ≤ 1 meter, it shall raise an error indicating insufficient spatial extent.

    If the Z dimension range is < -100 or > 5500 meters, it shall raise an error indicating an implausible Z range.

3. Point Density Check

    When calculating point density, it shall compute the convex hull of X and Y coordinates.

    If the point density is ≤ 1000 points/m², it shall raise an error and recommend checking for ALS data or multiple plots.

4. Extra Byte Data Type Check

    When checking extra byte attributes, it shall map each attribute’s data type to the expected R class.

    If any mismatch is found, it shall raise a warning and attempt to correct the attribute definition.

5. CRS Validation

    If the LAS file lacks a CRS, it shall raise an error.

    If the CRS exists, it shall attempt to transform the file to EPSG:3857.

    If the transformation fails, it shall raise an error.

6. Point Data Format Check

    If the LAS Point Data Format ID > 7, it shall:

        downgrade the format to 7,

        rename the NIR channel to "NIRe",

        raise a warning about potential loss of full waveform data.

7. Land Coverage Validation

    If the LAS file has a valid and transformable CRS, it shall:

        calculate the centroid of the bounding box,

        check whether the centroid lies within a known land area using a shapefile.

    If the centroid is not on land, it shall raise an error.

8. Output Writing

    When all checks are complete, it shall:

        write warnings and errors to a JSON file.

        If no errors are present, it shall write the LAZ file to the output folder.

        If any errors are present, it shall not write the LAZ file to the output folder.