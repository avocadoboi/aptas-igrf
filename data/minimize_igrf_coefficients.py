# Takes the IGRF coefficient text file downloaded from here: https://www.ncei.noaa.gov/products/international-geomagnetic-reference-field
# and extracts only the two columns that we need

HEADER_LINE_COUNT = 4

with open("igrf14coeffs.txt") as file_in, open("IGRF14_minified", "w") as file_out:
    # Skip the header
    for _ in range(HEADER_LINE_COUNT):
        next(file_in)

    for line in file_in:
        *_, coeffficient_2025, secular_variation = line.split()
        # Write only the 2025 coefficients and the secular variation 2025-2030.
        # The g/h coefficient type and n and m indices are ordered as in a nested sum with n going from 1 to 13 and m going from 0 to n 
        # and the g and h coefficients are alternating with g before h for each pair of indices. The h coefficients with m = 0 are skipped
        # because they would be multiplying sin(0) = 0 in the sum. (see Equation (1) in https://doi.org/10.1186/s40623-020-01288-x)
        print(coeffficient_2025, secular_variation, file=file_out)

