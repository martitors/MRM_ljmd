#include <gtest/gtest.h>
#include <stdio.h>
#include <stdlib.h>
#include "types.h" // Assuming mdsys_t and other necessary types are declared here
#include "input.h"

// Test case to check if read function correctly parses input data
TEST(InputReaderTest, MinimalInputFileTest){
    char line[200], restfile[200], trajfile[200], ergfile[200];
    mdsys_t sys;
    int nprint;

    const char* inputFilePath = "../tests/minimal_input.txt"; // Replace with the path to your minimal input file

    // Redirect stdin to read from the minimal input file
    FILE* original_stdin = freopen(inputFilePath, "r", stdin);

    // Call the read function with input from the file
    read_input(line, restfile, trajfile, ergfile, &sys, &nprint);

    // Assertions to check if the read function parsed the input correctly
    ASSERT_EQ(sys.natoms, 4);
    ASSERT_DOUBLE_EQ(sys.mass, 39.948);
    // Assert other values from the input file

    // Clean up: Restore original stdin
    fclose(stdin);
    stdin = original_stdin;
}
