

## Test if type is passed correctly in brquasiFitControl
types <- c("M", "iRBM", "eRBM", "MPQL_trace")
for (type in types) {
    brquasi_control <- brquasiFitControl(type = type)
    expect_equal(brquasi_control$type, type)
}
expect_error(brquasiFitControl(type = "magic_bias_reduction_method"), pattern = "'arg' should be one of")


## Test that error is produced if epsilon is <= 0
expect_error(brquasiFitControl(epsilon = -1), pattern = "'epsilon' must be > 0")
expect_error(brquasiFitControl(epsilon = 0), pattern = "'epsilon' must be > 0")
expect_error(brquasiFitControl(lambda = -1), pattern = "'lambda' must be >= 0")

## Test disp_factor
expect_error(brquasiFitControl(disp_factor = "n-p-1000"), pattern = "'arg' should be one of")
expect_identical(brquasiFitControl(disp_factor = "n-p")$disp_factor, "n-p")
expect_identical(brquasiFitControl(disp_factor = "n")$disp_factor, "n")
