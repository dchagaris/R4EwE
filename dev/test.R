

# Run unit tests with testthat
local({
  message("→ Running tests...")
  devtools::test()            # or testthat::test_dir("tests/testthat")
  message("✓ Tests completed.")
})


