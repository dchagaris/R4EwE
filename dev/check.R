
# Full R CMD check
local({
  message("→ R CMD check (CRAN-like)...")
  devtools::check(cran = TRUE)
  message("✓ Check completed.")
})



