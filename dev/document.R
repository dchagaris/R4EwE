
# Run roxygen to update NAMESPACE and man/
local({
  message("→ Documenting...")
  usethis::use_roxygen_md()   # safe to call; sets roxygen to markdown if not already
  devtools::document()
  message("✓ Documentation updated.")
})

