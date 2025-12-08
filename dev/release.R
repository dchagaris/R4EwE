
# Prep a GitHub release: version bump, NEWS, tag
local({
  # Update package version
  usethis::use_version(which = "minor")   # or "patch" / "major"
  # Ensure NEWS.md exists and is updated
  if (!file.exists("NEWS.md")) usethis::use_news_md()
  # Document and check before release
  devtools::document()
  devtools::check(cran = TRUE)
  # Create and push Git tag
  usethis::use_github_release()           # interactive; or git tag/push manually
})