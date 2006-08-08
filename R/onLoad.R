".onLoad" <- function(lib, pkg) {
  mylib <- dirname(system.file(package = "experiment"))
  title <- packageDescription("experiment", lib = mylib)$Title
  ver <- packageDescription("experiment", lib = mylib)$Version
  url <- packageDescription("experiment", lib = mylib)$URL
  cat(title, "\nVersion:", ver, "\nURL:", url, "\n")
}

