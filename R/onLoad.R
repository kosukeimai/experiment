".onLoad" <- function(lib, pkg) {
  mylib <- dirname(system.file(package = "are"))
  title <- packageDescription("are", lib = mylib)$Title
  ver <- packageDescription("are", lib = mylib)$Version
  url <- packageDescription("are", lib = mylib)$URL
  cat(title, "\nVersion:", ver, "\nURL:", url, "\n")
}

