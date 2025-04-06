expect_cli_warning <- function(object, index, result){
  warns <- capture_messages(object)
  warn <- cli::ansi_strip(warns[index])
  warn <- sub("\n$", "", warn)
  expect_equal(warn, result)
}
