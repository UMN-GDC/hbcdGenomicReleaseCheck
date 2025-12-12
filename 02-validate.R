library(testthat)
#library(xml2)

# junit_reporter = JunitReporter$new(file= "test-results.xml")
# 
# with_reporter(
#   code = test_dir("tests/testthat"),
#   reporter = JunitReporter$new(file= "test-results.xml")
# )

sink("test-log.txt")
test_dir("tests/testthat")

sink()
