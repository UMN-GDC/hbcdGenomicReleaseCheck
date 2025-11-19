library(testthat)
#library(xml2)

#junit_reporter = JunitReporter$new(file= "test-results.xml")
#
#with_reporter(
#  code = test_dir("tests/testthat"),
#  reporter = JunitReporter$new(file= "test-results.xml")
#)

sink("test-log.txt")

# 2. Run your tests (the results will be written to the file, not the console)
test_dir("tests/testthat", reporter = "summary")

# 3. Stop capturing (CRITICAL: closes the file and ensures it's fully written)
sink()
