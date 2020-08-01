
test_that("str_length", {
  expect_equal(str_length("a"),1)
})

check_api <- function() {
  #if (not_working()) {
    skip("API not available")
  #}
}

test_that("foo api returns bar when given baz", {
  check_api()
  ...
})
