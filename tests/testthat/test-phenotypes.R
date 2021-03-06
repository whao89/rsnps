# tests for phenotypes fxn in ropensnp
context("phenotypes")

sm <- suppressMessages
df <- sm(phenotypes(userid='1,6,8', df=TRUE))

test_that("phenotypes returns the correct class", {
	expect_that(df, is_a("list"))
	expect_that(df[[1]], is_a("data.frame"))
})

test_that("phenotypes returns the correct value", {
	expect_that(as.character(sm(phenotypes(userid=1))$user$name), equals("Bastian Greshake"))
})
