
# build package -----------------------------------------------------------

usethis::use_mit_license(copyright_holder = "Akira Terui")
usethis::use_roxygen_md()
devtools::document()
devtools::load_all()
devtools::check(vignettes=FALSE)


# check syntax ------------------------------------------------------------

lintr::lint_package()

#usethis::use_coverage()
# covr::package_coverage()


# test run ----------------------------------------------------------------

df0 <- data.frame(t = rep(1:10, 3),
                  species = rep(1:3, each = 10),
                  n = runif(30))

df1 <- lag_block(df0,
                 index = "species") |>
  dplyr::mutate(log_r = log(n) - log(n_lag)) |>
  na.omit()

loocv(log_r ~ n_lag + nt_lag + (1|species),
      data = df1,
      theta = 0.1,
      model = "glmmTMB")

