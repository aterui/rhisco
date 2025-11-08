
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

# a <- inla_lmer(log_r ~ n_lag + nt_lag,
#                data = df1,
#                control.compute = list(config = TRUE,
#                                       cpo = TRUE))

sapply(5:10*0.1, function(x) {
  loocv(log_r ~ n_lag + nt_lag + (1|species),
        data = df1,
        theta = x,
        group = "species",
        model = "lm",
        method = "max")
}) |>
  plot()

x_star <- xeq(log_r ~ nt_lag,
              df1 %>% filter(species == 1),
              theta = c(0.5, 1, 2, 4, 8))

get_psi(log_r ~ n_lag + nt_lag + (1 | species),
        data = df1,
        x_star = x_star,
        theta = c(0, 0.5, 1, 2, 4),
        model = "lm")
