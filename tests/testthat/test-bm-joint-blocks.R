
test_that("a block never contains derived nodes, even if all fixed effects are shared", {
  d <- bm_nb_data()
  both <- c("biomass", "richness")
  fx <- list(biomass = c("x1", "x2"), richness = c("x1", "x2"))
  fam <- c(biomass = "truncnorm", richness = "negbinomial")
  sh <- bm_resolve_sharing(both, fx, c("year", "replication"), "none", "none", TRUE, FALSE,
                           FALSE, FALSE, FALSE, TRUE)
  spec <- list(responses = both, family = fam, fixed_effects = fx,
               random_effects = c("year", "replication"), temporal_effect = "none",
               spatial_effect = "none", sharing = sh)
  prep <- bm_prepare_data(d, both, fx, spec$random_effects, fam, bm_default_columns())
  b <- bm_build_model(prep, spec, bm_default_priors(prep, fam))
  expect_false(any(grepl("\\[\\]", unlist(b$blocks))))
  # richness: every fixed effect is derived, but its year effects keep a block
  expect_true(any(vapply(b$blocks, function(x) any(grepl("^u_year_richness", x)), logical(1))))
  expect_false(any(grepl("^beta_richness", unlist(b$blocks))))
})
