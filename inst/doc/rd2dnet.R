## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7, fig.height = 4)

## ----setup--------------------------------------------------------------------
library(rd2dnet)
city <- load_mock_midtown()
obs <- city$obs
y <- obs$outcome
head(obs)

## ----map, fig.height = 6------------------------------------------------------
nd <- city$nodes
e <- city$edges
plot(nd$x, nd$y, type = "n", asp = 1, xlab = "x (m)", ylab = "y (m)",
     main = "mock city: treated (red) / control (blue)")
segments(nd$x[match(e$from, nd$name)], nd$y[match(e$from, nd$name)],
         nd$x[match(e$to, nd$name)], nd$y[match(e$to, nd$name)], col = "grey75", lwd = 0.6)
polygon(city$park, col = adjustcolor("forestgreen", 0.3), border = NA)
points(obs$x, obs$y, pch = 16, cex = 0.25, col = ifelse(obs$treated, "firebrick", "steelblue"))
lines(city$boundary$line, lwd = 2)
city$boundary

## ----metrics------------------------------------------------------------------
metrics <- list(
  "length (m)" = "length",
  "travel time (s)" = travel_time(),
  # 50 extra units per unit change in floor-area ratio between adjacent nodes
  "length + land-use contrast" = node_contrast("far", 50),
  # switching road class costs 500 m; ramps are transparent connectors, so
  # highway -> ramp -> highway is cheap but highway -> ramp -> street is not
  "length + class switches" = transition_cost("length", class_attr = "highway",
                                              default_penalty = 500, connectors = "motorway_link")
)
pts <- obs[1:5, c("x", "y")]
round(sapply(metrics, function(m) {
  network_coordinates(city$network(m), city$boundary, pts, n_eval = 5)$s
}), 1)

## ----precomputed--------------------------------------------------------------
base <- city$network("length")
D <- distances_from(base, seq_len(base$n_nodes)) # or your own matrix, in node order
pre <- city$network(precomputed(sqrt(D))) # e.g. a concave transform
round(network_coordinates(pre, city$boundary, pts, n_eval = 5)$s, 2)

## ----coords-------------------------------------------------------------------
net <- city$network("length")
coords <- network_coordinates(net, city$boundary, obs, n_eval = 21)
coords
head(as.data.frame(coords, j = 11))
head(eval_table(coords), 8)

## ----bandwidths---------------------------------------------------------------
bw <- rdbwnet(y, coords)
head(bw[, c("arc", "h0", "h1", "n0", "n1", "status")])

## ----estimate-----------------------------------------------------------------
res <- rdnet(y, coords, seed = 0)
res

## ----plot---------------------------------------------------------------------
plot(res, true_tau = city$true_tau(res$table$arc),
     main = "BATEC along the boundary: estimate, robust CIs, uniform band, truth")

## ----fixed--------------------------------------------------------------------
res_fixed <- rdnet(y, coords, h = 600, weights = obs$far, seed = 0)
head(res_fixed$table[, c("arc", "tau_bc", "se_rb", "status")])

## ----aggregate----------------------------------------------------------------
ok <- res$table$status == "ok"
summary(res, wbate = as.numeric(ok), lbate = TRUE)$aggregate
# Weight by frontage instead: here, by the number of parcels near each point.
wbate(res, ifelse(ok, res$table$n0 + res$table$n1, 0))

## ----correlated---------------------------------------------------------------
se_naive <- sqrt(sum((ok / sum(ok))^2 * ifelse(ok, res$table$se_rb, 0)^2))
c(with_covariance = wbate(res, as.numeric(ok))[["se_rb"]], ignoring_it = se_naive)

## ----covariates---------------------------------------------------------------
res_cov <- rdnet(y, coords, h = 600, covs = log(obs$far), seed = 0)
res_raw <- rdnet(y, coords, h = 600, seed = 0)
cmp <- data.frame(tau_bc = res_raw$table$tau_bc, tau_bc_covs = res_cov$table$tau_bc,
                  se_rb = res_raw$table$se_rb, se_rb_covs = res_cov$table$se_rb)
head(round(na.omit(cmp), 3))

## ----fuzzy--------------------------------------------------------------------
big <- make_mock_midtown(n_obs = 12000, seed = 11)
big_coords <- network_coordinates(big$network(), big$boundary, big$obs, n_eval = 21)
set.seed(2)
T <- big_coords$treated
D_takeup <- as.numeric(ifelse(T, runif(length(T)) < 0.75, runif(length(T)) < 0.25))
y_fuzzy <- 0.0008 * (big$obs$y - 1600) + 2 * D_takeup + rnorm(length(T), 0, 0.5)
res_fz <- rdnet(y_fuzzy, big_coords, fuzzy = D_takeup, seed = 0)
head(round(na.omit(res_fz$table[, c("arc", "fs", "fs_se", "tau", "ci_lo", "ci_hi")]), 3))

## ----kink---------------------------------------------------------------------
kink <- rdnet(y, coords, h = 800, deriv = c(1, 0), seed = 0)
head(signif(na.omit(kink$table[, c("arc", "tau_bc", "ci_lo", "ci_hi")]), 3))

## ----vce----------------------------------------------------------------------
fits <- list(
  HC1 = rdnet(y, coords, h = 600, seed = 0),
  cluster_block = rdnet(y, coords, h = 600, cluster = obs$block, seed = 0),
  conley_network = rdnet(y, coords, h = 600, conley_cutoff = 150, seed = 0),
  conley_euclidean = rdnet(y, coords, h = 600, conley_cutoff = 150, conley_distance = "euclidean", seed = 0)
)
head(round(na.omit(sapply(fits, function(r) r$table$se_rb)), 3))

## ----density------------------------------------------------------------------
density(coords)
density(coords, normalize = "none")

## ----donut--------------------------------------------------------------------
y_donut <- ifelse(abs(coords$s) < 40, NA, y)
donut <- rdnet(y_donut, coords, h = 700, seed = 0)
head(round(na.omit(data.frame(full = res_fixed$table$tau_bc, donut_40m = donut$table$tau_bc)), 3))

## ----balance------------------------------------------------------------------
bal <- rdnet(log(obs$far), coords, seed = 0)
summary(bal, wbate = as.numeric(bal$table$status == "ok"), lbate = TRUE)$aggregate

## ----placebo------------------------------------------------------------------
fake_y <- city$boundary_y + 8 * 80
north <- obs$y > city$boundary_y
fake <- boundary(rbind(c(-400, fake_y), c(2600, fake_y)),
                 treated_region = rect_polygon(-1e5, -1e5, 1e5, fake_y), mask = list(city$park))
pc <- network_coordinates(city$network(), fake, obs[north, ], n_eval = 15)
placebo <- rdnet(y[north], pc, seed = 0)
summary(placebo, wbate = as.numeric(placebo$table$status == "ok"))$aggregate

## ----did----------------------------------------------------------------------
set.seed(5)
f <- y - city$true_tau(coords$crossing_arc) * coords$treated # outcome without the effect
y_pre <- f + rnorm(length(y), 0, 0.5) # same parcels, earlier period
y_post <- y + rnorm(length(y), 0, 0.5)
did_panel <- rdnet(y_post - y_pre, coords, h = 700, seed = 0)

c_pre <- make_mock_midtown(n_obs = 3000, seed = 21)
c_post <- make_mock_midtown(n_obs = 3000, seed = 22)
co_pre <- network_coordinates(c_pre$network(), c_pre$boundary, c_pre$obs, n_eval = 21)
co_post <- network_coordinates(c_post$network(), c_post$boundary, c_post$obs, n_eval = 21)
y0 <- c_pre$obs$outcome - c_pre$true_tau(co_pre$crossing_arc) * co_pre$treated
r_pre <- rdnet(y0, co_pre, h = 700, seed = 0)$table
r_post <- rdnet(c_post$obs$outcome, co_post, h = 700, seed = 0)$table
did <- data.frame(arc = r_post$arc, panel = did_panel$table$tau_bc,
                  independent = r_post$tau_bc - r_pre$tau_bc,
                  se_independent = sqrt(r_pre$se_rb^2 + r_post$se_rb^2),
                  truth = city$true_tau(r_post$arc))
head(round(na.omit(did), 3))

## ----panel--------------------------------------------------------------------
panel_xy <- as.matrix(obs[rep(1:500, each = 20), c("x", "y")]) # 500 sensors x 20 readings
pcoords <- network_coordinates(net, city$boundary, panel_xy, n_eval = 11)
c(n_obs = pcoords$n_obs, n_locations = pcoords$n_locations)

