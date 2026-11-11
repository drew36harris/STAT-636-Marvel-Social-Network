# Libraries
library(igraph)
library(Matrix)
library(dplyr)
library(stringr)
library(ggplot2)
library(tibble)
library(tidyr)

# Added: Unweighted projections, Bipartite centralities, Projection centralities, and a quick comparison

# Read and Parse
path <-  "MarvelSocialNetwork.txt"
lines <- readLines(path, encoding="latin1") # File is Pajek-style. Vertices first, then Edges

### Get Counts
hdr <- as.integer(str_extract_all(lines[1], "\\d+")[[1]]) #Gets numbers out of the header line
n_total <- hdr[1] # Total Rows
n_char <- hdr[2] # Number of character rows
n_comic <- n_total - n_char # Number of comic rows

### Find start of edgelist
edge_start <- which(str_detect(lines, "^\\*Edgeslist"))

### Vertex labels
vlines <- lines[2:(edge_start - 1)]
vertices <- tibble(
  id = as.integer(str_extract(vlines, "^\\d+")),
  label = str_match(vlines, '^\\d+ "(.*)"$')[, 2]
) %>% 
  mutate(type = id > n_char, # Character = FALSE, Comic = TRUE
         kind = ifelse(type, "comic", "character"))

### Edge table
elines <- lines[(edge_start + 1):length(lines)]
elines <- elines[nzchar(trimws(elines))] # Drops blank lines
tok <- strsplit(trimws(elines), "\\s+") # split lines at white space
edges <- bind_rows(lapply(tok, function(t) {
  tibble(character = as.integer(t[1]), comic = as.integer(t[-1])) # one line = one character followed by its comics -> one row per (character, comic) pair
}))

# Validate
stopifnot(nrow(vertices) == n_total)

cat("Raw edges:", nrow(edges), "\n")
cat("Duplicate character-comic pairs:", sum(duplicated(edges)), "\n")
edges <- distinct(edges)

cat("Edges with character ID outside character range:", sum(edges$character > n_char), "\n")
cat("Edges with comic ID inside character range:", sum(edges$comic <= n_char), "\n")
cat("Unique characters in edges:", n_distinct(edges$character), "of", n_char, "\n")
cat("Unique comics in edges:", n_distinct(edges$comic), "of", n_comic, "\n")
cat("Isolated characters:", n_char - n_distinct(edges$character), "\n")
cat("Isolated comics:", n_comic - n_distinct(edges$comic), "\n")

# Cleaning
### Duplicates (exact)
dup_labels <- vertices %>%
  group_by(kind, label) %>% 
  filter(n() > 1) %>% 
  arrange(kind, label)

print(dup_labels %>% ungroup() %>% distinct(kind, label) %>% count(kind))

### Doppelgangers (shared start of name)
doppel <- vertices %>% 
  filter(!type) %>%
  mutate(stem = str_trim(str_split_fixed(label, "[/|]", 2)[, 1])) %>%
  group_by(stem) %>% 
  filter(n() > 1) %>% 
  arrange(stem)

cat("Stems shared by 2+ character nodes:", n_distinct(doppel$stem), "\n")

# Bipartite igraph
g_bip <- graph_from_data_frame(
  d = edges %>% transmute(from = as.character(character), to = as.character(comic)),
  directed = FALSE,
  vertices = vertices %>% transmute(name = as.character(id), label, type)
)

# Incidence matrix B: rows = chars, cols = comics
B <- sparseMatrix(
  i = edges$character,
  j = edges$comic - n_char,
  x = 1,
  dims = c(n_char, n_comic),
  dimnames = list(vertices$label[1:n_char], vertices$label[(n_char + 1):n_total])
)

# Bipartite descriptive stats

char_labels <- vertices$label[1:n_char]
app_count <- rowSums(B) # appearances per character
cast_size <- colSums(B) # characters per comic
m <- sum(B) # number of edges
stopifnot(m == sum(app_count), m == sum(cast_size))

char_df <- tibble(id = 1:n_char, label = char_labels, appearances = as.numeric(app_count))
comic_df <- tibble(cast_size = as.numeric(cast_size))

### Summary table
describe_deg <- function(x, what) {
  tibble(measure = what, n = length(x), mean = mean(x), median = median(x),
         sd = sd(x), max = max(x), share_degree_1 = mean(x == 1),
         p90 = quantile(x, 0.90), p99 = quantile(x, 0.99))
}
deg_summary <- bind_rows(
  describe_deg(char_df$appearances, "Comics per character"),
  describe_deg(comic_df$cast_size,  "Characters per comic")
)
print(deg_summary)

### Bipartite density: fraction of all possible character-comic pairs that are edges.
bip_density <- m / (n_char * n_comic)
cat("Bipartite density:", signif(bip_density, 3), "\n")

### Most-appearing characters (your original plot)
p_top <- char_df %>%
  slice_max(appearances, n = 15) %>%
  ggplot(aes(reorder(label, appearances), appearances)) +
  geom_col(fill = "darkred") +
  coord_flip() +
  labs(title = "Top 15 characters by comic appearances", x = NULL, y = "Comics")
print(p_top)

### Degree distributions on log-log axes via the CCDF. CCDF because heavy tails. "What proportion of the vertices have degree at least k"
ccdf <- function(x) {
  tab <- as.data.frame(table(x), stringsAsFactors = FALSE)
  tab$x <- as.numeric(tab$x)
  tab$p <- rev(cumsum(rev(tab$Freq))) / sum(tab$Freq)
  tab
}
both <- bind_rows(
  ccdf(char_df$appearances) %>% mutate(dist = "Comics per character"),
  ccdf(comic_df$cast_size)  %>% mutate(dist = "Characters per comic")
)
p_ccdf <- ggplot(both, aes(x, p, colour = dist)) +
  geom_point(alpha = 0.7) +
  scale_x_log10() + scale_y_log10() +
  labs(title = "Degree distributions of the bipartite network",
       x = "Degree k", y = "P(degree >= k)", colour = NULL)
print(p_ccdf)

fit_pl <- function(x) tryCatch({ # Checking power law
  f <- fit_power_law(x[x > 0]); tibble(alpha = f$alpha, xmin = f$xmin, KS_p = f$KS.p)
}, error = function(e) tibble(alpha = NA, xmin = NA, KS_p = NA))
pl_fits <- bind_rows(
  fit_pl(char_df$appearances) %>% mutate(measure = "Comics per character"),
  fit_pl(comic_df$cast_size) %>% mutate(measure = "Characters per comic")
)
print(pl_fits)

# Projection --> From character-comic to character-character
### Take cross prod to get number of times each character appears with each other character
W <- tcrossprod(B) # B * t(B), W symmetric is and diagonal is appearance counts
stopifnot(all(abs(diag(W) - app_count) < 1e-9))

### turns the upper triangle into igraph graph
proj_to_igraph <- function(M, labels = char_labels) {
  T_ <- as(triu(M, k = 1), "TsparseMatrix")
  g <- make_empty_graph(n = nrow(M), directed = FALSE)
  g <- add_edges(g, as.vector(rbind(T_@i + 1L, T_@j + 1L)))
  E(g)$weight <- T_@x
  V(g)$label <- labels
  g
}

g_proj <- proj_to_igraph(W) # Weighted

inv_c <- ifelse(cast_size >= 2, 1 / (cast_size - 1), 0) # Penalizes large cast issues, solo issues get weight 0
W_newman <- tcrossprod(B %*% Diagonal(x = inv_c), B)   # = B diag(inv_c) B'
g_proj_nw <- proj_to_igraph(W_newman)

# Projection descriptive stats
net_summary <- function(g, name, distances = TRUE) { # Summary function for one-mode graphs
  comp <- components(g)
  giant <- which(comp$membership == which.max(comp$csize))
  gg <- induced_subgraph(g, giant)
  tibble(
    network = name,
    nodes_non_isol = sum(degree(g) > 0),
    edges = ecount(g),
    density = edge_density(g),
    mean_degree = mean(degree(g)),
    n_components = comp$no,
    giant_share = max(comp$csize) / vcount(g), # Share of ALL nodes in the largest component
    transitivity = transitivity(g, type = "global"), # P(two neighbours of a node are linked)
    assortativity = assortativity_degree(g), # <0: hubs link to low-degree nodes
    mean_distance = if (distances) mean_distance(gg, directed = FALSE) else NA_real_
  )
}
proj_summary <- net_summary(g_proj, "Character projection")
print(proj_summary)

# Unweighted projection (proposal: P_ik > 0)
g_proj_unw <- proj_to_igraph(W, labels = char_labels)
E(g_proj_unw)$weight <- 1  # drop shared-comic counts; keep topology only

# Safety for centrality calls that may be slow or fail. Returns NA vector of the right length.
safely <- function(expr, n) {
  tryCatch({
    out <- expr
    out[!is.finite(out)] <- NA_real_
    out
  }, error = function(e) rep(NA_real_, n))
}

# Bipartite centralities
# Exact degree on the full bipartite graph, but restrict to CHARACTER nodes
# (we only rank characters; comics are the affiliation side).
n_bip <- vcount(g_bip)
char_vids_bip <- seq_len(n_char)  # vertices 1..n_char are characters (by construction)

bip_deg   <- degree(g_bip)

# Closeness/betweenness/eigenvector: estimate, using only character nodes as
# sources so the result is comparable to the projection-side measures.
# estimate* returns a vector for the *sources* only.
bip_close <- safely(estimate_closeness(g_bip, vids = char_vids_bip,
                                       cutoff = 4), n_char)
bip_betw  <- safely(estimate_betweenness(g_bip, vids = char_vids_bip,
                                         cutoff = 4), n_char)
bip_eig   <- safely({
  ev <- eigen_centrality(g_bip, directed = FALSE, scale = TRUE)
  ev$vector[char_vids_bip]
}, n_char)

# Projection centralities
# Given a projected graph g (n = n_char nodes, same order as char_labels),
# return a tibble of the four centralities plus weighted strength.
proj_centralities <- function(g, tag, weighted = TRUE) {
  n <- vcount(g)
  w_raw <- if (weighted) E(g)$weight else NULL
  # distance-based measures need inverse weights (strong tie = short path)
  w_dist <- if (weighted) {
    w <- E(g)$weight
    # zero-weight edges (cast-size-1 comics) become a large finite distance
    ifelse(w > 0, 1 / w, 1e6)
  } else NULL
  
  tibble(
    id        = seq_len(n),
    !!paste0(tag, "_degree")      := degree(g),
    !!paste0(tag, "_strength")    := strength(g, weights = E(g)$weight),
    !!paste0(tag, "_closeness")   := safely(
      estimate_closeness(g, weights = w_dist, cutoff = 4), n),
    !!paste0(tag, "_betweenness") := safely(
      estimate_betweenness(g, weights = w_dist, cutoff = 4), n),
    !!paste0(tag, "_eigen")       := safely(
      eigen_centrality(g, directed = FALSE, weights = w_raw, scale = TRUE)$vector, n)
  )
}

cent_unw <- proj_centralities(g_proj_unw, "unw")
cent_w   <- proj_centralities(g_proj,     "w")
cent_nw  <- proj_centralities(g_proj_nw,  "nw")

# Assemble one character-level table
char_cent <- tibble(
  id     = seq_len(n_char),
  label  = char_labels,
  appearances = as.numeric(app_count),
  # bipartite-side
  bip_degree    = bip_deg[char_vids_bip],
  bip_closeness = bip_close,
  bip_betweenness = bip_betw,
  bip_eigen     = bip_eig
) %>%
  left_join(cent_unw, by = "id") %>%
  left_join(cent_w,   by = "id") %>%
  left_join(cent_nw,  by = "id")

# Sanity checks
stopifnot(nrow(char_cent) == n_char)
cat("NAs per centrality column:\n")
print(colSums(is.na(char_cent[, -c(1, 2)])))

# Top 10 by each centrality on the unweighted projection 
top_by <- function(df, col, k = 10) {
  df %>% arrange(desc(.data[[col]])) %>%
    slice_head(n = k) %>%
    transmute(label, value = .data[[col]], measure = col)
}
bind_rows(
  top_by(char_cent, "bip_betweenness"),
  top_by(char_cent, "unw_betweenness"),
  top_by(char_cent, "w_betweenness"),
  top_by(char_cent, "nw_betweenness")
) %>% print(n = 40)

# Strength
### Raw -> sum over comics of cast-1 because each comic links every character cast-1 others
str_raw <- as.numeric(strength(g_proj, weights=E(g_proj)$weight))

### Newman -> number of appearances in comics with 2+ characters because (cast-1)* (1/(cast-1)) = 1
str_nw <- as.numeric(strength(g_proj_nw, weights = E(g_proj_nw)$weight))

# Distance -> Bipartite distances is double unweighted projection because bipartite routest character - comic - character for projection's character - character.
## Character specific carries the same info, so differences for closeness would come from other choices rather than the structure
samp <- sample(which(degree(g_proj_unw) > 0), 150)
d_bip <- distances(g_bip, v=samp, to=samp)
d_proj <- distances(g_proj_unw, v=samp, to=samp, weights=NA)
fin <- is.finite(d_bip) & is.finite(d_proj)

# Longer Centrality Comparisons
RUN_THE_LONG_VERSION <- TRUE # Toggle on/off because this is gonna take a while
char_vids <- seq_len(n_char)

proj_cent <- function(g, tag, weighted) {
  w_len <- if (weighted) 1 / E(g)$weight else NA # NA = ignore weights
  w_str <- if (weighted) E(g)$weight else NA
  out <- tibble(id = seq_len(vcount(g)))
  if (!weighted) out[[paste0(tag, "_degree")]] <- as.numeric(degree(g)) # All three projections have the same edges, so their degrees are the same
  if (weighted) out[[paste0(tag, "_strength")]] <- as.numeric(strength(g, weights = w_str))
  out[[paste0(tag, "_harmonic")]] <- as.numeric(harmonic_centrality(g, weights = w_len, normalized = TRUE)) # Closeness only averages reachable nodes, harmonic includes unreachable nodes
  out[[paste0(tag, "_betweenness")]] <- if (RUN_THE_LONG_VERSION) {
    as.numeric(betweenness(g, directed = FALSE, weights = w_len))
  } else NA_real_
  out[[paste0(tag, "_eigen")]] <- as.numeric(eigen_centrality(g, directed = FALSE, weights = w_str, scale = TRUE)$vector)
  out
}

cent_unw2 <- proj_cent(g_proj_unw, "unw", weighted = FALSE)
cent_w2 <- proj_cent(g_proj, "w", weighted = TRUE)
cent_nw2 <- proj_cent(g_proj_nw, "nw", weighted = TRUE)

bip_degree_v <- as.numeric(degree(g_bip)[char_vids])
bip_harm <- as.numeric(harmonic_centrality(g_bip, vids = char_vids, normalized = TRUE))
bip_betw2 <- if (RUN_THE_LONG_VERSION) { # Bipartite betweenness counts parallel routes so there could be multiple shortest paths
  as.numeric(betweenness(g_bip, v = char_vids, directed = FALSE))
} else rep(NA_real_, n_char)
bip_eig2     <- as.numeric(eigen_centrality(g_bip, directed = FALSE, scale = TRUE)$vector[char_vids])

### Bringing it all together
cent_all <- tibble(
  id = char_vids,
  label = char_labels,
  bip_degree = bip_degree_v, # = appearances
  bip_harmonic = bip_harm,
  bip_betweenness = bip_betw2,
  bip_eigen = bip_eig2
) %>%
  left_join(cent_unw2, by = "id") %>%
  left_join(cent_w2,   by = "id") %>%
  left_join(cent_nw2,  by = "id")

measures <- setdiff(names(cent_all), c("id", "label"))
cat("Measures:", paste(measures, collapse = ", "), "\n")

if (exists("char_cent") && RUN_THE_LONG_VERSION) {
  cutoff_check <- tibble(
    measure = c("bip_betweenness", "unw_betweenness", "w_betweenness", "nw_betweenness"),
    spearman_cutoff_vs_exact = sapply(measure, function(ms)
      cor(char_cent[[ms]], cent_all[[ms]], method = "spearman", use = "complete.obs"))
  )
  print(cutoff_check)
}

# Comparisons
### Rank Correlation
rho <- cor(cent_all[, measures], method = "spearman", use = "pairwise.complete.obs")
print(round(rho, 2))

rho_long <- as.data.frame(as.table(rho)) %>%
  rename(m1 = Var1, m2 = Var2, rho = Freq) %>%
  mutate(m1 = factor(m1, levels = measures), m2 = factor(m2, levels = rev(measures)))

p_rho <- ggplot(rho_long, aes(m1, m2, fill = rho)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 2.3) +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred", midpoint = 0.5, limits = c(0, 1), name = "Spearman") +
  labs(title = "Do the measures rank characters the same way?", x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_rho)

### Top overlap
TOP_K <- 20 # How many to check for overlap

topk <- function(x, k) head(order(x, decreasing = TRUE), k)

overlap <- outer(measures, measures, Vectorize(function(a, b)
  length(intersect(topk(cent_all[[a]], TOP_K), topk(cent_all[[b]], TOP_K))) / TOP_K))

dimnames(overlap) <- list(measures, measures)
print(round(overlap, 2))

top_lists <- sapply(measures, function(ms) cent_all$label[topk(cent_all[[ms]], 10)])
print(as_tibble(top_lists), width = Inf)

### Movement in the top
rank_change <- function(df, a, b, min_app = 20, k = 10) {
  d <- df %>%
    mutate(rank_a = rank(-.data[[a]], ties.method = "average"),
           rank_b = rank(-.data[[b]], ties.method = "average"),
           shift = rank_a - rank_b) %>% # > 0 -> ranked higher
    filter(bip_degree >= min_app)
  list(better_under_b = d %>% slice_max(shift, n = k, with_ties = FALSE) %>%
         select(label, bip_degree, rank_a, rank_b, shift),
       worse_under_b = d %>% slice_min(shift, n = k, with_ties = FALSE) %>%
         select(label, bip_degree, rank_a, rank_b, shift))
}
change_degree <- rank_change(cent_all, "bip_degree", "unw_degree")
change_betw <- rank_change(cent_all, "bip_degree", "unw_betweenness")
print(change_degree)
print(change_betw)

p_scatter <- ggplot(filter(cent_all, unw_degree > 0), aes(bip_degree, unw_degree)) +
  geom_point(alpha = 0.3) + scale_x_log10() + scale_y_log10() +
  labs(title = "Appearances vs. distinct co-stars", x = "Comics per character (bipartite degree)", y = "Distinct co-stars (unweighted projected degree)")
print(p_scatter)

### Clustering
w_up <- as(triu(W, k = 1), "TsparseMatrix")@x # number of shared comics
C4   <- sum(choose(w_up, 2)) # 4-cycles
L3   <- sum((app_count[edges$character] - 1) * (cast_size[edges$comic - n_char] - 1)) # 3-edge paths
RA   <- 4 * C4 / L3 # Robins-Alexander

proj_sum <- net_summary(g_proj_unw, "Unweighted projection")
bip_comp <- components(g_bip)
struct <- tibble(
  quantity = c("Nodes", 
                 "Edges", 
                 "Mean degree of characters", 
                 "Median degree of characters",
                 "Max degree of characters", 
                 "Share of nodes in largest component",
                 "Clustering", 
                 "Mean character-character distance (hops of the projection)"),
  bipartite = c(vcount(g_bip), # Nodes
                 ecount(g_bip), # Edges
                 mean(app_count), # Mean
                 median(app_count), # Median
                 max(app_count), # Max
                 max(bip_comp$csize) / vcount(g_bip), # Share
                 RA, # Clustering
                 NA), # This is just 2x the projection
  projection = c(vcount(g_proj_unw), # Nodes
                 ecount(g_proj_unw), # Edges
                 mean(cent_all$unw_degree), # Mean
                 median(cent_all$unw_degree), # Median
                 max(cent_all$unw_degree), # Max
                 proj_sum$giant_share, # Share
                 proj_sum$transitivity, # Clustering
                 proj_sum$mean_distance) # Mean distance
)
print(struct, n = Inf, width = Inf)

distort <- bind_rows(
  ccdf(cent_all$bip_degree) %>% mutate(dist = "Comics per character (bipartite)"),
  ccdf(as.numeric(cast_size)) %>% mutate(dist = "Characters per comic (bipartite)"),
  ccdf(cent_all$unw_degree[cent_all$unw_degree > 0]) %>% mutate(dist = "Distinct co-stars (projection)")
)
p_ccdf_distort <- ggplot(distort, aes(x, p, colour = dist)) +
  geom_point(alpha = 0.7) + scale_x_log10() + scale_y_log10() +
  labs(title = "Degree distributions: bipartite vs. projected",
       x = "Degree k", y = "P(degree >= k)", colour = NULL)
print(p_ccdf_distort)
