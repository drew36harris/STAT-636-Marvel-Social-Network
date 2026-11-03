# Libraries
library(igraph)
library(Matrix)
library(dplyr)
library(stringr)
library(ggplot2)

# Read and Parse
path <-  "MarvelSocialNetwork.txt"
lines <- readLines(path, encoding="latin1")

### Get Counts
hdr <- as.integer(str_extract_all(lines[1], "\\d+")[[1]])
n_total <- hdr[1]
n_char <- hdr[2]
n_comic <- n_total - n_char

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

### Edge labels
elines <- lines[(edge_start + 1):length(lines)]
elines <- elines[nzchar(trimws(elines))]
tok <- strsplit(trimws(elines), "\\s+")
edges <- bind_rows(lapply(tok, function(t) {
  tibble(character = as.integer(t[1]), comic = as.integer(t[-1]))
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
cat("Isolated characters:", n_char  - n_distinct(edges$character), "\n")
cat("Isolated comics:", n_comic - n_distinct(edges$comic), "\n")

# Cleaning
### Duplicates (exact)
dup_labels <- vertices %>%
  group_by(kind, label) %>% 
  filter(n() > 1) %>% 
  arrange(kind, label)

dup_labels %>% ungroup() %>% distinct(kind, label) %>% count(kind)

### Doppelgangers (shared start of name)
doppel <- vertices %>% 
  filter(!type) %>%
  mutate(stem = str_trim(str_split_fixed(label, "[/|]", 2)[, 1])) %>%
  group_by(stem) %>% 
  filter(n() > 1) %>% 
  arrange(stem)

n_distinct(doppel$stem)

# Bipartite igraph
g_bip <- graph_from_data_frame(
  d = edges %>% transmute(from = as.character(character), to = as.character(comic)),
  directed = FALSE,
  vertices = vertices %>% transmute(name = as.character(id), label, type)
)

# Incidence matrix B
B <- sparseMatrix(
  i = edges$character,
  j = edges$comic - n_char,
  x = 1,
  dims = c(n_char, n_comic),
  dimnames = list(vertices$label[1:n_char], vertices$label[(n_char + 1):n_total])
)

# Graphs
### Support
char_labels <- vertices$label[1:n_char]

app_count <- rowSums(B)
cast_size <- colSums(B)

char_df <- tibble(
  id          = 1:n_char,
  label       = char_labels,
  appearances = as.numeric(app_count)
)
comic_df <- tibble(cast_size = as.numeric(cast_size))

### Most appearing characters
char_df %>%
  slice_max(appearances, n = 15) %>%
  ggplot(aes(reorder(label, appearances), appearances)) +
  geom_col(fill = "darkred") +
  coord_flip() +
  labs(title = "Top 15 characters by comic appearances", x = NULL, y = "Comics")

### Degree distributions
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
ggplot(both, aes(x, p, colour = dist)) +
  geom_point(alpha = 0.7) +
  scale_x_log10() + scale_y_log10() +
  labs(title = "Degree distributions of the bipartite network",
       x = "Degree k", y = "P(degree >= k)", colour = NULL)

