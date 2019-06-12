# A <- shortest_path(graph, src, dest)
# 
# a <- igraph::shortest_paths(graph, src, dest,  mode = "out")
# 
# a %>%
#   map("vpath")
# 
# library(tidygraph)
# tbl_g <- tidygraph::as_tbl_graph(graph)
# 
# 
# edges <- tbl_g %>% tidygraph::activate("edges") %>% as_tibble()
# 
# start_index <- tbl_g %>%
#   activate("edges") %>%
#   mutate(from_name = .N()$name[from]) %>%
#   filter(from_name == src) %>%
#   as_tibble() %>%
#   pull(from) %>%
#   unique()
# 
# end_index <- map(dest, ~{
#   tbl_g %>%
#     activate("edges") %>%
#     mutate(to_name = .N()$name[to]) %>%
#     filter(to_name == .x) %>%
#     as_tibble() %>%
#     pull(to) %>%
#     unique()
#   })
# 
# pb <- progress::progress_bar$new(total = length(end_index))
# paths <- end_index %>%
#   map(~{
#   pb$tick()
#   yenpathy::k_shortest_paths(edges, start_vertex = start_index, end_vertex = .x, k = 2, verbose = FALSE)
# })
# 
# 
# bench <- bench::mark(
#   end_index %>%
#     map(~{
#       # print(.x)
#       yenpathy::k_shortest_paths(edges, start_vertex = start_index, end_vertex = .x, k = 1, verbose = FALSE)
#     })
# )
# 
# 
# 
# bench::mark(
#   for (i in end_index) {
#     # print(i)
#     yenpathy::k_shortest_paths(edges, start_vertex = start_index, end_vertex = i, k = 2, verbose = FALSE)
#   }
# )
# 
# yenpathy::k_shortest_paths(edges, start_vertex = start_index, end_vertex = end_index, k = 2)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# str(a)
# first_paths <- a$vpath %>% map(names)
# 
# first_paths %>%
#   map(~{
# 
#   })
# 
# B <- list()
# for (k_i in 2:k){
#   prev_path <- A[[k_i-1]] # Taake the path found
#   num_nodes_to_loop <- length(prev_path)-1
#   for(i in 1:num_nodes_to_loop){ # Loop over each of the nodes, deleting and finding a new path  on the total graph
#     spurNode <- prev_path[i] # spur might mean current node of inquiry???
#     rootPath <- prev_path[1:i]
# 
#     edgesToDelete <- find_edges_to_delete(A, i,rootPath) # A  is overall path list, i is the index of the current node, root_paath is the nodes  from the  staart up unti the index
#     t_g <- delete.edges(graph, edgesToDelete)
#     # for (edge in edgesToDelete) t_g <- delete.edges(t_g, edge)
# 
#     spurPath <- shortest_path(t_g,spurNode, dest)
# 
#     if (!is.null(spurPath)){
#       total_path <- unlist(c(rootPath[-i], spurPath))
#       # if (!list(total_path) %in% B) B[length(B)+1] <- total_path
#       # if (!purrr::has_element(B, total_path)) B[length(B)+1] <- total_path
#       if (!purrr::has_element(B, total_path)){
#         B <-  append(B, list(total_path))
#       }
#     }
#   }
#   if (length(B) == 0) break
#   B <- sort_paths(graph, B)
#   A[k_i] <- B[1]
#   B <- B[-1]
# }
