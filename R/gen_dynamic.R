#' Generate dynamic networks and simulate games
#'
#' This function generates dynamic networks based on the specified type
#' and simulates various types of games on the network.
#'
#' @param N Number of nodes in the network
#' @param L Number of simulation steps
#' @param aver_k Average degree of nodes (default is 6)
#' @param net_type Type of network to generate ("BA", "ER", "SW")
#' @param game_type Type of game to simulate ("PDG", "SDG", "DG", "SUG", etc.)
#' @return A list containing simulation results including accumulated payoff matrix,
#'         initial and final adjacency matrices, and payoff matrices for each node.
#' @export
#' @importFrom stats runif
#' @import igraph

gen_dynamic <- function(N, L, aver_k = 6, net_type = "BA", game_type = "PDG") {
  # Function to calculate payoff for SUG game
  payoff_sug <- function(s_i, s_j, S) {
    p_i <- S[1, s_i]
    q_j <- S[2, s_j]
    p_j <- S[1, s_j]
    q_i <- S[2, s_i]
    p_sug <- 0 # Initialize p_sug
    if ((p_i >= q_j) & (p_j >= q_i)) {
      p_sug <- p_j + 1 - p_i
    } else if ((p_i >= q_j) & (p_j < q_i)) {
      p_sug <- 1 - p_i
    } else if ((p_i < q_j) & (p_j >= q_i)) {
      p_sug <- p_j
    }
    return(p_sug)
  }
  
  # Parameter validation
  if (!is.numeric(N) || N <= 0) {
    stop("N must be a positive numeric value.")
  }
  if (!is.numeric(L) || L <= 0) {
    stop("L must be a positive numeric value.")
  }
  if (!is.numeric(aver_k) || aver_k <= 0) {
    stop("aver_k must be a positive numeric value.")
  }
  valid_net_types <- c("BA", "ER", "SW")
  if (!(net_type %in% valid_net_types)) {
    stop("net_type must be one of 'BA', 'ER', or 'SW'.")
  }
  valid_game_types <- c("PDG", "SDG", "DG", "SUG", "CG", "SG", "EG", "SeqG", "RG", "ZZG", "NCG", "IIG")
  if (!(game_type %in% valid_game_types)) {
    stop("game_type must be one of 'PDG', 'SDG', 'DG', 'SUG', 'CG', 'SG', 'EG', 'SeqG', 'RG', 'ZZG', 'NCG', 'IIG'.")
  }
  
  # Function to calculate payoff using matrix multiplication
  payoff <- function(si, sj, M) {
    return(t(as.matrix(si)) %*% M %*% as.matrix(sj))
  }
  
  # Generate network based on specified type
  generate_network <- function(N, aver_k, net_type) {
    switch(net_type,
           "BA" = igraph::sample_pa(N, power = 1, m = aver_k, directed = FALSE),
           "ER" = {
             p <- aver_k / (N - 1)
             igraph::erdos.renyi.game(n = N, p = p, type = "gnp", directed = FALSE)
           },
           "SW" = igraph::sample_smallworld(dim = 1, size = N, nei = aver_k/2, p = .0001),
           stop("Invalid net_type. Choose one of 'BA', 'ER', or 'SW'.")
    )
  }
  
  # Generate network
  G <- generate_network(N, aver_k, net_type)
  
  # Initialize Payoff Matrix
  payoff_matrices <- list(
    PDG = matrix(c(1, 0, 1.2, 0), 2, 2, byrow = TRUE),
    SDG = matrix(c(1, 1 - 1.2, 1 + 1.2, 0), 2, 2, byrow = TRUE),
    DG = matrix(c(2 - 1, -1, 2, 0), 2, 2, byrow = TRUE),
    CG = matrix(c(2, 2, 2, 2), 2, 2, byrow = TRUE),
    SG = matrix(c(1, -1, 0, 1), 2, 2, byrow = TRUE),
    EG = matrix(c(3, 0, 5, 1), 2, 2, byrow = TRUE),
    SeqG = matrix(c(2, 1, 0, 1), 2, 2, byrow = TRUE),
    RG = matrix(c(1, -1, -1, 1), 2, 2, byrow = TRUE),
    ZZG = matrix(c(1, -1, -1, 1), 2, 2, byrow = TRUE),
    NCG = matrix(c(1, 0, 0, 1), 2, 2, byrow = TRUE),
    IIG = matrix(c(1, -1, -1, 1), 2, 2, byrow = TRUE)
  )
  
  # Handle initialization of game matrix for SUG game
  M <- ifelse(game_type == "SUG", matrix(runif(4), nrow = 2), payoff_matrices[[game_type]])
  
  # Initialize Strategy
  S <- matrix(ifelse(game_type == "SUG", runif(N * 2), sample(c(0, 1), N * 2, replace = TRUE)), nrow = 2)
  
  # Record Initial Adjacency Matrix
  initial_adj <- as_adjacency_matrix(G)
  
  # Init X, Y
  Y <- matrix(NA, nrow = L, ncol = N) # Each column represents accumulated payoff at each time step
  X <- vector("list", N) # List containing payoff matrices for each node
  for (node_index in 1:N) {
    X[[node_index]] <- matrix(NA, nrow = L, ncol = N)
  }
  
  # Simulation process
  for (step in 1:L) {
    F_temp <- matrix(NA, nrow = N, ncol = N)
    if (game_type == "SUG") {
      for (u in 1:N) {
        for (v in 1:N) {
          F_temp[u, v] <- payoff_sug(u, v, S)
        }
      }
    } else {
      for (u in 1:N) {
        for (v in 1:N) {
          F_temp[u, v] <- payoff(S[, u], S[, v], M)
        }
      }
    }
    
    # Record Y and X
    F_cumulate <- rowSums(initial_adj * F_temp)
    Y[step, ] <- F_cumulate
    for (node_index in 1:N) {
      X[[node_index]][step, ] <- F_temp[node_index, ]
    }
    
    # Fermi Updating (Synchronization Update)
    for (i in 1:N) {
      neighbor_indices <- igraph::neighbors(G, i, "all")
      if (length(neighbor_indices) > 0) {
        update_index <- sample(as.vector(neighbor_indices), 1)
        w <- 1 / (1 + exp((F_cumulate[i] - F_cumulate[update_index]) / .1)) # Fermi parameter
        if (runif(1) >= w) {
          if (game_type == "SUG") {
            S[, i] <- S[, update_index] + 0.05
          } else {
            S[, i] <- S[, update_index]
          }
        }
      }
    }
  }
  
  # Record Final Adjacency Matrix
  final_adj <- as_adjacency_matrix(G)
  
  return(list(Y = Y, X = X, initial_adj = initial_adj, final_adj = final_adj))
}