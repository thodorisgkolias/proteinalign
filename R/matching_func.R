Norm <- function(X) {
  Y <- sqrt(sum(diag(X %*% t(X))))
  Y
}

SvdNew <- function(X) {
  SS <- svd(X)
  SS$d <- SS$d * c(1, 1, sign(det(SS$v %*% t(SS$u))))
  SS$u <- SS$u %*% diag(c(1, 1, sign(det(t(SS$u)))))
  SS$v <- SS$v %*% diag(c(1, 1, sign(det(t(SS$v)))))
  SS
}


IndRows <- function(comb2, mat) {
  ss1 <- which(rownames(mat) == comb2[1])
  return(ss1)
}

IndCols <- function(comb2, mat) {
  ss1 <- which(colnames(mat) == comb2[2])
  return(ss1)
}

VolFunc <- function(data) {
  vol <- prod(apply(data, 2, max))
  return(vol)
}


PAMLogLik <- function(seq1, seq2, p, PAM) {
  ## As in fallaize
  score <- 0
  if (!is.null(PAM)) {
    for (i in 1:p) {
      amino_row <- which(seq1[i] == rownames(PAM))
      amino_col <- which(seq2[i] == rownames(PAM))
      score <- score + log(PAM[amino_row, amino_col])
    }
  }
  return(score)
}

LogLikMatch <- function(pair = "NULL", data, seq, matched, vol, PAM,
                        gap_open, gap_ext) {
  ## Computes the likelihood of two matched molecules If pair != NULL the computes
  ## the likelihood of two matched molecules with the extra pair included
  M1 <- sapply(data, dim)[1, ]
  N <- length(M1)
  m <- dim(data[[1]])[2]
  hist <- matched
  hist <- rbind(hist, pair)
  pop <- sapply(M1, seq_len)
  Xnew <- list()
  p <- dim(hist)[1]
  seq_new <- list()
  for (i in 1:N) {
    matched_Xrows <- hist[, i]
    unmatched_Xrows <- which(!pop[[i]] %in% matched_Xrows)
    Xnew[[i]] <- data[[i]][c(matched_Xrows, unmatched_Xrows), ]
    seq_new[[i]] <- seq[[i]][c(matched_Xrows, unmatched_Xrows)]
    
  }
  EM_data <- lapply(Xnew, t)
  EM_data2 <- array(NA, dim = c(3, p, N))
  for (i in 1:N) {
    EM_data2[, , i] <- EM_data[[i]][, 1:p]
  }
  EM_data2 <- aperm(EM_data2, c(2, 1, 3))
  res <- EMSizeShape(EM_data2)
  sigma <- res$sigma
  lik <- NULL
  for (i in 1:N) {
    lik[i] <- LogLikSizeShape(res$rot.data[, , i], res$mean, sigma, p, M1[i], 
                              A = vol)
  }
  lik_str <- sum(lik, na.rm = TRUE)
  lik_seq <- PAMLogLik(seq_new[[1]], seq_new[[2]], p, PAM)
  gap <- GapPen(hist, gap_open, gap_ext, M1)
  lik0 <- lik_str + lik_seq + gap
  return(lik0)
}


HungarianStep <- function(data, seq, lik, lik0, matched, rows2, vol, PAM,
                          gap_open, gap_ext, comb2 = NULL) {
  ## Applies the Hungarian method for the matrix lik and returns the final alignment
  ## after testing if a pair can be added
  row_len <- sapply(rows2, length)
  mat <- matrix(lik, nrow = row_len[1], ncol = prod(row_len[-1]))
  rownames(mat) <- rows2[[1]]
  temp <- expand.grid(rows2[-1])
  colnames(mat) <- 1:prod(row_len[-1])
  if (all(mat < 0)) {
    y <- solve_LSAP(-mat)
  } else {
    y <- solve_LSAP(mat, maximum = TRUE)
  }
  # print(mem_used())
  ind <- cbind(1:dim(mat)[1], y)
  ali <- cbind(as.numeric(rownames(mat)), as.numeric(colnames(mat)[y]))
  liks <- mat[ind]
  tt <- cbind(ali, liks)
  tt <- tt[order(-tt[, 3]), ]
  if (!is.numeric(dim(tt))) {
    pair <- as.numeric(c(tt[1], temp[tt[2], ]))
    liknew <- LogLikMatch(pair, data, seq = seq, matched, vol, PAM = PAM,
                          gap_open = gap_open, gap_ext = gap_ext)
    if (liknew > lik0) {
      lik0 <- liknew
      matched <- rbind(matched, pair)
    }
  } else {
    k <- 1
    liknew <- NULL
    while (k <= dim(tt)[1]) {
      pair <- c(tt[k, 1], as.numeric(temp[tt[k, 2], ]))
      liknew <- LogLikMatch(pair, data, seq = seq, matched, vol, PAM = PAM, 
                            gap_open = gap_open, gap_ext = gap_ext)
      if (liknew > lik0) {
        lik0 <- liknew
        matched <- rbind(matched, pair)
        k <- k + 1
      } else {
        break
      }
    }
  }
  return(list(lik0 = lik0, hist = matched))
}

RemoveStep <- function(data, seq, matched, vol, lik0, PAM, gap_open, gap_ext) {
  remove <- TRUE
  hist <- matched
  while (remove) {
    lik_remove <- NULL
    for (i in 1:dim(hist)[1]) {
      remove_pair <- as.vector(hist[i, ])
      hist2 <- hist[-i, ]
      lik_remove[i] <- LogLikMatch(pair = NULL, data = data, seq = seq,
                                   matched = hist2, vol = vol, PAM = PAM,
                                   gap_open = gap_open, gap_ext = gap_ext)
    }
    if (max(lik_remove) > lik0) {
      print("Remove")
      pos_remove <- which.max(lik_remove)
      lik0 <- max(lik_remove)
      hist <- hist[-pos_remove, ]
      pnew <- dim(hist)[1]
      remove <- TRUE
    } else {
      print("No more to remove")
      remove <- FALSE
      pnew <- dim(hist)[1]
    }
  }
  
  return(list(lik0 = lik0, hist = hist))
}


GapPen <- function(seqs, gap_open, gap_ext, len) {
  seq <- seqs[order(seqs[, 1]), ]
  seq1 <- c(0, seq[, 1], len[1] + 1)
  seq2 <- c(0, seq[, 2], len[2] + 1)
  seq1 <- diff(seq1)
  seq2 <- diff(seq2)
  ff <- abs(c(seq1, seq2))
  gap <- 0
  for (i in 1:length(ff)) {
    if (ff[i] != 1) {
      gap <- gap + gap_open + (ff[i] - 2) * gap_ext
    }
  }
  return(-gap)
}

HungarianStep1 <- function(data, seq, lik, lik0, matched, comb2, rows2, vol,
                           PAM, gap_open, gap_ext) {
  ## Applies the Hungarian method for the matrix lik and returns the final alignment
  ## after testing if a pair can be added, when restrict search is TRUE
  
  row_len <- sapply(rows2, length)
  if (all(lik < 0)) {
    mat <- matrix(-10^5, nrow = row_len[1], ncol = prod(row_len[-1]))
  } else {
    mat <- matrix(10^5, nrow = row_len[1], ncol = prod(row_len[-1]))
    
  }
  rownames(mat) <- rows2[[1]]
  temp <- expand.grid(rows2[-1])
  colnames(mat) <- rows2[[2]]
  p1 <- apply(comb2, 1, IndRows, mat = mat)
  p2 <- apply(comb2, 1, IndCols, mat = mat)
  for (i in 1:dim(comb2)[1]) {
    mat[p1[i], p2[i]] <- lik[i]
  }
  if (all(mat < 0)) {
    y <- solve_LSAP(-mat)
  } else {
    y <- solve_LSAP(mat, maximum = TRUE)
  }
  ind <- cbind(1:dim(mat)[1], y)
  ali <- cbind(as.numeric(rownames(mat)), as.numeric(colnames(mat)[y]))
  liks <- mat[ind]
  tt <- cbind(ali, liks)
  tt <- tt[order(-tt[, 3]), ]
  if (!is.numeric(dim(tt))) {
    pair <- as.numeric(c(tt[1], temp[tt[2], ]))
    liknew <- LogLikMatch(pair, data, seq = seq, matched, vol, PAM = PAM,
                          gap_open = gap_open, gap_ext = gap_ext)
    if (liknew > lik0) {
      lik0 <- liknew
      matched <- rbind(matched, pair)
    }
  } else {
    k <- 1
    liknew <- NULL
    while (k <= dim(tt)[1]) {
      pair <- c(tt[k, 1], tt[k, 2])
      temp <- rbind(matched, pair)
      temp <- temp[order(temp[, 1]), ]
      if (all(diff(temp[, 1]) > 0) & all(diff(temp[, 2]) > 0)) {
        liknew <- LogLikMatch(pair, data, seq = seq, matched, vol, PAM = PAM, 
                              gap_open = gap_open, gap_ext = gap_ext)
        if (liknew > lik0) {
          lik0 <- liknew
          matched <- rbind(matched, pair)
          k <- k + 1
        } else {
          break
        }
      } else {
        k <- k + 1
      }
    }
  }
  return(list(lik0 = lik0, hist = matched))
}
