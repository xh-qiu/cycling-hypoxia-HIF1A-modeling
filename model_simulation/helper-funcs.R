# helper functions
# Function to compare two values with a tolerance level
compareWithTolerance <- function(x, y, epsilon, f) {
  if (abs(x - y) <= epsilon || abs(x - y) <= f * max(abs(x), abs(y))) {
    return("~")  # Considered equal
  } else if (x < y) {
    return("<")  # x is less than y
  } else {
    return(">")  # x is greater than y
  }
}

# Function to determine the order of H, N, O based on their values
getHNOorder <- function(H, N, O, epsilon=1e-4, fraction=1e-2) {
  HN = compareWithTolerance(H, N, epsilon, fraction)
  HO = compareWithTolerance(H, O, epsilon, fraction)
  NO = compareWithTolerance(N, O, epsilon, fraction)
  
  if (HN == "<") {
    if (NO == "<") {
      return("H<N<O")
    } else if (NO == "~") {
      return("H<N~O")
    } else {
      if (HO == "<") {
        return("H<O<N")
      } else if (HO == "~") {
        return("O~H<N")
      } else {
        return("O<H<N")
      }
    }
  } else if (HN == "~") {
    if (HO == "<") {
      return("H~N<O")
    } else if (HO == "~") {
      return("H~N~O")
    } else {
      return("O<H~N")
    }
  } else { # HN == ">"
    if (HO == "<") {
      return("N<H<O")
    } else if (HO == "~") {
      return("N<H~O")
    } else {
      if (NO == "<") {
        return("N<O<H")
      } else if (NO == "~") {
        return("O~N<H")
      } else {
        return("O<N<H")
      }
    }
  }
}

#Vectorize the function
getHNOorder <- Vectorize(getHNOorder)