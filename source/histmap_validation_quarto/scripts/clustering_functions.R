outersect <- function(x, y) {
  sort(c(
    x[!x %in% y],
    y[!y %in% x]
  ))
}


mask_follow_up <- function(variable,
                           clustering_position,
                           reference_position,
                           reverse = FALSE,
                           ...) {
  if (reverse) {
    variable <- stringi::stri_reverse(variable)


    position <- unique(nchar(variable) - clustering_position + 1)
  } else {
    position <- clustering_position
  }


  variable_char <- as.character(variable)


  if (reverse) {
    clustered <- sapply(variable_char, function(x) {
      first_LU <- substr(x, start = position, stop = position)


      chars <- strsplit(x, "")[[1]]


      for (i in c((position + 1):length(chars))) {
        if (chars[i] == chars[position]) {
          chars[i:length(chars)] <- "X"


          break
        }
      }


      stringi::stri_reverse(paste(chars, collapse = ""))
    })


    clustered %>% unname()
  } else {
    clustered <- sapply(variable_char, function(x) {
      first_LU <- substr(x, start = position, stop = position)


      chars <- strsplit(x, "")[[1]]


      for (i in c((position + 1):length(chars))) {
        if (chars[i] == chars[position]) {
          chars[i:length(chars)] <- "X"


          break
        }
      }


      paste(chars, collapse = "")
    })


    clustered %>% unname()
  }
}


past_follow_up <- function(variable,
                           position,
                           LU = "1",
                           akker = 1, bebouwing = 2, bos = 3,
                           grasland = 4, heide = 5, strand = 6,
                           moeras = 7, boomgaard = 8, water = 9,
                           years = c(2022, 1969, 1879, 1774),
                           ...) {
  if (akker == 1 & bebouwing == 2 & bos == 3 & grasland == 4 &


    heide == 5 & strand == 6 & moeras == 7 & boomgaard == 8 &


    water == 9) {
    message("Message: Consistent codeling of LU was used")
  } else {
    warning("Warning: There is a change in LU coding")
  }


  sapply(variable, past_follow_up_helper,
    position = position,
    LU = LU, years = years
  )
}


past_follow_up_helper <- function(x, position, LU,
                                  years = c(2022, 1969, 1879, 1774)) {
  if (str_sub(x, position, position) != LU) {
    warning(paste0(
      "Message: That LU (", LU, ") was not registered in ",
      years[position], "."
    ))


    x <- x
  } else {
    while (position > 0) {
      if (str_sub(x, (position - 1), (position - 1)) != LU) {
        substr(x, position - 1, position - 1) <- "x"
      }


      position <- position - 1
    }
  }


  x
}


mask_follow_up_intensity <- function(variable,
                                     LU = 2,
                                     akker = 1, bebouwing = 2, bos = 3,
                                     grasland = 4, heide = 5, strand = 6,
                                     moeras = 7, boomgaard = 8, water = 9,
                                     clustering_position,
                                     reference_position, ...) {
  position <- clustering_position


  if (akker == 1 & bebouwing == 2 & bos == 3 & grasland == 4 &


    heide == 5 & strand == 6 & moeras == 7 & boomgaard == 8 &


    water == 9) {
    message("Message: Consistent codeling of LU was used")
  } else {
    warning("Warning: There is a change in LU coding")
  }


  extens_1 <- as.character(c(
    bos, grasland, heide, strand, moeras,
    boomgaard, water
  ))


  extens_2 <- as.character(c(bos, heide, strand, moeras, water))


  intens_1 <- as.character(c(akker, bebouwing))


  intens_2 <- as.character(c(grasland, boomgaard, bebouwing))


  positions_1 <- c(3, 4)


  positions_2 <- c(
    2


    # ,1
  )


  diff_extens <- setdiff(extens_1, extens_2)


  variable_char <- as.character(variable)


  variable_clustered <- sapply(variable_char, function(x) {
    first_LU <- substr(x, start = position, stop = position)


    chars <- strsplit(x, "")[[1]]


    for (i in c((position + 1):length(chars))) {
      if (chars[i] != chars[position]) {
        chars[i:length(chars)] <- "X"


        break
      }
    }


    paste(chars, collapse = "")
  })


  data <- data.frame(
    original = variable_char,
    clustered = variable_clustered
  )


  data <- data %>%
    rowwise() %>%
    mutate(
      location = regexpr("X", clustered)[1],
      location = ifelse(location == -1, 0, location),
      all_locations = list(c((location):nchar(original))),
      intens_1_position = any(positions_1 %in% unlist(all_locations)),
      intens_2_position = any(positions_2 %in% unlist(all_locations)),
      intens_1_present = any(
        strsplit(original, "")[[1]][positions_1] %in% intens_1
      ),
      intens_2_present = any(
        strsplit(original, "")[[1]][positions_2] %in% intens_2
      ),
      intensity = ifelse(intens_1_present | intens_2_present,
        paste0(clustered, "_I"),
        paste0(clustered, "_E")
      )
    ) %>%
    select(-c(location, clustered)) %>%
    rename(clustered = intensity)


  return(data$clustered)
}


make_bin_change <- function(variable,
                            clustering_position,
                            reference_position, ...) {
  ## check if LU changed based on a between the given reference year and


  ## clustering year


  variable_char <- as.character(variable)


  clustering <- sapply(variable_char, function(x) {
    first <- strsplit(x, "")[[1]][clustering_position]


    last <- strsplit(x, "")[[1]][reference_position]


    ifelse(first == last, "Not changed", "Changed")
  })


  clustering %>% unname()
}


change_frequency <- function(variable,
                             reference_position,
                             clustering_position) {
  ### making different groups based on reference year and a clustering year


  variable_char <- as.character(variable)


  clustering <- sapply(variable, function(x) {
    chars <- str_split(x, "")[[1]]


    unique_char <- length(unique(
      chars[clustering_position:reference_position]
    ))


    if (unique_char == 1) {
      value <- "Stable (No Change)"
    } else if (unique_char == 2 & chars[clustering_position] ==


      chars[reference_position]) {
      value <- "Reverted (1 LU in between)"
    } else if (unique_char == 2 & chars[clustering_position] !=


      chars[reference_position]) {
      value <- "Shifted (1 Change)"
    } else if (unique_char == 3 & chars[clustering_position] ==


      chars[reference_position]) {
      value <- "Reverted (2 LUs in between)"
    } else if (unique_char == 3 & chars[clustering_position] !=


      chars[reference_position]) {
      value <- "Shifted (2 LU Changes)"
    } else if (unique_char == length(chars)) {
      value <- "Fully Dynamic (All Different)"
    } else {
      value <- "Multiple changes"
    }


    value
  })


  clustering
}


# presence <- function(variable,


#                      clustering_position,


#                      reference_position,


#                      LU = LU, ...){


#


#   ## checking if a certain LU is present between the reference year and


#   ## clustering year


#   ## check if during this period it is present or is gained or lost


#


#   if(reference_position > length(chars)) stop("reference_position out of bounds")


#   if(clustering_position > reference_position) stop("invalid positions"


#


#


#   variable_char <- as.character(variable)


#   LU <- as.character(LU)


#


#   clustering <- sapply(variable, function(x){


#


#     chars <- str_split(x, "")[[1]]


#     chars <- chars[clustering_position:reference_position]


#     n <- length(chars)


#     count <- sum(LU == chars)


#     upper <- ceiling(n/2)


#


#     if(all(LU == chars)){


#       value <- "Stable presence"


#


#     }else if(all(LU != chars)){


#       value <- "Stable absence"


#


#     }else if(count <= upper & chars[n] == LU & chars[1] != LU){


#       value <- "Eventual gain"


#


#     }else if(count >= upper & chars[n] != LU & chars[1] == LU){


#       value <- "Eventual loss"


#


#     }else if(count < upper & chars[n] == LU){


#       value <- "Eventual gain"


#


#     }else if(count >= upper & chars[n] == LU){


#       value <- "Stable gain"


#


#     }else {


#       value <- "Intermittent Presence"


#       }


#


#     value


#   })


#


#   clustering


#


# }


presence <- function(variable, clustering_position, reference_position, LU) {
  LU <- as.character(LU)


  sapply(as.character(variable), function(x) {
    chars <- strsplit(x, "")[[1]]


    if (reference_position > length(chars)) {
      stop("reference_position out of bounds")
    }


    chars <- chars[clustering_position:reference_position]


    present <- chars == LU


    n <- length(chars)


    first <- present[1]


    last <- present[n]


    count <- sum(present)


    upper <- ceiling(n / 2)


    if (all(present)) {
      return("Stable presence")
    } else if (!any(present)) {
      return("Stable absence")
    } else if (first && !last && count <= upper) {
      return("Stable loss")
    } else if (first && !last && count > upper) {
      return("Eventual loss")
    } else if (!first && last && count <= upper) {
      return("Eventual gain")
    } else if (!first && last && count > upper) {
      return("Stable gain")
    } else {
      return("Intermittent Presence")
    }
  })
}


transition_present <- function(variable,
                               transition,
                               clustering_position,
                               reference_position, ...) {
  ## Check if a transition is present between a reference year and a clustering


  ## year


  ## Count how frequent the transition is present


  variable <- as.character(variable)


  transition <- as.character(transition)


  first_char_transition <- substr(transition, 1, 1)


  sapply(variable, function(x) {
    # Subsetten van periode obv gegeven referentie en cluster jaar


    x <- str_sub(x, start = (clustering_position - 1), end = reference_position)


    # Aantal keer dat de transitie voorkomt (inclusief overlappende matches)


    count <- str_count(x, paste0("(?=(", transition, "))"))


    # Posities van de transitie in x


    transition_positions <- str_locate_all(
      x,
      paste0("(?=(", transition, "))")
    )[[1]][, 1]


    if (length(transition_positions) == 0) {
      last_transition_position <- NA


      first_transition_position <- NA
    } else {
      last_transition_position <- max(transition_positions) # Laatste keer transitie start


      first_transition_position <- min(transition_positions) # Eerste keer transitie start
    }


    # Substring vóór de laatste transitie


    if (!is.na(first_transition_position) & first_transition_position > 1) {
      LUs_after_trans <- substr(x, 1, first_transition_position - 1)
    } else {
      LUs_after_trans <- ""
    }


    # Aantal keer dat first_char_transition voorkomt in het stuk vóór de laatste transitie


    nr_LU_after_trans <- sum(first_char_transition == str_split(LUs_after_trans, "")[[1]])


    # Controleren of de transitie aanwezig is in x


    x_contains_trans <- grepl(transition, x)


    # Logische condities voor correcte uitspraken


    if (x_contains_trans & count == 1 & nr_LU_after_trans == 0 &


      nr_LU_after_trans < nchar(LUs_after_trans)) {
      return("Transitie 1 keer aanwezig en laatste LU van de transitie duikt niet opnieuw op")
    } else if (x_contains_trans & count == 1 & nr_LU_after_trans > 0 &


      nr_LU_after_trans < nchar(LUs_after_trans)) {
      return("Transitie 1 keer aanwezig en laatste LU van de transitie duikt opnieuw op")
    } else if (x_contains_trans & count == 1 & nr_LU_after_trans > 0 &


      nr_LU_after_trans == nchar(LUs_after_trans)) {
      return("Transitie 1 keer aanwezig en laatste LU van de blijft angehouden")
    } else if (x_contains_trans & count > 1 & nr_LU_after_trans == 0 &


      nr_LU_after_trans < nchar(LUs_after_trans)) {
      return("Transitie minstens 1 keer aanwezig en laatste LU van laatste transitie duikt niet opnieuw op")
    } else if (x_contains_trans & count > 1 & nr_LU_after_trans > 0 &


      nr_LU_after_trans < nchar(LUs_after_trans)) {
      return("Transitie minstens 1 keer aanwezig en laatste LU van laatste transitie duikt opnieuw op")
    } else if (x_contains_trans & count > 1 & nr_LU_after_trans > 0 &


      nr_LU_after_trans == nchar(LUs_after_trans)) {
      return("Transitie minstens 1 keer aanwezig en laatste LU van laatste blijft aanhouden")
    } else if (x_contains_trans & first_transition_position == 1) {
      return("Aanwezig en laatste LU van transitie blijft aangehouden")
    } else if (x_contains_trans & nr_LU_after_trans > 0) {
      return("Aanwezig maar laatste LU wel aanwezig")
    } else if (x_contains_trans & first_transition_position != 1 & nr_LU_after_trans == 0) {
      return("Aanwezig en laatste LU niet aanwezig")
    } else if (!x_contains_trans & grepl(first_char_transition, x)) {
      return("Niet aanwezig maar laatste LU wel aanwezig")
    } else if (!x_contains_trans & !grepl(first_char_transition, x)) {
      return("Niet aanwezig en laatste LU niet aanwezig")
    } else {
      return("Niet aanwezig en laatste LU niet aanwezig")
    }
  })
}


clustering <- function(years = c(2022, 1969, 1873, 1774),
                       prediction = data$observed,
                       reference = data$reference,
                       LU = "1",
                       transition = "5555",
                       clustering_year = 2022,
                       clustering_method = "mask_follow_up",
                       akker = 1, bebouwing = 2, bos = 3,
                       grasland = 4, heide = 5, strand = 6,
                       moeras = 7, boomgaard = 8, water = 9,
                       reference_year = 1774) {
  ### check properties of the clustering method


  # assert_that(length(clustering_method) == 1,


  #             msg = "Please provide 1 clustering method")


  assert_that(
    all(clustering_method %in% c(
      "mask_follow_up",
      "past_follow_up",
      "mask_follow_up_intensity",
      "make_bin_change",
      "change_frequency",
      "transition",
      "presence"
    )),
    msg = "clustering_method should be either mask_follow_up, mask_follow_up_intensity, make_bin_change, past_follow_up, change_frequency, transition, presence"
  )


  ### Check properties of coding LU


  assert_that(
    is.numeric(c(
      akker, bebouwing, bos, grasland,
      heide, strand, moeras, boomgaard,
      water
    )),
    msg = "The provided LU should be in numerical form"
  )


  if (akker == 1 & bebouwing == 2 & bos == 3 & grasland == 4 &


    heide == 5 & strand == 6 & moeras == 7 & boomgaard == 8 &


    water == 9) {
    message("Consistent coding of LU was used")
  } else {
    warning("There is a change in LU coding")
  }


  ### Check properties of supplied years


  assert_that(identical(years, sort(years, decreasing = TRUE)),
    msg = "Supplied years are not in descending order"
  )


  assert_that(length(years) == 4,
    msg = "The number of provided years is different from 4"
  )


  assert_that(all(years %in% c(2022, 1969, 1873, 1774)),
    msg = "One ore more of the provided year(s) is/are not 2022, 1969, 1873 or 1774"
  )


  assert_that(is.numeric(years),
    msg = "The supplied years are not in numerical form"
  )


  assert_that(is.numeric(reference_year) & is.numeric(clustering_year),
    msg = "Reference year and clustering year should be numeric"
  )


  assert_that(clustering_year %in% years,
    msg = "The clustering year is not in the provided years"
  )


  assert_that(length(clustering_year) == 1,
    msg = "Provide 1 clusterig year"
  )


  assert_that(reference_year %in% years,
    msg = "The reference year is not in the provided years"
  )


  assert_that(length(reference_year) == 1,
    msg = "Provide 1 reference year"
  )


  assert_that(reference_year != clustering_year,
    msg = "Clustering year and reference year are the same"
  )


  assert_that(clustering_year > reference_year,
    msg = "Clustering year is prior to reference year"
  )


  ### Checking properties of the supplied data

  assert_that((any(str_detect(reference, "x", negate = TRUE)) &
                 any(str_detect(reference, "X", negate = TRUE)) &
                 any(str_detect(prediction, "x", negate = TRUE)) &
                 any(str_detect(prediction, "X", negate = TRUE))),
              msg = "the input already has been clustered (contains an x or an X")


  assert_that(length(prediction) == length(reference),
    msg = "Prediction and reference are not of the same length"
  )


  assert_that(all(levels(prediction) == levels(reference)),
    msg = "Prediction and reference have not the same levels"
  )


  assert_that(
    all(nchar(as.character(prediction)) ==


      nchar(as.character(prediction[1]))),
    msg = "Some predictions do not have the same length"
  )


  assert_that(
    all(nchar(as.character(reference)) ==


      nchar(as.character(reference[1]))),
    msg = "Some references do not have the same length"
  )


  assert_that(
    nchar(as.character(reference[1])) ==


      nchar(as.character(prediction[1])) &


      nchar(as.character(reference[1])) == length(years),
    msg = "The reference and prediction differ in length of characters"
  )


  ### Check properties of supplied LU


  assert_that(is.character(LU) & nchar(LU) == 1,
    msg = "Provided LU not in right format"
  )


  ### Check properties of transition


  assert_that(is.character(transition) & nchar(transition) > 1,
    msg = "Provided transition is in wrong format"
  )


  ### Converting supplied years to the position in the given years


  clustering_position <- which(clustering_year == years)[[1]]


  reference_position <- which(reference_year == years)[[1]]


  assert_that(clustering_position < reference_position,
    msg = "Clustering year is prior to the reference year"
  )


  tempdata <- data.frame(
    observed = prediction,
    reference = reference
  )


  for (method in clustering_method) {
    if (method == "mask_follow_up") {
      tempdata[[paste0("prediction_", method)]] <- mask_follow_up(prediction,
        clustering_position =


          clustering_position,
        reference_position =


          reference_position
      )


      tempdata[[paste0("reference_", method)]] <- mask_follow_up(reference,
        clustering_position =


          clustering_position,
        reference_position =


          reference_position
      )


      ## set them at the same level


      all_levels <- lubridate::union(
        levels(factor(
          tempdata[[paste0("prediction_", method)]]
        )),
        levels(factor(
          tempdata[[paste0("reference_", method)]]
        ))
      )


      tempdata[[paste0("reference_", method)]] <- factor(tempdata[[paste0("reference_", method)]],
        levels = all_levels
      )


      tempdata[[paste0("prediction_", method)]] <- factor(tempdata[[paste0("prediction_", method)]],
        levels = all_levels
      )
    }


    if (method == "past_follow_up") {
      tempdata[[paste0("prediction_", method)]] <- past_follow_up(variable = as.character(prediction),
                                                                  position =
                                                                    clustering_position,
                                                                  LU = LU,
                                                                  years = years)


      tempdata[[paste0("reference_", method)]] <- past_follow_up(variable = as.character(reference),
                                                                 position =
                                                                   reference_position,
                                                                 LU = LU,
                                                                 years = years,
      )


      ## set them at the same level


      all_levels <- lubridate::union(
        levels(factor(
          tempdata[[paste0("prediction_", method)]]
        )),
        levels(factor(
          tempdata[[paste0("reference_", method)]]
        ))
      )


      tempdata[[paste0("reference_", method)]] <- factor(tempdata[[paste0("reference_", method)]],
                                                         levels = all_levels
      )


      tempdata[[paste0("prediction_", method)]] <- factor(tempdata[[paste0("prediction_", method)]],
                                                          levels = all_levels
      )
    }


    if (method == "mask_follow_up_intensity") {
      tempdata[[paste0("prediction_", method)]] <- mask_follow_up_intensity(prediction,
        clustering_position =


          clustering_position,
        reference_position =


          reference_position
      )


      tempdata[[paste0("reference_", method)]] <- mask_follow_up_intensity(reference,
        clustering_position =


          clustering_position,
        reference_position =


          reference_position
      )


      ## set them at the same level


      all_levels <- lubridate::union(
        levels(factor(
          tempdata[[paste0("prediction_", method)]]
        )),
        levels(factor(
          tempdata[[paste0("reference_", method)]]
        ))
      )


      tempdata[[paste0("reference_", method)]] <- factor(tempdata[[paste0("reference_", method)]],
        levels = all_levels
      )


      tempdata[[paste0("prediction_", method)]] <- factor(tempdata[[paste0("prediction_", method)]],
        levels = all_levels
      )
    }


    if (method == "make_bin_change") {
      tempdata[[paste0("prediction_", method)]] <- make_bin_change(prediction,
        clustering_position =


          clustering_position,
        reference_position =


          reference_position
      )


      tempdata[[paste0("reference_", method)]] <- make_bin_change(reference,
        clustering_position =


          clustering_position,
        reference_position =


          reference_position
      )


      ## set them at the same level


      all_levels <- lubridate::union(
        levels(factor(
          tempdata[[paste0("prediction_", method)]]
        )),
        levels(factor(
          tempdata[[paste0("reference_", method)]]
        ))
      )


      tempdata[[paste0("reference_", method)]] <- factor(tempdata[[paste0("reference_", method)]],
        levels = all_levels
      )


      tempdata[[paste0("prediction_", method)]] <- factor(tempdata[[paste0("prediction_", method)]],
        levels = all_levels
      )
    }


    if (method == "change_frequency") {
      tempdata[[paste0("prediction_", method)]] <- change_frequency(prediction,
        clustering_position =


          clustering_position,
        reference_position =


          reference_position
      )


      tempdata[[paste0("reference_", method)]] <- change_frequency(reference,
        clustering_position =


          clustering_position,
        reference_position =


          reference_position
      )


      ## set them at the same level


      all_levels <- lubridate::union(
        levels(factor(
          tempdata[[paste0("prediction_", method)]]
        )),
        levels(factor(
          tempdata[[paste0("reference_", method)]]
        ))
      )


      tempdata[[paste0("reference_", method)]] <- factor(tempdata[[paste0("reference_", method)]],
        levels = all_levels
      )


      tempdata[[paste0("prediction_", method)]] <- factor(tempdata[[paste0("prediction_", method)]],
        levels = all_levels
      )
    }


    if (method == "presence") {
      tempdata[[paste0("prediction_", method)]] <- presence(prediction,
        clustering_position =


          clustering_position,
        reference_position =


          reference_position,
        LU = LU
      )


      tempdata[[paste0("reference_", method)]] <- presence(reference,
        clustering_position =


          clustering_position,
        reference_position =


          reference_position,
        LU = LU
      )


      ## set them at the same level


      all_levels <- lubridate::union(
        levels(factor(
          tempdata[[paste0("prediction_", method)]]
        )),
        levels(factor(
          tempdata[[paste0("reference_", method)]]
        ))
      )


      tempdata[[paste0("reference_", method)]] <- factor(tempdata[[paste0("reference_", method)]],
        levels = all_levels
      )


      tempdata[[paste0("prediction_", method)]] <- factor(tempdata[[paste0("prediction_", method)]],
        levels = all_levels
      )
    }


    if (method == "transition") {
      tempdata[[paste0("prediction_", method)]] <- transition_present(prediction,
        transition = transition,
        clustering_position =


          clustering_position,
        reference_position =


          reference_position
      )


      tempdata[[paste0("reference_", method)]] <- transition_present(reference,
        transition = transition,
        clustering_position =


          clustering_position,
        reference_position =


          reference_position
      )


      ## set them at the same level


      all_levels <- lubridate::union(
        levels(factor(
          tempdata[[paste0("prediction_", method)]]
        )),
        levels(factor(
          tempdata[[paste0("reference_", method)]]
        ))
      )


      tempdata[[paste0("reference_", method)]] <- factor(tempdata[[paste0("reference_", method)]],
        levels = all_levels
      )


      tempdata[[paste0("prediction_", method)]] <- factor(tempdata[[paste0("prediction_", method)]],
        levels = all_levels
      )
    }
  }


  return(tempdata)
}


clustering_area <- function(variable,
                            clustering_method = "change_frequency") {
  tempdata <- clustering(
    reference = variable,
    prediction = variable,
    clustering_method = clustering_method
  )


  tempdata <- tempdata %>%
    select(reference, paste0("reference_", clustering_method))
}
