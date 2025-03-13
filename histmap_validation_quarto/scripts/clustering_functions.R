outersect <- function(x, y) {
  sort(c(x[!x%in%y],
         y[!y%in%x]))
}

mask_follow_up <- function(variable){

  position <- clustering_position

  variable_char <- as.character(variable)

  clustered <- sapply(variable_char, function(x){
        first_LU <- substr(x, start = position, stop = position)
    chars <- strsplit(x, "")[[1]]
    for(i in c((position + 1):length(chars))) {

      if(chars[i] != chars[position]){
        chars[i:length(chars)] <- "X"
        break
      }

    }
    paste(chars, collapse = "")

  })

  clustered %>% unname()

}

mask_follow_up_intensity <- function(variable,
                                     akker = 1, bebouwing = 2, bos = 3,
                                     grasland = 4, heide = 5, strand = 6,
                                     moeras = 7, boomgaard = 8, water = 9){

  position <-  clustering_position

  if(akker == 1 & bebouwing == 2 & bos == 3 & grasland == 4 &
     heide == 5 & strand == 6 & moeras == 7 & boomgaard == 8 &
     water == 9){

    message("Message: Consistent codeling of LU was used")

  } else{
    warning("Warning: There is a change in LU coding")
  }

  extens_1 <- as.character(c(bos, grasland, heide, strand, moeras,
                             boomgaard, water))
  extens_2 <- as.character(c(bos, heide, strand, moeras, water))
  intens <- as.character(c(akker, bebouwing))
  diff_extens <- setdiff(extens_1, extens_2)

  variable_char <- as.character(variable)

  variable_clustered <- sapply(variable_char, function(x){
    first_LU <- substr(x, start = position, stop = position)
    chars <- strsplit(x, "")[[1]]
    for(i in c((position + 1):length(chars))) {

      if(chars[i] != chars[position]){
        chars[i:length(chars)] <- "X"
        break
      }

    }
    paste(chars, collapse = "")


  })

  data <- data.frame(original = variable_char,
                     clustered = variable_clustered)

  data <- data %>% rowwise() %>%
    mutate(location = regexpr("X", clustered)[1],
           location = ifelse(location == -1, 0, location),
           intensity = ifelse(any(intens %in%
                                    strsplit(original, "")[[1]][c((location):
                                                                    nchar(original)
                                    )]) &
                                location != -1,
                              paste0(clustered, "_I"),
                              paste0(clustered, "_E")),
           intensity = ifelse(any(location == 2 &
                                    strsplit(original, "")[[1]][location] %in%
                                    diff_extens),
                              paste0(clustered, "_I"),
                              paste0(intensity, ""))
    ) %>% select(-c(location, clustered)) %>%
    rename(clustered = intensity)

  return(data$clustered)


}

make_bin_change <- function(variable){

  variable_char <- as.character(variable)

  clustering <- sapply(variable_char, function(x){
    first <- strsplit(x, "")[[1]][clustering_position]
    last <- strsplit(x, "")[[1]][reference_position]
    ifelse(first == last, "Not changed", "Changed")
  })

  clustering %>% unname()

}

change_frequency <- function(variable){

### making different groups based on reference year and a clustering year

  variable_char <- as.character(variable)

  clustering <- sapply(variable, function(x){

    chars <- str_split(x, "")[[1]]
    unique_char <- length(unique(
      chars[clustering_position:reference_position]))
    if(unique_char == 1){
      value <- "Stable (No Change)"
    }else if(unique_char == 2 & chars[clustering_position] ==
              chars[reference_position]){
      value <- "Reverted (1 LU in between)"
    }else if(unique_char == 2 & chars[clustering_position] !=
             chars[reference_position]){
      value <- "Shifted (1 Change)"
    } else if(unique_char == 3 & chars[clustering_position] ==
              chars[reference_position]){
      value <- "Reverted (2 LUs in between)"
    } else if(unique_char == 3 & chars[clustering_position] !=
              chars[reference_position]){
      value <- "Shifted (2 LU Changes)"
    } else if(unique_char == length(chars)){
      value <- "Fully Dynamic (All Different)"
    } else{
      value <- "Multiple changes"
    }

    value

  })

}

presence <- function(variable){


  variable_char <- as.character(variable)

  clustering <- sapply(variable, function(x){

    chars <- str_split(x, "")[[1]]
    chars <- chars[clustering_position:reference_position]
    n <- length(chars)
    count <- sum(LU == chars)
    upper <- ceiling(n/2)

    if(all(LU == chars)){
      value <- "Stable presence"

    }else if(all(LU != chars)){
      value <- "Stable absence"

    }else if(count <= upper & chars[n] == LU){
      value <- "Stable loss"

    }else if(count >= upper & chars[n] == LU){
      value <- "Eventual loss"

    }else if(count <= upper & chars[1] == LU){
      value <- "Eventual gain"

    }else if(count >= upper & chars[1] == LU){
      value <- "Stable gain"

    }else {
      value <- "Intermittent Presence"
      }


  })

}


transition_present <- function(variable) {

  variable <- as.character(variable)
  transition <- as.character(transition)

  first_char_transition <- substr(transition, 1, 1)

  sapply(variable, function(x) {

    # Subsetten van periode obv gegeven referentie en cluster jaar
    x <- str_sub(x, start = (clustering_position-1), end = reference_position)

    # Aantal keer dat de transitie voorkomt (inclusief overlappende matches)
    count <- str_count(x, paste0("(?=(", transition, "))"))

    # Posities van de transitie in x
    transition_positions <- str_locate_all(x,
                                           paste0("(?=(", transition, "))")
                                           )[[1]][,1]

    if (length(transition_positions) == 0) {
      last_transition_position <- NA
      first_transition_position <- NA
    } else {
      last_transition_position <- max(transition_positions)  # Laatste keer transitie start
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
                       reference_year = 1774){

  ### check properties of the clustering method

  # assert_that(length(clustering_method) == 1,
  #             msg = "Please provide 1 clustering method")

  assert_that(all(clustering_method %in% c("mask_follow_up",
                                           "mask_follow_up_intensity",
                                           "make_bin_change",
                                           "change_frequency",
                                           "transition",
                                           "presence")),
              msg = "clustering_method should be either
              mask_follow_up, mask_follow_up_intensity, make_bin_change,
              change_frequency, transition, presence")

  ### Check properties of coding LU

  assert_that(is.numeric(c(akker, bebouwing, bos, grasland,
                           heide, strand, moeras, boomgaard,
                           water)),
              msg = "The provided LU should be in numerical form")

  if(akker == 1 & bebouwing == 2 & bos == 3 & grasland == 4 &
     heide == 5 & strand == 6 & moeras == 7 & boomgaard == 8 &
     water == 9){

    message("Consistent coding of LU was used")

  } else{
    warning("There is a change in LU coding")
  }

  ### Check properties of supplied years

  assert_that(identical(years, sort(years, decreasing = TRUE)),
              msg = "Supplied years are not in descending order")

  assert_that(length(years) == 4,
              msg = "The number of provided years is different from 4")

  assert_that(all(years %in%  c(2022, 1969, 1873, 1774)),
              msg = "One ore more of the provided year(s) is/are not
              2022, 1969, 1873 or 1774")

  assert_that(is.numeric(years),
              msg = "The supplied years are not in numerical form")

  assert_that(is.numeric(reference_year) & is.numeric(clustering_year),
              msg = "Reference year and clustering year should be numeric")

  assert_that(clustering_year %in% years,
              msg = "The clustering year is not in the provided years")

  assert_that(length(clustering_year) == 1,
              msg = "Provide 1 clusterig year")

  assert_that(reference_year %in% years,
              msg = "The reference year is not in the provided years")

  assert_that(length(reference_year) == 1,
              msg = "Provide 1 reference year")

  assert_that(reference_year != clustering_year,
              msg = "Clustering year and reference year are the same")

  assert_that(clustering_year > reference_year,
              msg = "Clustering year is prior to reference year")

  ### Checking properties of the supplied data

  assert_that(length(prediction) == length(reference),
              msg = "Prediction and reference are not of the same length")

  assert_that(all(levels(prediction) == levels(reference)),
              msg = "Prediction and reference have not the same levels")

  assert_that(all(nchar(as.character(prediction)) ==
                    nchar(as.character(prediction[1]))),
              msg = "Some predictions do not have the same length")

  assert_that(all(nchar(as.character(reference)) ==
                    nchar(as.character(reference[1]))),
              msg = "Some references do not have the same length")

  assert_that(nchar(as.character(reference[1])) ==
                nchar(as.character(prediction[1])) &
                nchar(as.character(reference[1])) == length(years),
              msg = "The reference and prediction differ in length of characters")

  ### Check properties of supplied LU

  assert_that(is.character(LU) & nchar(LU) == 1,
              msg = "Provided LU not in right format")

  ### Check properties of transition

  assert_that(is.character(transition) & nchar(transition) > 1,
              msg = "Provided transition is in wrong format")


  ### Converting supplied years to the position in the given years

  clustering_position <- which(clustering_year == years)[[1]]
  reference_position <- which(reference_year == years)[[1]]

  assert_that(clustering_position < reference_position,
              msg = "Clustering year is prior to the reference year")

  tempdata <- data.frame(observed = prediction,
                         reference = reference)


  for (method in clustering_method) {
    if (method == "mask_follow_up") {
      tempdata[[paste0("prediction_", method)]] <- mask_follow_up(prediction)
      tempdata[[paste0("reference_", method)]] <- mask_follow_up(reference)

      ## set them at the same level

      all_levels <- lubridate::union(levels(factor(
                                       tempdata[[paste0("prediction_", method)]])),
                                     levels(factor(
                                       tempdata[[paste0("reference_", method)]])))

      tempdata[[paste0("reference_", method)]] <- factor(tempdata[[paste0("reference_", method)]],
                                                         levels = all_levels)

      tempdata[[paste0("prediction_", method)]] <- factor(tempdata[[paste0("prediction_", method)]],
                                                         levels = all_levels)

    }

    if (method == "mask_follow_up_intensity") {
      tempdata[[paste0("prediction_", method)]] <- mask_follow_up_intensity(prediction)
      tempdata[[paste0("reference_", method)]] <- mask_follow_up_intensity(reference)

      ## set them at the same level

      all_levels <- lubridate::union(levels(factor(
        tempdata[[paste0("prediction_", method)]])),
        levels(factor(
          tempdata[[paste0("reference_", method)]])))

      tempdata[[paste0("reference_", method)]] <- factor(tempdata[[paste0("reference_", method)]],
                                                         levels = all_levels)

      tempdata[[paste0("prediction_", method)]] <- factor(tempdata[[paste0("prediction_", method)]],
                                                         levels = all_levels)

    }

    if (method == "make_bin_change") {
      tempdata[[paste0("prediction_", method)]] <- make_bin_change(prediction)
      tempdata[[paste0("reference_", method)]] <- make_bin_change(reference)

      ## set them at the same level

      all_levels <- lubridate::union(levels(factor(
        tempdata[[paste0("prediction_", method)]])),
        levels(factor(
          tempdata[[paste0("reference_", method)]])))

      tempdata[[paste0("reference_", method)]] <- factor(tempdata[[paste0("reference_", method)]],
                                                         levels = all_levels)

      tempdata[[paste0("prediction_", method)]] <- factor(tempdata[[paste0("prediction_", method)]],
                                                         levels = all_levels)
    }

    if (method == "change_frequency") {
      tempdata[[paste0("prediction_", method)]] <- change_frequency(prediction)
      tempdata[[paste0("reference_", method)]] <- change_frequency(reference)

      ## set them at the same level

      all_levels <- lubridate::union(levels(factor(
        tempdata[[paste0("prediction_", method)]])),
        levels(factor(
          tempdata[[paste0("reference_", method)]])))

      tempdata[[paste0("reference_", method)]] <- factor(tempdata[[paste0("reference_", method)]],
                                                         levels = all_levels)

      tempdata[[paste0("prediction_", method)]] <- factor(tempdata[[paste0("prediction_", method)]],
                                                         levels = all_levels)
    }

    if (method == "presence") {
      tempdata[[paste0("prediction_", method)]] <- presence(prediction)
      tempdata[[paste0("reference_", method)]] <- presence(reference)

      ## set them at the same level

      all_levels <- lubridate::union(levels(factor(
        tempdata[[paste0("prediction_", method)]])),
        levels(factor(
          tempdata[[paste0("reference_", method)]])))

      tempdata[[paste0("reference_", method)]] <- factor(tempdata[[paste0("reference_", method)]],
                                                         levels = all_levels)

      tempdata[[paste0("prediction_", method)]] <- factor(tempdata[[paste0("prediction_", method)]],
                                                         levels = all_levels)
    }

    if (method == "transition") {
      tempdata[[paste0("prediction_", method)]] <- transition_present(prediction)
      tempdata[[paste0("reference_", method)]] <- transition_present(reference)

      ## set them at the same level

      all_levels <- lubridate::union(levels(factor(
        tempdata[[paste0("prediction_", method)]])),
        levels(factor(
          tempdata[[paste0("reference_", method)]])))

      tempdata[[paste0("reference_", method)]] <- factor(tempdata[[paste0("reference_", method)]],
                                                         levels = all_levels)

      tempdata[[paste0("prediction_", method)]] <- factor(tempdata[[paste0("prediction_", method)]],
                                                         levels = all_levels)
    }
  }

  return(tempdata)
}



