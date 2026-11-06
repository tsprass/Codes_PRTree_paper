# ============================================================================
# FUNCTIONS TO EXTRACT UCI METADATA VIA API
# ============================================================================

# Check and create directory if needed
if (!dir.exists("data")) { dir.create("data") }

#-----------------------------------------------------------------------------
# Function to fetch metadata for a dataset by ID
#-----------------------------------------------------------------------------
fetch_uci_metadata <- function(id) {

  cat("  ID", id, "- Buscando...")

  # API URL (no trailing slash)
  api_url <- paste0("https://archive.ics.uci.edu/api/dataset?id=", id)

  # Make request with error handling
  response <- tryCatch({
    jsonlite::fromJSON(api_url)
  }, error = function(e) {
    cat(" Erro de conexão\n")
    return(NULL)
  })

  # Check whether the response is valid
  if (is.null(response) || response$status != 200) {
    cat(" Não encontrado\n")
    return(NULL)
  }

  # Extract data
  d <- response$data

  # Extract name (if available)
  name <- ifelse(is.null(d$name), NA, d$name)

  if (is.na(name)) {
    cat(" Sem nome\n")
    return(NULL)
  }

  cat(" Encontrado:", name, "\n")

  # Process tasks (may be a list or a vector)
  tasks <- if (!is.null(d$tasks)) {
    paste(d$tasks, collapse = ", ")
  } else {
    NA
  }

  # Process area
  area <- ifelse(is.null(d$area), NA, d$area)

  # Process characteristics
  characteristics <- if (!is.null(d$characteristics)) {
    paste(d$characteristics, collapse = ", ")
  } else {
    NA
  }

  # Process feature_types
  feature_types <- if (!is.null(d$feature_types)) {
    paste(d$feature_types, collapse = ", ")
  } else {
    NA
  }

  # Process has_missing_values
  has_missing <- ifelse(is.null(d$has_missing_values), "Unknown", d$has_missing_values)

  # Process instances and features
  instances <- ifelse(is.null(d$num_instances), NA, d$num_instances)
  features <- ifelse(is.null(d$num_features), NA, d$num_features)

  # !!! IMPORTANT: Target is now stored in target_col !!!
  target_names <- if (!is.null(d$target_col)) {
    paste(d$target_col, collapse = ", ")
  } else {
    NA
  }

  target_count <- if (!is.null(d$target_col)) {
    length(d$target_col)
  } else {
    NA
  }

  # Index/ID columns
  index_col <- if (!is.null(d$index_col)) {
    paste(d$index_col, collapse = ", ")
  } else {
    NA
  }

  # Create a single-row data frame
  df <- data.frame(
    id = id,
    name = name,
    instances = as.numeric(instances),
    features = as.numeric(features),
    tasks = tasks,
    area = area,
    characteristics = characteristics,
    feature_types = feature_types,
    has_missing = has_missing,
    target_count = target_count,
    target_names = target_names,
    index_col = index_col,
    data_url = ifelse(is.null(d$data_url), NA, d$data_url),
    stringsAsFactors = FALSE
  )

  return(df)
}

#-----------------------------------------------------------------------------
# Function to scan multiple IDs
#-----------------------------------------------------------------------------
scan_uci_datasets <- function(max_id = 700, delay = 0.2) {

  all_results <- list()
  found_count <- 0

  cat("\n🔍 Escaneando IDs de 1 até", max_id, "...\n\n")

  for (id in 1:max_id) {

    # Show progress every 50 IDs
    if (id %% 50 == 0) {
      cat("\n  Progresso: ID", id, "/", max_id, "- Encontrados:", found_count, "datasets\n")
    }

    # Fetch metadata
    result <- fetch_uci_metadata(id)

    # If found, add to the list
    if (!is.null(result)) {
      all_results[[as.character(id)]] <- result
      found_count <- found_count + 1

      # Save progress every 20 datasets found
      if (found_count %% 20 == 0) {
        partial_df <- dplyr::bind_rows(all_results)
        saveRDS(partial_df, "data/uci_metadata_partial.rds")
        write.csv(partial_df, "data/uci_metadata_partial.csv", row.names = FALSE)
        cat("    💾 Progresso salvo -", found_count, "datasets encontrados\n")
      }
    }

    # Short delay to avoid overloading the API
    Sys.sleep(delay)
  }

  cat("\n✅ Escaneamento completo!\n")
  cat("📊 Total de datasets encontrados:", found_count, "de", max_id, "\n")

  # Combine all results
  if (found_count > 0) {
    all_datasets <- dplyr::bind_rows(all_results)
    return(all_datasets)
  } else {
    return(NULL)
  }
}

#-----------------------------------------------------------------------------
# RUN THE SCAN
#-----------------------------------------------------------------------------

# Scan IDs from 1 to 700
resultados <- scan_uci_datasets(max_id = 700, delay = 0.2)

#-----------------------------------------------------------------------------
# PROCESS AND FILTER RESULTS
#-----------------------------------------------------------------------------

if (!is.null(resultados)) {

  # Clean data
  resultados <- resultados %>%
    dplyr::mutate(
      instances = as.numeric(instances),
      features = as.numeric(features),
      target_count = as.numeric(target_count),
      has_missing = dplyr::case_when(
        tolower(has_missing) %in% c("yes", "true", "y", "t") ~ "Yes",
        tolower(has_missing) %in% c("no", "false", "f", "n") ~ "No",
        TRUE ~ has_missing
      )
    )

  #-----------------------------------------------------------------------------
  # SAVE ALL RESULTS
  #-----------------------------------------------------------------------------

  # Save full results
  saveRDS(resultados, "data/uci_all_metadata.rds")
  write.csv(resultados, "data/uci_all_metadata.csv", row.names = FALSE)

  cat("\n💾 Todos os dados salvos na pasta 'data/'\n")
  cat("   - uci_all_metadata.csv (todos os datasets)\n")

} else {
  cat("\n❌ Nenhum dataset encontrado. Verifique sua conexão.\n")
}