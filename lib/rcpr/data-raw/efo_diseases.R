PATH_TO_OPEN_TARGETS_PARQUESTS <- "<<ENTER_PATH>>"

efo_diseases <- "<<LOAD OPENTARGETS DISEASE TABLE>>" %>%
  dplyr::select(id, name) %>%
  dplyr::as_tibble() %>%
  dplyr::arrange(name)

usethis::use_data(efo_diseases, overwrite = TRUE)
