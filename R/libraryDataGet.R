#' Retrieve a reference dataset from ExperimentHub
#'
#' @description
#' Retrieves one ExperimentHub resource matching the supplied title.
#'
#' @param title A single, non-empty character string identifying the
#'   ExperimentHub resource, for example `"FlowSorted.Blood.EPIC"`.
#'
#' @return The ExperimentHub resource matching `title`.
#'
#' @examples
#' FlowSorted.Blood.EPIC <- libraryDataGet("FlowSorted.Blood.EPIC")
#' FlowSorted.Blood.EPIC
#'
#' @export
libraryDataGet <- function(title) {
  if (
    !is.character(title) ||
    length(title) != 1L ||
    is.na(title) ||
    !nzchar(title)
  ) {
    stop(
      "`title` must be a single, non-empty character string.",
      call. = FALSE
    )
  }

  hub <- ExperimentHub::ExperimentHub()
  matches <- AnnotationHub::query(hub, title)
  resource_ids <- unique(matches$ah_id)

  if (length(resource_ids) == 0L) {
    stop(
      sprintf(
        "No ExperimentHub resource matched title '%s'.",
        title
      ),
      call. = FALSE
    )
  }

  if (length(resource_ids) > 1L) {
    stop(
      sprintf(
        paste0(
          "Multiple ExperimentHub resources matched title '%s': ",
          "%s."
        ),
        title,
        paste(resource_ids, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  hub[[resource_ids[[1L]]]]
}
