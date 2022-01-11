#' libraryDataGet
#' @description
#' Function to load the library data from ExperimentHub
#' @import ExperimentHub
#' @importFrom AnnotationHub query
#' @param
#' title        title of the data, e.g., 'FlowSorted.Blood.EPIC'
#' @return
#' The function will look for the dataset in ExperimentHub and load the object
#' @examples
#' FlowSorted.Blood.EPIC <-
#'     libraryDataGet("FlowSorted.Blood.EPIC")
#' FlowSorted.Blood.EPIC
#' @return
#' This function will return an object matching the title of the ExperimenHub
#' @export
libraryDataGet <- function(title) {
    assign(title, ExperimentHub()[[query(
        ExperimentHub(),
        title
    )$ah_id]])
}
