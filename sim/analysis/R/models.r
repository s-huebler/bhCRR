# sim/analysis/R/models.r -- the model registry
#
# OWNS
#   default_models()          the registry: one row per model id with display
#                             label, colour, role and plotting order
#   select_models(models, ids)
#                             keep only `ids`, in that order; the result will
#                             NOT pick up unregistered ids from the data
#   resolve_models(models, present)
#                             reconcile a registry with the model ids actually
#                             present in a results table; drops registry rows
#                             with no data, appends unregistered ids with a
#                             fallback colour (and a message), returns the
#                             registry in plotting order
#   label_models(df, models)  turn df$Model into an ordered factor of display
#                             labels and drop rows whose Model is not in the
#                             registry
#   scale_colour_models(models), scale_fill_models(models)
#                             ggplot2 manual scales keyed on display label
#
# WHY THIS FILE EXISTS
#   Every other analysis function takes a `models` registry argument and never
#   names a model. To rename a comparator, change its colour, drop it from the
#   report, or add a new one, edit default_models() (or pass a modified copy
#   from the .qmd). Nothing in bias.r or selection.r needs to change.
#
# REGISTRY COLUMNS
#   model   character  id exactly as it appears in the parsed Model column
#   label   character  display name used in legends and tables
#   colour  character  hex colour
#   role    character  "proposed" | "oracle" | "comparator"; free text, used
#                      only for the reader (and available for filtering)

default_models <- function() {
  data.frame(
    model  = c("bhcrr",    "fg_oracle", "lasso",   "alasso",  "scad",    "mcp"),
    label  = c("bhCRR",    "Oracle FG", "LASSO",   "aLASSO",  "SCAD",    "MCP"),
    colour = c("#E31A1C",  "#4D4D4D",   "#1F78B4", "#33A02C", "#6A3D9A", "#FF7F00"),
    role   = c("proposed", "oracle",    "comparator", "comparator",
               "comparator", "comparator"),
    stringsAsFactors = FALSE
  )
}

# Colours handed out to model ids that are in the data but not the registry.
.fallback_colours <- c("#A6761D", "#E7298A", "#66A61E", "#1B9E77", "#D95F02",
                       "#7570B3", "#E6AB02", "#666666")

resolve_models <- function(models, present) {
  stopifnot(is.data.frame(models),
            all(c("model", "label", "colour") %in% names(models)))
  if (anyDuplicated(models$model)) stop("Duplicate model ids in registry.")
  if (anyDuplicated(models$label)) stop("Duplicate labels in registry.")
  present <- unique(as.character(present))
  add_unreg <- !identical(attr(models, "add_unregistered"), FALSE)
  if (is.null(models$role)) models$role <- NA_character_
  cols <- c("model", "label", "colour", "role")

  keep <- models[models$model %in% present, cols, drop = FALSE]

  extra <- setdiff(present, models$model)
  if (length(extra) && add_unreg) {
    # Only add unregistered ids when the registry was not deliberately
    # filtered with select_models(), which sets add_unregistered = FALSE.
    message("Models in the data but not in the registry (given fallback ",
            "colours): ", paste(extra, collapse = ", "))
    add <- data.frame(
      model  = extra,
      label  = extra,
      colour = rep_len(.fallback_colours, length(extra)),
      role   = "unregistered",
      stringsAsFactors = FALSE
    )
    keep <- rbind(keep, add[, cols])
  }
  rownames(keep) <- NULL
  keep
}

# Select a subset of the registry for a report. Unregistered ids in the data
# are then ignored rather than added back.
select_models <- function(models, ids) {
  bad <- setdiff(ids, models$model)
  if (length(bad)) stop("Not in the registry: ", paste(bad, collapse = ", "))
  out <- models[match(ids, models$model), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "add_unregistered") <- FALSE
  out
}

label_models <- function(df, models) {
  df <- df[df$Model %in% models$model, , drop = FALSE]
  df$Model <- factor(models$label[match(df$Model, models$model)],
                     levels = models$label)
  df
}

scale_colour_models <- function(models, ...) {
  ggplot2::scale_colour_manual(values = stats::setNames(models$colour,
                                                        models$label),
                               drop = TRUE, ...)
}

scale_fill_models <- function(models, ...) {
  ggplot2::scale_fill_manual(values = stats::setNames(models$colour,
                                                      models$label),
                             drop = TRUE, ...)
}
