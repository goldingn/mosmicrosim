# make some figures to check and display the climatic abundance outputs

library(terra)
library(tidyverse)

# load all the rasters into one object
files <- list.files("~/Dropbox/github/mosmicrosim/processing/vector_rasters/",
                    pattern = "^an_gambiae",
                    full.names = TRUE)

names <- files |>
  basename() |>
  stringr::str_remove("^an_gambiae.") |>
  stringr::str_remove(".adult.tif$")

ag <- terra::rast(files)
names(ag) <- names

# fill in bad NAs (all well outside the area of interest) with 0s
mask <- terra::rast("temp/raster_mask.tif")
ag[is.na(ag)] <- 0
ag <- mask(ag, mask)

# find years and months for plotting and summarising
slice_year <- names |>
  stringr::str_split_i("\\.", 1) |>
  as.numeric()
slice_month <- names |>
  stringr::str_split_i("\\.", 2) |>
  as.numeric()
years <- sort(unique(slice_year))
months <- sort(unique(slice_month))

# create monthly synoptic version
monthly_list <- list()
for(month in months) {
  monthly_list[[month]] <- mean(ag[[slice_month == month]])
}
synoptic <- terra::rast(monthly_list)
names(synoptic) <- month.name[months]

# create an overall average
average <- mean(synoptic)

# load some things for plotting
borders <- readRDS("~/Dropbox/github/ir_cube/data/clean/gadm_polys.RDS")
pf_water_mask <- rast("~/Dropbox/github/ir_cube/data/clean/pfpr_water_mask.tif")

# # mask out open water and arid areas?
# ag_plot <- terra::mask(ag, pf_water_mask)
ag_plot <- ag

# plot all versions with the same max
max <- max(terra::global(ag_plot, "max", na.rm = TRUE))

library(tidyterra)

theme_maps <- function() {
  theme_minimal() +
    theme(axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_blank())
}

border_col <- grey(0.4)

country_borders <- geom_sf(
  data = borders,
  col = border_col,
  linewidth = 0.1,
  fill = "transparent"
)

for (year in years) {

  print(year)

  sub <- which(slice_year == year)

  data_sub <- ag_plot[[sub]]
  names(data_sub) <- month.name[slice_month[sub]]

  year_plot <- ggplot() +
    geom_spatraster(
      data = data_sub
    ) +
    facet_wrap(
      ~lyr
    ) +
    country_borders +
    scale_fill_distiller(
      guide = "none",
      palette = "Purples",
      direction = 1,
      limits = c(0, max),
      transform = "log1p",
      na.value = "transparent"
    ) +
    theme_maps() +
    ggtitle(
      sprintf("Climatic abundance in %i", year),
      "Vector abundance per unit of larval habitat"
    )

  year_plot_fpath <- file.path(
    "figures",
    sprintf(
      "monthly_year_%s.png",
      year)
  )

  ggsave(
    filename = year_plot_fpath,
    plot = year_plot,
    bg = "white",
    width = 6.5,
    height = 6.5
  )

}


# plot synoptic monthly version
synoptic_max <- max(terra::global(synoptic, "max", na.rm = TRUE))
synoptic_plot <- ggplot() +
  geom_spatraster(
    data = synoptic
  ) +
  facet_wrap(
    ~lyr
  ) +
  country_borders +
  scale_fill_distiller(
    guide = "none",
    palette = "Purples",
    direction = 1,
    limits = c(0, synoptic_max),
    transform = "log1p",
    na.value = "transparent"
  ) +
  theme_maps() +
  ggtitle(
    sprintf("Average monthly climatic abundance 2000-2024"),
    "Vector abundance per unit of larval habitat"
  )

ggsave(
  filename = file.path(
    "figures",
    "synoptic.png"
  ),
  plot = synoptic_plot,
  bg = "white",
  width = 6.5,
  height = 6.5
)



# plot overall average version
average_max <- max(terra::global(average, "max", na.rm = TRUE))
average_plot <- ggplot() +
  geom_spatraster(
    data = average
  ) +
  country_borders +
  scale_fill_distiller(
    guide = "none",
    palette = "Purples",
    direction = 1,
    limits = c(0, average_max),
    transform = "log1p",
    na.value = "transparent"
  ) +
  theme_maps() +
  ggtitle(
    sprintf("Average climatic abundance 2000-2024"),
    "Vector abundance per unit of larval habitat"
  )

ggsave(
  filename = file.path(
    "figures",
    "average.png"
  ),
  plot = average_plot,
  bg = "white",
  width = 6.5,
  height = 6.5
)

# plot overall average version

# convert to probability of detection scale, under an arbitrary surveillance
# effort that looks nice
average_suitability <- 1 - exp(-5 * average)
average_suitability_max <- max(terra::global(average_suitability, "max", na.rm = TRUE))
average_suitability_plot <- ggplot() +
  geom_spatraster(
    data = average_suitability
  ) +
  country_borders +
  scale_fill_distiller(
    guide = "none",
    palette = "Purples",
    direction = 1,
    limits = c(0, average_suitability_max),
    transform = "log1p",
    na.value = "transparent"
  ) +
  theme_maps() +
  ggtitle(
    sprintf("Average climatic suitability 2000-2024"),
    "Relative probability of detection per unit of larval habitat"
  )

ggsave(
  filename = file.path(
    "figures",
    "average_suitability.png"
  ),
  plot = average_suitability_plot,
  bg = "white",
  width = 6.5,
  height = 6.5
)


# animate a gif of all the slice data

ag_animate <- ag
animate_names <- paste(slice_year, month.name[slice_month])
names(ag_animate) <- animate_names

library(gganimate)
slice_max <- max(terra::global(ag_animate, "max", na.rm = TRUE))

slice_plot <- ggplot() +
  geom_spatraster(
    data = ag_animate[[1:3]]
  ) +
  facet_wrap(
    ~lyr
  ) +
  country_borders +
  scale_fill_distiller(
    guide = "none",
    palette = "Purples",
    direction = 1,
    limits = c(0, slice_max),
    transform = "log1p",
    na.value = "transparent"
  ) +
  theme_maps()

slice_plot_animated <- slice_plot +
  transition_states(
    lyr,
    transition_length = 1,
    state_length = 0
  )

anim <- gganimate::animate(slice_plot_animated, fps = 2)
gganimate::anim_save("~/Desktop/mosmicrosim_animation.gif",
                     anim)

magick::image_write(anim,
                    path = "~/Desktop/mosmicrosim_animation.gif")




# create monthly anomaly maps

idx_2020 <- which(slice_year == 2020)
anomaly_2020 <- synoptic - ag[[idx_2020]]


names(anomaly_2020) <- month.name[slice_month[idx_2020]]

anomaly_max <- max(terra::global(anomaly_2020, "max", na.rm = TRUE))
anomaly_min <- min(terra::global(anomaly_2020, "min", na.rm = TRUE))
anomaly_extreme <- max(abs(anomaly_min), abs(anomaly_max))
anomaly_range <- c(-anomaly_extreme, anomaly_extreme)

anomaly_2020_plot <- ggplot() +
  geom_spatraster(
    data = anomaly_2020
  ) +
  facet_wrap(
    ~lyr
  ) +
  country_borders +
  scale_fill_distiller(
    guide = "none",
    palette = "RdBu",
    direction = 1,
    limits = anomaly_range,
    na.value = "transparent"
  ) +
  theme_maps() +
  ggtitle(
    "Climatic abundance anomaly in 2020",
    "Difference from 2000-2024 mean vector abundance per unit of larval habitat"
  )

anomaly_2020_plot_fpath <- file.path(
  "figures",
  "anomaly_2020.png"
)

ggsave(
  filename = anomaly_2020_plot_fpath,
  plot = anomaly_2020_plot,
  bg = "white",
  width = 6.5,
  height = 6.5
)


# zoom in on anomaly in May 2020

anomaly_may_2020_plot <- ggplot() +
  geom_spatraster(
    data = anomaly_2020[["May"]]
  ) +
  country_borders +
  scale_fill_distiller(
    guide = "none",
    palette = "RdBu",
    direction = 1,
    limits = anomaly_range,
    na.value = "transparent"
  ) +
  coord_sf(
    xlim = c(-1, 50),
    ylim = c(-17, 15)
  ) +
  theme_maps()

anomaly_may_2020_plot_fpath <- file.path(
  "figures",
  "anomaly_may_2020.png"
)

ggsave(
  filename = anomaly_may_2020_plot_fpath,
  plot = anomaly_may_2020_plot,
  bg = "white",
  width = 6.5,
  height = 6.5
)
