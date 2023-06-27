library(tidyverse)
library(sf)

data_path <- "./data/monthly/"
file_names <- list.files(data_path)

raw_data <- file_names |>
    map(~ readLines(paste0(data_path, .x), encoding = "latin1"))

# get a vector of correctly formatted station names
station_names <- file_names |>
    str_replace(".csv$", "") |>
    str_split("_") |>
    map(~ paste(.x, collapse = " ")) |>
    unlist()
names(raw_data) <- station_names

# get station coordinates from each file
stations_sf <- raw_data |>
    map(
        function(file) {
            x_line <- grep("Coordinata X", file)[1]

            # get Gauss-Boaga coordinates
            gb_coords <- file[c(x_line, x_line + 1)] |>
                str_extract("[[:digit:]]+") |>
                as.integer() |>
                t()

            # convert to WGS84
            sf_project(
                from = st_crs(3003),
                to = st_crs(4326),
                pts = gb_coords
            ) |> as.vector()
        }
    ) |>
    enframe() |>
    unnest_wider(value, names_sep = "_") |>
    set_names("station", "lon", "lat") |>
    st_as_sf(coords = c("lon", "lat"), crs = st_crs(4326))

maps_folder <- "./data/maps/province"
select_veneto <- sprintf("
    SELECT shape_area
    FROM %s
    WHERE DEN_PROV IN ('Belluno', 'Padova', 'Rovigo',
                       'Treviso', 'Verona', 'Vicenza')
       OR DEN_CM = 'Venezia'",
    st_layers(maps_folder)$name[1]
)
veneto_sf <- st_read(maps_folder, query = select_veneto)

ggplot() +
    geom_sf(data = veneto_sf) +
    geom_sf(data = stations_sf) +
    coord_sf()
