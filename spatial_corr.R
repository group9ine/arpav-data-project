library(tidyverse)
library(sf)
theme_set(theme_minimal(base_size = 16, base_family = "Roboto Condensed"))

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
    SELECT DEN_PROV + DEN_CM AS province
    FROM %s
    WHERE DEN_PROV IN ('Belluno', 'Padova', 'Rovigo',
                       'Treviso', 'Verona', 'Vicenza')
       OR DEN_CM = 'Venezia'",
    st_layers(maps_folder)$name[1]
)
veneto_sf <- st_read(maps_folder, query = select_veneto) |>
    mutate(province = factor(str_replace(province, "-", "")))

ggplot() +
    geom_sf(aes(fill = province), data = veneto_sf) +
    geom_sf(data = stations_sf) +
    coord_sf() +
    labs(x = "Longitude", y = "Latitude", fill = "Province")

full_data <- enframe(raw_data, name = "station", value = "dump") |>
    unnest(dump) |>
    group_by(station) |>
    group_modify(
        function(x, y) {
            mins <- x$dump[seq(
                grep("^Parametro Temperatura.*minime", x$dump) + 5,
                grep("valore medio delle minime", x$dump) - 3
            )]
            mins <- read.csv(
                text = mins,
                sep = ";", dec = ".",
                header = FALSE, na.strings = ">>",
                col.names = c("year", month.abb, "yavg")
            ) |>
                select(-yavg) |>
                pivot_longer(-year, names_to = "month", values_to = "min")

            avgs <- x$dump[seq(
                grep("^Parametro Temperatura.*medie", x$dump) + 5,
                grep("valore medio delle medie", x$dump) - 3
            )]
            avgs <- read.csv(
                text = avgs,
                sep = ";", dec = ".",
                header = FALSE, na.strings = ">>",
                col.names = c("year", month.abb, "yavg")
            ) |>
                select(-yavg) |>
                pivot_longer(-year, names_to = "month", values_to = "avg")

            maxs <- x$dump[seq(
                grep("^Parametro Temperatura.*massime", x$dump) + 5,
                grep("valore medio delle massime", x$dump) - 3
            )]
            maxs <- read.csv(
                text = maxs,
                sep = ";", dec = ".",
                header = FALSE, na.strings = ">>",
                col.names = c("year", month.abb, "yavg")
            ) |>
                select(-yavg) |>
                pivot_longer(-year, names_to = "month", values_to = "max")

            right_join(mins, avgs, by = c("year", "month")) |>
                left_join(maxs, by = c("year", "month"))
        }
    ) |>
    drop_na()

station_corr <- function(stat) {
    require(dplyr)
    require(purrr)

    # keep only full years
    full_data <- full_data |>
        group_by(.data$station, .data$year) |>
        filter(n() == 12) |>
        ungroup()

    stat_data <- full_data |>
        filter(.data$station == stat) |>
        group_by(.data$year) |>
        summarize(
            min = mean(.data$min),
            avg = mean(.data$avg),
            max = mean(.data$max)
        )

    # get years when the chosen station was always operational
    stat_years <- full_data |>
        filter(.data$station == stat) |>
        distinct(.data$year) |>
        unlist()

    # compute correlation with compatible stations
    compatibles <- full_data |>
        group_by(.data$station) |>
        distinct(.data$year) |>
        # operational for the same number of years
        filter(n() == length(stat_years)) |>
        # and the years should match too
        filter(.data$year == stat_years) |>
        distinct(.data$station) |>
        unlist()

    map_dfr(
        compatibles,
        function(s) {
            s_data <- full_data |>
                filter(.data$station == s) |>
                group_by(.data$year) |>
                summarize(
                    min = mean(.data$min),
                    avg = mean(.data$avg),
                    max = mean(.data$max)
                )

            tibble(
                min = cor(stat_data$min, s_data$min),
                avg = cor(stat_data$avg, s_data$avg),
                max = cor(stat_data$max, s_data$max)
            )
        }
    ) |>
        mutate(station = compatibles, .before = 1)
}