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
veneto_sf <- st_transform(veneto_sf, 4326)

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
        ) |>
        # compute differences between years
        mutate(across(min:max, ~ .x - lag(.x))) |>
        na.omit()

    stat_years <- full_data |>
        filter(.data$station == stat) |>
        pull(.data$year) |>
        unique()

    compatibles <- full_data |>
        # filter out selected station
        filter(.data$station != stat) |>
        group_by(.data$station) |>
        distinct(.data$year) |>
        # operational for the same number of years
        filter(n() == length(stat_years)) |>
        # and the years should match too
        filter(.data$year == stat_years) |>
        pull(.data$station) |>
        unique()

    # compute correlation with compatible stations
    map_dfr(
        compatibles,
        function(s) {
            s_data <- full_data |>
                filter(.data$station == s) |>
                # compute yearly data
                group_by(.data$year) |>
                summarize(
                    min = mean(.data$min),
                    avg = mean(.data$avg),
                    max = mean(.data$max)
                ) |>
                # compute differences between years
                mutate(across(min:max, ~ .x - lag(.x))) |>
                na.omit()

            tibble(
                min = cor(stat_data$min, s_data$min),
                avg = cor(stat_data$avg, s_data$avg),
                max = cor(stat_data$max, s_data$max)
            )
        }
    ) |>
        # add station coordinates and convert to sf object
        mutate(
            station = compatibles,
            coords = stations_sf |>
                filter(.data$station %in% compatibles) |>
                pull(.data$geometry),
            .before = 1
        ) |>
        st_as_sf(sf_column_name = "coords")
}

ggplot() +
    geom_sf(data = veneto_sf) +
    geom_sf(
        aes(colour = avg),
        data = station_corr("Malo") |> select(avg, coords),
        size = 5
    ) +
    scale_colour_viridis_c() +
    coord_sf() +
    labs(x = "Longitude", y = "Latitude")

# IDW
library(gstat)

corr_sf <- station_corr("Malo") |> select(avg, coords)

veneto_grid <- st_make_grid(veneto_sf, n = 100, what = "corners")
veneto_grid <- veneto_grid[veneto_sf]

corr_idw <- idw(
    formula = avg ~ 1,
    locations = corr_sf,
    newdata = veneto_grid,
    idp = 2
)

corr_idw |>
    mutate(
        lon = unlist(map(geometry, 1)),
        lat = unlist(map(geometry, 2))
    ) |>
    ggplot() +
        geom_raster(aes(x = lon, y = lat, fill = var1.pred)) +
        geom_sf(data = veneto_sf, fill = NA, colour = "black") +
        geom_sf_text(
            aes(label = station),
            data = stations_sf |>
                filter(station == "Malo"),
            colour = "firebrick",
            nudge_x = 0.125,
            nudge_y = 0.03
        ) +
        geom_sf(
            data = stations_sf |>
                filter(station == "Malo"),
            colour = "firebrick",
            size = 3
        ) +
        geom_sf(
            data = stations_sf |>
                filter(station != "Malo"),
            size = 0.5
        ) +
        scale_fill_viridis_c() +
        coord_sf() +
        labs(
            x = "Longitude", y = "Latitude",
            fill = "Correlation",
            title = paste(
                "Correlation with the temperature records in",
                "Malo", "across Veneto"
            )
        )

yearly_fits <- full_data |>
    group_by(station, year) |>
    # filter only full years
    filter(n() == 12) |>
    summarize(min = mean(min), avg = mean(avg), max = mean(max)) |>
    # keep only stations that were always operational
    filter(n() == length(1994:2022)) %>%
    # perform the linear fits
    do(
        min = lm(min ~ year, data = .) |>
            summary() |>
            coefficients() %>%
            .[2, 1],
        avg = lm(avg ~ year, data = .) |>
            summary() |>
            coefficients() %>%
            .[2, 1],
        max = lm(max ~ year, data = .) |>
            summary() |>
            coefficients() %>%
            .[2, 1]
    ) |>
    unnest(min:max)

# add station coordinates
yearly_fits <- yearly_fits |>
    mutate(
        coords = stations_sf |>
            filter(station %in% unique(yearly_fits$station)) |>
            pull(geometry),
        .after = 1
    ) |>
    st_as_sf(sf_column_name = "coords")

idw(
    formula = avg ~ 1,
    locations = yearly_fits,
    newdata = veneto_grid,
    idp = 2
) |>
    mutate(
        lon = unlist(map(geometry, 1)),
        lat = unlist(map(geometry, 2))
    ) |>
    ggplot() +
        geom_raster(aes(x = lon, y = lat, fill = var1.pred)) +
        geom_sf(data = veneto_sf, fill = NA, colour = "black") +
        geom_sf(data = stations_sf, size = 0.5) +
        scale_fill_viridis_c() +
        coord_sf() +
        labs(
            x = "Longitude", y = "Latitude",
            fill = "Increase (°C/year)",
            title = paste(
                "Yearly increase in the average temperatures",
                "across Veneto"
            ),
            subtitle = "in the 1994–2022 period"
        )

# KRIGING
library(gstat)
library(automap)

corr_sf <- station_corr("Auronzo") |> select(avg, coords)

# fit the variogram (variance over distance)
corr_vgm <- variogram(avg ~ 1, corr_sf, cutoff = 120, width = 4)
corr_fit <- fit.variogram(
    corr_vgm,
    model = vgm(psill = 0.005, "Lin", nugget = 0.002, range = 0),
    fit.method = 6
)
plot(corr_vgm, corr_fit)

veneto_grid <- st_make_grid(veneto_sf, n = 40, what = "corners")
veneto_grid <- veneto_grid[veneto_sf]

corr_krig <- krige(
    formula = avg ~ 1,
    locations = corr_sf,
    newdata = veneto_grid,
    model = corr_fit
)

tmp_sf <- st_set_precision(corr_krig, precision = 10^3)
st_write(tmp_sf, "./data/tmp.shp", append = FALSE)
corr_krig <- st_read("./data/tmp.shp")

corr_krig |>
    mutate(
        lon = unlist(map(corr_krig$geometry, 1)),
        lat = unlist(map(corr_krig$geometry, 2))
    ) |>
    ggplot() +
        geom_raster(aes(x = lon, y = lat, fill = var1_pred)) +
        coord_sf()
