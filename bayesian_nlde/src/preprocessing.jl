

function import_edcd_data(filepath = "", url = "https://opendata.ecdc.europa.eu/covid19/casedistribution/csv")
    if filepath == "" && startswith(url, "https://")
        filepath = string(@__DIR__) * "/data/raw/latest_edcd_covid.csv" 
        download(url, filepath)
    end

    raw_dat = DataFrame(CSV.File(filepath))

    raw_dat.dateRep = Date.(raw_dat[:, :dateRep], "dd/mm/yyyy")
    sort!(raw_dat, (:countriesAndTerritories, :dateRep));
    filter!(row -> !ismissing(row.popData2018), raw_dat);

    function add_cumsum_time(df)
        cum_case_prop = cumsum(df.cases ./ df.popData2018)
        t = df.dateRep - minimum(df.dateRep)
        t = map(x -> convert(Float64, x.value), t)
        return (cum_case_prop = cum_case_prop, time_since_first_case = t, popData2018 = df.popData2018)
    end

    dat = @pipe (groupby(raw_dat, :countriesAndTerritories) |>
        map(add_cumsum_time, _) |>
        DataFrame(_));
    filter!(row -> row.cum_case_prop > 0.0, dat);

    return dat
end

function import_hopkins_data(filepath, value_col_name)

    function aggregate_across_regions(df)
        colnames = (Symbol(value_col_name), :date)
        colvalues = (sum(df[:, Symbol(value_col_name)]), df.date[1])
        res = namedtuple(colnames, colvalues)
        return res
    end

    dat = @pipe (DataFrame(CSV.File(filepath)) |>
        stack(_, 5:size(_)[2]) |>
        rename(_, Dict( "Country/Region" => "Country_Region",
            "Province/State" => "Province_State",
            "variable" => "date",
            "value" => value_col_name)) |>
            groupby(_, [:Country_Region, :date]) |>
            map(aggregate_across_regions, _) |>
            DataFrame(_))
    return dat

end

function import_data(drop_threshold = 100; download_data = false)

    confirmed_fn = string(@__DIR__) * "/data/raw/time_series_covid19_confirmed_global.csv"
    deaths_fn = string(@__DIR__) * "/data/raw/time_series_covid19_deaths_global.csv"
    recovered_fn = string(@__DIR__) * "/data/raw/time_series_covid19_recovered_global.csv"

    if download_data
        download("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv",
            confirmed_fn)
        download("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv",
            deaths_fn)
        download("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_recovered_global.csv",
            recovered_fn)
    end

    hopkins_case_dat = import_hopkins_data(confirmed_fn, "confirmed_cases")
    hopkins_death_dat = import_hopkins_data(deaths_fn, "deaths")
    hopkins_recov_dat = import_hopkins_data(recovered_fn, "recovered")

    hopkins_dat = @pipe (join(hopkins_case_dat, hopkins_recov_dat, on = [:date, :Country_Region]) |>
        join(_, hopkins_death_dat, on = [:date, :Country_Region]))

    edcd_dat = @pipe(import_edcd_data() |>
        select(_, :countriesAndTerritories, :popData2018) |>
        groupby(_, :countriesAndTerritories) |>
        map(df -> (countriesAndTerritories = df.countriesAndTerritories[1], popData2018 = df.popData2018[1]), _) |>
        DataFrame(_))

    replace!(edcd_dat.countriesAndTerritories, "United_States_of_America" => "US")

    hopkins_ecdc = join(hopkins_dat, edcd_dat, on = [:Country_Region => :countriesAndTerritories])

    # drop observations bewlow 100 cases
    filter!(row -> row.confirmed_cases > drop_threshold, hopkins_ecdc);

    function add_SIR_time(df)
        case_prop = df.confirmed_cases ./ df.popData2018
        cum_case_prop = cumsum(df.confirmed_cases ./ df.popData2018)
        death_prop = df.deaths ./ df.popData2018
        recov_prop = df.recovered ./ df.popData2018
        susc_prop = ((df.popData2018 - df.confirmed_cases - df.deaths - df.recovered) ./ df.popData2018)

        dates = Date.(string.(df.date), "mm/dd/yy") + Dates.Year("2000")
        t = dates - minimum(dates)
        res = (time_since_case_thresh = map(x -> convert(Float64, x.value), t),
            case_prop = case_prop,
            cum_case_prop = cum_case_prop,
            death_prop = death_prop,
            recov_prop = recov_prop,
            susc_prop = susc_prop,
            pop = df.popData2018)
        return res
    end

    dat = @pipe (groupby(hopkins_ecdc, :Country_Region) |>
        map(add_SIR_time, _) |>
        DataFrame(_));

    return dat
end
