
using DataFrames
using CSV
import HTTP.request
using Pipe
using JSON
import NamedTupleTools.namedtuple


function import_hopkins_data(filepath, value_col_name)

    function aggregate_across_regions(df)
        colnames = (Symbol(value_col_name), :date)
        colvalues = (sum(df[:, Symbol(value_col_name)]), df.date[1])
        res = namedtuple(colnames, colvalues)
        return res
    end

    dat = @pipe (DataFrame(CSV.File(filepath)) |>
        DataFrames.stack(_, 5:size(_)[2]) |>
        rename(_, Dict( "Country/Region" => "Country_Region",
            "Province/State" => "Province_State",
            "variable" => "date",
            "value" => value_col_name)) |>
            groupby(_, [:Province_State, :Country_Region, :date]) |>
            map(aggregate_across_regions, _) |>
            DataFrame(_))
    return dat

end

function get_region_pop(region)
    region_param = replace(region, " " => "%20")
    url = "https://query.wikidata.org/sparql?format=json&query=SELECT%20%3Fpopulation%20WHERE%20%7B%0A%20%20SERVICE%20wikibase%3Amwapi%20%7B%0A%20%20%20%20%20%20bd%3AserviceParam%20mwapi%3Asearch%20%22$region_param%22%20.%20%20%20%20%0A%20%20%20%20%20%20bd%3AserviceParam%20mwapi%3Alanguage%20%22en%22%20.%20%20%20%20%0A%20%20%20%20%20%20bd%3AserviceParam%20wikibase%3Aapi%20%22EntitySearch%22%20.%0A%20%20%20%20%20%20bd%3AserviceParam%20wikibase%3Aendpoint%20%22www.wikidata.org%22%20.%0A%20%20%20%20%20%20bd%3AserviceParam%20wikibase%3Alimit%201%20.%0A%20%20%20%20%20%20%3Fitem%20wikibase%3AapiOutputItem%20mwapi%3Aitem%20.%0A%20%20%7D%0A%20%20%3Fitem%20wdt%3AP1082%20%3Fpopulation%0A%7D"
    r = request("GET", url)
    result = missing
    try
        results = JSON.parse((String(r.body)))
        result = results["results"]["bindings"][1]["population"]["value"]
    catch e
        println("No population found for $region")
    end
    return result
end

dat = import_hopkins_data("data/raw/time_series_covid19_confirmed_global.csv", "confirmed_cases")
# regions = [("region" => x[1] + " " + x[2], population => get_region_pop()) for province_state in unique(dat.Province_State)]

regions_dat = @pipe (dat |>
    groupby(_, [:Province_State, :Country_Region]) |>
    map(df -> (df.Country_Region[1], df.Province_State[1]), _) |>
    DataFrame(_))


population_data = DataFrame(Province_State = Union{Missing, String}[], Country_Region = String[], population = Union{Missing, String}[])


for i in 1:size(regions_dat, 1)
    lookup_name = ismissing(regions_dat[i, :Province_State]) ? regions_dat[i, :Country_Region] : regions_dat[i, :Province_State]
    push!(population_data, (regions_dat[i, :Province_State], regions_dat[i, :Country_Region], get_region_pop(lookup_name)))
    println(population_data[i, :])
end

CSV.write("data/raw/wikidata_region_populations.csv", population_data)