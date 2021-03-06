# Metadata 
This file specifies metadata for generated outputs found in `./generated_data`. Metadata follows the Frictionless Standards template in https://specs.frictionlessdata.io/table-schema/ 

## generated_data/best_seasonality_estimates_schema.json

 ```json 
{
  "fields": [
    {
      "name": ["country"],
      "type": ["string"],
      "format": ["default"],
      "title": ["Country ISOcode"],
      "description": ["Country ISO 3 code"]
    },
    {
      "name": ["gid"],
      "type": ["string"],
      "format": ["default"],
      "title": ["GADM v3.6 id"],
      "description": ["Unique id in the GADM dataset"]
    },
    {
      "name": ["name"],
      "type": ["string"],
      "format": ["default"],
      "title": ["GADM name"],
      "description": ["Name in the GADM dataset"]
    },
    {
      "name": ["run_level"],
      "type": ["integer"],
      "format": ["default"],
      "title": ["Run level"],
      "description": ["Administrive level at which the model was run depending on data availability (either 1 or 2)"]
    },
    {
      "name": ["model"],
      "type": ["string"],
      "format": ["default"],
      "title": ["Best model"],
      "description": ["Model variant which was retained by model comparison (one of null, base, offset, mixture or mixture_offset). See Supplementary material for descriptions."]
    },
    {
      "name": ["month"],
      "type": ["integer"],
      "format": ["default"],
      "title": ["Month"],
      "description": ["Month of the year"]
    },
    {
      "name": ["mean"],
      "type": ["number"],
      "format": ["default"],
      "title": ["Mean estimate"],
      "description": ["Posterior mean of the seasonality coefficient (odds ratio of observing excess cholera in a given month. For the null model values were set to NA. Estimates are GADM unit-specific for the mixture and mixture_offset models, and country-specific for the base and offset models."]
    },
    {
      "name": ["q025"],
      "type": ["number"],
      "format": ["default"],
      "title": ["2.5% quantile"],
      "description": ["Posterior 2.5% quantile of the seasonality coefficient"]
    },
    {
      "name": ["q975"],
      "type": ["number"],
      "format": ["default"],
      "title": ["97.5% quantile"],
      "description": ["Posterior 2.5% quantile of the seasonality coefficient"]
    }
  ],
  "missingValues": ["NA"]
} 
``` 

## generated_data/run_data/all_run_data_schema.json

```json 
{
    "fields": [
        {
            "name": "country",
            "type": "string",
            "format": "default",
            "title": "Country ISOcode",
            "description": "Country ISO 3 code"
        },
        {
            "name": "identifier",
            "type": "string",
            "format": "default",
            "title": "Cholera presence identifier",
            "description": "Threshold used to define cholera persence/absence"
        },
        {
            "name": "gid",
            "type": "string",
            "format": "default",
            "title": "GADM v3.6 id",
            "description": "Unique id in the GADM dataset"
        },
        {
            "name": "gadm_lev",
            "type": "integer",
            "format": "default",
            "title": "GADM level",
            "description": "Administrive level"
        },
        {
            "name": "year",
            "type": "integer",
            "format": "default",
            "title": "Year",
            "description": "Year of the data"
        },
        {
            "name": "month_left",
            "type": "date",
            "format": "default",
            "title": "Month left",
            "description": "Month of the start of the observations"
        },
        {
            "name": "month_right",
            "type": "date",
            "format": "default",
            "title": "Month right",
            "description": "Month of the end of the observations"
        },
        {
            "name": "n_obs",
            "type": "integer",
            "format": "default",
            "title": "Number of observations",
            "description": "Number of distintct observations covering the given gid and time span."
        },
         {
            "name": "cholera",
            "type": "integer",
            "format": "default",
            "title": "Number of observations of cholera presence",
            "description": "Number of obserations for which incidence data meets the threshold for cholera presence"
        },
    ],
    "missingValues": "NA"
}
``` 
