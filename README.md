# Reconnect: physical vs conversational contacts

Analysing age- and setting-specific physical and non-physical (conversational) social contacts using data from the Reconnect social contact survey. 

## Data

[Social contact patterns in the United Kingdom following the COVID-19 pandemic: The Reconnect cross-sectional survey](https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1005038), L Goodfellow, B J Quilty, K van Zandvoort, W J Edmunds. PLOS Medicine. 2026.

[Zenodo dataset](https://doi.org/10.5281/zenodo.21218377) (10.5281/zenodo.16845074).

## Running the analysis

Run `phys_contacts.R` to produce this analysis.

Functions for analysis are defined in `functions.R`, e.g. fitting a negative binomial distribution to contact counts, enforcing reciprocity of social contacts, weighting participants, and calculating age assortativity of mixing. Some are taken directly from the Reconnect [GitHub repo](https://github.com/cmmid/reconnect_uk_social_contact_survey), e.g. `weight_participants()`.

## Limitations 

This analysis does not include large group contacts, as there is no information on the physical/non-physical nature of these contacts. Due to an error in data collection, data on physical interaction was missing for 4% of individually-reported contacts.
