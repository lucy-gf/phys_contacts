# Reconnect: physical vs conversational contacts

Analysing age- and setting-specific physical and non-physical (conversational) social contacts using data from the Reconnect social contact survey. 

Run `phys_contacts.R` to produce this analysis.

Functions for analysis are defined in `functions.R`, such as for fitting a negative binomial distribution to contact counts, enforcing reciprocity of social contacts, weighting participants, and calculating age assortativity of mixing.

**Limitations:** This analysis does not include large group contacts, as there is no information on the physical/non-physical nature of these contacts. Due to an error in data collection, data on physical interaction was missing for 4% of individually-reported contacts.
