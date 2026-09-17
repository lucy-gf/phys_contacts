## RECONNECT PHYSICAL CONTACTS ##

library(ggplot2)
library(data.table)
library(contactsurveys)
library(socialmixr)
library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(viridis)
library(patchwork)

age_limits <- seq(0, 80, 5) # for matrices 

# for age-specific contacts
age_breaks <- c(-Inf, 5*1:(75/5), Inf)
age_vals <- age_breaks[is.finite(age_breaks)]
age_labels <- c(paste0(c(0, age_vals[1:length(age_vals)-1]), '-', c(age_vals-1)), paste0(age_vals[length(age_vals)], '+'))

source('functions.R')

## load reconnect data

ZENODO_WORKING <- F # ZENODO NOT WORKING FOR SOME REASON

if(ZENODO_WORKING){
  
  c_survey <- contactsurveys::download_survey("https://doi.org/10.5281/zenodo.16845074",
                                              overwrite = T)
  reconnect <- load_survey(c_survey, participant_key = c("part_id", "cont_id"))
  
}else{
  
  # ELSE MANUAL DOWNLOAD
  rcp1 <- data.table(read_csv(file.path('data','zenodo','reconnect_participant_common.csv'),show_col_types = F))
  rcp2 <- data.table(read_csv(file.path('data','zenodo','reconnect_participant_extra.csv'),show_col_types = F))
  rcc1 <- data.table(read_csv(file.path('data','zenodo','reconnect_contact_common.csv'),show_col_types = F))
  rcc2 <- data.table(read_csv(file.path('data','zenodo','reconnect_contact_extra.csv'),show_col_types = F))
  
  reconnect <- as_contact_survey(
    x = list(
      participants = rcp1[rcp2, on = 'part_id'],
      contacts = rcc1[rcc2, on = 'cont_id']
    )
  )
  
}

## age structure

age_dat <- read_xlsx(file.path('data','mye24tablesuk.xlsx'),
                     sheet = 7, skip = 7)

age_formatted <- age_dat[1, 4:ncol(age_dat)] %>% 
  pivot_longer(!`All ages`) %>% select(!`All ages`) %>% 
  rename(age = name, population = value) %>%
  mutate(age = as.numeric(gsub('[+]', '', age)))

age_population <- group_ages(age_formatted, 
                             age_limits)

polymod_dist <- socialmixr::contact_age_distribution(polymod)

## adding weights

part_data <- reconnect$participants
colnames(part_data) <- gsub('part_', 'p_', colnames(part_data))

reconnect$participants <- reconnect$participants %>% left_join(weight_participants(part_data, 
                                                         weighting = c('p_age_group','p_gender','p_ethnicity','day_week'),
                                                         group_vars = 'p_age_group') %>% 
                                       select(p_id, post_strat_weight) %>% rename(part_id = p_id),
                                        by = 'part_id')

## adding contacts

contacts_data <- reconnect$contacts %>% 
  mutate(phys_contact = case_when(
    is.na(phys_contact) ~ 'large_group', 
    phys_contact==1 ~ 'physical',
    phys_contact==2 ~ 'conversational'
  )) %>% 
  group_by(part_id, phys_contact) %>% 
  count() %>% pivot_wider(names_from = phys_contact, values_from = n) %>% 
  mutate(large_group = case_when(is.na(large_group) ~ 0, T ~ large_group),
         physical = case_when(is.na(physical) ~ 0, T ~ physical),
         conversational = case_when(is.na(conversational) ~ 0, T ~ conversational))

contacts_data <- rbind(contacts_data, 
                       data.table(
                         part_id = reconnect$participants$part_id[reconnect$participants$part_id %notin% contacts_data$part_id],
                         large_group = 0,
                         physical = 0, 
                         conversational = 0
                       ))
           
reconnect$participants <- reconnect$participants %>% left_join(contacts_data,
                                               by = 'part_id')

#### BOOTSTRAPPING #### 

sampled_ids <- sample(reconnect$participants$part_id,
                      size = 1000*nrow(reconnect$participants),
                      replace = T,
                      prob = reconnect$participants$post_strat_weight)

bootstrapped_data <- data.table(
  bootstrap_index = rep(1:1000, each = nrow(reconnect$participants)),
  part_id = sampled_ids
) %>% left_join(reconnect$participants %>% select(part_id, part_age_group, large_group, physical, conversational),
                by = 'part_id') %>% select(!part_id)

mean_contacts <- bootstrapped_data[, lapply(.SD, neg_bin_fcn), by = c('bootstrap_index','part_age_group')]

write_csv(mean_contacts, file.path('outputs', 'mean_contacts.xlsx'))

mean_contacts$part_age_group <- factor(mean_contacts$part_age_group,
                                       levels = age_labels)

mean_contacts_summ <- mean_contacts %>% 
  pivot_longer(!c(bootstrap_index, part_age_group)) %>%
  summarise(mean = mean(value), l95 = quantile(value, 0.025), u95 = quantile(value, 0.975),
            .by = c(part_age_group, name)) 

mean_contacts_summ %>%
  filter(name != 'large_group') %>% 
  ggplot() + 
  geom_errorbar(aes(x = part_age_group, ymin = l95, ymax = u95, col = name),
                width = 0.3, position = position_dodge(width = 0.4)) +
  geom_point(aes(x = part_age_group, y = mean, col = name),
             position = position_dodge(width = 0.4), size = 2) +
  labs(x = 'Age group', y = 'Mean daily contacts', col = '') +
  scale_colour_manual(values = c('#F17105', '#1A8FE3')) +
  ylim(c(0,NA)) +
  theme_bw()
ggsave(file.path('figures','mean_contacts.png'),
       width = 10, height = 4)

mean_contacts_summ_neat <- mean_contacts_summ %>% 
  mutate(mean_contacts = paste0(round(mean,2), ' (', round(l95,2), ', ', round(u95,2), ')')) %>% 
  select(part_age_group, name, mean_contacts) %>% 
  pivot_wider(names_from = name, values_from = mean_contacts) %>% 
  arrange(part_age_group)
write_csv(mean_contacts_summ, file.path('outputs', 'mean_contacts_summ.xlsx'))

#### SETTING SPECIFIC CONTACTS ####

sett_contacts_data <- reconnect$contacts %>% 
  mutate(phys_contact = case_when(
    phys_contact==1 ~ 'physical',
    phys_contact==2 ~ 'conversational'
  )) %>% filter(!is.na(phys_contact)) %>% 
  group_by(part_id, phys_contact, cnt_location) %>% 
  count() %>% pivot_wider(names_from = phys_contact, values_from = n) %>% 
  mutate(physical = case_when(is.na(physical) ~ 0, T ~ physical),
         conversational = case_when(is.na(conversational) ~ 0, T ~ conversational)) %>% 
  group_by(part_id) %>% complete(cnt_location = c('Home', 'School', 'Work', 'Other'), fill = list(conversational = 0, physical = 0))

sett_contacts_data <- rbind(sett_contacts_data, 
                       data.table(
                         CJ(part_id = reconnect$participants$part_id[reconnect$participants$part_id %notin% sett_contacts_data$part_id],
                            cnt_location = c('Home', 'School', 'Work', 'Other')),
                         physical = 0, 
                         conversational = 0
                       ))

reconnect_part_sett <- reconnect$participants %>% 
  select(!any_of(c('conversational','physical'))) %>% 
  right_join(sett_contacts_data,
             by = 'part_id') %>% 
  pivot_wider(names_from = cnt_location, values_from = c('physical','conversational'))

bootstrapped_data <- data.table(
  bootstrap_index = rep(1:1000, each = nrow(reconnect$participants)),
  part_id = sampled_ids
) %>% left_join(reconnect_part_sett %>% select(part_id, part_age_group, contains('physical'), contains('conversational')),
                by = 'part_id') %>% select(!part_id)

mean_contacts_setting <- bootstrapped_data[, lapply(.SD, neg_bin_fcn), by = c('bootstrap_index','part_age_group')]

write_csv(mean_contacts_setting, file.path('outputs', 'mean_contacts_setting.xlsx'))

mean_contacts_setting$part_age_group <- factor(mean_contacts_setting$part_age_group,
                                               levels = age_labels)

mean_contacts_setting_summ <- mean_contacts_setting %>% 
  pivot_longer(!c(bootstrap_index, part_age_group)) %>%
  summarise(mean = mean(value), l95 = quantile(value, 0.025), u95 = quantile(value, 0.975),
            .by = c(part_age_group, name)) %>% 
  mutate(location = sub("^[^_]*_", "", name),
         name = gsub("_Home|_Other|_Work|_School", "", name))

mean_contacts_setting_summ %>%
  ggplot() + 
  geom_errorbar(aes(x = part_age_group, ymin = l95, ymax = u95, col = name),
                width = 0.3, position = position_dodge(width = 0.4)) +
  geom_point(aes(x = part_age_group, y = mean, col = name),
             position = position_dodge(width = 0.4), size = 2) +
  labs(x = 'Age group', y = 'Mean daily contacts', col = '') +
  scale_colour_manual(values = c('#F17105', '#1A8FE3')) +
  ylim(c(0,NA)) +
  facet_grid(location ~ ., scales = 'free') + 
  theme_bw()
ggsave(file.path('figures','mean_contacts_setting.png'),
       width = 10, height = 10)

mean_contacts_setting_summ_neat <- mean_contacts_setting_summ %>% 
  mutate(mean_contacts = paste0(round(mean,2), ' (', round(l95,2), ', ', round(u95,2), ')')) %>% 
  select(part_age_group, location, name, mean_contacts) %>% 
  pivot_wider(names_from = name, values_from = mean_contacts) %>% 
  arrange(part_age_group)
write_csv(mean_contacts_setting_summ_neat, file.path('outputs', 'mean_contacts_setting_summ_neat.xlsx'))

#### CONTACT MATRICES ####

## plot contact matrix 

ggplot_matrix(reconnect)

## FILTER TO INDIVIDUAL

reconnect_individual <- copy(reconnect)
reconnect_individual$contacts <- reconnect_individual$contacts[reconnect_individual$contacts$phys_contact %in% 1:2, ]

ggplot_matrix(reconnect_individual)

## FILTER TO PHYSICAL

reconnect_physical <- copy(reconnect)
reconnect_physical$contacts <- reconnect_physical$contacts[reconnect_physical$contacts$phys_contact == 1, ]

ggplot_matrix(reconnect_physical)

## FILTER TO CONVERSATIONAL

reconnect_conversational <- copy(reconnect)
reconnect_conversational$contacts <- reconnect_conversational$contacts[reconnect_conversational$contacts$phys_contact == 2, ]

ggplot_matrix(reconnect_conversational)

## PATCHWORK

ggplot_matrix(reconnect) + ggplot_matrix(reconnect_individual) +
  ggplot_matrix(reconnect_physical) + ggplot_matrix(reconnect_conversational) 

physical_plot <- ggplot_matrix(reconnect_physical) + 
  scale_fill_viridis(limits = c(0,1.3)) + ggtitle('Physical contacts') +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))
conversational_plot <- ggplot_matrix(reconnect_conversational) + 
  scale_fill_viridis(limits = c(0,1.3)) + ggtitle('Conversational contacts') +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))

physical_plot + conversational_plot + plot_layout(guides = 'collect')
ggsave(file.path('figures','contact_matrices.png'),
       width = 14, height = 7)

## PLOTTING WITH RECIPROCITY

ggplot_reciprocal_matrix(reconnect) + ggplot_reciprocal_matrix(reconnect_individual) +
  ggplot_reciprocal_matrix(reconnect_physical) + ggplot_reciprocal_matrix(reconnect_conversational) 

physical_plot <- ggplot_reciprocal_matrix(reconnect_physical) + 
  scale_fill_viridis(limits = c(0,1.3)) + ggtitle('Physical contacts') +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))
conversational_plot <- ggplot_reciprocal_matrix(reconnect_conversational) + 
  scale_fill_viridis(limits = c(0,1.3)) + ggtitle('Conversational contacts') +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))

physical_plot + conversational_plot + plot_layout(guides = 'collect')
ggsave(file.path('figures','contact_matrices_reciprocal.png'),
       width = 14, height = 7)


## MORE FILTERS

matrix_data <- data.table()

for(cont_type in 1:2){
  for(location in c('Home','Work','School','Other')){
    
    reconnect_in <- copy(reconnect)
    reconnect_in$contacts <- reconnect_in$contacts[reconnect_in$contacts$phys_contact == cont_type & 
                                                     reconnect_in$contacts$cnt_location == location, ]
    
    matrix_data <- rbind(matrix_data,
                         calc_matrix(reconnect_in) %>% 
                           mutate(cnt_location = location,
                                  type = ifelse(cont_type==1, "Physical", "Conversational")))
    
  }
}

## PLOT

plot_cm_setting <- function(location,
                            m_data = matrix_data){
  
  m_data %>% 
    filter(cnt_location == location) %>% 
    ggplot() +
    geom_tile(aes(x = part_age, y = cont_age, 
                  fill = value)) + 
    labs(x = 'Participant age group', y = 'Contact age group', fill = 'Mean\ncontacts') +
    theme_bw() + scale_fill_viridis(limits = c(0,NA)) +
    facet_grid(type ~ cnt_location) + 
    theme(strip.background = element_rect(fill = 'white',
                                          color = 'white'),
          strip.text = element_text(face = 'bold', size = 12)) +
    coord_fixed() +
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))

}

plot_cm_setting('Home') +
  plot_cm_setting('Work') + 
  plot_cm_setting('School') + 
  plot_cm_setting('Other') + 
  plot_layout(nrow = 1)

ggsave(file.path('figures','contact_matrices_setting.png'),
       width = 22, height = 9)

## with reciprocity

matrix_data_reciprocal <- data.table()

for(cont_type in 1:2){
  for(location in c('Home','Work','School','Other')){
    
    reconnect_in <- copy(reconnect)
    reconnect_in$contacts <- reconnect_in$contacts[reconnect_in$contacts$phys_contact == cont_type & 
                                                     reconnect_in$contacts$cnt_location == location, ]
    
    matrix_data_reciprocal <- rbind(matrix_data_reciprocal,
                                    calc_matrix(reconnect_in, reciprocal = T) %>% 
                                      mutate(cnt_location = location,
                                             type = ifelse(cont_type==1, "Physical", "Conversational")))
    
  }
}

plot_cm_setting('Home', matrix_data_reciprocal) +
  plot_cm_setting('Work', matrix_data_reciprocal) + 
  plot_cm_setting('School', matrix_data_reciprocal) + 
  plot_cm_setting('Other', matrix_data_reciprocal) + 
  plot_layout(nrow = 1)

ggsave(file.path('figures','contact_matrices_setting_recip.png'),
       width = 22, height = 9)


#### ASSORTATIVITY ####

assortativity <- data.table()

for(cont_type in unique(matrix_data$type)){
  for(location in c('Home','Work','School','Other')){
    
    asst <- fcn_assortativity('age_group',
                              matrix_data %>% filter(type == cont_type,
                                                                cnt_location == location),
                              age_population)
    
    assortativity <- rbind(assortativity,
                                      data.table(
                                        cnt_location = location,
                                        type = cont_type,
                                        Q = asst$mean
                                      ))
    
  }
}

write_csv(assortativity,
          file.path('outputs', 'assortativity.csv'))

assortativity %>% 
  ggplot() +
  geom_point(aes(x = cnt_location, y = Q, col = type),
             size = 3) +
  geom_text(aes(x = cnt_location, y = Q + 0.012, label = round(Q, 3)), size = 3) +
  theme_bw() +
  scale_colour_manual(values = c('#F17105', '#1A8FE3')) +
  ylim(c(0,0.4)) +
  labs(x = 'Setting', y = 'Age assortativity (Q)', color = '')

ggsave(file.path('figures','age_assortativity.png'),
       width = 7, height = 5)

## on reciprocal matrices

assortativity_reciprocal <- data.table()

for(cont_type in unique(matrix_data_reciprocal$type)){
  for(location in c('Home','Work','School','Other')){
    
    asst <- fcn_assortativity('age_group',
                      matrix_data_reciprocal %>% filter(type == cont_type,
                                                        cnt_location == location),
                      age_population)
    
    assortativity_reciprocal <- rbind(assortativity_reciprocal,
                                      data.table(
                                        cnt_location = location,
                                        type = cont_type,
                                        Q = asst$mean
                                      ))
    
  }
}

write_csv(assortativity_reciprocal,
          file.path('outputs', 'assortativity_reciprocal.csv'))

assortativity_reciprocal %>% 
  ggplot() +
  geom_point(aes(x = cnt_location, y = Q, col = type),
             size = 3) +
  geom_text(aes(x = cnt_location, y = Q + 0.015, label = round(Q, 3)), size = 3) +
  theme_bw() +
  scale_colour_manual(values = c('#F17105', '#1A8FE3')) +
  ylim(c(0,0.4)) +
  labs(x = 'Setting', y = 'Age assortativity (Q)', color = '')

ggsave(file.path('figures','age_assortativity_recip.png'),
       width = 7, height = 5)





