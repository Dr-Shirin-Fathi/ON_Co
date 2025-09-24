setwd('C:\Users\Administrator\Desktop\maybe_accept')
setwd('C:/Users/Administrator/Desktop/maybe_accept')


library(readxl)    
library(dplyr)     
library(stringr)   
library(ggplot2)   
library(tidyr)    


dep_full  <- read_excel("deprescription_impact_FULL_with_risk.xlsx")
dep_x     <- read_excel("deprescription_impact_X_only.xlsx")
edges_x   <- read_excel("high_risk_interactions_X_only.xlsx")
top10_xf  <- read_excel("top10_riskX_fragmentation_drugs.xlsx")
View(dep_full)
View(dep_x)
View(edges_x)
View(top10_xf)

# ATC Level-1 (1st letter) and Level-3/4 prefixes if needed
dep_full  <- dep_full %>%
  mutate(ATC_L1 = substr(`ATC Code`, 1, 1),
         ATC_L3 = substr(`ATC Code`, 1, 3),
         ATC_L4 = substr(`ATC Code`, 1, 4))
dep_x <- dep_x %>%
  mutate(ATC_L1 = substr(`ATC Code`, 1, 1),
         ATC_L3 = substr(`ATC Code`, 1, 3),
         ATC_L4 = substr(`ATC Code`, 1, 4))
top10_xf <- top10_xf %>%
  mutate(ATC_L1 = substr(`ATC Code`, 1, 1),
         ATC_L3 = substr(`ATC Code`, 1, 3),
         ATC_L4 = substr(`ATC Code`, 1, 4))

#“Deprescription fragmentation” signature (drug-level)

frag_full <- dep_full %>%
  mutate(
    # components growth beyond a connected baseline (1)
    frag_components = pmax(`Num Components After Removal` - 1, 0),
    # smaller giant component => more fragmentation (inverse-scaled by a typical total size;
    # if you don't know the original N, you can treat LargeComp drop as the signal itself)
    frag_giant_drop = max(`Largest Component Size`, na.rm = TRUE) - `Largest Component Size`,
    # simple normalized signature 
    frag_components_z = (frag_components - min(frag_components, na.rm=TRUE)) /
      (max(frag_components, na.rm=TRUE) - min(frag_components, na.rm=TRUE)),
    frag_giant_drop_z = (frag_giant_drop - min(frag_giant_drop, na.rm=TRUE)) /
      (max(frag_giant_drop, na.rm=TRUE) - min(frag_giant_drop, na.rm=TRUE)),
    # equally weighted composite fragmentation index
    frag_index = 0.5*frag_components_z + 0.5*frag_giant_drop_z
  ) %>%
  arrange(desc(frag_index))

# Do the same for X-only 
frag_x <- dep_x %>%
  mutate(
    frag_components = pmax(`Num Components After Removal` - 1, 0),
    frag_giant_drop = max(`Largest Component Size`, na.rm = TRUE) - `Largest Component Size`,
    frag_components_z = (frag_components - min(frag_components, na.rm=TRUE)) /
      (max(frag_components, na.rm=TRUE) - min(frag_components, na.rm=TRUE)),
    frag_giant_drop_z = (frag_giant_drop - min(frag_giant_drop, na.rm=TRUE)) /
      (max(frag_giant_drop, na.rm=TRUE) - min(frag_giant_drop, na.rm=TRUE)),
    frag_index = 0.5*frag_components_z + 0.5*frag_giant_drop_z
  ) %>%
  arrange(desc(frag_index))

#  Oncology focus (ATC_L1 == 'L')
frag_full_onc <- frag_full %>% filter(ATC_L1 == "L")
frag_x_onc    <- frag_x    %>% filter(ATC_L1 == "L")
# Top "fragmenters" overall and within oncology
top_frag_overall <- frag_full %>% select(`ATC Code`,`Generic Name`, frag_index, `Original Degree`) %>% head(20)
top_frag_onc     <- frag_full_onc %>% select(`ATC Code`,`Generic Name`, frag_index, `Original Degree`) %>% head(20)
print(top_frag_overall)
print(top_frag_onc)
write.csv(top_frag_overall, "top_frag_overall.csv")
write.csv(top_frag_onc, "top_frag_onc.csv")
#A ranked list of drugs whose removal most fragments the network (network signature of centrality-to-robustness).
ggplot(top_frag_onc, aes(x=reorder(`Generic Name`, frag_index), y=frag_index)) +
  geom_col() +
  coord_flip() +
  labs(title = "Oncology (ATC L*) drugs with strongest fragmentation on removal",
       x = "", y = "Fragmentation Index")
# Build a long node list (degree) 
edges_x_clean <- edges_x %>%
  rename(d1 = `Drug 1 (Generic)`,
         d2 = `Drug 2 (Generic)`,
         w  = `Interaction Weight`)
#  Degree (unique partners) and weighted degree (sum of weights)
deg_df <- bind_rows(
  edges_x_clean %>% select(drug = d1, partner = d2, w),
  edges_x_clean %>% select(drug = d2, partner = d1, w)
) %>%
  group_by(drug) %>%
  summarise(
    degree = n_distinct(partner),
    wdegree = sum(w, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(wdegree), desc(degree))
# Top X-only hubs
top_hubs_x <- deg_df %>% head(20)
print(top_hubs_x)
# If ATC codes are needed for these drugs, join from dep_full (best map via Generic Name)
deg_with_atc <- deg_df %>%
  left_join(dep_full %>% select(`Generic Name`, `ATC Code`, ATC_L1, ATC_L3, ATC_L4) %>% distinct(),
            by = c("drug" = "Generic Name"))
top_hubs_x_onc <- deg_with_atc %>% filter(ATC_L1 == "L") %>% arrange(desc(wdegree)) %>% head(20)
print(top_hubs_x_onc)
summary_depresc <- frag_x %>%
  mutate(is_oncology = if_else(ATC_L1 == "L", "Oncology (L*)", "Non-Oncology")) %>%
  group_by(is_oncology) %>%
  summarise(
    med_frag_components = median(frag_components, na.rm=TRUE),
    med_frag_index      = median(frag_index, na.rm=TRUE),
    n = n(),
    .groups = "drop"
  )
print(summary_depresc)
View(summary_depresc)
write.csv(summary_depresc, "summary_depresc.csv")
