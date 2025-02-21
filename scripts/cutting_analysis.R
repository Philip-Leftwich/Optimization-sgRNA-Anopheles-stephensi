# Cutting and homing when crossed to white eyes ====
## Multiple Cas9 lines tested in a split-drive

font_add_google("Open Sans", "Sans")

source("scripts/functions.R")
source("scripts/custom_theme.R")


cutting_tbl <- tbl <-read_data(4)


cutting_data <- cutting_tbl %>% 
  rename(`F1 cross`=1,
         genotype_A=4,
         genotype_AB=5,
         genotype_B=6,
         genotype_WT=7,
         mosaic_A=15,
         mosaic_AB=16,
         mosaic_B=17,
         mosaic_WT=18) %>% 
  select(c(1:7, 15:18, 19)) %>% 
  fill(`F1 cross`) %>% 
  drop_na(genotype_A:mosaic_WT) %>% 
  mutate(cross= str_remove_all(`F1 cross`, "[()Xx]")) %>% ### keep F1 cross variable
  mutate(cross = str_remove_all(cross, "Cd ")) %>% 
  filter(!str_detect(cross,"2072")) %>% 
  separate(cross, into =c("male_parent", "female_parent"), 
           sep=" ", extra = "merge", fill = "left") %>% 
  mutate(Cas9_parent=if_else(male_parent=="KO", "Female", "Male")) %>% 
  separate(male_parent, into=c("m_male_grandparent", "m_female_grandparent"), sep=":", remove=FALSE) %>% 
  separate(female_parent, into=c("f_male_grandparent", "f_female_grandparent"), sep=":", remove=FALSE) %>% 
  mutate(m_male_grandparent=str_trim(m_male_grandparent, side = "both")) %>% # extra white spaces causing issues with ifelse statement
  mutate(f_male_grandparent=str_trim(f_male_grandparent, side = "both")) %>% 
  mutate(Cas9_grandparent=ifelse(m_male_grandparent=="1590B", "Male",
                                 ifelse(f_male_grandparent=="1590B", "Male", "Female"))) %>% 
  rename("Line"=filename)

## Models====

cutting_data_model <- cutting_data %>% 
  mutate(win = (genotype_A+genotype_AB),
         loss = (genotype_B+genotype_WT),
         win1=(mosaic_A+mosaic_AB),
         loss1=((win+loss)- win1)) %>% 
  mutate(id=row_number())

# use the same parameters as the models for crosses to WT (homing_mosaic_analysis.R)

cut.model.1 <- glmmTMB((cbind(win,loss))~ Cas9_parent+Cas9_grandparent+Line+Cas9_parent:Cas9_grandparent+Cas9_parent:Line+Cas9_grandparent:Line+(1|Line/Cas9_grandparent/Cas9_parent/id), family=binomial, data=cutting_data_model, REML=FALSE)


mosaic.model.1 <- glmmTMB((cbind(win1,loss1))~ Cas9_parent+Cas9_grandparent+Line+Cas9_parent:Cas9_grandparent+Cas9_parent:Line+Cas9_grandparent:Line+(1|Line/Cas9_grandparent/Cas9_parent/id), family=binomial, data=cutting_data_model, REML=FALSE)





## Figure====

cutting_summary <- emmeans::emmeans(cut.model.1, specs=pairwise~Cas9_parent+Cas9_grandparent+Line+Cas9_parent:Cas9_grandparent+Cas9_parent:Line+Cas9_grandparent:Line, type="response") %>% .$emmeans %>% as_tibble() %>% bind_cols(homing_or_cutting="cutting")

homing_summary <- emmeans::emmeans(mosaic.model.1, specs=pairwise~Cas9_parent+Cas9_grandparent+Line+Cas9_parent:Cas9_grandparent+Cas9_parent:Line+Cas9_grandparent:Line, type="response") %>% .$emmeans %>% as_tibble() %>% bind_cols(homing_or_cutting="homing")

total_CI <- rbind(cutting_summary, homing_summary) %>% 
  mutate(Cas9_parent = factor(Cas9_parent, levels=c("Female", "Male"), labels=c("\u2640", "\u2642"))) %>%
  mutate(Cas9_grandparent = factor(Cas9_grandparent, levels=c("Female", "Male"), labels=c("\u2640", "\u2642"))) %>% 
  mutate(Line = factor(Line, levels=c("1626", "1757", "1758", "1759"), labels=c("cd^U6A", "cd^U6B", "cd^U6C", "cd^{7~SK}"))) 


custom_labeller <- labeller(
  Line = as_labeller(label_parsed),
  Cas9_grandparent = label_value, # keep the default labeller
  Cas9_parent = label_value,
  homing_or_cutting = (label_parsed)# keep the default labeller
)

cutting_data_model %>% 
  mutate(cutting=((win1/(win+loss))*100)) %>% 
  mutate(homing=((win/(win+loss)*100))) %>% 
  mutate(Cas9_parent = factor(Cas9_parent, levels=c("Female", "Male"), labels=c("\u2640", "\u2642"))) %>%
  mutate(Cas9_grandparent = factor(Cas9_grandparent, levels=c("Female", "Male"), labels=c("\u2640", "\u2642"))) %>%
  mutate(Line = factor(Line, levels=c("1626", "1757", "1758", "1759"), labels=c("cd^U6A", "cd^U6B", "cd^U6C", "cd^{7~SK}"))) %>% 
  pivot_longer(cols=homing:cutting, names_to="homing_or_cutting", values_to="count") %>% 
  ggplot(aes(x=homing_or_cutting, y=count, group=id, colour=Cas9_parent, shape=homing_or_cutting))+
  geom_point(aes(size=win+loss), fill="white", alpha=0.6)+
  geom_line(alpha=0.4, size=0.1)+
  geom_errorbar(data=total_CI, aes(min=(asymp.LCL*100), max=(asymp.UCL*100), y=response, group=NA),width=0,  size=1.2, position=position_nudge(x=0.4))+
  geom_point(data=total_CI, aes(y=response*100, group=NA, fill=after_scale(desaturate(lighten(colour, .6), .6))), size=2, position=position_nudge(x=0.4), stroke=0.6)+
  scale_shape_manual(values=c(21,24), labels=c("pink eyes",bquote(cd^sgRNAs)))+
  scale_size(range=c(0,6),
             breaks=c(50,100,150))+
  scale_color_manual(values=c("#a86826", "#006c89"))+
  guides(fill="none", colour="none")+
  labs(x="", 
       y="Percentage of individuals scored",
       size="Number of offspring",
       shape="")+
  scale_y_continuous(limits=c(50,100),
    labels=scales::percent_format(scale=1) # automatic percentages
  )+
  facet_nested_wrap(~Line+Cas9_grandparent+Cas9_parent, 
                    nrow=1,
                    strip.position="bottom",
                    scales="free_x",
                    labeller = custom_labeller)+ # reduce lines
  labs(tag = bquote(atop(Cas9[F1], Cas9[F0]))) +
  theme_custom2() +
  theme(axis.text.x=element_blank(),
        plot.tag.position = c(.02,.25),
        legend.position="bottom",
        legend.direction = "vertical",
        strip.background = element_blank(),
        strip.placement="outside",
        panel.spacing = unit(.5, "lines"),
        strip.text.x = element_text(margin = margin(5,0,5,0),
                                    hjust = 0.5),
        plot.margin = margin(1,1,1,1.2, "cm"),
        ggh4x.facet.nestline = element_line(colour = "black"))  # center strip text



                  