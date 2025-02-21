# Cutting and homing when crossed to wildtype ====
## Multiple Cas9 lines tested in a split-drive


source("scripts/functions.R")
source("scripts/custom_theme.R")


tbl <-read_data(2)


data <- tbl %>% 
  rename(genotype_A=4,
         genotype_AB=5,
         genotype_B=6,
         genotype_WT=7,
         mosaic_A=15,
         mosaic_AB=16,
         mosaic_B=17,
         mosaic_WT=18) %>% 
  select(c(1:7,12, 15:18, 23)) %>% 
  fill(`F1 cross`) %>% 
  mutate(cross = str_remove_all(`F1 cross`, "[()Xx]")) %>% 
  filter(!str_detect(cross,"2072")) %>% 
  separate(cross, into =c("male_parent", "female_parent"), 
           sep=" ", extra = "merge", fill = "left") %>% 
  mutate(Cas9_parent=if_else(male_parent=="WT", "Female", "Male")) %>% 
  separate(male_parent, into=c("m_male_grandparent", "m_female_grandparent"), sep=":", remove=FALSE) %>% 
  separate(female_parent, into=c("f_male_grandparent", "f_female_grandparent"), sep=":", remove=FALSE) %>% 
  mutate(m_male_grandparent=str_trim(m_male_grandparent, side = "both")) %>% # extra white spaces causing issues with ifelse statement
  mutate(f_male_grandparent=str_trim(f_male_grandparent, side = "both")) %>% 
  mutate(Cas9_grandparent=ifelse(m_male_grandparent=="1590B", "Male",
                                 ifelse(f_male_grandparent=="1590B", "Male", "Female"))) %>% 
  rename("Line"=filename)


## Homing model====


homing_data <- data %>% 
  mutate(win = (genotype_A+genotype_AB),
         loss = (genotype_B+genotype_WT),
         unhatched=`No. of embryos`-Total) %>% 
  drop_na(genotype_A:genotype_WT) %>% 
  mutate(id = row_number())

homing.model <- glmmTMB((cbind(win,loss))~ Cas9_parent+Cas9_grandparent+Line+Cas9_parent:Cas9_grandparent+Cas9_parent:Line+Cas9_grandparent:Line+(1|Line/Cas9_grandparent/Cas9_parent/id), family=binomial, data=homing_data)

DHARMa_check(homing.model)

total_CI <- emmeans::emmeans(homing.model, specs=pairwise~Cas9_parent+Cas9_grandparent+Line+Cas9_parent:Cas9_grandparent+Cas9_parent:Line+Cas9_grandparent:Line, type="response") %>% 
  .$emmeans %>% as_tibble() %>% 
  mutate(Cas9_parent = factor(Cas9_parent, levels=c("Female", "Male"), labels=c("♀", "♂"))) %>%
  mutate(Cas9_grandparent = factor(Cas9_grandparent, levels=c("Female", "Male"), labels=c("♀", "♂"))) %>% 
  mutate(Line = factor(Line, levels=c("1626", "1757", "1758", "1759"), labels=c("cd^U6A", "cd^U6B", "cd^U6C", "cd^{7~SK}"))) 
# This message reflects a small issue in that it gets confused due to the formula not having a left-hand side, and thinks it detects a response transformation. But we can ignore it, or use update() to get rid of it

custom_labeller <- labeller(
  Line = as_labeller(label_parsed),
  Cas9_grandparent = label_value, # keep the default labeller
  Cas9_parent = label_value       # keep the default labeller
)

homing_plot <- homing_data %>% 
  mutate(homing=((win/(win+loss)*100))) %>% 
  mutate(Cas9_parent = factor(Cas9_parent, levels=c("Female", "Male"), labels=c("♀", "♂"))) %>%
  mutate(Cas9_grandparent = factor(Cas9_grandparent, levels=c("Female", "Male"), labels=c("♀", "♂"))) %>%
  mutate(Line = factor(Line, levels=c("1626", "1757", "1758", "1759"), labels=c("cd^U6A", "cd^U6B", "cd^U6C", "cd^{7~SK}"))) %>% 
  ggplot(aes(x=Cas9_parent, y=homing, group=id, colour=Cas9_parent))+
  geom_point(aes(size=win+loss), fill="white", alpha=0.6, position=position_nudge(x=-0.3), shape=21)+
  geom_errorbar(data=total_CI, aes(min=(asymp.LCL*100), max=(asymp.UCL*100), y=response, group=NA),width=0,  linewidth=1.2, position=position_nudge(x=0.3))+
  geom_point(data=total_CI, shape=21, aes(y=response*100, group=NA, fill=after_scale(desaturate(lighten(colour, .6), .6))), size=3, position=position_nudge(x=0.3), stroke=0.6)+
  scale_size(range=c(0,6),
             breaks=c(50,100,150))+
  scale_color_manual(values=c("#a86826", "#006c89"))+
  guides(fill=FALSE)+
  labs(x="", 
       y="Percentage of individuals scored",
       size="Number of offspring",
       shape="",
       colour = "Cas9 parent")+
  scale_y_continuous(limits=c(50,100),
                     labels=scales::percent_format(scale=1) # automatic percentages
  )+
  facet_nested_wrap(~Line+Cas9_grandparent+Cas9_parent, 
                    nrow=1,
                    strip.position="bottom",
                    scales="free_x",
                    labeller = custom_labeller) # reduce lines


 homing_plot <- homing_plot +  theme_custom()+
  theme(axis.text.x=element_blank(),
        legend.position="bottom",
        legend.direction = "vertical",
        strip.background = element_blank(),
        strip.placement="outside",
        panel.spacing = unit(.5, "lines"),
        strip.text.x = element_text(margin = margin(5,0,5,0),
                                    hjust = 0.5))# center strip text



homing_plot <- homing_plot + 
  labs(tag = "Cas9 F1\n\nCas9 F0") +
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
        ggh4x.facet.nestline = element_line(colour = "black")) + # center strip text
 
  
  coord_cartesian(clip="off")


# Mosaicism with each genotype====

counts <- data %>% 
  mutate(genotype_A= (genotype_A-mosaic_A),
         genotype_B= (genotype_B-mosaic_B),
         genotype_AB=(genotype_AB-mosaic_AB),
         genotype_WT=(genotype_WT-mosaic_WT)
          ) %>% 
  pivot_longer(cols=genotype_A:mosaic_WT, names_to="genotype", values_to="count" ) %>% 
  separate(genotype, into=c("mosaic", "genotype"), sep= "_") %>% 
  filter(!mosaic %in%  "Total") # remove total column which throws off numbers

count_type <- counts %>% 
  group_by(Cas9_parent,Cas9_grandparent, Line, mosaic, genotype) %>% 
  summarise(count= sum(count, na.rm=TRUE))

total <- counts %>% 
  group_by(Cas9_parent, Cas9_grandparent, Line) %>% 
  summarise(total= sum(count, na.rm=TRUE))
  
full_join(count_type, total) %>% 
  drop_na() %>% 
  mutate(percentage=(count/total)*100) %>% 
  
  ggplot(aes(x=interaction(Cas9_parent, Cas9_grandparent, Line), y=percentage, fill=mosaic))+
  geom_bar(position=position_stack(reverse=TRUE), 
           stat="identity",
           colour="black",
           alpha=0.8)+
  facet_wrap(~genotype, nrow=1)+
  coord_flip()+
  scale_fill_manual(values=c("#C64974", "#6EB067"),
                    labels=c("wildtype", "mosaic"))+
  scale_y_continuous(limits=c(0,
                              55),
                     labels=scales::percent_format(scale=1) # automatic percentages
                     )+
  theme_custom()+
  labs(x="",
       y= "",
       fill="Phenotype")

custom_labeller <- labeller(
  Line = as_labeller(label_parsed),
  Cas9_grandparent = label_value, # keep the default labeller
  Cas9_parent = label_value,
  genotype = as_labeller(label_parsed)# keep the default labeller
)

  mosaic_plot <-  
    full_join(count_type, total) %>% 
      drop_na() %>% 
      mutate(percentage=(count/total)*100) %>% 
      mutate(Cas9_parent = factor(Cas9_parent, levels=c("Female", "Male"), labels=c("\u2640", "\u2642"))) %>%
      mutate(Cas9_grandparent = factor(Cas9_grandparent, levels=c("Female", "Male"), labels=c("\u2640", "\u2642"))) %>%
      mutate(Line = factor(Line, levels=c("1626", "1757", "1758", "1759"), labels=c("cd^U6A", "cd^U6B", "cd^U6C", "cd^{7~SK}"))) %>% 
      mutate(genotype = factor(genotype, levels = c("A", "AB", "B", "WT"), labels = c("cd^sgRNAs", " ", "zpg^Cas9", "WT"))) %>% 
      ggplot(aes(x=interaction(Cas9_parent, Cas9_grandparent, Line), y=percentage, fill=mosaic))+
      geom_bar(position=position_stack(reverse=TRUE), 
               stat="identity",
               colour="black",
               alpha=0.8)+
      facet_nested(Line+Cas9_grandparent+Cas9_parent~genotype, scales="free_y", switch="both",
                   labeller = custom_labeller)+
      
      scale_fill_manual(values=c("#C64974", "#6EB067"),
                        labels=c("wildtype", "mosaic"))+
      scale_y_continuous(limits=c(0,
                                  55),
                         labels=scales::percent_format(scale=1) # automatic percentages
      )+
      theme_custom2()+
      theme(axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            strip.text.y.left = element_text(angle = 0),
            strip.background = element_blank(),
            panel.spacing.y = unit(.5, "lines"),
            panel.grid = element_blank(),
            legend.position = "bottom",
            legend.direction = "vertical",
            plot.tag.position = c(.2,1.05),
            ggh4x.facet.nestline = element_line(colour = "black"),
            plot.margin = margin(3,1,1,1, "cm"))+
      labs(x="",
           y= "",
           fill="Phenotype"
           )+
    labs(tag = " Cas9    Cas9  \nF0      F1") + coord_flip(clip = "off")
  
  
  homing_plot + mosaic_plot
  

# Stats for mosaicism requested by reviewer====
  
  wide_mosaic <- count_type %>% 
    mutate(mosaic = if_else(mosaic == "genotype", "wildtype", "mosaic")) %>% 
    pivot_wider(names_from = "mosaic", values_from = "count") %>% 
    mutate(proportion = mosaic/(mosaic+wildtype),
           total = mosaic + wildtype)

null_model <- glm(cbind(mosaic,wildtype) ~ 1, family = binomial, data = wide_mosaic)  
  
mosaic_model <- glm(cbind(mosaic,wildtype) ~ Cas9_parent + Cas9_grandparent + Line + genotype, family = binomial, data = wide_mosaic)  

mosaic_model2 <- glm(cbind(mosaic,wildtype) ~ Cas9_grandparent+ Line + genotype, family = binomial, data = wide_mosaic)  


drop1(mosaic_model, test = "Chi")

# Extract Deviance
dev_null <- deviance(null_model)
dev_full <- deviance(mosaic_model)
dev_partial <- deviance(mosaic_model2)

# Compute R-squared
r2_full <- 1 - (dev_full / dev_null)   # Full model R^2
r2_partial <- 1 - (dev_partial / dev_null)  # Partial model R^2
r2_x2 <- r2_full - r2_partial  # Contribution of x2

# Print Results
cat("Full Model R^2:", r2_full, "\n")
cat("Partial Model (x1) R^2:", r2_partial, "\n")
cat("Additional Contribution of x2:", r2_x2, "\n")
  