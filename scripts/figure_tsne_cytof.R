plot_group <- ggplot(tsne_2plot, aes(x = X1, y = X2, col = as.factor(GROUP),
                       label = Id.Cryostem.R)) +
  geom_point() +
  geom_text(aes(label=Id.Cryostem.R),hjust=0, vjust=0) +
  ggtitle("Tolerance group of the recipients") +
  theme_minimal() +
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank()) +
  theme(legend.text=element_text(size=12)) + labs(color='GROUP')

plot_cohort <- ggplot(tsne_2plot, aes(x = X1, y = X2, col = as.factor(cohort),
                                      label = Id.Cryostem.R)) +
  geom_point() +
  geom_text(aes(label=Id.Cryostem.R),hjust=0, vjust=0) +
  ggtitle("Origin cohort of the recipients") +
  theme_minimal() +
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank()) +
  theme(legend.text=element_text(size=12)) + labs(color='COHORT')

plot_date <- ggplot(tsne_2plot, aes(x = X1, y = X2, col = as.factor(DATEOFCYTOFEXPERIMENT),
                                      label = Id.Cryostem.R)) +
  geom_point() +
  geom_text(aes(label=Id.Cryostem.R),hjust=0, vjust=0) +
  ggtitle("Date of experiment") +
  theme_minimal() +
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank()) +
  theme(legend.text=element_text(size=12)) + labs(color='DATE')

plot_cmv <- ggplot(tsne_2plot, aes(x = X1, y = X2, col = as.factor(CMVStatus),
                                      label = Id.Cryostem.R)) +
  geom_point() +
  geom_text(aes(label=Id.Cryostem.R),hjust=0, vjust=0) +
  scale_color_manual(values=c("grey", "black", "white")) +
  ggtitle("CMV status of the recipients") +
  theme_minimal() +
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank()) +
  theme(legend.text=element_text(size=12)) + labs(color='CMV')

plot_gender <- ggplot(new_tsne_2plot, aes(x = X1, y = X2, col = as.factor(gender_comp),
                           label = Id.Cryostem.R)) +
  geom_point() +
  geom_text(aes(label=Id.Cryostem.R),hjust=0, vjust=0) +
  ggtitle("Gender compatibility between donors and recipients") +
  #scale_color_manual(values=c("brown1", "orange", "lightgreen", "royalblue")) +
  theme_minimal() +
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank()) +
  theme(legend.text=element_text(size=12)) + labs(color='GENDER COMPATIBILITY')

plot_age <- ggplot(new_tsne_2plot, aes(x = X1, y = X2, col = age_recip,
                           label = Id.Cryostem.R)) +
  geom_point() +
  geom_text(aes(label=Id.Cryostem.R),hjust=0, vjust=0) +
  ggtitle("Age of the recipients") +
  #scale_color_manual(values=c("brown1", "orange", "lightgreen", "royalblue")) +
  theme_minimal() +
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank()) +
  theme(legend.text=element_text(size=12)) + labs(color='AGE')

library(patchwork)

{plot_group + plot_cohort + plot_date} / {plot_gender + plot_age + plot_cmv}
