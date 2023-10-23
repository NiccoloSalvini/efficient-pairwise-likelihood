source(here("notebooks", "11.2_simulation_functions.R"))

n <- 200

set.seed(1234)

coord_1 <- runif(n, min = 400, max = 5000)
coord_2 <- runif(n, min = 400, max = 5000)

data <- data.frame(
  coord_1 = coord_1,
  coord_2 = coord_2
)

couplets  <- GenerateSpatialCoupletsBuffered(data,  radius_filter = "custom", buffer = 500)

viz_couples(data = data, coppiette = couplets[[1]])

library(patchwork)
library(latex2exp)
couplets  <- GenerateSpatialCoupletsBuffered(data,  radius_filter = "mean")
gg_mean = viz_couples(data = data, coppiette = couplets[[1]]) +
  labs(title = TeX(r'($\textbf{r_{mean}} = $ 2371.868)'))


couplets  <- GenerateSpatialCoupletsBuffered(data,  radius_filter = "custom", buffer = 50)
gg_50 = viz_couples(data = data, coppiette = couplets[[1]]) +
  labs(title = TeX(r'($\textbf{r_{buff50}}$)'))

couplets  <- GenerateSpatialCoupletsBuffered(data,  radius_filter = "custom", buffer = 200)
gg_200 = viz_couples(data = data, coppiette = couplets[[1]]) +
  labs(title = TeX(r'($\textbf{r_{buff200}}$)'))

couplets  <- GenerateSpatialCoupletsBuffered(data,  radius_filter = "custom", buffer = 500)
gg_500 = viz_couples(data = data, coppiette = couplets[[1]])+
  labs(title = TeX(r'($\textbf{r_{buff500}}$)'))

couplets  <- GenerateSpatialCoupletsBuffered(data,  radius_filter = "custom", buffer = 700)
gg_700 = viz_couples(data = data, coppiette = couplets[[1]])+
  labs(title = TeX(r'($\textbf{r_{buff700}}$)'))

couplets  <- GenerateSpatialCoupletsBuffered(data,  radius_filter = "custom", buffer = 1000)
gg_1000 = viz_couples(data = data, coppiette = couplets[[1]])+
  labs(title = TeX(r'($\textbf{r_{buff1000}}$)'))

couplets  <- GenerateSpatialCoupletsBuffered(data,  radius_filter = "max")
gg_max = viz_couples(data = data, coppiette = couplets[[1]])+
  labs(title = TeX(r'($\textbf{r_{max}}  = $ 6053.066)'))

bufferd_couples = (gg_mean + gg_max)/(gg_50 + gg_200)
ggsave(plot = bufferd_couples,filename = here("img", "plots", "bufferd_couples.png"))
