##################################################
# UPSET PLOT MATRIZ DE GENES
##################################################

library(tidyverse)
library(ComplexUpset)
library(openxlsx)
library(ggplot2)

##################################################
# DIRECTORIO DE TRABAJO
##################################################

setwd("C:/Users/alema/Desktop/ESTUDIOS/UCO/TFG/Actividades del TFG/Filtrado Matrices/.CSV/TFG_UpSet_R")

##################################################
# ARCHIVOS ORIGINALES
##################################################

archivos <- c(
  "data/PEG_UP_HOJA.csv",
  "data/PEG_UP_RAIZ.csv",
  "data/PEG_DOWN_HOJA.csv",
  "data/PEG_DOWN_RAIZ.csv",
  "data/SAL_UP_HOJA.csv",
  "data/SAL_UP_RAIZ.csv",
  "data/SAL_DOWN_HOJA.csv",
  "data/SAL_DOWN_RAIZ.csv"
)

##################################################
# CREAR TABLA BINARIA
##################################################

lista_genes <- list()

for (archivo in archivos) {
  
  datos <- read.csv(archivo)
  
  nombre <- tools::file_path_sans_ext(basename(archivo))
  
  lista_genes[[nombre]] <- unique(datos$Gene_id)
  
}

todos_los_genes <- unique(unlist(lista_genes))

tabla_binaria <- data.frame(
  Gene_id = todos_los_genes
)

for (condicion in names(lista_genes)) {
  
  tabla_binaria[[condicion]] <- ifelse(
    tabla_binaria$Gene_id %in% lista_genes[[condicion]],
    1,
    0
  )
  
}

##################################################
# GUARDAR TABLA BINARIA
##################################################

write.csv2(
  tabla_binaria,
  "results/tabla_binaria_upset.csv",
  row.names = FALSE
)

##################################################
# ORDEN DE CONDICIONES
##################################################

condiciones <- c(
  "PEG_UP_HOJA",
  "PEG_UP_RAIZ",
  "SAL_UP_HOJA",
  "SAL_UP_RAIZ",
  "PEG_DOWN_HOJA",
  "PEG_DOWN_RAIZ",
  "SAL_DOWN_HOJA",
  "SAL_DOWN_RAIZ"
)

##################################################
# TABLA DE NETWORKS
##################################################

tabla_networks <- tabla_binaria %>%
  mutate(
    vector_binario = apply(
      select(., all_of(condiciones)),
      1,
      paste0,
      collapse = ""
    ),
    
    n_condiciones = rowSums(
      select(., all_of(condiciones))
    ),
    
    condiciones_presentes = apply(
      select(., all_of(condiciones)),
      1,
      function(x){
        
        presentes <- condiciones[
          as.numeric(x) == 1
        ]
        
        paste(
          presentes,
          collapse = " + "
        )
        
      }
    )
  ) %>%
  group_by(
    vector_binario,
    n_condiciones,
    condiciones_presentes
  ) %>%
  summarise(
    n_genes = n(),
    Gene_id = paste(
      Gene_id,
      collapse = "; "
    ),
    .groups = "drop"
  ) %>%
  arrange(
    desc(n_condiciones),
    desc(n_genes)
  ) %>%
  mutate(
    Network = paste0(
      "Network_",
      row_number()
    )
  )

##################################################
# GUARDAR TABLA NETWORKS
##################################################

write.xlsx(
  tabla_networks,
  "results/tabla_networks.xlsx",
  rowNames = FALSE
)

##################################################
# PREPARAR UPSET
##################################################

tabla_binaria <- tabla_binaria %>%
  mutate(
    n_condiciones = rowSums(
      select(., all_of(condiciones))
    ),
    
    Condiciones = paste0(
      "Condiciones ",
      n_condiciones
    )
  )

tabla_binaria[condiciones] <- lapply(
  tabla_binaria[condiciones],
  as.logical
)

##################################################
# UPSET PLOT
##################################################

plot_upset <- upset(
  tabla_binaria,
  condiciones,
  width_ratio = 0.3,
  encode_sets = FALSE,
  sort_sets = FALSE,
  base_annotations = list(
    "Genes por intersección" =
      intersection_size(
        mapping = aes(
          fill = Condiciones
        ),
        text = list(
          size = 3,
          color = "black",
          vjust = -0.5
        )
      )
  ),
  set_sizes = upset_set_size()
) +
  scale_fill_brewer(
    palette = "Set2",
    name = "Condiciones"
  ) +
  theme_minimal(
    base_size = 12
  ) +
  theme(
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      size = 7
    ),
    axis.text.y = element_text(
      size = 10
    ),
    legend.position = "right"
  )

plot_upset

##################################################
# GUARDAR FIGURAS
##################################################

ggsave(
  "figures/upset_plot.pdf",
  plot_upset,
  width = 30,
  height = 14,
  units = "in"
)

ggsave(
  "figures/upset_plot.png",
  plot_upset,
  width = 30,
  height = 14,
  units = "in",
  dpi = 600
)

##################################################
# FIN
##################################################
