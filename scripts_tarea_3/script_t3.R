########################### TAREA 3: PHYLOSEQ ##################################

# Carga de librerías requeridas para realizar la actividad

library(phyloseq)
library(microbiome)
library(ggplot2)
library(vegan)
library(dplyr)

# Obtención del objeto phyloseq 

data("dietswap", package = "microbiome")
obj_psq <- dietswap
obj_psq 

## Identificación de variables

sample_variables(obj_psq)

########################### Curvas de rarefracción ############################# 

## Genarar la matriz de abundancias: obtenido de "https://github.com/joey711/phyloseq/issues/1641"

taxa_are_rows(obj_psq)

mat <- as(t(otu_table(obj_psq)), "matrix")

class(mat)

raremax <- min(rowSums(mat))

system.time(rarecurve(mat, step = 100, sample = raremax, col = "blue", label = FALSE))

## Como la curva de rarefracción de ve muy saturada podemos usar esta alternativa para ver muestras específicas:

system.time(rarecurve(mat [1:2, ], step = 100, sample = raremax, col = "blue", label = FALSE))


############################## DIVERSIDAD ALFA ################################

div_alfa_psq <- plot_richness(obj_psq, measures = c("Observed", "Shannon", "Simpson"),
              color = "group") +theme(
                axis.text.x = element_text(angle = 90, hjust = 1, size = 2)
              )

div_alfa_psq 

# Filtrado de datos

psq_filtrado <- obj_psq %>%
  transform_sample_counts(function(x) x/sum(x)) %>% # Normalizar una vez
  filter_taxa(function(x) sum(x > 0.001) >= (0.1 * nsamples(obj_psq)), TRUE)
psq_filtrado

############################## Diversidad Beta #################################

# Calcular la ordenación PCoA usando distancia Bray-Curtis

distancia_bray <- distance(obj_psq, method = "bray")

orden_pcoa <- ordinate(obj_psq, method = "PCoA", distance = distancia_bray)

# Visualizar la ordenación (en este caso para nationality, sex y group que fueron las que elegí)

ordn_nat <- plot_ordination(obj_psq, orden_pcoa, color = "nationality") +
  geom_point(size = 2) 
ordn_nat

ordn_sex <- plot_ordination(obj_psq, ord_bc, color = "sex") +
  geom_point(size = 2) 

ordn_group <- plot_ordination(obj_psq, ord_bc, color = "group") +
  geom_point(size = 2) 
ordn_group

################# Grafico de rango abundancia ##################################

# Calculamos la abundancia total de cada taxón y ordenamos de mayor a menor
abundancia_total <- taxa_sums(obj_psq)
rank_abund <- sort(abundancia_total, decreasing = TRUE)

# Creamos un data frame para poder usarlo en ggplot
df_rank <- data.frame(
  Rank = 1:length(rank_abund),
  Abundance = rank_abund
)

# Generamos el gráfico en ggplot2 

rank_psq <- ggplot(df_rank, aes(x = Rank, y = Abundance)) +
  geom_line(color = "darkred") +
  scale_y_log10() +
  labs(title = "Curva de Rank-Abundance", y = "Abundancia (escala log10)", x = "Rango") 

rank_psq

##################### Gráficas apiladas de abundancia por taxón ################

graf_apilada <- plot_bar(psq_filtrado, fill = "Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") +  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 5)
  )

graf_apilada
####################### Exportación de resultados ##############################



########################### GLOBALPATTERNS #####################################

data("GlobalPatterns")
gp <- GlobalPatterns
gp  

# Primero realizamos el procesamiento de datos indicado ()

#1. Filtrar taxa con menos de 5 lecturas en al menos 20% de las muestras

gp_filtrado <- filter_taxa(gp, function(x) sum(x >= 5) >= 0.2 * nsamples(gp), prune = TRUE)

# 2. Aglomerar a nivel de Familia

gp_familia <- tax_glom(gp_filtrado, taxrank = "Family")

# Subset para incluir solo muestras de: Soil, Feces, Skin

gp_subset <- subset_samples(gp_familia, SampleType %in% c("Soil", "Feces", "Skin"))

gp_subset

########################## DIVERSIDAD ALFA ####################################

# Calculo de los 3 indices de diverssidad alfa

diversidad_alfa <- estimate_richness(gp_subset, measures = c("Shannon", "Simpson", "Observed"))
diversidad_alfa

diversidad_alfa$SampleType <- sample_data(gp_subset)$SampleType
head(diversidad_alfa)

# Plots comparativos

p_gp <- plot_richness(gp_subset, x = "SampleType", measures = c("Observed", "Shannon", "Simpson"), color = "SampleType")
p_gp

# Pruebas estadísticas 

kruskal_observados <- kruskal.test(Observed ~ SampleType, data = diversidad_alfa)
kruskal_shannon <- kruskal.test(Shannon ~ SampleType, data = diversidad_alfa)
kruskal_simpson <- kruskal.test(Simpson ~ SampleType, data = diversidad_alfa)


kruskal_observados
kruskal_shannon
kruskal_simpson

######################### Curvas de rango abundanncia ##########################

# Calcular la abundancia relativa
gp_abund_rel <- transform_sample_counts(gp_subset, function(x) x / sum(x))

# Convertir el objeto phyloseq a un data frame
df_abund <- psmelt(gp_abund_rel)

# Ordenar los taxones por abundancia 
df_abund <- df_abund %>%
  group_by(SampleType, OTU) %>%
  summarise(Abundance = mean(Abundance)) %>%
  arrange(SampleType, desc(Abundance)) %>%
  ungroup()

# Asignar rangos
df_abund <- df_abund %>%
  group_by(SampleType) %>%
  mutate(Rank = row_number()) %>%
  ungroup()

ggplot(df_abund, aes(x = Rank, y = Abundance, color = SampleType)) +
  geom_line() +
  scale_y_log10() +
  labs(title = "Curvas de Rango-Abundancia por Tipo de Muestra",
       x = "Rango",
       y = "Abundancia Relativa (escala log10)") +
  theme_minimal()
