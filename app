

library(shiny)
library(AlphaSimR)
library(dplyr)
library(ggplot2)
library(ggridges)
library(tidyr)
library(data.table)
library(stringr)

pob <- load("cruda188.RData")
pop <- get(pob[1])

SP <- SimParam$new(pop) # ge

# UI ----
ui <- fluidPage(
  titlePanel("Simulación de Mejora Genética - AlphaSimR"),
  sidebarLayout(
    sidebarPanel(
      tabsetPanel(
        type = "tabs",
        tabPanel("Hembras",
                 numericInput("n_hembras", "Número de hembras", 
                              value = 500, min = 2, 
                              max = 2000),
                # numericInput("n_partos", "N° de partos por hembra",                             value = 2, min = 1, max = 10),
                 sliderInput("edad_servicio_hembra", 
                             "Edad 1° servicio hembras (meses)", 
                             min = 18, max = 72, value = 18),
                 sliderInput("porc_preñez", "% preñez", value = 90, 
                             min = 60, max = 99),
                 sliderInput("porc_destete", "% destete", 
                             value = 90, min = 50, max = 99),
                 numericInput("repo_hembras", "% reposición hembras", 
                              value = 20, min = 0, max = 50),
                 selectInput("sel_hembras", "Criterio selección hembras",
                             c("rand", "gv", "pheno"))
        ),
        tabPanel("Machos",
                 sliderInput("porc_machos", "% machos", 
                             min = 0.1, max = 5, value = 3),
                 sliderInput("edad_servicio_macho", 
                             "Edad 1° servicio machos (meses)",
                             min = 18, max = 72, value = 22),
                 numericInput("repo_machos", "% reposición machos",
                              value = 25, min = 0, max = 100),
                 selectInput("sel_machos", "Criterio selección machos",
                             choices = c("rand", "gv", "pheno"))
        ),
        tabPanel("Simulación",
                 sliderInput("n_anos", "N° de años de simulación", 
                             min = 1, max = 15, value = 10),
                 numericInput("h2", "Heredabilidad del carácter", 
                              value = 0.3, min = 0.05, max = 1,
                              step = 0.05),
                 actionButton("simular", "Simular", class = "btn-primary")
        )
      )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Generación y Reposición",
                 plotOutput("plot_intervalo"),
                 plotOutput("plot_repo")
        ),
        tabPanel("Distribuciones",
                 plotOutput("plot_feno"),
                 plotOutput("plot_media"),
                 plotOutput("plot_promedios") # nueva pestaña para promedios
        )
      )
    )
  )
)

# Server ----
server <- function(input, output, session) {
  sim_result <- eventReactive(input$simular, {
    
    tryCatch({
    #pob <- load(file.path("H:/Mi unidad/garra/fundadoras/cruda188.RData"))
    #pop <- get(pob[1])
    
    SP <- SimParam$new(pop) # get(pob[2])
    SP$setSexes("yes_sys")
    SP$setTrackPed(TRUE)

    SP$resetPed()
    SP$addTraitA(3500, mean = 0, var = 1)
    pop <- newPop(pop, simParam = SP)
    pop <- setPheno( pop, h2 = input$h2, simParam = SP)
   
    # Inicializar
    df_summary <- tibble()
    df_medias <- tibble()
    df_medias_sexo <- tibble()
    # Definir tamaños
    npadres <- floor((input$repo_machos / 100) * input$n_hembras)
    nmadres <- input$n_hembras

    naños_macho <- floor(100/input$repo_machos)
    naños_hembras <- floor(100/input$repo_hembras)
    
    min_machos <- ceiling(input$edad_servicio_macho / 12)
    min_hembras <- ceiling(input$edad_servicio_hembra / 12)
    
    sex <- c(rep("M", floor((npadres * 2 * min_machos))),
             rep("F", (nInd(pop) - floor((npadres * 2 * min_machos))) ))
    pop@sex <- sample(sex, nInd(pop), replace = F)
    
    candidatos <- vector("list", length = input$n_anos + 1)
    ig <- 0.5 * (((input$edad_servicio_macho/12) + naños_macho) +
                   ((input$edad_servicio_hembra/12) + naños_hembras) ) 
    
    # Loop por años (generaciones)
    for (year in 1:input$n_anos) {
      # Selección padres según edad mínima
      if(year < min_machos) {
        
        if (year == 1){
          macho <-  pop[pop@sex == "M"]
          padre <- selectInd(macho, nInd = npadres , 
                             use = input$sel_machos,
                             sex = "M", selectTop = T)
          padre@misc$edad <- (-1) * sample(min_machos:(naños_macho+min_machos),
                                           size = nInd(padre),
                                           replace = TRUE)  
         
        } else { 
          
          machoid <- setdiff( macho@id, padre@id)
          
          papa <- padre[padre@misc$edad != min((padre@misc$edad))]
          length(papa)
          padrecito <- selectInd(macho, nInd = npadres - nInd(papa) , 
                                 use = input$sel_machos,
                                 sex = "M", selectTop = T,
                                 candidates = machoid)
          padrecito@misc$edad <- rep((-1) * year, nInd(padrecito)) 
          
          padre <- c(papa, padrecito)
          

        }
      } else{
        
        macho <- candidatos[[year]]
        papa <- padre[padre@misc$edad != min((padre@misc$edad))]
        
        padrecito <- selectInd(macho, nInd = npadres - nInd(papa),
                           use = input$sel_machos,
                           sex = "M", selectTop = T)
        
        padre <- c(papa, padrecito)
      }
      
      if(year < min_hembras) {
        
        hembra <- pop[pop@sex == "F"]
       
        if (year == 1){
          
          madre <- selectInd(hembra, nInd = input$n_hembras,
                             use = input$sel_hembras, 
                             sex = "F", selectTop = T)
          madre@misc$edad <- (-1) * sample(min_hembras:(naños_hembras + min_hembras),
                                           size = nInd(madre),
                                           replace = TRUE)  
        } else { 
          
          
          ########
          hembraid <- setdiff( hembra@id, madre@id)
          
          mama <- madre[madre@misc$edad != min((madre@misc$edad))]
          
          mamacita <- selectInd(hembra, 
                                nInd = input$n_hembras - nInd(mama), 
                                use = input$sel_hembras,
                                sex = "F", selectTop = T,
                                candidates = hembraid)
          mamacita@misc$edad <- rep((-1) * year, nInd(mamacita)) 
          
          madre <- c(mama, mamacita)
         
        }
        
      } else{
        
        hembra <- candidatos[[year]]
        mama <- madre[madre@misc$edad != min((madre@misc$edad))]
        
        mamacita <- selectInd(hembra, 
                           nInd = input$n_hembras - nInd(mama), 
                           use = input$sel_hembras,
                           sex = "F", selectTop = T)
        
        madre <- c(mama, mamacita)
      }
      # Cruzamientos
      ncruzas <- ceiling((input$n_hembras * input$porc_preñez / 100) * (input$porc_destete / 100))

      base <- randCross2(madre, padre, nCrosses = ncruzas, nProgeny = 1,
                         simParam = SP)
      base <- setPheno(base, h2 = input$h2, simParam = SP)
      base@misc$edad <- rep( year, nInd(base)) 

      candidatos[[year + 1]] <- base
      
      # Resúmenes
      genera <- ifelse(year == 0, 0, ceiling(ig / year))
      
      df_summary <- bind_rows(df_summary,
                              tibble(
                                sex = base@sex,
                                año = year,
                                generacion = rep(genera,nInd(base) ),
                                pheno = base@pheno,
                                media = meanP(base)
                              )
      )
      

      df_medias <- bind_rows(df_medias,
                              tibble(
                                mediaP = meanP(base),
                                mediaG = meanG(base),
                                varP = varP(base),
                                varG = varG(base),
                                generacion = genera,
                                año = year
                              )
      )
      
      df_medias_sexo <- bind_rows(df_medias_sexo,
                             tibble(
                               mediaPh = meanP(madre),
                               mediaPm = meanP(padre),
                               mediaGh = meanG(madre),
                               mediaGm = meanG(padre),
                               varPh = varP(madre),
                               varPm = varP(padre),
                               varGh = varG(madre),
                               varGm = varG(padre),
                               año = year,
                               generacion = genera
                             )
      )
     
      
    }
    list(df = df_summary, 
         medias = df_medias,
         media_sexo = df_medias_sexo)
    
    }, error = function(e) {
      # En caso de error, imprimo mensaje y hago return de NULL
      message("Error dentro de sim_result(): ", e$message)
      NULL
    })
  })
  
  # Histograma de intervalo generacional
  output$plot_intervalo <- renderPlot({
    #req(sim_result())
    naños_macho <- floor(100/input$repo_machos)
    naños_hembras <- floor(100/input$repo_hembras)
    
    ig_m <-((input$edad_servicio_macho/12) + naños_macho) 
    ig_h <-((input$edad_servicio_hembra/12) + naños_hembras) 
    ig <- 0.5 * (((input$edad_servicio_macho/12) + naños_macho) +
                 ((input$edad_servicio_hembra/12) + naños_hembras) ) 
    
    df_int <-   tibble(  intervalo = c(ig_h, ig, ig_m),
                         sex = factor(c("Hembra", "Promedio", "Macho"),
                                      levels = c("Macho", "Promedio", "Hembra"))
                         
    )
    
    ggplot(df_int, aes(x = sex, y = intervalo, 
                       fill = sex)) +
      geom_col(alpha = 0.6) +
      geom_text(aes(label = round(intervalo, 2)), vjust = -0.5) +
      scale_fill_manual(values = c(Macho = "#1f77b4", Promedio = "red",
                                   Hembra = "#ff7f0e")) +
      labs(
        x = "Sexo",
        y = "Intervalo (años)",
        fill = "Sexo",
        title = "Intervalo Generacional en machos, hembras y promedio"
      ) +
      theme_minimal()
   
  })
  
  # Disponibles vs Reposición
  output$plot_repo <- renderPlot({
    
    candidates <- input$n_hembras * (input$porc_preñez / 100) * (input$porc_destete / 100)
    df_repo <- tibble(
      grupo = c("Hembras", "Machos"),
      disponibles = c(
        #input$n_hembras * (input$porc_preñez / 100) * (input$porc_destete / 100),
        candidates * 0.5, 
        #input$n_hembras * (input$porc_machos / 100) * (input$porc_preñez / 100 * input$porc_destete / 100)
        candidates * 0.5
      ),
      reposicion = c(
        input$n_hembras * (input$repo_hembras / 100),
        input$n_hembras * (input$porc_machos / 100) * (input$repo_machos / 100)
      )
    ) %>% pivot_longer(c(disponibles, reposicion), names_to = "tipo", values_to = "cantidad")
    
    ggplot(df_repo, aes(x = grupo, y = cantidad, fill = tipo)) +
      geom_col(position = "dodge") +
      geom_text(aes(label = round(cantidad, 2)), vjust = -0.5) +
      labs(y = "Cantidad", title = "Disponibles vs Reposición") +
      theme_minimal()
  })
  
  # Densidades con color por sexo
  output$plot_feno <- renderPlot({
    req(sim_result())
    
    df <- sim_result()$df
    
    
    ggplot(df, aes(x = pheno, y = factor(año), fill = sex)) +
      geom_density_ridges(alpha = 0.6, scale = 1.2) +
      geom_vline( aes(xintercept = media),
                  data = df %>% filter(!is.na(año)) %>% slice(1),
                  linetype = "dashed", color = "red"
      ) +
      labs(x = "Valor fenotípico", y = "Año", 
           title = "Distribución fenotípica de la progenie por año y sexo") +
      theme_minimal()
  })
  
  # Promedio por generación
  output$plot_promedios <- renderPlot({
    req(sim_result())
    df <- sim_result()$medias
    
    ggplot(df, aes(x = año)) +
      ######
    geom_vline( aes(xintercept = generacion),
      data = df %>% filter(!is.na(generacion)) %>% slice(1),
      linetype = "dashed", color = "red"
    ) +
    ########
      #geom_vline(aes(xintercept = generacion, linetype = "dashed"), color = "red") +
      geom_line(aes(y = mediaP, color = "Media Fenotípica"), 
                linewidth = 2) +
      geom_point(aes(y = mediaP, color = "Media Fenotípica"),
                 linewidth = 2) +

      geom_line(aes(y = mediaG, color = "Media Genotípica"), 
                linewidth = 2) +
      geom_point(aes(y = mediaG, color = "Media Genotípica"), 
                 linewidth = 2) +
      
      scale_color_manual(
        name   = "", 
        values = c(
          "Media Fenotípica" = "#2ca02c", 
          "Media Genotípica" = "#d62728"
        )
      ) +
 
      labs(x = "Año", 
           y = "Valor medio", 
           title = "Promedio fenotípico y genotípico por año en progenie") +
      theme_minimal()
  })
 
  output$plot_media <- renderPlot({
    req(sim_result())
    df <- sim_result()$media_sexo
    
    ggplot(df, aes(x = año)) +
      geom_line(aes(y = mediaPh, color = "Fenotípica Hembra"), size = 1) +
      geom_point(aes(y = mediaPh, color = "Fenotípica Hembra"), size = 2) +
      geom_line(aes(y = mediaPm, color = "Fenotípica Macho"), size = 1) +
      geom_point(aes(y = mediaPm, color = "Fenotípica Macho"), size = 2) +
      geom_line(aes(y = mediaGh, color = "Genotípica Hembra"), size = 1) +
      geom_point(aes(y = mediaGh, color = "Genotípica Hembra"), size = 2) +
      geom_line(aes(y = mediaGm, color = "Genotípica Macho"), size = 1) +
      geom_point(aes(y = mediaGm, color = "Genotípica Macho"), size = 2) +
      geom_vline(
        aes(xintercept = generacion),
        data = df %>% filter(!is.na(generacion)) %>% slice(1),
        linetype = "dashed", color = "red"
      ) +
      scale_color_manual(
        name = "",
        values = c(
          "Fenotípica Hembra" = "#1f77b4",
          "Fenotípica Macho"  = "#ff7f0e",
          "Genotípica Hembra" = "purple",
          "Genotípica Macho"  = "maroon"
          
        )
      ) +
      labs(
        x     = "Generación (Año)",
        y     = "Media Fenotípica",
        title = "Fenotípica en progenitores (hembra y Macho) por generación"
      ) +
      theme_minimal() +
      theme(
        legend.position = "top",
        legend.title    = element_blank()
      )
  })
  
}

shinyApp(ui = ui, server = server)
