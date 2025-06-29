# ---------------------------------------------------
#  Dockerfile para desplegar tu Shiny App en Railway
# ---------------------------------------------------

# 1) Partimos de una imagen oficial de R + Shiny Server con R 4.4.3
FROM rocker/shiny:4.4.3

# 2) Instalamos aquí los paquetes de R que necesite tu app.
#    Ajusta esta lista (o usa renv/packrat si ya los gestionas de otro modo).
RUN R -e "install.packages(c('shiny','dplyr','ggplot2', 'AlphaSimR', 'ggridges','tidyr','stringr'), repos='https://cloud.r-project.org')"

# 3) Copiamos todo el código de tu carpeta Shiny (asumiendo que tu app 
#    está en la raíz del repo o en carpeta llamada 'masal'; ajústalo si 
#    cambia la ruta).
COPY . /srv/shiny-server/masal

# 4) Damos permisos al usuario shiny para que ejecute la app
RUN chown -R shiny:shiny /srv/shiny-server/masal

# 5) Exponemos el puerto 3838 (puerto por defecto de Shiny Server)
EXPOSE 3838

# 6) Arrancamos Shiny Server al iniciar el contenedor
CMD ["/usr/bin/shiny-server"]
