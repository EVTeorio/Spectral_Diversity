

library(leaflet)
library(sf)
points_circles <- st_transform(points_circles, 4326)
# create a color palette based on species
pal <- colorFactor(palette = "Set2", domain = points_circles$sp)

leaflet(points_circles) %>%
  addTiles() %>%  # base map
  addPolygons(
    fillColor = ~pal(sp),
    fillOpacity = 0.7,
    color = NA,
    popup = ~paste(
      "Species:", sp,
      "<br>DBH:", DBH.2024,
      "<br>Crown Position:", crown.position
    )
  ) %>%
  addLegend(
    "topright",
    pal = pal,
    values = ~sp,
    title = "Species",
    opacity = 1
  )
