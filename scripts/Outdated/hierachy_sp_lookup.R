

library(dplyr)

# 1. Create FIA species code lookup table
fia_lookup <- tibble::tribble(
  ~sp_code, ~Species,
  "ACSA3C", "Acer saccharum",
  "ACRU",   "Acer rubrum",
  "ACNE2",  "Acer negundo",
  "AEFL",   "Aesculus flava",
  "AIAL",   "Ailanthus altissima",
  "ALJU",   "Alnus jorullensis",
  "AMAR3",  "Acer macrophyllum",
  "ARSP2",  "Aronia sp.",
  "ASTR",   "Asimina triloba",
  "CAAM2",  "Carya amara",
  "CACA18", "Carya carolinae-septentrionalis",
  "CACA38", "Carya cahawbaensis",
  "CACO15", "Carya cordiformis",
  "CAGL8",  "Carya glabra",
  "CALA21", "Carya ovata",
  "CARYA",  "Carya <d7> hybrid",
  "CAOV2",  "Carya ovata",
  "CAOV3",  "Carya ovata var. ovata",
  "CAPA24", "Capsicum annum",
  "CATO6",  "Carya tomentosa",
  "CECA4",  "Celtis laevigata",
  "CELA",   "Celtis laevigata var. laevigata",
  "COFL2",  "Cornus florida",
  "COFUZ",  "Cornus florida <d7> nuttallii",
  "COOB2",  "Cornus obliqua",
  "CRCR2",  "Calocedrus decurrens",
  "DIVI5",  "Diospyros virginiana",
  "EUAM9",  "Euonymus americana",
  "EUAT5",  "Euonymus atropurpureus",
  "FAGR",   "Fagus grandifolia",
  "FOLI",   "Fothergilla major",
  "FRAMCO", "Fraxinus americana / pennsylvanica",
  "FRCA13", "Fraxinus caroliniana",
  "FRQU",   "Fraxinus quadrangulata",
  "GLTR",   "Gleditsia triacanthos",
  "HAVI4",  "Hamamelis virginiana",
  "HYFR",   "Hydrangea quercifolia",
  "ILDE",   "Ilex decidua",
  "ILEX",   "Ilex spp.",
  "ILLO",   "Ilex opaca",
  "ILOP",   "Ilex opaca",
  "JUNI",   "Juniperus virginiana",
  "JUVI",   "Juniperus virginiana",
  "LIBE3",  "Liquidambar styraciflua",
  "LISI",   "Liquidambar styraciflua",
  "LIST2",  "Liquidambar styraciflua <d7> hybrid",
  "LITU",   "Liriodendron tulipifera",
  "MAAC",   "Magnolia acuminata",
  "MORU2",  "Morus rubra",
  "NYSY",   "Nyssa sylvatica",
  "OSVI",   "Ostrya virginiana",
  "OXAR",   "Oxydendrum arboreum",
  "PATO2",  "Platanus occidentalis",
  "PIEC2",  "Pinus echinata",
  "PIVI2",  "Pinus virginiana",
  "PLOC",   "Prunus serotina",
  "PRME",   "Prunus mexicana",
  "PRSES",  "Prunus serotina",
  "QUAL",   "Quercus alba",
  "QUMU",   "Quercus muehlenbergii",
  "QUFA",   "Quercus falcata",
  "QUMO4",  "Quercus montana",
  "QUNI",   "Quercus nigra",
  "QUPA5",  "Quercus pagoda",
  "QURU",   "Quercus rubra",
  "QUSH",   "Quercus shumardii",
  "QUVE",   "Quercus velutina",
  "RHCO",   "Rhus copallinum",
  "ROPS",   "Robinia pseudoacacia",
  "SAAL5",  "Salix alba",
  "SANIC4", "Sapindus saponaria",
  "SILY",   "Salix lyallii",
  "STTR",   "Styrax grandifolius",
  "TIAM",   "Tilia americana",
  "ULAL",   "Ulmus alata",
  "ULAM",   "Ulmus americana",
  "ULMUS",  "Ulmus sp.",
  "ULRU",   "Ulmus rubra",
  "ULSE",   "Ulmus serotina",
  "VIRU",   "Vitis rotundifolia"
)

# 2. Keep only FIA codes that actually appear in census$sp
fia_lookup <- fia_lookup %>%
  filter(sp_code %in% census$sp)

# 3. Join FIA codes onto hierarchy table
heirachy_fia <- heirachy %>%
  left_join(fia_lookup, by = "Species")

# 4. Remove rows with no matched FIA code
heirachy_fia <- heirachy_fia %>%
  filter(!is.na(sp_code))

# Result
write.csv(heirachy_fia, "SOFOR/big_data/heirchy table prfdp.csv")


