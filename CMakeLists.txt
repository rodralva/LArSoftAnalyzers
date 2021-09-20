art_make(
   LIB_LIBRARIES
    ${FHICLCPP}
 
  MODULE_LIBRARIES
    ${ART_FRAMEWORK_CORE}
    ${ART_FRAMEWORK_PRINCIPAL}
    ${ART_ROOT_IO_TFILE_SUPPORT}
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_ROOT_IO_TFILESERVICE_SERVICE}
    ${FHICLCPP}
    ${ROOT_BASIC_LIB_LIST}
    canvas
    art_Persistency_Common
    art_Persistency_Provenance
    art_Utilities
    nusimdata_SimulationBase
    sbndcode_OpDetSim
)
install_fhicl()
install_source()
FILE(GLOB fcl_files *.fcl)

