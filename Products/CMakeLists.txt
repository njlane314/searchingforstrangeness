art_make(
  NO_PLUGINS
  DICT_LIBRARIES
    canvas::canvas
    art::Persistency_Common
    art::Persistency_Provenance
    cetlib::cetlib
    cetlib_except::cetlib_except
)

install_headers()
install_source()
install_fhicl()
