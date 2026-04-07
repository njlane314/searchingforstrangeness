install_fhicl()

install(DIRECTORY data/
        DESTINATION job/reinteractions
        FILES_MATCHING
        PATTERN "*.root")
