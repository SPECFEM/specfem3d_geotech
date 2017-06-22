# build documentation
# REVISION:
#   HNG, Sep 06,2011;
IF(ISDOCUMENTATION)

  # copy all LaTeX source and required files
  FILE(COPY ${SDIR}/doc DESTINATION ${CMAKE_BINARY_DIR})

  # compile documenent
  # pdflatex1
  ADD_CUSTOM_COMMAND(
    OUTPUT    manual_SPECFEM3D_GEOTECH.aux
    DEPENDS   manual_SPECFEM3D_GEOTECH.tex
    COMMAND   ${PDFLATEX_COMPILER}
    ARGS      -interaction=batchmode manual_SPECFEM3D_GEOTECH
    COMMENT   "pdflatex1"
  )

  # bibtex
  ADD_CUSTOM_COMMAND(
  OUTPUT    manual_SPECFEM3D_GEOTECH.bbl
  DEPENDS   manual_SPECFEM3D_GEOTECH.aux
  COMMAND   ${BIBTEX_COMPILER}
  ARGS      -terse manual_SPECFEM3D_GEOTECH
  COMMENT   "bibtex"
  )

  # pdflatex2
  ADD_CUSTOM_COMMAND(
    OUTPUT    manual_SPECFEM3D_GEOTECH.pdf
    DEPENDS   manual_SPECFEM3D_GEOTECH.bbl
    COMMAND   ${PDFLATEX_COMPILER}
    ARGS      -interaction=batchmode manual_SPECFEM3D_GEOTECH

    COMMENT   "pdflatex2"
  )

  # pdflatex3
  ADD_CUSTOM_COMMAND(
    OUTPUT    manual_SPECFEM3D_GEOTECH.log
    DEPENDS   manual_SPECFEM3D_GEOTECH.bbl
              manual_SPECFEM3D_GEOTECH.pdf
    COMMAND   ${PDFLATEX_COMPILER}
    ARGS      -interaction=batchmode manual_SPECFEM3D_GEOTECH

    COMMENT   "pdflatex3"
  )

  ADD_CUSTOM_TARGET(Documentation ALL echo
    DEPENDS   manual_SPECFEM3D_GEOTECH.log
  )
ENDIF(ISDOCUMENTATION)

