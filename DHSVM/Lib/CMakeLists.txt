# -------------------------------------------------------------
# file: CMakeLists.txt
# -------------------------------------------------------------
# -------------------------------------------------------------
# Battelle Memorial Institute
# Pacific Northwest Laboratory
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created February  1, 2019 by Perkins
# Last Change: 2019-02-01 12:38:29 d3g096
# -------------------------------------------------------------


set(CSSrcDir "${CMAKE_CURRENT_SOURCE_DIR}/../Source")
get_filename_component(CSSrcDir ${CSSrcDir} ABSOLUTE)
set(CSIncDir "${CMAKE_CURRENT_SOURCE_DIR}/../Include")
get_filename_component(CSIncDir ${CSIncDir} ABSOLUTE)
set(CSSrc0rig
  cs_add.c cs_amd.c cs_chol.c cs_cholsol.c cs_counts.c cs_cumsum.c 
  cs_droptol.c cs_dropzeros.c cs_dupl.c cs_entry.c 
  cs_etree.c cs_fkeep.c cs_gaxpy.c cs_happly.c cs_house.c cs_ipvec.c 
  cs_lsolve.c cs_ltsolve.c cs_lu.c cs_lusol.c cs_util.c cs_multiply.c 
  cs_permute.c cs_pinv.c cs_post.c cs_pvec.c cs_qr.c cs_qrsol.c 
  cs_scatter.c cs_schol.c cs_sqr.c cs_symperm.c cs_tdfs.c cs_malloc.c 
  cs_transpose.c cs_compress.c cs_usolve.c cs_utsolve.c cs_scc.c 
  cs_maxtrans.c cs_dmperm.c cs_updown.c cs_print.c cs_norm.c cs_load.c 
  cs_dfs.c cs_reach.c cs_spsolve.c cs_leaf.c cs_ereach.c cs_randperm.c
  )

set(CSSrc "")
foreach (f ${CSSrc0rig})
  list(APPEND CSSrc "${CSSrcDir}/${f}")
endforeach (f)

include_directories(BEFORE ${CSIncDir} )


add_library(cs_di_lib OBJECT ${CSSrc})

add_library(cs_dl_lib OBJECT ${CSSrc})
target_compile_definitions(cs_dl_lib PUBLIC CS_LONG=1)

add_library(cs_ci_lib OBJECT ${CSSrc})
target_compile_definitions(cs_ci_lib PUBLIC CS_COMPLEX=1)

add_library(cs_cl_lib OBJECT ${CSSrc})
target_compile_definitions(cs_cl_lib PUBLIC CS_LONG=1 CS_COMPLEX=1)

add_library(cxsparse
  $<TARGET_OBJECTS:cs_di_lib>
  $<TARGET_OBJECTS:cs_dl_lib>
  $<TARGET_OBJECTS:cs_ci_lib>
  $<TARGET_OBJECTS:cs_cl_lib>
)