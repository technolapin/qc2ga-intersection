# Find Libqc2ga
#
#  QC2GA_FOUND        - True if Libqc2ga found.
#  QC2GA_INCLUDE_DIRS - where to find Libqc2ga include files
#  QC2GA_LIBRARIES    - where to find Libqc2ga binary libraries


# find include path
find_path (QC2GA_INCLUDE_DIR 
		NAMES qc2ga
		HINTS /usr/local/include)

# find library file
find_library (QC2GA_LIBRARY 
		NAMES qc2ga
		HINTS /usr/local/lib/ )


set(QC2GA_LIBRARIES ${QC2GA_LIBRARY})
set(QC2GA_INCLUDE_DIRS ${QC2GA_INCLUDE_DIR})

# handle the QUIETLY and REQUIRED arguments and set QC2GA_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (QC2GA DEFAULT_MSG 
                                   QC2GA_LIBRARIES QC2GA_INCLUDE_DIRS)

mark_as_advanced(QC2GA_INCLUDE_DIR QC2GA_LIBRARIE)





