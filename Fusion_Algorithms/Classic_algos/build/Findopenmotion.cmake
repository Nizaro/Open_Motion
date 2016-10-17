# ===================================================================================
#  openmotion CMake configuration file
#
#             ** File generated automatically, do not modify **
#
#  Usage from an external project:
#    In your CMakeLists.txt, add these lines:
#
#    FIND_PACKAGE(openmotion REQUIRED )
#    TARGET_LINK_LIBRARIES(MY_TARGET_NAME )
#
#    This file will define the following variables:
#      - openmotion_LIBS          : The list of libraries to links against.
#      - openmotion_LIB_DIR       : The directory where lib files are. Calling LINK_DIRECTORIES
#                                with this path is NOT needed.
#      - openmotion_VERSION       : The  version of this PROJECT_NAME build. Example: "1.2.0"
#      - openmotion_VERSION_MAJOR : Major version part of VERSION. Example: "1"
#      - openmotion_VERSION_MINOR : Minor version part of VERSION. Example: "2"
#      - openmotion_VERSION_PATCH : Patch version part of VERSION. Example: "0"
#
# ===================================================================================
INCLUDE_DIRECTORIES("/usr/local/include")
SET(openmotion_INCLUDE_DIRS "/usr/local/include")

LINK_DIRECTORIES("/usr/local/lib")
SET(openmotion_LIB_DIR "/usr/local/lib")

SET(openmotion_LIBS  openmotion)

SET(openmotion_FOUND 1)
SET(openmotion_VERSION        1.0.0)
SET(openmotion_VERSION_MAJOR  1)
SET(openmotion_VERSION_MINOR  0)
SET(openmotion_VERSION_PATCH  0)
