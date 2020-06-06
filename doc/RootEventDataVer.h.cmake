#ifndef ROOTEVENTDATAVER_H
#define ROOTEVENTDATAVER_H

// This is automatically generated file.
// Edit "cmake/RootEventDataVer.h.cmake" file.

#define ROOTLIBDIR "${ROOT_LIBRARY_DIR}"

#define ROOTEVENTDATA "${ROOTEVENTDATA}"

#define ROOTEVENTDATA_VERSION "${ROOTEVENTDATA_VERSION}"

#define ROOTEVENTDATA_VER(a,b,c) ((a*100) + (b*10) + (c))

#define ROOTEVENTDATA_VERSION_NUMERIC \
        ROOTEVENTDATA_VER(${ROOTEVENTDATA_MAJOR_VERS},\
                          ${ROOTEVENTDATA_MINOR_VERS},\
                          ${ROOTEVENTDATA_PATCH_VERS})

#endif // ROOTEVENTDATAVER_H
