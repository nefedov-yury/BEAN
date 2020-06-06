#ifndef _Shared_DLLDEFINES_H_
#define _Shared_DLLDEFINES_H_

/*
  Cmake defines the variable "Library_EXPORTS" on Windows
  when it configures to build a shared library with name "Library".
*/

#if defined (_WIN32) 
    // We are using the Visual Studio Compiler and building Shared libraries
    #if defined(BeanUser_EXPORTS)
        #define  BeanUserShared_EXPORT __declspec(dllexport)
    #else
        #define  BeanUserShared_EXPORT __declspec(dllimport)
    #endif /* BeanUser_EXPORTS */
#else  /* UNIX */
    #define BeanUserShared_EXPORT
#endif /* _WIN32 */

#endif /* _Shared_DLLDEFINES_H_ */
