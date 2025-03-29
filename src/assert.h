
#include <cstdlib>

// Error by the programmer
// Ex: Shader compilation error
#define LGL_PRGRM_ASSERT(check, ...) if(!(check)){ __VA_ARGS__; exit(EXIT_FAILURE); }

// Error by the user's system
// Ex: OpenGL version not supported
#define LGL_SYS_ASSERT(check, ...) if(!(check)){ __VA_ARGS__; exit(EXIT_FAILURE); }
