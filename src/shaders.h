#pragma once

#include <glad/glad.h>

namespace lgl{

void compileShader(GLuint shader, const GLchar* source);

void debugShader(GLuint shader, const char* fileLoc);	

void debugShaderProgram(GLuint shaderProgram);


[[nodiscard]]
GLuint shaderFromVertFrag(const GLchar* vertSrc, const GLchar* fragSrc);

}
