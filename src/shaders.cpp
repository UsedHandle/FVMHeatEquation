#include "shaders.h"

#include <cstdio>
#include <cstdlib>

#include "assert.h"

namespace lgl{

void compileShader(GLuint shader, const GLchar* source){
  glShaderSource(shader, 1, &source, nullptr);
  glCompileShader(shader);
}

void debugShader(GLuint shader, const char* fileLoc){	
	int  success;
	char infoLog[512];
	glGetShaderiv(shader, GL_COMPILE_STATUS, &success);

	LGL_PRGRM_ASSERT(success,
		glGetShaderInfoLog(shader, 512, nullptr, infoLog),
		std::fprintf(stderr, "%s\n%s\n\n", fileLoc, infoLog)
	)
}

void debugShaderProgram(GLuint shaderProgram){
	int  success;
	char infoLog[512];

	glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
	LGL_PRGRM_ASSERT(success,
		glGetProgramInfoLog(shaderProgram, 512, nullptr, infoLog),
		std::fprintf(stderr, "Shader Program Failed to link/create:\n%s\n", infoLog)
	)
}


[[nodiscard]]
GLuint shaderFromVertFrag(const GLchar* vertSrc, const GLchar* fragSrc){
	GLuint shaderProgram = glCreateProgram();
  GLuint vertShader = glCreateShader(GL_VERTEX_SHADER);
	compileShader(vertShader, vertSrc);
	debugShader(vertShader, vertSrc);
	
	GLuint fragShader = glCreateShader(GL_FRAGMENT_SHADER);
	compileShader(fragShader, fragSrc);
	debugShader(fragShader, fragSrc);

	glAttachShader(shaderProgram, vertShader);
	glAttachShader(shaderProgram, fragShader);

	glLinkProgram(shaderProgram);
	debugShaderProgram(shaderProgram);

	glDeleteShader(vertShader);
	glDeleteShader(fragShader);
	
	return shaderProgram;
}

}
