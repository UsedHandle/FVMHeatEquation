#version 400 core
layout (location = 0) in vec2 inPos;
layout (location = 1) in float inT;

out float T;
/* uniform mat4 model; */
/* uniform mat4 view; */
/* uniform mat4 proj; */

void main(){
  gl_Position = vec4(inPos, 0.0, 1.0);
  T = inT;
}
