#pragma once

#include <string>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <array>

#include "shaders.h"
#include "window.h"

#include <glm/glm.hpp>


class App{
  lgl::GLFW glfwInit;
  lgl::Window m_window;
  
  const int m_width = 800;
  const int m_height = 600;
  
  static constexpr std::size_t nx = 100;
  static constexpr std::size_t ny = 100;
 
  static constexpr std::size_t nt = 100;

  static constexpr float dx = 2.f/((float)nx - 1.f);
  static constexpr float dy = 2.f/((float)ny - 1.f);
  
  static constexpr float dt  = 0.0001f;
  
  std::array<std::array<float, nx>, ny> T{};
 

  struct Vertex{
    glm::vec2 coord;
    float T;
  };

  // width, height, triangles per cell,
  // vertices per triangle, coord and U per vertice
  std::array<std::array<Vertex, nx>, ny> vertices;
  std::array<unsigned int, (nx-1)*(ny-1)*2*3> indices;
  GLuint m_VAO, m_VBO, m_EBO;
  GLuint m_rectangleShader;

public:
 ~App();

  void setup();
  void run();

  void iterate(auto& T, float dt, float dx, float dy);

  void build_up_b(auto& b, const auto& U,
                   float rho, float dt,
                   float dx, float dy);

  void pressure_poisson(auto& p, const auto& b,
                        float dx, float dy,
                        std::size_t nit);

  void cavity_flow(auto& U, auto& p,
                   float dt, float dx, float dy,
                   float rho, float nu,
                   std::size_t nt, std::size_t nit);
};
