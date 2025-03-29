#include "app.h"

#include <glad/glad.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "assert.h"

#include <array>
#include <cmath>
#include <algorithm>

namespace {
std::string getFile(const std::string& fileLoc){
  std::fstream fileStream(fileLoc);
  LGL_PRGRM_ASSERT(fileStream,
    std::fprintf(stderr, "%s not found\n", fileLoc.c_str())
  )

  std::stringstream stringStream;
  stringStream << fileStream.rdbuf();
  return stringStream.str();
}

template<typename Func>
void iterate2D(auto& A,  Func f){
  for(std::size_t j = 0; j < A.size(); ++j){
    for(std::size_t i = 0; i < A[0].size(); ++i){
      A[j][i] = f(j, i); 
    }
  }
}
template<typename Func>
void iterate2DVal(auto& A,  Func f){
  for(std::size_t j = 0; j < A.size(); ++j){
    for(std::size_t i = 0; i < A[0].size(); ++i){
      A[j][i] = f(j, i, A[j][i]); 
    }
  }
}

void setColumn(auto& A, size_t ival, auto val){
  for(std::size_t j = 0; j < A.size(); ++j){
    A[j][ival] = val;
  }
}

void setRow(auto& A, size_t jval, auto val){
  std::fill(A[jval].begin(), A[jval].end(), val);
}
}

App::~App(){
  glDeleteVertexArrays(1, &m_VAO);
  glDeleteBuffers(1, &m_VBO);
  glDeleteProgram(m_rectangleShader);
}

void App::setup(){
  // Sets minimum required version but uses latest
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef __APPLE__
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

  glfwWindowHint(GLFW_RESIZABLE, GLFW_FALSE);
  
  m_window = lgl::Window(lgl::WinProp{m_width, m_height, "Engine Test", nullptr});

  printf("%s\n", glGetString(GL_VERSION));

  iterate2D(T, [](std::size_t j, std::size_t i){
    const float x = ((float)i - 1.f)*dx;
    const float y = ((float)j - 1.f)*dy;
    if(x < 1.25f && x > 0.75f && y < 1.75f)
      return 273.f + 100.f;
    else
      return 273.f + 25.f;
    
  });
  
  for(int i = 0; i < 10; ++i)
    iterate(T, dt, dx, dy);




  
  for(std::size_t j = 0; j < ny; ++j){
    for(std::size_t i = 0; i < nx; ++i){
      using glm::vec2;
      const vec2 normalizedCoord = vec2(
        (float)i/(float)(nx-1), 
        (float)j/(float)(ny-1)
      );
      
      const vec2 coord = normalizedCoord*2.0f - vec2(1.f, 1.f);
      vertices[j][i] = Vertex{coord, T[j][i]};
 
    }
  }

  for(std::size_t j = 0; j < ny-1; ++j){
    for(std::size_t i = 0; i < nx-1; ++i){
      indices[6*(j*(nx-1)+i)  ]   = ((j)  *nx)+i;  
      indices[6*(j*(nx-1)+i)+1] = ((j+1)*nx)+i;  
      indices[6*(j*(nx-1)+i)+2] = ((j+1)*nx)+i+1;  
      
      indices[6*(j*(nx-1)+i)+3] = ((j)  *nx)+i;  
      indices[6*(j*(nx-1)+i)+4] = ((j)  *nx)+i+1;  
      indices[6*(j*(nx-1)+i)+5] = ((j+1)*nx)+i+1;  
    }
  }


  /* const float vertices[]{ */
  /*   -1.0f, -1.0f,  0.0f, */
  /*   -1.0f,  1.0f,  0.0f, */
  /*    1.0f, -1.0f,  0.0f, */

  /*   -1.0f,  1.0f,  0.0f, */
  /*    1.0f, -1.0f,  0.0f, */
  /*    1.0f,  1.0f,  0.0f */
  /* }; */

  glGenVertexArrays(1, &m_VAO);
  glGenBuffers(1, &m_VBO);
  glGenBuffers(1, &m_EBO);
   
  glBindVertexArray(m_VAO);
  glBindBuffer(GL_ARRAY_BUFFER, m_VBO);
  glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices.begin(), GL_STREAM_DRAW);
  
  glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), nullptr);
  glEnableVertexAttribArray(0);

  glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE,
      sizeof(Vertex), (void*)offsetof(Vertex, T));
  glEnableVertexAttribArray(1);

  
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_EBO);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices.begin(), GL_STREAM_DRAW);
  {
    
    std::string vertSrc = getFile("shaders/vertexshader.glsl");
    std::string fragSrc = getFile("shaders/fragmentshader.glsl");

    m_rectangleShader =
      lgl::shaderFromVertFrag(vertSrc.c_str(), fragSrc.c_str());

  }
}

void App::run(){
  /* glm::mat4 proj, view, model; */
  /* { */
  /*   float aspectrat = (float)m_width/(float)m_height; */
  /*   proj = glm::perspective(90.0f, aspectrat, 0.1f, 100.0f); */
  /* } */

  /* glm::vec3 pos = glm::vec3(0.0, 0.0, -5.0f); */
  /* view = glm::lookAt(pos, pos + glm::vec3(0.0f, 0.0f, 1.0f), glm::vec3(0.0f, 1.0f, 0.0f)); */
  /* model = glm::mat4(1.0f); */

  glUseProgram(m_rectangleShader);
  glBindVertexArray(m_VAO);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_EBO);

  /* GLint projLoc  =  glGetUniformLocation(m_rectangleShader, "proj"); */
  /* GLint viewLoc  =  glGetUniformLocation(m_rectangleShader, "view"); */
  /* GLint modelLoc =  glGetUniformLocation(m_rectangleShader, "model"); */
  
  /* std::array<std::array<float, nx>, ny> b; */
  // problem is nt
  //
  /* for(int it = 0; it < 400; ++it){ */
  /*   cavity_flow(U, p, dt, dx, dy, rho, nu, nt, nit); */
    
  /*   for(std::size_t j = 0; j < U.size(); ++j){ */
  /*     for(std::size_t i = 0; i < U[0].size(); ++i){ */
  /*       if(U[j][i].x > 0.5) */
  /*         printf("#"); */
  /*       else if(std::isnan(U[j][i].x)) */
  /*         printf("@"); */
  /*       else if(U[j][i].x < -20.0f) */
  /*         printf("*"); */
  /*       else */
  /*         printf("."); */
  /*     } */
  /*     printf("\n"); */
  /*   } */

  /*   printf("\n\n"); */
  /* } */
  /* build_up_b(p, U, rho, dt, dx, dy); */
  /* for(std::size_t j = 0; j < p.size(); ++j){ */
  /*   for(std::size_t i = 0; i < p[0].size(); ++i){ */
  /*     printf("%.1f ", p[j][i]); */
  /*   } */
  /*   printf("\n"); */
  /* } */
  
  glClearColor(0.0, 0.0, 0.0, 1.0);
  while(!m_window.shouldClose()){
    m_window.onUpdate();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    iterate(T, dt, dx, dy);
    glBindBuffer(GL_ARRAY_BUFFER, m_VBO);
    float* data = (float*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    size_t t = 0;
    for(std::size_t j = 0; j < ny; ++j){
      for(std::size_t i = 0; i < nx; ++i){
        using glm::vec2;
        const vec2 normalizedCoord = vec2(
          (float)i/(float)(nx-1), 
          (float)j/(float)(ny-1)
        );
          
        const vec2 coord = normalizedCoord*2.0f - vec2(1.f, 1.f);
        data[t] = coord.x;
        data[++t] = coord.y;
        data[++t] = T[j][i];
        ++t;
       }
    }
    glUnmapBuffer(GL_ARRAY_BUFFER);
    /* model = glm::rotate(model, 0.1f, glm::vec3(0.0f, 1.0f, 0.0f)); */
    /* glUniformMatrix4fv(projLoc, 1, GL_FALSE, glm::value_ptr(proj)); */
    /* glUniformMatrix4fv(viewLoc, 1, GL_FALSE, glm::value_ptr(view)); */
    /* glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model)); */
    glDrawElements(GL_TRIANGLES, 3*2*(nx-1)*(ny-1), GL_UNSIGNED_INT, 0);

  }
}

void App::iterate(auto& T, float dt, float dx, float dy){
  auto Tn = T;
  for(std::size_t j = 1; j < T.size()-1; ++j){
    for(std::size_t i = 1; i < T[0].size()-1; ++i){
      const float faceN = (Tn[j+1][i]-Tn[j][i])/dy * dx; 
      const float faceS = (Tn[j-1][i]-Tn[j][i])/dy * dx; 
      const float faceE = (Tn[j][i+1]-Tn[j][i])/dx * dy; 
      const float faceW = (Tn[j][i-1]-Tn[j][i])/dx * dy;
      
      const float dTdt = (faceN + faceS + faceE + faceW)/(dx*dy);
      T[j][i] = dt * dTdt + Tn[j][i]; 
      if(std::abs(T[j][i]) > 373.){
        printf("%f %f %f %f: ",
      (Tn[j+1][i]-Tn[j][i])/dy * dx, 
      (Tn[j-1][i]-Tn[j][i])/dy * dx, 
      (Tn[j][i+1]-Tn[j][i])/dx * dy, 
      (Tn[j][i-1]-Tn[j][i])/dx * dy);
        printf("%f\n", T[j][i]);
      }
      /* const float faceN =  (Tn[j][i]+Tn[j+1][i])/2.0f * dx; */
      /* const float faceS = -(Tn[j][i]+Tn[j-1][i])/2.0f * dx; */
      /* const float faceE =  (Tn[j][i]+Tn[j][i+1])/2.0f * dy; */
      /* const float faceW = -(Tn[j][i]+Tn[j][i-1])/2.0f * dy; */
      /* const float dudt = -(faceN+faceS+faceE+faceW)/(dx*dy); */
      /* T[j][i] = dt*dudt + Tn[j][i]; */
    }
  }
  setColumn(T, 0, 273.f + 25.f);
  setColumn(T, T.size()-1, 273.f + 25.f);
  setRow(T, 0, 273.f + 25.f);
  setRow(T, T[0].size()-1, 273.f + 25.f);
}

void App::build_up_b(auto& b, const auto& U,
                      float rho,
                      float dt, float dx, float dy)
{
  for(std::size_t j = 1; j < b.size()-1; ++j){
    for(std::size_t i = 1; i < b[0].size()-1; ++i){
      const float UuDx = (U[j][i+1].x - U[j][i-1].x)/(2.f*dx);
      const float UuDy = (U[j+1][i].x - U[j-1][i].x)/(2.f*dy);
      const float UvDy = (U[j+1][i].y - U[j-1][i].y)/(2.f*dy);
      const float UvDx = (U[j][i+1].y - U[j][i-1].y)/(2.f*dx);
      b[j][i] = rho * ((1.f/dt) * (UuDx + UvDy)
                       - (UuDx*UuDx)
                       - 2.f*UuDy*UvDx
                       - (UvDy*UvDy));
      
    }
  }
}

void App::pressure_poisson(auto& p, const auto& b,
                           float dx, float dy,
                           std::size_t nit)
{
  const float dx2 = dx*dx;
  const float dy2 = dy*dy;
  auto pn = p;
  for(std::size_t iter = 0; iter < nit; ++iter){
    pn = p;
    for(std::size_t j = 1; j < p.size()-1; ++j){
      for(std::size_t i = 1; i < p[0].size()-1; ++i){
        p[j][i] = ((pn[j][i+1]+pn[j][i-1])*dy2 + (pn[j+1][i]+pn[j-1][i])*dx2
                   -(dx2*dy2*b[j][i]))
                  /(2.f*(dx2+dy2));
                  
      }
    }
    p[0] = p[1];
    std::fill_n(p[p.size()-1].begin(), p[0].size(), 0.0f);
    for(std::size_t j = 1; j < p.size()-1; ++j){
      p[j][0] = p[j][1];
      p[j][p[0].size()-1] = p[j][p[0].size()-2];
    }
  }
}

void App::cavity_flow(auto& U, auto& p,
                      float dt, float dx, float dy,
                      float rho, float nu,
                      std::size_t nt, std::size_t nit)
{
  std::remove_reference_t<decltype(U)> Un; 
  std::remove_reference_t<decltype(p)> b;
  
  for(std::size_t n = 0; n < nt; ++n){
    Un = U;
    build_up_b(b, U, rho, dt, dx, dy);
    pressure_poisson(p, b, dx, dy, nit);
    for(std::size_t j = 0; j < U.size(); ++j){
      if(j == 0){
        std::fill_n(U[j].begin(), U[j].size(), glm::vec2(0.0f, 0.0f));
        continue;
      } else if(j == U.size()-1){
        std::fill_n(U[j].begin(), U[j].size(), glm::vec2(1.0f, 0.0f));
        continue;
      }
      for(std::size_t i = 0; i < U[0].size(); ++i){
        if(i == 0 || i == U[0].size()-1){
          U[j][i] = glm::vec2(0.0f, 0.0f);
          continue;
        }
        const float Unx = Un[j][i].x;
        const float UnxW = Un[j][i-1].x;
        const float UnxE = Un[j][i+1].x;
        const float UnxS = Un[j-1][i].x;
        const float UnxN = Un[j+1][i].x;

        const float Uny = Un[j][i].y;
        const float UnyW = Un[j][i-1].y;
        const float UnyE = Un[j][i+1].y;
        const float UnyS = Un[j-1][i].y;
        const float UnyN = Un[j+1][i].y;

        const float pW = p[j][i-1];
        const float pE = p[j][i+1];
        const float pN = p[j+1][i];
        const float pS = p[j-1][i];
        
        const float dx2 = dx*dx;
        const float dy2 = dy*dy;

        U[j][i].x = Unx + (-Unx*(Unx-UnxW)/dx - Uny*(Unx-UnxS)/dy
                        - (pE-pW)/(2.f*rho*dx)
                        + nu*((UnxE-(2.f*Unx)+UnxW)/dx2 +
                              (UnxN-(2.f*Unx)+UnxS)/dy2))*dt;
        U[j][i].y = Uny + (Unx*(Uny-UnyW)/dx - Uny*(Uny-UnyS)/dy
                        - (pN-pS)/(2.f*rho*dy)
                        + nu*((UnyE-(2.f*Uny)+UnyW)/dx2 +
                              (UnyN-(2.f*Uny)+UnyS)/dy2))*dt;
      }
    }
  }
}


