#pragma once

#include <string>
#include <cstdio>
#include <cstdlib>

#include <glad/glad.h>
#include <GLFW/glfw3.h>


namespace lgl{

struct GLFW{
  GLFW();
  ~GLFW();

  GLFW(const GLFW&) = delete;
  GLFW& operator=(const GLFW) = delete;
};

struct WinProp{
  int width, height;
  std::string title;
  GLFWmonitor* monitor;
};

class Window{
  WinProp m_prop;

public:
  GLFWwindow* m_window; 

  Window();
  Window(const WinProp& new_prop);
  Window(const Window&) = delete;
  Window& operator=(const Window&) = delete;
  Window(Window&& b);
  Window& operator=(Window&& b);
  ~Window();

  void makeContextCurrent() const;
  bool shouldClose() const;
  void onUpdate() const;
  void setSize(int width, int height);
  void setTitle(const std::string& title);
  
  // When GLFWmonitor* is nullptr it is windowed
  // refreshRate can be GLFW_DONT CARE
  void setMonitor( GLFWmonitor* monitor,
                   int xpos, int ypos,
                   int width, int height,
                   int refreshRate );

};
}
