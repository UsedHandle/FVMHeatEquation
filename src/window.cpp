#include "window.h"
#include <GLFW/glfw3.h>

#include "assert.h"

namespace lgl{

GLFW::GLFW(){
  LGL_SYS_ASSERT(glfwInit(),
    std::fprintf(stderr, "Could not initialize GLFW\n")
  )
}

GLFW::~GLFW(){
  glfwTerminate();
}

Window::Window() : m_window(nullptr){}

Window::Window(const WinProp& new_prop){
  m_window = glfwCreateWindow( new_prop.width, new_prop.height,
                               new_prop.title.c_str(),
                               new_prop.monitor,
                               nullptr );

  LGL_SYS_ASSERT(m_window,
    glfwTerminate(),
    std::fprintf(stderr, "Failed to create a window with GLFW\n")
  )

  glfwSetWindowUserPointer(m_window, reinterpret_cast<void*>(&m_prop));

  glfwSetErrorCallback([](int error, const char* description){
    LGL_SYS_ASSERT(false,
      std::fprintf(stderr, "GLFW ERROR (%d):\n%s\n", error, description);
    )
  });

  glfwSetWindowSizeCallback(m_window, [](GLFWwindow* window, int width, int height){
    WinProp* win_prop =
      reinterpret_cast<WinProp*>(glfwGetWindowUserPointer(window));
    win_prop->width = width;
    win_prop->height = height;
  });

  glfwMakeContextCurrent(m_window);

  LGL_SYS_ASSERT(
    gladLoadGLLoader(reinterpret_cast<GLADloadproc>(glfwGetProcAddress)),
    std::fprintf(stderr, "Failed to initialize GLAD")
  )

  glfwSetFramebufferSizeCallback(m_window, [](GLFWwindow* window, int width, int height){
    glViewport(0, 0, width, height);
  });
 
  // On MacOS the framebuffer resolution is not the same as the window size
  {
    int fbWidth, fbHeight;
    glfwGetFramebufferSize(m_window, &fbWidth, &fbHeight);
    glViewport(0, 0, fbWidth, fbHeight);
  }

}

Window::Window(Window&& b) :
  m_prop(std::move(b.m_prop)),
  m_window(b.m_window)
{
  b.m_window = nullptr;  
}

Window& Window::operator=(Window&& b){
  if(this == &b)
    return *this;

  glfwDestroyWindow(m_window); 
  std::swap(m_prop, b.m_prop);
  std::swap(m_window, b.m_window);

  return *this;
}

Window::~Window(){
  glfwDestroyWindow(m_window);
}

void Window::makeContextCurrent() const{
  glfwMakeContextCurrent(m_window);
}

bool Window::shouldClose() const{
  return glfwWindowShouldClose(m_window); 
}

void Window::onUpdate() const{
    glfwPollEvents();
    glfwSwapBuffers(m_window);
}

void Window::setSize(int width, int height){
  m_prop.width = width;
  m_prop.height = height;
  glfwSetWindowSize(m_window, width, height);
}

void Window::setTitle(const std::string& title){
  m_prop.title = title;
  glfwSetWindowTitle(m_window, title.c_str());
}

void Window::setMonitor( GLFWmonitor* monitor,
                         int xpos, int ypos,
                         int width, int height,
                         int refreshRate )
{
  m_prop.monitor = monitor;
  glfwSetWindowMonitor( m_window,
                        monitor,
                        xpos, ypos,
                        width, height,
                        refreshRate );
}

}
