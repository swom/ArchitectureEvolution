# Entry point for user

HEADERS +=  \
  Stopwatch.hpp \
  environment.h \
  individual.h \
  network.h \
  observer.h \
  parser.h \
  population.h \
  rndutils.hpp \
  simulation.h \
  utilities.h


SOURCES +=  \
  environment.cpp \
  individual.cpp \
  main.cpp \
  network.cpp \
  observer.cpp \
  parser.cpp \
  population.cpp \
  simulation.cpp \
  utilities.cpp


CONFIG += c++17
QMAKE_CXXFLAGS += -std=c++17
CONFIG += resources_big

# High warning levels
# SFML goes bad with -Weffc++
QMAKE_CXXFLAGS += -Wall -Wextra -Wshadow -Wnon-virtual-dtor -pedantic

# A warning is an error
QMAKE_CXXFLAGS += -Werror

# Debug and release settings
CONFIG += debug_and_release
CONFIG(release, debug|release) {
  DEFINES += NDEBUG
}

# Qt5
QT += core gui

# SFML, default compiling
# GNU/Linux
unix:!macx {

  CONFIG(debug, debug|release) {
    # gcov
    QMAKE_CXXFLAGS += -fprofile-arcs -ftest-coverage
    LIBS += -lgcov
  }
}

win32{
  #LIBS += -lopenal32              #Dependency
  #LIBS += -lfreetype              #Dependency
  LIBS += -lopengl32              #Dependency
  LIBS += -lgdi32                 #Dependency
  LIBS += -lwinmm                 #Dependency
}




