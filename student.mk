# Editable build settings. Uncomment or add assignments below.
# Standard options: CC, CPPFLAGS, CFLAGS, LDFLAGS, LDLIBS.
# Use += to extend defaults, or = to replace them.
# CFLAGS += -O2
# CPPFLAGS += -DMY_OPTION=1
# LDFLAGS += -L/path/to/libraries
# LDLIBS += -lmy_library
# Command-line assignments take precedence, e.g. make CFLAGS="-O2 -Wall".
CFLAGS += -O3 -flto -march=sapphirerapids -mprefer-vector-width=512 -ffast-math -funroll-loops