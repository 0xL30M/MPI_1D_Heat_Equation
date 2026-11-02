.DEFAULT_GOAL := all

EXEC = ex01.app
CC = mpicc
CPPFLAGS =
CFLAGS = -Wall -O2 -MMD -MP
LDFLAGS =
LDLIBS = -lm
OBJ_DIR = obj

include student.mk

OBJS = $(OBJ_DIR)/main.o $(OBJ_DIR)/student.o
DEPS = $(OBJS:.o=.d)

all: $(EXEC)

$(EXEC): $(OBJS) Makefile student.mk
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(LDLIBS)

$(OBJ_DIR)/main.o: src/main.c src/utils_1d.h Makefile student.mk
	@mkdir -p $(OBJ_DIR)
	$(CC) $(CPPFLAGS) -DUSE_MPI -Isrc $(CFLAGS) -c $< -o $@

$(OBJ_DIR)/student.o: student.c src/utils_1d.h Makefile student.mk
	@mkdir -p $(OBJ_DIR)
	$(CC) $(CPPFLAGS) -DUSE_MPI -Isrc $(CFLAGS) -c $< -o $@

-include $(DEPS)

# FULL=1 enables all fixtures; MPIEXEC may override the MPI launcher.
test: ex01.app
	FULL="$(FULL)" python3 test_smoke.py

clean:
	rm -rf $(OBJ_DIR) $(EXEC)

.PHONY: all clean test
