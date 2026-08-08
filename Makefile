# ProteinGo: HP Lattice Protein Folding Prototyping Toolkit
# Copyright (c) 2026 ProteinGo. Licensed under the MIT License.

CC = g++
CFLAGS = -O3 -Wall -Wextra -Iinclude
SRC = src/core/protein.c \
      src/analytics/energy.c \
      src/analytics/metrics.c \
      src/solvers/branch_bound.c \
      src/solvers/simulated_annealing.c \
      src/io/ascii_render.c \
      src/io/pdb_exporter.c \
      src/main.c

OBJ = $(SRC:.c=.o)
TARGET = hp_protein_solver

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(SRC) -o $(TARGET)

test: tests/test_energy.c $(SRC)
	$(CC) $(CFLAGS) tests/test_energy.c src/core/protein.c src/analytics/energy.c src/analytics/metrics.c -o test_runner
	./test_runner

clean:
	rm -f $(TARGET) test_runner *.pdb *.o
