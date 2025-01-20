# Schrodinger Equation in a Well

This project simulates the Schrodinger equation in a potential well using the Crank-Nicolson method. The simulation is implemented in C and visualized using Python.

## Files

- `mides.txt`: Contains the dimensions of the simulation grid.
- `gaussT.txt`: Output file for the simulation results.
- `README.md`: This file.
- `cranck_nicolson.c`: Main C program implementing the Crank-Nicolson method for solving the Schrodinger equation.
- `plot CN.py`: Python script for visualizing the simulation results.

## Usage

1. **Compile the C Program**:
    Use the provided VSCode task to compile the `cranck_nicolson.c` file.

2. **Run the Simulation**:
    Execute the compiled program to generate the `mides.txt` and `gaussT.txt` files.

3. **Visualize the Results**:
    Run the `plot CN.py` script to generate an animation of the simulation results.

## Requirements

- GCC compiler
- Python 3.x
- NumPy
- Matplotlib

## How to Run

1. **Compile and Run the C Program**:
    ```sh
    gcc -o cranck_nicolson cranck_nicolson.c -lm
    ./cranck_nicolson
    ```

2. **Run the Python Visualization**:
    ```sh
    python plot\ CN.py
    ```

## License

This project is licensed under the MIT License.
