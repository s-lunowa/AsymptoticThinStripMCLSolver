# AsymptoticThinStripMCLSolver

Solver for the asymptotic model of immiscible two-phase flow with moving contact line in a thin strip.

The solver is writen in MATLAB. The model is derived in the article *On an averaged model for immiscible two-phase
flow with surface tension and dynamic contact angle in a thin strip* by Stephan B. Lunowa, Carina Bringedal and
Iuliu Sorin Pop (2021) ([UHasselt CMAT Preprint](http://www.uhasselt.be/Documents/CMAT/Preprints/2020/UP2006.pdf)).
<!-- TODO extend description -->

## Installation

No installation is necessary.
You can just download or pull the software.
The only requirement is that you have a current version of MATLAB&reg; installed (at least R2019a).
MATLAB&reg; is a product and trademark of MathWorks, Inc., which is available at <http://www.mathworks.com/products/matlab.html>.

## Usage

The file [example.m](example.m) shows in detail the usage of the different program components,
such as solving the ode and dae system, plotting the solutions and saving them.

The file [Model.m](Model.m) contains the model as a class with all routines necessary to simulate the process.
The scripts [solver_constant.m](solver_constant.m) and [solver_constricted.m](solver_constricted.m), [solver_dynamic.m](solver_dynamic.m),
and [solver_hysteretic.m](solver_hysteretic.m) contain examples of usage for uniform and constricted thin strips with
constant, dynamic and hysteretic contact angle model, respectively.

## License

Copyright &copy; 2021 Stephan B. Lunowa.

This software is licensed under the Creative Commons Attribution 4.0 International License.
You should have obtained a [LICENSE](LICENSE) file alongside this file.
To view a copy of this license, visit <http://creativecommons.org/licenses/by/4.0/>.
This licence allows for usage, modification and redistribution of the software, but EXCLUDES ANY WARRANTY AND LIABILITY.
