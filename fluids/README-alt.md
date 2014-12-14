3D Fluid Solver
Kevin Yeh (kky226)
====================

This program is equipped with basic fluid interactions. The implementation is based off of the Stable Fluids paper by Jos Stam.

The fluids are rendered using GL_POINTS at a large radius, so it will only be visible at certain angles. No true voxel rendering is provided at this time.

On startup, the user can click and drag on the screen to generate a very dense fluid particle at (x,y,gridSize/2). This will diffuse rapidly when the simulation is running to create a large stream of fluid, so very small user interactions go a long way. When the fluids drag force checkbox is enabled, mouse dragging will instead apply external forces on the system, again at (x,y,gridSize/2). This is much more noticeable if gravity is turned off.

The bounds will keep fluid density in, as long as the gravity is very small.

Tunable Parameters: Diffusion Rate, Timestep, Gravity, User Drag Force
