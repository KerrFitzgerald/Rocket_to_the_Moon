# Rocket_to_the_Moon
Calculate the positions of a rocket travelling between the Earth and Moon using the Runge-Kutta method for solving differential equations of motion.
Trajectories are defined in list format.
The first list item is a string defining the starting location of the rocket. Current options are:

 - 'L1' - Langrangian Point 1

 - 'L2' - Langrangian Point 2

 - 'Earth_Orbit' - Orbit above the Earth

The second list item is a list defining time array index values at which to apply a boost.

The third list item is a list defining boost magnitudes.

Example trajectories include:

'Approximate Free Return'

![Free_Return](https://user-images.githubusercontent.com/60627318/116542399-85755080-a8e4-11eb-9aba-47af7e7b4500.png)

'Lagrangian Point 1'

![Approximate L1](https://user-images.githubusercontent.com/60627318/116543063-4e536f00-a8e5-11eb-8618-49333afada4f.png)

'Lagrangian Point 2'

![Approximate L2](https://user-images.githubusercontent.com/60627318/116543453-c9b52080-a8e5-11eb-9d10-f06895828947.png)
