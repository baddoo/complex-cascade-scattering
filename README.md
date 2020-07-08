# Complex cascade scattering

This repository contains the codes used to produce the results in the paper
["Acoustic scattering by cascades with complex boundary conditions: compliance, porosity and impedance"](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/acoustic-scattering-by-cascades-with-complex-boundary-conditions-compliance-porosity-and-impedance/51EEDE5E2F82DEC45F38A4A85B37B8AB), Journal of Fluid Mechanics, 2020.

The paper is concerned with calculating the scattering by an acoustic or vortical wave interacting with a cascade of blades with complex boundary conditions. The boundary conditions can model a range of physical phenomena, including compliance, porosity, and impedance. 
In particular, we derive exact solutions using the Wiener--Hopf method.
Equipped with these solutions, it is now straightforward to calculate important aeroacoustic quantities, as illustrated below:

<center>

Pressure field             |  Chord-wise velocity field
:-------------------------:|:-------------------------:
![](animations/totalPressure.gif?raw=true)  |  ![](animations/totalHVelocity.gif?raw=true)

</center>

Note the vortex shedding in the chord-wise velocity field.

In order to use the code, simply clone the repository onto your machine with:

```
git clone https://github.com/baddoo/complex-cascade-scattering.git
```

Make sure that all the subfolders are in the MATLAB search path by using, for example, the genpath command.

To get started try playing around with the `testFields.m` file inside the `tests` folder. This should give you an idea of what the code can do.

Feel free to get in touch if you have any questions or identify any bugs.


