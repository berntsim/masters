This is the code used for a masters thesis trying to simulate the mixing of
black carbon aerosols and mineral dust particles. The code uses diffusion-limited cluster aggregation and multifractal velocity fields to accomplish this. It is licensed under the MiT licence (https://opensource.org/licenses/MIT). In the event that someone finds this code useful, please contact me and let me know, and I'll be happy to help in any way possible!

The code is structured relatively poorly, just to make uploading it to the cluster on which simulations were ran easy. This is also partly due to the inexperience of the author. The important aspects are however:

main:
This is where all the parameters are defined, and is in general a bit messy. This is due to the fact that the code was among other things tested visually, by utilising the SFML library to see that the behaviour of the particles were satisfactory. Using SFML necessitates that everything runs in one thread, which is called from main. No further alterations were made after the tests were done.

routines:
This is simply where most of the routines are defined. It is also a bit messy, and may contain some out-of-date routines which are not actually in use in the final code. Caution is advised.

containers:
This is where all the structures and containers used in the code are defined. This includes member functions.

To run the code:
to actually run the code should be pretty straight forward. Just set the parameters are you wish in the main.cpp file, and compile using whatever compiler you would like. It should not require any non-standard libraries to run this code. 