/* GLOBALS CONSTANTS */

//* QDot parameters
#define Q_Shape "Lens"  // available: Cyl, Trap, Cone & Lens
#define P_pos 2         // QDot Position at the structure
#define ho 70.0         // QD height
#define Ro 150.0        // QD Radious >> Major radious
#define ro 150.0        // QD Radious >> Minor radious
#define Rs (8.0*Ro)     // Estructure Radious

//* ....
#define HOLE  "si" 

//* maps:  to plot maps
//* energ: to plot energ. diag.
//* diag:  to compute energies

#define what2do "diag"
#define nthread 24

//* Substrate infos
#define a_subs a_GaAs

//* Diagonalization Parameters 
#define Nmax 20          //* Number modes in z-direction
#define Mmax 10          //* Number modes for Bessel's zeros
#define Lmax 3           //* Angular momenta L  
#define Nlevels DIM       //* Number of nivels desired

#define DIM Nmax*Mmax    //* Base dimension taken
#define K_coup 2         //* K states by L to compute OS
#define E_min 0.01       //* Energy range to
#define E_max 0.4        //* compute the OS
//****************************************************************
// Physics Constants

#define Ry 13.6          // Energy unit (eV)
#define ao 0.529         // Length unit (Anstroms)
#define hb 6.58211e-16   // Planck's constant
#define me 9.1009e-31    // Free electron mass (Kg)
#define c_luz 3.0e14     // Light velocity (um/s)
#define depl_L 5.547e-9
#define depl_e -0.6 * 0
//****************************************************************

//****************************************************************
//* BIAS

#define BIAS 1                 //* Turn On=1 (Off=0) the electric field
#define E_field (-1.6/230.0 * BIAS)   //* Electric field strength

//--- Potentials modified

#define type_V 0               //* type 1-> Linear 2-> Cuadratic

//--- Potential V1

#define V1_pos 20              //* potential position in the structure
#define x_con_a1 0.11
#define x_con_b1 0.20

//--- Potential V2

#define V2_pos 20             //* potential position in the structure
#define x_con_a2 0.11
#define x_con_b2 0.198

//**************************************************************************


/***********************
  Material parameters 
************************/

// Binaries...

#define Eg_GaAs 1.519
#define m_GaAs 0.067
#define mhh_100_GaAs 1.0/(g1_GaAs-2.0*g2_GaAs)
#define mhh_110_GaAs 2.0/(2.0*g1_GaAs-g2_GaAs-3.0*g3_GaAs)
#define g1_GaAs 6.98
#define g2_GaAs 2.06
#define g3_GaAs 2.93
#define VBO_GaAs -0.80
#define a_GaAs 5.65325
#define C11_GaAs 1221.0
#define C12_GaAs 566.0
#define C44_GaAs 600.0
#define ac_GaAs -7.17
#define av_GaAs -1.16

#define Eg_InAs 0.417
#define m_InAs 0.026
#define mhh_100_InAs 1.0/(g1_InAs-2.0*g2_InAs)
#define mhh_110_InAs 2.0/(2.0*g1_InAs-g2_InAs-3.0*g3_InAs)
#define g1_InAs 20.0
#define g2_InAs 8.5
#define g3_InAs 9.2
#define VBO_InAs -0.59
#define a_InAs 6.0583
#define C11_InAs 832.9
#define C12_InAs 452.6
#define C44_InAs 395.9
#define ac_InAs -5.08
#define av_InAs -1.00

#define Eg_AlAs 3.099
#define m_AlAs 0.15
#define mhh_100_AlAs 1.0/(g1_AlAs-2.0*g2_AlAs)
#define mhh_110_AlAs 2.0/(2.0*g1_AlAs-g2_AlAs-3.0*g3_AlAs)
#define g1_AlAs 3.76
#define g2_AlAs 0.82
#define g3_AlAs 1.42
#define VBO_AlAs -1.33
#define a_AlAs 5.6611
#define C11_AlAs 1250.0
#define C12_AlAs 534.0
#define C44_AlAs 542.0
#define ac_AlAs -5.64
#define av_AlAs -2.47

#define Eg_InP 1.4236
#define m_InP 0.0795
#define mhh_100_InP 1.0/(g1_InP-2.0*g2_InP)
#define mhh_110_InP 2.0/(2.0*g1_InP-g2_InP-3.0*g3_InP)
#define g1_InP 5.08
#define g2_InP 1.6
#define g3_InP 2.10
#define VBO_InP -0.94
#define a_InP 5.8697
#define C11_InP 1011.0
#define C12_InP 561.0
#define C44_InP 456.0
#define ac_InP -6.00
#define av_InP -0.60

// Ternaries... A_(1-x)B_xC

#define bow_m_AlGaAs 0.0
#define bow_VBO_AlGaAs 0.0
#define bow_a_AlGaAs 0.0
#define bow_C11_AlGaAs 0.0
#define bow_C12_AlGaAs 0.0
#define bow_C44_AlGaAs 0.0
#define bow_ac_AlGaAs 0.0
#define bow_av_AlGaAs 0.0

#define bow_Eg_GaInAs 0.477
#define bow_m_GaInAs 0.0091
#define bow_VBO_GaInAs -0.38
#define bow_a_GaInAs 0.0
#define bow_C11_GaInAs 0.0
#define bow_C12_GaInAs 0.0
#define bow_C44_GaInAs 0.0
#define bow_ac_GaInAs 2.61
#define bow_av_GaInAs 0.0

#define bow_Eg_AlInAs 0.70
#define bow_m_AlInAs 0.049
#define bow_VBO_AlInAs -0.64
#define bow_a_AlInAs 0.0
#define bow_C11_AlInAs 0.0
#define bow_C12_AlInAs 0.0
#define bow_C44_AlInAs 0.0
#define bow_ac_AlInAs -1.4
#define bow_av_AlInAs 0.0
