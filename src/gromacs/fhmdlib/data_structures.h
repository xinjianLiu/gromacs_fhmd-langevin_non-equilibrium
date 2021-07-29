#ifndef FHMD_DATA_STRUCTURES_H_
#define FHMD_DATA_STRUCTURES_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "params.h"
#include "gromacs/topology/topology.h"      /* GROMACS gmx_mtop_t  structure definition */
#include "gromacs/domdec/domdec_struct.h"   /* GROMACS gmx_domdec_t structure definition */
#include "gromacs/mdtypes/commrec.h"        /* GROMACS MPI definitions, 't_commrec', MASTER(), PAR(), etc. */
#include "gromacs/gmxlib/network.h"         /* GROMACS MPI functions, gmx_bcast(), gmx_sumf(), etc. */
#include "gromacs/math/vectypes.h"          /* GROMACS vector types: rvec, dvec, ivec, etc. */
#include "gromacs/math/vec.h"               /* GROMACS vector operations: copy_ivec, dvec_add, etc. */


typedef struct FH_arrays                    /* FH/MD arrays */
{
    double      ro_md, ro_fh;               /* densities */
    double      inv_ro;                     /* inverse density: 1/ro_md */
    dvec        u_md, u_fh, u_fh_flow, u_md_cell,u_avg_P,u_sum_P;      /* velocities */
    dvec        uro_md, uro_md_flow;        /* momentum */
    dvec        f_fh;                       /* FH force */
    dvec        alpha_term;                 /* alpha term for du/dt equation */
    dvec        beta_term;                  /* beta term for du/dt equation */
    dvec        Uflow;               /* Non-equilibrium flow velocity distribution*/

    double      delta_ro;                   /* delta of MD and FH densities */
    dvec        grad_ro;                    /* grad of density */
    matrix      alpha_u_grad;               /* preliminary alpha-term [u-index][grad-index] */

    double      S;                          /* S parameter in the FH cell centres */
    dvec        Sf;                         /* S parameter in the FH cell faces */

    double      p, pn;                      /* FH pressure */
    dvec        rof, rofn, pf, pfn;         /* FH flux variables */
    matrix      uf, ufn;                    /* FH flux velocities */
    double      ro_fh_n;                    /* FH density (new time layer) */
    dvec        u_fh_n;                     /* FH velocity (new time layer) */
    matrix      rans;                       /* FH random stress */

    double      ro_md_prime;                /* MD density from layer n */
    dvec        uro_md_prime;               /* MD momentum from layer n */

    double      ro_prime, ron_prime, ronn_prime;    /* density prime */
    double      ro_prime_b;                         /* density prime on the boundary */
    double      ro_star, ron_star;                  /* density star */
    dvec        m_prime, mn_prime, mnn_prime;       /* momentum prime */
    dvec        m_prime_b;                          /* momentum prime on the boundary */
    dvec        m_star, mn_star;                    /* momentum star */

    /* sources of MD density and MD momentum */
    double      ro_md_s, MDS_RHO_N;
    dvec        uro_md_s, f_md_s, UURHO_S, MDS_MOM_N, F_MD_S_AVA, UURHO_S_AVA;
    matrix      uuro_md_s;

    double     AVT_RHO_P, AVT_RHO_C;       /* new artificial viscocity term for balancing the MD density source*/
    dvec       AVT_MOM_P, AVT_MOM_C;	   /* new artificial viscocity term for balancing the MD momentum source*/

    dvec        Ke_cell;                            /* kinetic energy component of every cell */
    double      T_cell;                             /* Temperature of every cell */
    double      Natom_cell;                         /* atom number of every cell*/
    double      lambda_cell;                        /* kinetic energy component of every cell */
    double      gamma_cell;                         /* cell_based gamma parameter for FH in Langevin-type thermostat */

    dvec        q_prime;                            /* q' of particle zone (for testing) */
} FH_arrays;


typedef struct FH_grid          /* Computational grid */
{
    dvec       *c;              /* FH cell centres coordinates */
    dvec       *n;              /* FH cell nodes coordinates */
    dvec       *h;              /* FH cell steps */
    double     *vol;            /* FH cell volume */
    double     *ivol;           /* 1/cellVolume */
    FHMD_CELL  *md;             /* The cell is in the FH_zone, hybrid_zone or boundary */
} FH_grid;


typedef struct MD_stat          /* Particle statistics */
{
    int         N, n, n_flow, n_flow_p;
    double      invN;

    double     *avg_rho_md_cell, *avg_rho_fh_cell;
    double      avg_rho_md_ccell, avg_rho_fh_ccell; //central cell

    double      davg_rho_md,  davg_rho_fh;
    double      davg_rho2_md, davg_rho2_fh, davg_rho2_md_ccell, davg_rho2_fh_ccell;
    dvec        avg_u_md_ccell, avg_u_fh_ccell;
    dvec        davg_u_md,    davg_u_fh;
    dvec        davg_u2_md,   davg_u2_fh,   davg_u2_md_ccell,   davg_u2_fh_ccell;

    double      avg_rho_md,   avg_rho_fh;
    double      std_rho_md,   std_rho_fh, std_rho_md_ccell, std_rho_fh_ccell;
    dvec        avg_u_md,     avg_u_fh;
    dvec        std_u_md,     std_u_fh,   std_u_md_ccell,   std_u_fh_ccell;
} MD_stat;


typedef struct FHMD
{
    FH_arrays  *arr;            /* FH/MD arrays */
    FH_grid     grid;           /* FH grid */
    MD_stat     stat;           /* Particle statistics */
    int        *ind;            /* FH cell number for each atom */
    ivec       *indv;           /* 3-component FH cell number for each atom (vector) */
    double     *mpi_linear;     /* Linear array to summarise MDFH arrays */
    dvec       *vel;            /* Real velocity of each particle (dx/dt) */

    FHMD_S      S_function;     /* S = const or S = S(x,y,z) - fixed or moving */
    FHMD_NF     NF_type;        /* Non-equilibrium flow selection (0- Const velocity: Uflow_x=0, 1-poiseuille, 2-couette)*/
    int         scheme;         /* 0 - Pure MD, 1 - One-way coupling, 2 - Two-way coupling */
    double      S;              /* Parameter S (-1 - fixed sphere, -2 - moving sphere) */
    double      NF;             /* Non-equilibrium flow selection (0- Const velocity: Uflow_x=0, 1-poiseuille, 2-couette)*/
    double      R1;             /* MD sphere radius for variable S, [0..1] */
    double      R2;             /* FH sphere radius for variable S, [0..1] */
    double      Smin;           /* Minimum S for variable S */
    double      Smax;           /* Maximum S for variable S */
    double      R12, R22, RS;   /* Derived variables from R1, R2, Smin, Smax */
    double      alpha;          /* Alpha parameter for dx/dt and du/dt equations, nm^2/ps */
    double      beta;           /* Beta parameter, nm^2/ps or ps^-1 depending on the scheme */
    double      gamma_x;        /* Gamma_x parameter (density fluctuations dissipator), ps^-1 */
    double      gamma_u;        /* Gamma_u parameter (velocity fluctuations dissipator), ps^-1 */
    double      eps_rho;        /* Eps_rho parameter for ro_prime FH equation (dissipator factor, 0 <= eps_rho <= 1) */
    double      eps_mom;        /* Eps_mom parameter for m_prime FH equation (dissipator factor, 0 <= eps_mom <= 1) */
    double      eps_hy_rho;      /* epsilon-hybrid scheme for rho' (0 <= eps_rho < 1)*/
    double      eps_hy_mom;      /* epsilon-hybrid scheme for mom' (0 <= eps_rho < 1)*/
    double      eps_av_rho;      /* artificial viscocity term dissipator for rho* (0 <= eps_rho <= 2/3)*/
    double      eps_av_mom;      /* artificial viscocity term dissipator for mom* (0 <= eps_rho <= 2/3)*/
    double      S_berendsen;    /* If S_berendsen >= 0, Berendsen thermostat works for S <= S_berendsen, otherwise factor (1-S^(-S_berendsen)) is applied */
    double      tau;            /* Langevin-type thermostat parameter, ps^-1 */

    int         T_MD_N;         /* Number of temperature values to be averaged (in MD zone) */
    double      T_MD;           /* Current temperature in the MD zone */
    double      gamma_MD;       /* Current gamma parameter for Langevin-type thermostat in pure MD zone */

    ivec        N, N_md;        /* Number of FH cells along each direction */
    ivec        N_shift;        /* N_shift = (N - N_md)/2 */
    dvec        box;            /* Box size */
    dvec        box05;          /* Half of box size */
    dvec        protein_com;    /* Protein COM coordinates */
    double      protein_mass;   /* Protein mass */
    int         protein_N;      /* Number of atoms in the protein */
    double      box_volume;     /* Volume of the box, nm^3 */
    double      total_density;  /* Total density of the box, a.m.u./nm^3 */
    int         Ntot, Ntot_md;  /* Total number of FH cells */
    int         step_MD;        /* Current MD time step */
    int         Noutput;        /* Write arrays to files every Noutput MD time steps (0 - do not write) */
    int         step_flow_output;  /*MD time step for flow output*/

    FHMD_EOS    eos;            /* Equation of state */
    int         FH_EOS;         /* EOS: 0 - Liquid Argon, 1 - SPC/E water */
    int         FH_step;        /* dt_FH = FH_step * dt_MD */
    int         FH_equil;       /* Number of time steps for the FH model equilibration */
    double      FH_dens;        /* FH mean density */
    double      FH_temp;        /* FH mean temperature */
    double      FH_blend;       /* FH Blending: -1 - dynamic, or define static blending parameter (0.0 = Central Diff., 1.0 = Classic CABARET) */
    double      dt_FH;          /* FH time step */
    double      std_rho;        /* Analytical STD of density */
    double      std_u;          /* Analytical STD of velocity */

    dvec        Ke_cell;                            /* kinetic energy component of every cell */
    double      T_cell;                             /* Temperature of every cell */
    double      Natom_cell;                         /* atom number of every cell*/
    double      lambda_cell;                        /* kinetic energy component of every cell */
    double      gamma_cell;                         /* cell_based gamma parameter for FH in Langevin-type thermostat */
    int         Ccell;                              /* central cell*/

    /*flow*/
    double      Uflow_x;

    /* Averaging (flow) */
    /*Flow_n, based on number of particles*/
    double      avg_vel_n[FHMD_FLOW_LAYERS],     avg_vel_S_n[FHMD_FLOW_LAYERS];
    int         avg_n_n[FHMD_FLOW_LAYERS],       avg_n_S_n[FHMD_FLOW_LAYERS];
    double      avg_vel_tot_n[FHMD_FLOW_LAYERS], avg_vel_S_tot_n[FHMD_FLOW_LAYERS];
    int         avg_n_tot_n[FHMD_FLOW_LAYERS];

    /*FLow_m, based on mass of particles*/
    double      avg_vel[FHMD_FLOW_LAYERS],     avg_vel_S[FHMD_FLOW_LAYERS];
    int         avg_n[FHMD_FLOW_LAYERS],       avg_n_S[FHMD_FLOW_LAYERS];
    double      avg_m[FHMD_FLOW_LAYERS],       avg_m_S[FHMD_FLOW_LAYERS];
    double      avg_vel_tot[FHMD_FLOW_LAYERS], avg_vel_S_tot[FHMD_FLOW_LAYERS];
    int         avg_n_tot[FHMD_FLOW_LAYERS];

    /*ufh_FLow*/
    double      *ufh_flow;            /* Sum of ufh[0] in the flow direction of the certain layer */
    double      *ufh_flow_tot;        /* Sum of ufh_flow with steps */
    double      *avg_ufh_flow_tot;    /* Time average */
} FHMD;

#endif /* FHMD_DATA_STRUCTURES_H_ */
