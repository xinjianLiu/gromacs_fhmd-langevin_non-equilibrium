#include "data_structures.h"
#include "sfunction.h"
#include "macro.h"
#include "interpolation.h"


void fhmd_update_MD_in_FH(rvec x[], rvec v[], real mass[], rvec f[], int N_atoms, FHMD *fh)
{
    FH_arrays *arr = fh->arr;
    dvec       xn;
    int        ind;
    double     S = fh->S;
    int        lr;
    double     h_half;
    const double lrs = (double)(FHMD_FLOW_LAYERS_Y);
    double     Uflow_lr[FHMD_FLOW_LAYERS_Y][DIM], layer_h[FHMD_FLOW_LAYERS_Y];
    dvec       U_flow; // A new scheme for non-equilibrium flow
    ASSIGN_DVEC(U_flow, fh->Uflow_x, 0, 0);

   	h_half = fh->box[1]/lrs/2;

    for(int k = 0; k < FHMD_FLOW_LAYERS_Y; k++)
        {
        	layer_h[k] = (2*k+1)*h_half;
        }

        for(int k = 0; k < FHMD_FLOW_LAYERS_Y; k++)
        {
        	for (int d = 0; d < DIM; d++)
        	{
        		Uflow_lr[k][d] = U_flow[d]*1/3*fh->N_md[1]*fh->N_md[1]/fh->N[1]/fh->N[1]-U_flow[d]*((layer_h[k]-fh->box05[1])*(layer_h[k]-fh->box05[1])/(fh->box[1]/fh->N_md[1]*fh->N[1]/2)/(fh->box[1]/fh->N_md[1]*fh->N[1]/2));
        	}
        }

    /* Reset statistics */
    for(int i = 0; i < fh->Ntot; i++)
    {
        arr[i].ro_md   = 0;
        arr[i].ro_md_s = 0;
        arr[i].Natom_cell = 0;
        arr[i].T_cell = 0;

        for(int d = 0; d < DIM; d++)
        {
        	arr[i].u_md_cell[d]   = 0;
        	arr[i].Ke_cell[d]   = 0;
            arr[i].f_md_s[d]	= 0;
        	arr[i].uro_md[d]    = 0;
            arr[i].uro_md_s[d]  = 0;
            arr[i].uro_md_flow[d]  = 0; // A new scheme for non-equilibrium flow
            
            for(int dim = 0; dim < DIM; dim++)
            {
            	arr[i].uuro_md_s[d][dim] = 0;
            }
        }
    }

    fh->T_MD = 0;
    fh->T_MD_N = 0;

    /* Collect statistics */
    for(int n = 0; n < N_atoms; n++)
    {
        PBC(xn, x[n], fh->box);
        for(int d = 0; d < DIM; d++)
        {
            fh->indv[n][d] = (int)(xn[d]/fh->box[d]*(double)(fh->N_md[d])) + fh->N_shift[d];
        }

        ind = I(fh->indv[n], fh->N);

        if(ind < 0 || ind >= fh->Ntot)      // This should never happen... only if the coordinates are NaN.
        {
            printf(MAKE_RED "\nFHMD: ERROR: Solution diverged. Atom #%d coordinates: (%g, %g, %g) nm\n" RESET_COLOR "\n", n, xn[0], xn[1], xn[2]);
            exit(21);
        }

        fh->ind[n] = ind;
        arr[ind].Natom_cell ++;

        for(int d = 0; d < DIM; d++)
        {
        	arr[ind].u_md_cell[d] += v[n][d];
        }
    }

    /* Update statistics */
    for(int i = 0; i < fh->Ntot; i++)
       {
           for(int d = 0; d < DIM; d++)
           {
               arr[i].u_md_cell[d] = arr[i].u_md_cell[d]/(arr[i].Natom_cell);
           }
       }

    /* Collect statistics */
    for(int n = 0; n < N_atoms; n++)
    {
        PBC(xn, x[n], fh->box);

        for(int d = 0; d < DIM; d++)
        {
            fh->indv[n][d] = (int)(xn[d]/fh->box[d]*(double)(fh->N_md[d])) + fh->N_shift[d];
        }

        ind = I(fh->indv[n], fh->N);

        if(ind < 0 || ind >= fh->Ntot)      // This should never happen... only if the coordinates are NaN.
        {
            printf(MAKE_RED "\nFHMD: ERROR: Solution diverged. Atom #%d coordinates: (%g, %g, %g) nm\n" RESET_COLOR "\n", n, xn[0], xn[1], xn[2]);
            exit(21);
        }

        fh->ind[n] = ind;

        if(fh->S_function == moving_sphere)
            S = fhmd_Sxyz_r(x[n], fh->protein_com, fh);     // MD/FH sphere follows protein
        else if(fh->S_function == fixed_sphere)
            S = fhmd_Sxyz_r(x[n], fh->box05, fh);           // Fixed MD/FH sphere

        lr = (int)(xn[YY]/fh->box[YY]*lrs);
        if(lr < 0 || lr >= FHMD_FLOW_LAYERS_Y)             // This should never happen... only if the coordinates are NaN.
        {
            printf(MAKE_RED "\nFHMD: ERROR: Solution diverged. Atom #%d coordinates: (%g, %g, %g) nm\n" RESET_COLOR "\n", n, xn[0], xn[1], xn[2]);
            exit(21);
        }

        arr[ind].ro_md   += mass[n];
        arr[ind].ro_md_s += (1 - S)*mass[n];
        //arr[ind].Natom_cell ++;
        if(S < 1e-8) fh->T_MD_N++;

        for(int d = 0; d < DIM; d++)
        {
            //arr[ind].Ke_cell[d]  += v[n][d]*v[n][d]*mass[n];
        	//arr[ind].Ke_cell[d] += (v[n][d]-arr[ind].Uflow[d])*(v[n][d]-arr[ind].Uflow[d])*mass[n];
        	//arr[ind].Ke_cell[d] += (v[n][d]-Uflow_lr[lr][d])*(v[n][d]-Uflow_lr[lr][d])*mass[n];
        	arr[ind].Ke_cell[d] += (v[n][d]-arr[ind].u_md_cell[d])*(v[n][d]-arr[ind].u_md_cell[d])*mass[n];
           	arr[ind].f_md_s[d]   += (1 - S)*f[n][d];

           	if (fh->eos == eos_spce && S < 1e-8)
           		//fh->T_MD += 1.0/FHMD_T_DOF/FHMD_kB*(v[n][d]-arr[ind].Uflow[d])*(v[n][d]-arr[ind].Uflow[d])*mass[n];
           		//fh->T_MD += 1.0/FHMD_T_DOF/FHMD_kB*(v[n][d]-Uflow_lr[lr][d])*(v[n][d]-Uflow_lr[lr][d])*mass[n];
           		fh->T_MD += 1.0/FHMD_T_DOF/FHMD_kB*(v[n][d]-arr[ind].u_md_cell[d])*(v[n][d]-arr[ind].u_md_cell[d])*mass[n];
           	else if (fh->eos == eos_argon && S < 1e-8 )
           		//fh->T_MD += 1.0/FHMD_T_DOF_Ar/FHMD_kB*(v[n][d]-arr[ind].Uflow[d])*(v[n][d]-arr[ind].Uflow[d])*mass[n];
           	    //fh->T_MD += 1.0/FHMD_T_DOF_Ar/FHMD_kB*(v[n][d]-Uflow_lr[lr][d])*(v[n][d]-Uflow_lr[lr][d])*mass[n];
           		fh->T_MD += 1.0/FHMD_T_DOF_Ar/FHMD_kB*(v[n][d]-arr[ind].u_md_cell[d])*(v[n][d]-arr[ind].u_md_cell[d])*mass[n];

           	// A new scheme for non-equilibrium flow
           	if (fh-> step_MD == 0)
           	{
            	arr[ind].uro_md[d]   += v[n][d]*mass[n];
                arr[ind].uro_md_s[d] += (1 - S)*v[n][d]*mass[n];
           	}
           	else
           	{
           		//arr[ind].uro_md[d]   += (v[n][d]-arr[ind].Uflow[d])*mass[n];
           		//arr[ind].uro_md_s[d] += (1 - S)*(v[n][d]-arr[ind].Uflow[d])*mass[n];
           		//arr[ind].uro_md[d]   += (v[n][d]-Uflow_lr[lr][d])*mass[n];
           		//arr[ind].uro_md_s[d] += (1 - S)*(v[n][d]-Uflow_lr[lr][d])*mass[n];.
           		arr[ind].uro_md[d]   += (v[n][d]-arr[ind].u_md_cell[d])*mass[n];
           		arr[ind].uro_md_s[d] += (1 - S)*(v[n][d]-arr[ind].u_md_cell[d])*mass[n];
           	}
           	arr[ind].uro_md_flow[d]   += v[n][d]*mass[n];

            for(int dim = 0; dim < DIM; dim++)
            {
               	if (fh-> step_MD == 0)
               	{
               		arr[ind].uuro_md_s[d][dim] += (1 - S)*mass[n]*v[n][d]*v[n][dim];
               	}
               	else
               	{
               		//arr[ind].uuro_md_s[d][dim] += (1 - S)*mass[n]*(v[n][d]-arr[ind].Uflow[d])*(v[n][dim]-arr[ind].Uflow[dim]);
               		//arr[ind].uuro_md_s[d][dim] += (1 - S)*mass[n]*(v[n][d]-Uflow_lr[lr][d])*(v[n][dim]-Uflow_lr[lr][dim]);
               		arr[ind].uuro_md_s[d][dim] += (1 - S)*mass[n]*(v[n][d]-arr[ind].u_md_cell[d])*(v[n][dim]-arr[ind].u_md_cell[dim]);
               	}
            }
        }
    }

    /* Update statistics */
    for(int i = 0; i < fh->Ntot; i++)
    {
        arr[i].ro_md   *= fh->grid.ivol[i];
        arr[i].ro_md_s *= fh->grid.ivol[i];

        for(int d = 0; d < DIM; d++)
        {
            arr[i].uro_md[d]   *= fh->grid.ivol[i];
            arr[i].uro_md_s[d] *= fh->grid.ivol[i];
            arr[i].uro_md_flow[d]   *= fh->grid.ivol[i];
            
            for(int dim = 0; dim < DIM; dim++)
            {
                arr[i].uuro_md_s[d][dim] *= fh->grid.ivol[i];
            }
        }
    }
}


void fhmd_sum_arrays(t_commrec *cr, FHMD *fh)
{
    FH_arrays *arr = fh->arr;

    /* Pack FHMD arrays to linear array */
    for(int i = 0; i < fh->Ntot; i++)
    {
        fh->mpi_linear[i]              = arr[i].ro_md;
        fh->mpi_linear[i + fh->Ntot]   = arr[i].uro_md[0];
        fh->mpi_linear[i + fh->Ntot*2] = arr[i].uro_md[1];
        fh->mpi_linear[i + fh->Ntot*3] = arr[i].uro_md[2];

        fh->mpi_linear[i + fh->Ntot*4] = arr[i].ro_md_s;
        fh->mpi_linear[i + fh->Ntot*5] = arr[i].uro_md_s[0];
        fh->mpi_linear[i + fh->Ntot*6] = arr[i].uro_md_s[1];
        fh->mpi_linear[i + fh->Ntot*7] = arr[i].uro_md_s[2];

        fh->mpi_linear[i + fh->Ntot*8] = arr[i].f_md_s[0];
        fh->mpi_linear[i + fh->Ntot*9] = arr[i].f_md_s[1];
        fh->mpi_linear[i + fh->Ntot*10] = arr[i].f_md_s[2];

        fh->mpi_linear[i + fh->Ntot*11] = arr[i].uuro_md_s[0][0];
        fh->mpi_linear[i + fh->Ntot*12] = arr[i].uuro_md_s[0][1];
        fh->mpi_linear[i + fh->Ntot*13] = arr[i].uuro_md_s[0][2];
        fh->mpi_linear[i + fh->Ntot*14] = arr[i].uuro_md_s[1][0];
        fh->mpi_linear[i + fh->Ntot*15] = arr[i].uuro_md_s[1][1];
        fh->mpi_linear[i + fh->Ntot*16] = arr[i].uuro_md_s[1][2];
        fh->mpi_linear[i + fh->Ntot*17] = arr[i].uuro_md_s[2][0];
        fh->mpi_linear[i + fh->Ntot*18] = arr[i].uuro_md_s[2][1];
        fh->mpi_linear[i + fh->Ntot*19] = arr[i].uuro_md_s[2][2];

        fh->mpi_linear[i + fh->Ntot*20] = arr[i].Natom_cell;
        fh->mpi_linear[i + fh->Ntot*21] = arr[i].Ke_cell[0];
        fh->mpi_linear[i + fh->Ntot*22] = arr[i].Ke_cell[1];
        fh->mpi_linear[i + fh->Ntot*23] = arr[i].Ke_cell[2];

        fh->mpi_linear[i + fh->Ntot*24] = arr[i].uro_md_flow[0];
        fh->mpi_linear[i + fh->Ntot*25] = arr[i].uro_md_flow[1];
        fh->mpi_linear[i + fh->Ntot*26] = arr[i].uro_md_flow[2];
    }

    /* Broadcast linear array */
    gmx_sumd(fh->Ntot*27, fh->mpi_linear, cr);
    gmx_sumd(1, &fh->T_MD, cr);
    gmx_sumi(1, &fh->T_MD_N, cr);

    /* Unpack linear array */
    for(int i = 0; i < fh->Ntot; i++)
    {
        arr[i].ro_md     = fh->mpi_linear[i];
        arr[i].uro_md[0] = fh->mpi_linear[i + fh->Ntot];
        arr[i].uro_md[1] = fh->mpi_linear[i + fh->Ntot*2];
        arr[i].uro_md[2] = fh->mpi_linear[i + fh->Ntot*3];

        arr[i].ro_md_s     = fh->mpi_linear[i + fh->Ntot*4];
        arr[i].uro_md_s[0] = fh->mpi_linear[i + fh->Ntot*5];
        arr[i].uro_md_s[1] = fh->mpi_linear[i + fh->Ntot*6];
        arr[i].uro_md_s[2] = fh->mpi_linear[i + fh->Ntot*7];

        arr[i].f_md_s[0] = fh->mpi_linear[i + fh->Ntot*8];
        arr[i].f_md_s[1] = fh->mpi_linear[i + fh->Ntot*9];
        arr[i].f_md_s[2] = fh->mpi_linear[i + fh->Ntot*10];

        arr[i].uuro_md_s[0][0] = fh->mpi_linear[i + fh->Ntot*11];
        arr[i].uuro_md_s[0][1] = fh->mpi_linear[i + fh->Ntot*12];
        arr[i].uuro_md_s[0][2] = fh->mpi_linear[i + fh->Ntot*13];
        arr[i].uuro_md_s[1][0] = fh->mpi_linear[i + fh->Ntot*14];
        arr[i].uuro_md_s[1][1] = fh->mpi_linear[i + fh->Ntot*15];
        arr[i].uuro_md_s[1][2] = fh->mpi_linear[i + fh->Ntot*16];
        arr[i].uuro_md_s[2][0] = fh->mpi_linear[i + fh->Ntot*17];
        arr[i].uuro_md_s[2][1] = fh->mpi_linear[i + fh->Ntot*18];
        arr[i].uuro_md_s[2][2] = fh->mpi_linear[i + fh->Ntot*19];

        arr[i].Natom_cell     = fh->mpi_linear[i + fh->Ntot*20];
        arr[i].Ke_cell[0]     = fh->mpi_linear[i + fh->Ntot*21];
        arr[i].Ke_cell[1]     = fh->mpi_linear[i + fh->Ntot*22];
        arr[i].Ke_cell[2]     = fh->mpi_linear[i + fh->Ntot*23];

        arr[i].uro_md_flow[0] = fh->mpi_linear[i + fh->Ntot*24];
        arr[i].uro_md_flow[1] = fh->mpi_linear[i + fh->Ntot*25];
        arr[i].uro_md_flow[2] = fh->mpi_linear[i + fh->Ntot*26];
    }
}


void fhmd_calculate_MDFH_terms(FHMD *fh)
{
    FH_arrays *arr = fh->arr;

    ivec ind, ind1; // inds for 27 cells average
    dvec alpha_term, f_md_term, MDS_RHO, MDS_MOM;
	matrix MDS_RHO_L, MDS_RHO_R, MDS_RHO_M, MDS_MOM_L, MDS_MOM_R,MDS_MOM_M;

    fh->T_MD /= (double)fh->T_MD_N;

	for(int k = 0; k < NZ; k++)
	{
	    for(int j = 0; j < NY; j++)
	    {
	        for(int i = 0; i < NX; i++)
	        {
	            ASSIGN_IND(ind, i, j, k);

	            if (arr[C].S < 1)
	            {
	            	if (fh->eos == eos_spce)
	            		arr[C].T_cell= (arr[C].Ke_cell[0]+arr[C].Ke_cell[1]+arr[C].Ke_cell[2])/(arr[C].Natom_cell)/FHMD_T_DOF*166.053886/1.380649;
	                else if (fh->eos == eos_argon)
	                	arr[C].T_cell= (arr[C].Ke_cell[0]+arr[C].Ke_cell[1]+arr[C].Ke_cell[2])/(arr[C].Natom_cell)/FHMD_T_DOF_Ar*166.053886/1.380649;
	            }
	            else if (arr[C].S == 1)
	            {
	                arr[C].T_cell = fh->FH_temp *(1-arr[C].S)/(1-arr[C].S/2);
	                //arr[C].T_cell = fh->FH_temp *sqrt(1-arr[C].S);
	            }

	            double T_ref = fh->FH_temp *(1-arr[C].S)/(1-arr[C].S/2); //local thermostat
	            //double T_ref = fh->FH_temp;  //constant thermostat
	            //double T_ref = fh->FH_temp *sqrt(1-arr[C].S);

	            if (arr[C].S < 1)
	            {
	   	             arr[C].lambda_cell = sqrt(1.0 + fh->dt_FH/(double)(fh->FH_step)/fh->tau*(T_ref/arr[C].T_cell-1.0));
	   	             arr[C].gamma_cell = (1.0-arr[C].lambda_cell)/fh->dt_FH*(double)(fh->FH_step);
	            }
	            else if (arr[C].S == 1)
	            {
	            	arr[C].gamma_cell = 0;
	            }
	        }
	    }
	}

    for(int k = fh->N_shift[2]; k < (fh->N_md[2] + fh->N_shift[2]); k++)
    {
        for(int j = fh->N_shift[1]; j < (fh->N_md[1] + fh->N_shift[1]); j++)
        {
            for(int i = fh->N_shift[0]; i < (fh->N_md[0] + fh->N_shift[0]); i++)
            {
                ASSIGN_IND(ind, i, j, k);

                if(arr[C].ro_md <= 0) {
                    printf(MAKE_RED "\nFHMD: ERROR: Zero or NaN MD density in the cell %d-%d-%d (ro_md = %g)\n" RESET_COLOR "\n", i, j, k, arr[C].ro_md);
                    exit(22);
                }

                arr[C].inv_ro = 1.0/arr[C].ro_md;

                for(int d = 0; d < DIM; d++)
                    arr[C].u_md[d] = arr[C].uro_md_flow[d]*arr[C].inv_ro; // A new scheme for non-equilibrium flow

                if(fh->scheme == One_Way)
                {
                    arr[C].delta_ro = arr[C].ro_fh - arr[C].ro_md;
                    for(int d = 0; d < DIM; d++)
                        arr[C].beta_term[d] = fh->beta*(arr[C].u_fh[d]*arr[C].ro_fh - arr[C].uro_md[d]);
                }
                else if(fh->scheme == Two_Way)
                {
//                  arr[C].delta_ro = arr[C].ron_prime;
//                  for(int d = 0; d < DIM; d++)
//                      arr[C].beta_term[d] = fh->beta*arr[C].mn_prime[d];
                    arr[C].delta_ro = arr[C].ro_fh - arr[C].ro_md;                                                    // Layer n may work better than n+1/2
                    for(int d = 0; d < DIM; d++)
                        //arr[C].beta_term[d] = fh->beta*(arr[C].u_fh[d]*arr[C].ro_fh - arr[C].uro_md[d]);            // Layer n may work better than n+1/2
                        arr[C].beta_term[d] = fh->beta*(arr[C].u_fh_flow[d]*arr[C].ro_fh - arr[C].uro_md_flow[d]);    // A new scheme for non-equilibrium flow
                    if(fh->grid.md[C] == FH_zone) arr[C].delta_ro = 0;
                }
            }
        }
    }


/*
 * Output the q' of particle zone
 */

    for(int k = fh->N_shift[2]; k < (fh->N_md[2] + fh->N_shift[2]); k++)
    {
        for(int j = fh->N_shift[1]; j < (fh->N_md[1] + fh->N_shift[1]); j++)
        {
            for(int i = fh->N_shift[0]; i < (fh->N_md[0] + fh->N_shift[0]); i++)
            {
                ASSIGN_IND(ind, i, j, k);

                for(int d = 0; d < DIM; d++)
                {
                   arr[Cm].q_prime[d] = arr[Cm].u_fh_flow[d]*arr[Cm].ro_fh - arr[Cm].uro_md_flow[d];
                }
            }
        }
    }
    
    for(int k = fh->N_shift[2]; k < (fh->N_md[2] + fh->N_shift[2]); k++)
    {
        for(int j = fh->N_shift[1]; j < (fh->N_md[1] + fh->N_shift[1]); j++)
        {
            for(int i = fh->N_shift[0]; i < (fh->N_md[0] + fh->N_shift[0]); i++)
            {
                ASSIGN_IND(ind, i, j, k);

                for(int d = 0; d < DIM; d++)
                {
                    arr[Cm].grad_ro[d] = fh->alpha*(arr[CRm].delta_ro - arr[CLm].delta_ro)/(0.5*(fh->grid.h[CLm][d] + 2.0*fh->grid.h[Cm][d] + fh->grid.h[CRm][d]));

                    for(int du = 0; du < DIM; du++)
                        arr[Cm].alpha_u_grad[du][d] = arr[Cm].grad_ro[d]*arr[Cm].S*(1 - arr[Cm].S)*arr[Cm].u_md[du];    // TODO: Fast but rough estimation!
                }
            }
        }
    }

    for(int k = fh->N_shift[2]; k < (fh->N_md[2] + fh->N_shift[2]); k++)
    {
        for(int j = fh->N_shift[1]; j < (fh->N_md[1] + fh->N_shift[1]); j++)
        {
            for(int i = fh->N_shift[0]; i < (fh->N_md[0] + fh->N_shift[0]); i++)
            {
                ASSIGN_IND(ind, i, j, k);

                if((fh->scheme == Two_Way) && (fh->grid.md[C] == boundary)) continue;

                for(int du = 0; du < DIM; du++)
                {
                    for(int d = 0; d < DIM; d++)
                    {
                        alpha_term[d] = (arr[CRm].alpha_u_grad[du][d] - arr[CLm].alpha_u_grad[du][d])
                                /(0.5*(fh->grid.h[CLm][d] + 2.0*fh->grid.h[Cm][d] + fh->grid.h[CRm][d]));
                    }

                    arr[Cm].alpha_term[du] = SUM(alpha_term);
                }
            }
        }
    }

    for(int k = fh->N_shift[2]; k < (fh->N_md[2] + fh->N_shift[2]); k++)
    {
        for(int j = fh->N_shift[1]; j < (fh->N_md[1] + fh->N_shift[1]); j++)
        {
            for(int i = fh->N_shift[0]; i < (fh->N_md[0] + fh->N_shift[0]); i++)
            {
                ASSIGN_IND(ind, i, j, k);

                for(int d = 0; d < DIM; d++)
                {
                	MDS_RHO[d]=0;
                	MDS_MOM[d]=0;
                	f_md_term[d]=0;
                }

                for(int d = 0; d < DIM; d++)
                {
                	MDS_RHO[d] = -0.5*(arr[CRm].uro_md_s[d]-arr[CLm].uro_md_s[d])/fh->grid.h[Cm][d];
                }

                arr[Cm].MDS_RHO_N  = SUM(MDS_RHO);

                for(int dim  = 0; dim < DIM; dim++)
                {
                	 for(int d = 0; d < DIM; d++)
                	 {
                         MDS_MOM[d] = -0.5*(arr[CRm].uuro_md_s[dim][d]-arr[CLm].uuro_md_s[dim][d])/fh->grid.h[Cm][d];
                	 }
                	 arr[Cm].UURHO_S_AVA[dim] = SUM(MDS_MOM);
                }

               /* for(int d = 0; d < DIM; d++)
                {
                 	arr[Cm].MDS_MOM_N[d]  = arr[Cm].f_md_s[d]+arr[Cm].UURHO_S_AVA[d];
                }*/

                for(int k1 = k-1; k1 <=k+1; k1++)
                {
                     for(int j1 = j-1; j1 <= j+1; j1++)
                     {
                         for(int i1 = i-1; i1 <= i+1; i1++)
                         {
                        	 copy_ivec(ind, ind1);        // Backup global (i,j,k) ind
                             ASSIGN_IND(ind, i1, j1, k1); // Create new 'local' ind
                             for(int d = 0; d < DIM; d++)
                             {
                            	 f_md_term[d]+=arr[Cm].f_md_s[d];
                             }
                         }
                     }
                 }
                copy_ivec(ind1, ind);

                for(int d = 0; d < DIM; d++)
                {
                	arr[Cm].F_MD_S_AVA[d] = f_md_term[d]/27;
                }

                for(int d = 0; d < DIM; d++)
                {
                	arr[Cm].MDS_MOM_N[d]  = arr[Cm].F_MD_S_AVA[d]+arr[Cm].UURHO_S_AVA[d];
                }
            }
        }
    }

    //27 cells for the average of the MD force in the MD_SOURCE term
  /*  for(int k = fh->N_shift[2]; k < (fh->N_md[2] + fh->N_shift[2]); k++)
     {
         for(int j = fh->N_shift[1]; j < (fh->N_md[1] + fh->N_shift[1]); j++)
         {
             for(int i = fh->N_shift[0]; i < (fh->N_md[0] + fh->N_shift[0]); i++)
             {
                 ASSIGN_IND(ind, i, j, k);

                 for(int dim = 0; dim < DIM; dim++)
                 {
                 	f_md_term[dim]=0;
                 }

                 for(int k1 = k-1; k1 <=k+1; k1++)
                 {
                      for(int j1 = j-1; j1 <= j+1; j1++)
                      {
                          for(int i1 = i-1; i1 <= i+1; i1++)
                          {
                         	 copy_ivec(ind, ind1);        // Backup global (i,j,k) ind
                              ASSIGN_IND(ind, i1, j1, k1); // Create new 'local' ind
                              for(int d = 0; d < DIM; d++)
                              {
                             	 f_md_term[d]+=arr[Cm].f_md_s[d];
                              }
                          }
                      }
                  }
                 copy_ivec(ind1, ind);

                 for(int d = 0; d < DIM; d++)
                 {
                 	arr[Cm].F_MD_S_AVA[d] = f_md_term[d]/27;
                 }

                 for(int d = 0; d < DIM; d++)
                 {
                 	arr[Cm].MDS_MOM_N[d]  = arr[Cm].F_MD_S_AVA[d]+arr[Cm].UURHO_S_AVA[d];
                 }
             }
         }
     }*/

    //27 cells for the average of MD_SOURCE term
   /*for(int k = fh->N_shift[2]; k < (fh->N_md[2] + fh->N_shift[2]); k++)
    {
        for(int j = fh->N_shift[1]; j < (fh->N_md[1] + fh->N_shift[1]); j++)
        {
            for(int i = fh->N_shift[0]; i < (fh->N_md[0] + fh->N_shift[0]); i++)
            {
                ASSIGN_IND(ind, i, j, k);

                for(int dim = 0; dim < DIM; dim++)
                {
                	f_md_term[dim]=0;
                	for(int d = 0; d < DIM; d++)
                    {
                    	MDS_RHO_L[dim][d]=0;
                    	MDS_RHO_R[dim][d]=0;
                    	MDS_RHO_M[dim][d]=0;
                    	MDS_MOM_L[dim][d]=0;
                    	MDS_MOM_R[dim][d]=0;
                    	MDS_MOM_M[dim][d]=0;
                    }
                }

                for(int k1 = k-1; k1 <=k+1; k1++)
                {
                     for(int j1 = j-1; j1 <= j+1; j1++)
                     {
                         for(int i1 = i-1; i1 <= i+1; i1++)
                         {
                        	 copy_ivec(ind, ind1);        // Backup global (i,j,k) ind
                             ASSIGN_IND(ind, i1, j1, k1); // Create new 'local' ind

                             if(i1 == i-1)
                             {
                            	 MDS_RHO_L[0][0] += arr[Cm].uro_md_s[0];  //X-direction
                            	 for(int d = 0; d < DIM; d++)
                            	 {
                            		 MDS_MOM_L[0][d] += arr[Cm].uuro_md_s[0][d];
                            	 }
                             }

                             if(i1 == i+1)
                             {
                        		 MDS_RHO_R[0][0] += arr[Cm].uro_md_s[0];
                            	 for(int d = 0; d < DIM; d++)
                            	 {
                            		 MDS_MOM_R[0][d] += arr[Cm].uuro_md_s[0][d];
                            	 }
                             }

                             if(j1 == j-1)
                             {
                        		 MDS_RHO_L[1][1] += arr[Cm].uro_md_s[1];  //Y-direction
                            	 for(int d = 0; d < DIM; d++)
                            	 {
                            		 MDS_MOM_L[1][d] += arr[Cm].uuro_md_s[1][d];
                            	 }
                             }

                             if(j1 == j+1)
                             {
                        		 MDS_RHO_R[1][1] += arr[Cm].uro_md_s[1];
                            	 for(int d = 0; d < DIM; d++)
                            	 {
                            		 MDS_MOM_R[1][d] += arr[Cm].uuro_md_s[1][d];
                            	 }
                             }

                             if(k1 == k-1)
                             {
                        		 MDS_RHO_L[2][2] += arr[Cm].uro_md_s[2];   //Z-direction
                            	 for(int d = 0; d < DIM; d++)
                            	 {
                            		 MDS_MOM_L[2][d] += arr[Cm].uuro_md_s[2][d];
                            	 }
                             }

                             if(k1 == k+1)
                             {
                        		 MDS_RHO_R[2][2] += arr[Cm].uro_md_s[2];
                            	 for(int d = 0; d < DIM; d++)
                            	 {
                            		 MDS_MOM_R[2][d] += arr[Cm].uuro_md_s[2][d];
                            	 }
                             }

                             for(int d = 0; d < DIM; d++)
                             {
                            	 f_md_term[d]+=arr[Cm].f_md_s[d];
                             }
                         }
                     }
                 }
                copy_ivec(ind1, ind);

                for(int dim = 0; dim < DIM; dim++)
                {
                    for(int d = 0; d < DIM; d++)
                    {
                    	if(dim == d)
                    	{
                    		MDS_RHO_M[dim][d] = -0.5*(MDS_RHO_R[dim][d]-MDS_RHO_L[dim][d])/fh->box[d]/(double)(fh->N_md[d])/9;
                    	}
                    	MDS_MOM_M[dim][d] = -0.5*(MDS_MOM_R[dim][d]-MDS_MOM_L[dim][d])/fh->box[d]/(double)(fh->N_md[d])/9;
                    }
                }

                arr[Cm].MDS_RHO_N = MDS_RHO_M[0][0]+MDS_RHO_M[1][1]+MDS_RHO_M[2][2];

                for(int d = 0; d < DIM; d++)
                {
                	arr[Cm].F_MD_S_AVA[d] = f_md_term[d]/27;
                }

                for(int d = 0; d < DIM; d++)
                {
                	arr[Cm].UURHO_S_AVA[d] = MDS_MOM_M[0][d]+MDS_MOM_M[1][d]+MDS_MOM_M[2][d];
                }

                for(int d = 0; d < DIM; d++)
                {
                	arr[Cm].MDS_MOM_N[d]  = arr[Cm].F_MD_S_AVA[d]+arr[Cm].UURHO_S_AVA[d];
                	//arr[Cm].MDS_MOM_N[d]  = arr[Cm].UURHO_S_AVA[d];
                }
            }
        }
    }*/
}

