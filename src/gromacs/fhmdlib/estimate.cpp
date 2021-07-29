#include "data_structures.h"
#include "sfunction.h"
#include "macro.h"


void fhmd_reset_statistics(FHMD *fh)
{
    MD_stat *st = &fh->stat;

    st->N            = 0;
    st->n            = 0;
    st->invN         = 0;
    st->n_flow       = 0;

    st->davg_rho_md  = 0;
    st->davg_rho_fh  = 0;
    st->davg_rho2_md = 0;
    st->davg_rho2_fh = 0;
    st->davg_rho2_md_ccell = 0;
    st->davg_rho2_fh_ccell = 0;
    st->avg_rho_md_ccell   = 0;
    st->avg_rho_fh_ccell   = 0;

    for(int d = 0; d < DIM; d++)
    {
        st->davg_u_md[d]  = 0;
        st->davg_u2_md[d] = 0;
        st->davg_u_fh[d]  = 0;
        st->davg_u2_fh[d] = 0;
        st->davg_u2_md_ccell[d] = 0;
        st->davg_u2_fh_ccell[d] = 0;
        st->avg_u_md_ccell[d]  = 0;
        st->avg_u_fh_ccell[d]  = 0;
    }
}


void fhmd_collect_statistics(FHMD *fh)
{
    MD_stat   *st  = &fh->stat;
    FH_arrays *arr =  fh->arr;

    st->n++;
    double invn = 1.0/(double)(st->n);

    for(int i = 0; i < fh->Ntot; i++)
    {
        if(fh->grid.md[i] == FH_zone) continue;     // Skip pure FH region

        st->N++;

        st->avg_rho_md_cell[i] += arr[i].ro_md;
        st->avg_rho_fh_cell[i] += arr[i].ro_fh;
        if (i == fh->Ccell)
        {
        	st->avg_rho_md_ccell += arr[i].ro_md;
        	st->avg_rho_fh_ccell += arr[i].ro_fh;
        	for(int d = 0; d < DIM; d++)
        	{
            	st->avg_u_md_ccell[d] += arr[i].u_md[d];
            	st->avg_u_fh_ccell[d] += arr[i].u_fh[d];
        	}
        }

        st->davg_rho_md  += arr[i].ro_md;
        st->davg_rho_fh  += arr[i].ro_fh;

        // Faster convergence but less accurate
        //st->davg_rho2_md += (arr[i].ro_md - fh->total_density)*(arr[i].ro_md - fh->total_density);
        //st->davg_rho2_fh += (arr[i].ro_fh - fh->FH_dens)*(arr[i].ro_fh - fh->FH_dens);

        // More accurate estimation but much longer convergence
        st->davg_rho2_md += (arr[i].ro_md - st->avg_rho_md_cell[i]*invn)*(arr[i].ro_md - st->avg_rho_md_cell[i]*invn);
        st->davg_rho2_fh += (arr[i].ro_fh - st->avg_rho_fh_cell[i]*invn)*(arr[i].ro_fh - st->avg_rho_fh_cell[i]*invn);
        if (i == fh->Ccell)
        {
        	st->davg_rho2_md_ccell = (arr[i].ro_md - st->avg_rho_md_ccell*invn)*(arr[i].ro_md - st->avg_rho_md_ccell*invn);
        	st->davg_rho2_fh_ccell = (arr[i].ro_fh - st->avg_rho_fh_ccell*invn)*(arr[i].ro_fh - st->avg_rho_fh_ccell*invn);
        	for(int d = 0; d < DIM; d++)
        	{
            	st->davg_u2_md_ccell[d] = (arr[i].u_md[d] - st->avg_u_md_ccell[d]*invn)*(arr[i].u_md[d] - st->avg_u_md_ccell[d]*invn);
            	st->davg_u2_fh_ccell[d] = (arr[i].u_fh[d] - st->avg_u_fh_ccell[d]*invn)*(arr[i].u_fh[d] - st->avg_u_fh_ccell[d]*invn);
        	}
        }

        for(int d = 0; d < DIM; d++)
        {
            st->davg_u_md[d]  += arr[i].u_md[d];
            st->davg_u_fh[d]  += arr[i].u_fh[d];
            st->davg_u2_md[d] += arr[i].u_md[d]*arr[i].u_md[d];
            st->davg_u2_fh[d] += arr[i].u_fh[d]*arr[i].u_fh[d];
        }
    }

    st->invN = 1.0/(double)(st->N);
}


void fhmd_update_statistics(FHMD *fh)
{
    MD_stat *st = &fh->stat;

    st->avg_rho_md = st->davg_rho_md*st->invN;
    st->avg_rho_fh = st->davg_rho_fh*st->invN;
    st->std_rho_md = sqrt(st->davg_rho2_md*st->invN);
    st->std_rho_fh = sqrt(st->davg_rho2_fh*st->invN);
    st->std_rho_md_ccell = sqrt(st->davg_rho2_md_ccell);
    st->std_rho_fh_ccell = sqrt(st->davg_rho2_fh_ccell);

    for(int d = 0; d < DIM; d++)
    {
        st->avg_u_md[d] = st->davg_u_md[d]*st->invN;
        st->avg_u_fh[d] = st->davg_u_fh[d]*st->invN;
        st->std_u_md[d] = sqrt(fabs(st->davg_u2_md[d]*st->invN - st->avg_u_md[d]*st->avg_u_md[d]));
        st->std_u_fh[d] = sqrt(fabs(st->davg_u2_fh[d]*st->invN - st->avg_u_fh[d]*st->avg_u_fh[d]));
        st->std_u_md_ccell[d] = sqrt(st->davg_u2_md_ccell[d]);
        st->std_u_fh_ccell[d] = sqrt(st->davg_u2_fh_ccell[d]);
    }
}


void fhmd_print_statistics(FHMD *fh, t_commrec *cr)
{
    MD_stat *st = &fh->stat;

    fhmd_collect_statistics(fh);
    fhmd_update_statistics(fh);     // Every MD time step -- for stochastic integration

    if(MASTER(cr))
    {
        if((!(fh->step_MD % (fh->FH_step))    && (fh->step_MD < (fh->FH_step*10)))  ||
           (!(fh->step_MD % (fh->FH_step*10)) && (fh->step_MD < (fh->FH_step*100))) ||
            !(fh->step_MD % (fh->FH_step*100)))
        {
            if(!fh->step_MD)
            {
                printf("%8s " MAKE_LIGHT_BLUE "%10s " MAKE_BLUE "%10s " MAKE_LIGHT_BLUE "%9s %9s %9s " MAKE_BLUE "%9s %9s %9s " MAKE_LIGHT_BLUE "%10s " MAKE_BLUE "%10s " MAKE_LIGHT_BLUE "%9s %9s %9s " MAKE_BLUE "%9s %9s %9s " MAKE_LIGHT_BLUE "%9s" MAKE_BLUE "%9s"
                       RESET_COLOR "\n",
                       "Step", "STD_rho_MD", "STD_rho_FH", "STD_Ux_MD", "STD_Uy_MD", "STD_Uz_MD", "STD_Ux_FH", "STD_Uy_FH", "STD_Uz_FH", "STD_rho_MD_ccell", "STD_rho_FH_ccell", "STD_Ux_MD_C", "STD_Uy_MD_C", "STD_Uz_MD_C", "STD_Ux_FH_C", "STD_Uy_FH_C", "STD_Uz_FH_C", "T_MD, K","rho(C)_T");
                       //"Step", "STD_rho_MD", "STD_rho_FH", "STD_Ux_MD", "STD_Uy_MD", "STD_Uz_MD", "STD_Ux_FH", "STD_Uy_FH", "STD_Uz_FH",
                       //"<Ux_MD>", "<Uy_MD>", "<Uz_MD>", "rho(C)_MD");
                printf("----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
            }

            printf("\r%8d " MAKE_LIGHT_BLUE "%10.4f " MAKE_BLUE "%10.4f " MAKE_LIGHT_BLUE "%9.5f %9.5f %9.5f " MAKE_BLUE "%9.5f %9.5f %9.5f "
            		MAKE_LIGHT_BLUE "%10.4f " MAKE_BLUE "%10.4f " MAKE_LIGHT_BLUE "%9.5f %9.5f %9.5f " MAKE_BLUE "%9.5f %9.5f %9.5f " MAKE_LIGHT_BLUE "%9.3f" MAKE_BLUE "%10.3f",
					fh->step_MD, st->std_rho_md, st->std_rho_fh, st->std_u_md[0], st->std_u_md[1], st->std_u_md[2], st->std_u_fh[0], st->std_u_fh[1], st->std_u_fh[2],
					st->std_rho_md_ccell, st->std_rho_fh_ccell, st->std_u_md_ccell[0], st->std_u_md_ccell[1], st->std_u_md_ccell[2], st->std_u_fh_ccell[0], st->std_u_fh_ccell[1], st->std_u_fh_ccell[2], fh->T_MD,fh->arr[fh->Ccell].T_cell);

            printf(RESET_COLOR "\n");

            fflush(stdout);
        }
    }
}

void fhmd_flow_n_avg(FHMD *fh, rvec x[], rvec v[], int N_atoms, t_commrec *cr)
{
    FH_arrays *arr = fh->arr;
    dvec       xn;
    double     S = fh->S;
    int        lr;

    const double eps = 1e-6;
    const double lrs = (double)(FHMD_FLOW_LAYERS);

    for(int k = 0; k < FHMD_FLOW_LAYERS; k++)
    {
        fh->avg_vel_n[k]   = 0;
        fh->avg_n_n[k]     = 0;
        fh->avg_vel_S_n[k] = 0;
        fh->avg_n_S_n[k]   = 0;
    }

    for(int n = 0; n < N_atoms; n++)
    {
        PBC(xn, x[n], fh->box);

        if(fh->S_function == moving_sphere)
            S = fhmd_Sxyz_r(x[n], fh->protein_com, fh);     // MD/FH sphere follows protein
        else if(fh->S_function == fixed_sphere)
            S = fhmd_Sxyz_r(x[n], fh->box05, fh);           // Fixed MD/FH sphere

        lr = (int)(xn[YY]/fh->box[YY]*lrs);

        if(lr < 0 || lr >= FHMD_FLOW_LAYERS)             // This should never happen... only if the coordinates are NaN.
        {
            printf(MAKE_RED "\nFHMD: ERROR: Solution diverged. Atom #%d coordinates: (%g, %g, %g) nm\n" RESET_COLOR "\n", n, xn[0], xn[1], xn[2]);
            exit(21);
        }

        fh->avg_n_n[lr]++;
        fh->avg_vel_n[lr] += fh->vel[n][XX];

        if(S < eps)
        {
            fh->avg_n_S_n[lr]++;
            fh->avg_vel_S_n[lr] += fh->vel[n][XX];
        }
    }

    if(PAR(cr))
    {
        gmx_sumi(FHMD_FLOW_LAYERS, fh->avg_n_n, cr);
        gmx_sumd(FHMD_FLOW_LAYERS, fh->avg_vel_n, cr);
        gmx_sumi(FHMD_FLOW_LAYERS, fh->avg_n_S_n, cr);
        gmx_sumd(FHMD_FLOW_LAYERS, fh->avg_vel_S_n, cr);
    }

    for(int k = 0; k < FHMD_FLOW_LAYERS; k++)
    {
        if(fh->avg_n_n[k])
            fh->avg_vel_tot_n[k]   += fh->avg_vel_n[k]/(double)(fh->avg_n_n[k]);
        if(fh->avg_n_S_n[k])
            fh->avg_vel_S_tot_n[k] += fh->avg_vel_S_n[k]/(double)(fh->avg_n_S_n[k]);
        fh->avg_n_tot_n[k]++;
    }
}

void fhmd_flow_m_avg(FHMD *fh, rvec x[], rvec v[],real mass[], int N_atoms, t_commrec *cr)
{
    FH_arrays *arr = fh->arr;
    MD_stat   *st  = &fh->stat;
    dvec       xn;
    double     S = fh->S;
    int        lr;
    st->n_flow_p++;

    const double eps = 1e-6;
    const double lrs = (double)(FHMD_FLOW_LAYERS);

    for(int k = 0; k < FHMD_FLOW_LAYERS; k++)
    {
        fh->avg_vel[k]   = 0;
        fh->avg_n[k]     = 0;
        fh->avg_m[k]     = 0;
        fh->avg_vel_S[k] = 0;
        fh->avg_n_S[k]   = 0;
        fh->avg_m_S[k]   = 0;
    }

    if (st->n_flow_p > fh->step_flow_output)
    {
        for(int n = 0; n < N_atoms; n++)
        {
            PBC(xn, x[n], fh->box);

            if(fh->S_function == moving_sphere)
                S = fhmd_Sxyz_r(x[n], fh->protein_com, fh);     // MD/FH sphere follows protein
            else if(fh->S_function == fixed_sphere)
                S = fhmd_Sxyz_r(x[n], fh->box05, fh);           // Fixed MD/FH sphere

            lr = (int)(xn[YY]/fh->box[YY]*lrs);

            if(lr < 0 || lr >= FHMD_FLOW_LAYERS)             // This should never happen... only if the coordinates are NaN.
            {
                printf(MAKE_RED "\nFHMD: ERROR: Solution diverged. Atom #%d coordinates: (%g, %g, %g) nm\n" RESET_COLOR "\n", n, xn[0], xn[1], xn[2]);
                exit(21);
            }

            fh->avg_n[lr]++;
            fh->avg_vel[lr] += mass[n]*fh->vel[n][XX];
            fh->avg_m[lr]+= mass[n];

            if(S < eps)
            {
                fh->avg_n_S[lr]++;
                fh->avg_vel_S[lr] += mass[n]*fh->vel[n][XX];
                fh->avg_m_S[lr]+= mass[n];
            }
        }

        if(PAR(cr))
        {
            gmx_sumi(FHMD_FLOW_LAYERS, fh->avg_n, cr);
            gmx_sumd(FHMD_FLOW_LAYERS, fh->avg_vel, cr);
            gmx_sumd(FHMD_FLOW_LAYERS, fh->avg_m, cr);
            gmx_sumi(FHMD_FLOW_LAYERS, fh->avg_n_S, cr);
            gmx_sumd(FHMD_FLOW_LAYERS, fh->avg_vel_S, cr);
            gmx_sumd(FHMD_FLOW_LAYERS, fh->avg_m_S, cr);
        }

        for(int k = 0; k < FHMD_FLOW_LAYERS; k++)
        {
            if(fh->avg_n[k])
                fh->avg_vel_tot[k]   += fh->avg_vel[k]/(double)(fh->avg_m[k]);
            if(fh->avg_n_S[k])
                fh->avg_vel_S_tot[k] += fh->avg_vel_S[k]/(double)(fh->avg_m_S[k]);
            fh->avg_n_tot[k]++;
        }
    }
}

void fhmd_flow_ufh_avg(FHMD *fh, t_commrec *cr)
{
    FH_arrays *arr = fh->arr;
    MD_stat   *st  = &fh->stat;
    ivec ind;
    int l = 0;
    st->n_flow++;
    double invn_flow = 1.0/(double)(st->n_flow);

    for(int l = 0; l < fh->N[1]; l++)
    {
        fh->ufh_flow[l]   = 0;
    }

    if (st->n_flow > fh->step_flow_output)
    {
        for(int k = 0; k < NZ; k++)
        {
        	for(int j = 0; j < NY; j++)
        	{
        		for(int i = 0; i < NX; i++)
        		{
        			ASSIGN_IND(ind, i, j, k);
        			for(int l = 0; l < fh->N[1]; l++)
        			{
        				if(l == j)
        				{
        					fh->ufh_flow[l] += arr[C].u_fh_flow[0];
        				}
        			}
        		}
        	}
        }

        if(PAR(cr))
        {
            gmx_sumd(fh->N[1], fh->ufh_flow, cr);
        }

        for(int l = 0; l < fh->N[1]; l++)
        {
            fh->ufh_flow_tot[l] += fh->ufh_flow[l]/fh->Ntot*fh->N[1]/4;
        }

        for(int l = 0; l < fh->N[1]; l++)
        {
            fh->avg_ufh_flow_tot[l]  = fh->ufh_flow_tot[l]*invn_flow;
        }
    }
}
