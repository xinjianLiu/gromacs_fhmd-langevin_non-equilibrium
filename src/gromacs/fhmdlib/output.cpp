#include "data_structures.h"
#include "macro.h"


#define write_dump(var, name) \
    sprintf(fname, "dump_%s.txt", name); \
    if(fh->step_MD == 0) { \
        fw = fopen(fname, "w"); \
        write_header(fw, fh); \
    } else { \
        fw = fopen(fname, "a"); \
    } \
    fprintf(fw, "\n%d\t", fh->step_MD); \
    for(int i = 0; i < fh->Ntot; i++) fprintf(fw, "%g\t", fh->arr[i].var); \
    fclose(fw);

#define write_dump_Ccell(var, name) \
    sprintf(fname, "dump_%s_Ccell.txt", name); \
    if(fh->step_MD == 0) { \
        fw = fopen(fname, "w"); \
        write_header_Ccell(fw, fh); \
    } else { \
        fw = fopen(fname, "a"); \
    } \
    fprintf(fw, "\n%d\t", fh->step_MD); \
    fprintf(fw, "%g\t", fh->arr[fh->Ccell].var); \
    fclose(fw);

void write_header(FILE *fw, FHMD *fh)
{
    fprintf(fw, "step\t");

    for(int k = 0; k < NZ; k++)
        for(int j = 0; j < NY; j++)
            for(int i = 0; i < NX; i++)
                fprintf(fw, "cell %d-%d-%d\t", i, j, k);
}

void write_header_Ccell(FILE *fw, FHMD *fh)
{
    fprintf(fw, "step\t");

    fprintf(fw, "cell %d\t", fh->Ccell);
}

void fhmd_dump_all(FHMD *fh)
{
    FILE *fw;
    char  fname[64];

    //write_dump(T_cell,      "T");
    write_dump(ro_prime,      "ro_prime");
    write_dump(m_prime[0],    "m_prime_X");
    write_dump(m_prime[1],    "m_prime_Y");
    write_dump(m_prime[2],    "m_prime_Z");
    write_dump(u_fh_flow[0], "u_fh_flow_X");
   /* write_dump(ro_md,   "ro_md");
    write_dump(ro_fh,   "ro_fh");
    write_dump(u_md[0], "u_md_X");
    write_dump(u_md[1], "u_md_Y");
    write_dump(u_md[2], "u_md_Z");
    write_dump(u_fh[0], "u_fh_X");
    write_dump(u_fh[1], "u_fh_Y");
    write_dump(u_fh[2], "u_fh_Z");
    write_dump(u_fh_flow[0], "u_fh_flow_X");
    write_dump(u_fh_flow[1], "u_fh_flow_Y");
    write_dump(u_fh_flow[2], "u_fh_flow_Z");
    write_dump(q_prime[0],  "q_prime_X");
    write_dump(q_prime[1],  "q_prime_Y");
    write_dump(q_prime[2],  "q_prime_Z");

    if (fh->step_MD == 0 )
    {
    	write_dump(Uflow[0],  "U_flow_Analytical");
    }*/

    /*
    write_dump(m_star[0],     "m_star_X");
    write_dump(m_star[1],     "m_star_Y");
    write_dump(m_star[2],     "m_star_Z");
    write_dump(MDS_RHO_N,   "MDS_RHO_N");
    write_dump(MDS_MOM_N[0], "MDS_MOM_N_X");
    write_dump(MDS_MOM_N[1], "MDS_MOM_N_Y");
    write_dump(MDS_MOM_N[2], "MDS_MOM_N_Z");
    write_dump(F_MD_S_AVA[0], "F_MD_S_AVA_X");
    write_dump(F_MD_S_AVA[1], "F_MD_S_AVA_Y");
    write_dump(F_MD_S_AVA[2], "F_MD_S_AVA_Z");
    write_dump(UURHO_S_AVA[0], "UURHO_S_AVA_X");
    write_dump(UURHO_S_AVA[1], "UURHO_S_AVA_Y");
    write_dump(UURHO_S_AVA[2], "UURHO_S_AVA_Z");*/

#ifdef FHMD_DEBUG
    write_dump(F_MD_S_AVA[0], "F_MD_S_AVA_X");
    write_dump(F_MD_S_AVA[1], "F_MD_S_AVA_Y");
    write_dump(F_MD_S_AVA[2], "F_MD_S_AVA_Z");
    write_dump(AVT_RHO_P,   "AVT_RHO_P");
    write_dump(AVT_RHO_C,   "AVT_RHO_C");
    write_dump(AVT_MOM_P[0], "AVT_MOM_P_X");
    write_dump(AVT_MOM_P[1], "AVT_MOM_P_Y");
    write_dump(AVT_MOM_P[2], "AVT_MOM_P_Z");
    write_dump(AVT_MOM_C[0], "AVT_MOM_C_X");
    write_dump(AVT_MOM_C[1], "AVT_MOM_C_Y");
    write_dump(AVT_MOM_C[2], "AVT_MOM_C_Z");
    write_dump(f_fh[0],       "f_fh_X");
    write_dump(f_fh[1],       "f_fh_Y");
    write_dump(f_fh[2],       "f_fh_Z");
    write_dump(alpha_term[0], "alpha_term_X");
    write_dump(alpha_term[1], "alpha_term_Y");
    write_dump(alpha_term[2], "alpha_term_Z");
    write_dump(beta_term[0],  "beta_term_X");
    write_dump(beta_term[1],  "beta_term_Y");
    write_dump(beta_term[2],  "beta_term_Z");
    write_dump(grad_ro[0],    "grad_ro_X");
    write_dump(grad_ro[1],    "grad_ro_Y");
    write_dump(grad_ro[2],    "grad_ro_Z");

    write_dump(ro_prime,      "ro_prime");
    write_dump(ro_star,       "ro_star");
    write_dump(m_prime[0],    "m_prime_X");
    write_dump(m_prime[1],    "m_prime_Y");
    write_dump(m_prime[2],    "m_prime_Z");
    write_dump(m_star[0],     "m_star_X");
    write_dump(m_star[1],     "m_star_Y");
    write_dump(m_star[2],     "m_star_Z");

    write_dump(S,             "S");
    write_dump(Sf[0],         "Sf_X");
    write_dump(Sf[1],         "Sf_Y");
    write_dump(Sf[2],         "Sf_Z");
#endif
}

void fhmd_dump_Ccell(FHMD *fh)
{
    FILE *fw;
    char  fname[64];
    write_dump_Ccell(ro_md,   "ro_md");
    write_dump_Ccell(ro_fh,   "ro_fh");
    write_dump_Ccell(u_md[0], "u_md_X");
    write_dump_Ccell(u_md[1], "u_md_Y");
    write_dump_Ccell(u_md[2], "u_md_Z");
    write_dump_Ccell(u_fh[0], "u_fh_X");
    write_dump_Ccell(u_fh[1], "u_fh_Y");
    write_dump_Ccell(u_fh[2], "u_fh_Z");
}

void fhmd_flow_n_avg_write(FHMD *fh)
{
    FILE *fw;

    const double eps = 1e-15;

    if(fh->step_MD == 0)
    {
        fw = fopen("flow_n.txt", "w");
    } else {
        fw = fopen("flow_n.txt", "a");
    }

    fprintf(fw, "MD step = %d\n", fh->step_MD);
    fprintf(fw, "layer\ty, nm\t<Ux>, nm/ps\t<Ux> (MD), nm/ps\tUx, nm/ps\tUx (MD), nm/ps\n");

    for(int k = 0; k < FHMD_FLOW_LAYERS; k++)
    {
        fprintf(fw, "%d\t%g\t%g\t%g\t%g\t%g\n", k, ((double)(k) + 0.5)/(double)(FHMD_FLOW_LAYERS)*fh->box[YY],
                fh->avg_vel_tot_n[k]/fh->avg_n_tot_n[k], fh->avg_vel_S_tot_n[k]/fh->avg_n_tot_n[k],
                fh->avg_vel_n[k]/((double)(fh->avg_n_n[k]) + eps), fh->avg_vel_S_n[k]/((double)(fh->avg_n_S_n[k]) + eps));
    }

    fprintf(fw, "\n");
    fclose(fw);
}

void fhmd_flow_m_avg_write(FHMD *fh)
{
    FILE *fw;

    const double eps = 1e-15;

    if(fh->step_MD == 0)
    {
        fw = fopen("flow_m.txt", "w");
    } else {
        fw = fopen("flow_m.txt", "a");
    }

    fprintf(fw, "MD step = %d\n", fh->step_MD);
    fprintf(fw, "layer\ty, nm\t<Ux>, nm/ps\t<Ux> (MD), nm/ps\tUx, nm/ps\tUx (MD), nm/ps\n");

    for(int k = 0; k < FHMD_FLOW_LAYERS; k++)
    {
        fprintf(fw, "%d\t%g\t%g\t%g\t%g\t%g\n", k, ((double)(k) + 0.5)/(double)(FHMD_FLOW_LAYERS)*fh->box[YY],
                fh->avg_vel_tot[k]/fh->avg_n_tot[k], fh->avg_vel_S_tot[k]/fh->avg_n_tot[k],
                fh->avg_vel[k]/((double)(fh->avg_m[k]) + eps), fh->avg_vel_S[k]/((double)(fh->avg_m_S[k]) + eps));
    }

    fprintf(fw, "\n");
    fclose(fw);
}

void fhmd_flow_ufh_avg_write(FHMD *fh)
{
    FILE *fw;

    const double eps = 1e-15;

    if(fh->step_MD == 0)
    {
        fw = fopen("flow_ufh.txt", "w");
    }
    else
    {
        fw = fopen("flow_ufh.txt", "a");
    }

    fprintf(fw, "MD step = %d\n", fh->step_MD);
    fprintf(fw, "layer\t Ufh_flow , nm/ps\n");

    for(int l = 0; l < fh->N[1]; l++)
    {
        fprintf(fw, "%d\t%g\n", l, fh->avg_ufh_flow_tot[l]);
    }

    fprintf(fw, "\n");
    fclose(fw);
}

void writePoint(FILE *fw,FHMD *fh, const char *line, char ch)
{
	int c = 0;
    ivec ind;

    fprintf(fw, "%s", line);

    for(int k = 0; k <= NZ; k++) {
        for(int j = 0; j <= NY; j++) {
            for(int i = 0; i <= NX; i++) {
                ASSIGN_IND(ind, i - (int)(i/NX), j - (int)(j/NY), k - (int)(k/NZ));
                     if(ch == 'X' && i < NX)  fprintf(fw, "%e ", fh->grid.n[C][0]);
                else if(ch == 'X' && i == NX) fprintf(fw, "%e ", fh->grid.n[C][0] + fh->grid.h[C][0]);
                else if(ch == 'Y' && j < NY)  fprintf(fw, "%e ", fh->grid.n[C][1]);
                else if(ch == 'Y' && j == NY) fprintf(fw, "%e ", fh->grid.n[C][1] + fh->grid.h[C][1]);
                else if(ch == 'Z' && k < NZ)  fprintf(fw, "%e ", fh->grid.n[C][2]);
                else if(ch == 'Z' && k == NZ) fprintf(fw, "%e ", fh->grid.n[C][2] + fh->grid.h[C][2]);
                if(c++ == 10) {fprintf(fw, "\n"); c = 0;}
            }
        }
    }
}

void writeCons(FILE *fw,FHMD *fh, const char *line, char ch)
{
	int c = 0;
    ivec ind;

    fprintf(fw, "%s", line);

    for(int k = 0; k < NZ; k++) {
        for(int j = 0; j < NY; j++) {
            for(int i = 0; i < NX; i++) {
                ASSIGN_IND(ind, i, j, k);
                     if(ch == 'U') fprintf(fw, "%e ", fh->arr[C].u_fh_flow[0]);
                else if(ch == 'V') fprintf(fw, "%e ", fh->arr[C].u_fh_flow[1]);
                else if(ch == 'W') fprintf(fw, "%e ", fh->arr[C].u_fh_flow[2]);
                else if(ch == 'R') fprintf(fw, "%e ", fh->arr[C].ro_fh);
                if(c++ == 10) {fprintf(fw, "\n"); c = 0;}
            }
        }
    }
}

// TODO: This function is a copy/paste from the previous code -- should be rewritten completely
void fhmd_write_tecplot_data(FHMD *fh, int step, double time)
{
    FILE *fw;

    static int zc = 0;

    char fname[64];

    sprintf(fname, "tecplot/data_%7.7d.dat", step);             // Tecplot data filename

    if((fw = fopen(fname, "w")) == NULL)                        // Open file
        printf("\n ERROR creating %s for output!\n", fname);

    zc++;

    fprintf(fw, "TITLE=\"OUT\"\nVARIABLES=\"X\",\"Y\",\"Z\",\"U_tilde\",\"V_tilde\",\"W_tilde\",\"RHO_tilde\"\n");
    fprintf(fw, "ZONE T=\"%f\", I=%d, J=%d, K=%d, DATAPACKING=BLOCK\nVARLOCATION=([4-7]=CELLCENTERED)\n", time, NX+1, NY+1, NZ+1);
    fprintf(fw, "STRANDID=%d\nSOLUTIONTIME=%g", zc, time);


    writePoint(fw,fh, "\n\n# X:\n", 'X');                        // Write X
    writePoint(fw,fh, "\n\n# Y:\n", 'Y');                        // Write Y
    writePoint(fw,fh, "\n\n# Z:\n", 'Z');                        // Write Z


    writeCons(fw,fh, "\n\n# U:\n",   'U');                       // Write U
    writeCons(fw,fh, "\n\n# V:\n",   'V');                       // Write V
    writeCons(fw,fh, "\n\n# W:\n",   'W');                       // Write W
    writeCons(fw,fh, "\n\n# RHO:\n", 'R');                       // Write RHO

   	fclose(fw);
}


#ifdef FHMD_PARAVIEW
void write_paraview_data(HHMD *hhmd) {

    const int ln = 5;       // Numbers per one line

    FILE *fout;
    char fname[32];
    static int zc = 0;
    int lc = 1;

    int nx = hhmd->cellNum.x;
    int ny = hhmd->cellNum.y;
    int nz = hhmd->cellNum.z;

    sprintf(fname, "paraview_%d.vtk", zc);

    if((fout = fopen(fname, "w")) == NULL)
        printf("\n ERROR creating %s for output!\n", fname);

    zc++;

    fprintf(fout, "# vtk DataFile Version 3.0\nData Output\nASCII\nDATASET STRUCTURED_GRID\nDIMENSIONS %d %d %d\n", nx + 1, ny + 1, nz + 1);
    fprintf(fout, "POINTS %d float\n", (nx + 1)*(ny + 1)*(nz + 1));

    for(int k = 0; k < nz + 1; k++) {
        for(int j = 0; j < ny + 1; j++) {
            for(int i = 0; i < nx + 1; i++) {
                fprintf(fout, "%f %f %f\n", hhmd->grid.n[i].x, hhmd->grid.n[j].y, hhmd->grid.n[k].z);
            }
        }
    }

    fprintf(fout, "CELL_DATA %d\n", nx*ny*nz);
    fprintf(fout, "FIELD FieldData %d\n", 1);
    fprintf(fout, "Density %d %d float\n", 1, nx*ny*nz);

    for(int k = 0; k < nz; k++) {
        for(int j = 0; j < ny; j++) {
            for(int i = 0; i < nx; i++) {
                fprintf(fout, "%f", hhmd->arr->ro_fh[i][j][k]);
                if(lc++ >= ln) {
                    fprintf(fout, "\n");
                    lc = 1;
                } else {
                    fprintf(fout, " ");
                }
            }
        }
    }

    fprintf(fout, "VECTORS Velocity float\n");

    lc = 1;

    for(int k = 0; k < nz; k++) {
        for(int j = 0; j < ny; j++) {
            for(int i = 0; i < nx; i++) {
                fprintf(fout, "%f %f %f", hhmd->arr->u_fh[i][j][k].x, hhmd->arr->u_fh[i][j][k].y, hhmd->arr->u_fh[i][j][k].z);
                if(lc++ >= ln) {
                    fprintf(fout, "\n");
                    lc = 1;
                } else {
                    fprintf(fout, "   ");
                }
            }
        }
    }

    fprintf(fout, "POINT_DATA %d\n", (nx + 1)*(ny + 1)*(nz + 1));

    fclose(fout);

}
#endif

