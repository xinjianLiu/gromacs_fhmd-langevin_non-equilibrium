#ifndef FHMD_OUTPUT_H_
#define FHMD_OUTPUT_H_

void fhmd_dump_all(FHMD *fh);
void fhmd_dump_Ccell(FHMD *fh);
void fhmd_flow_n_avg_write(FHMD *fh);
void fhmd_flow_m_avg_write(FHMD *fh);
void fhmd_flow_ufh_avg_write(FHMD *fh);
void fhmd_write_tecplot_data(FHMD *fh, int step, double time);

#endif /* FHMD_OUTPUT_H_ */
