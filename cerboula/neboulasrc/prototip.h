void read_param(void);
void read_flow_conditions(double *,double *,double *,double *,double *,double *,double *,double *,
                          double **,double *,double *,double **,double *,double *,
                          double *,int *,double *,double **,
                          double ***,double **,double *,double *,double *,
                          double *,int **);
void write_solution(double *,double *,double *,
                    double *,double *,double *,double **,double,double,
                    double,double,double,double,int,int,double,
                    double,double,double *,double *,double **,double *, int);
void grid_gen(double *,double *,double *,double *,double *,double *,double *,
              double *,double *,double *);
void spline_gen(double *,double *,double *,double *,double *,double *,
                double *,double **,double *,double *,double *,double *,
                double *,double *,double *,double *,double **,double *,double *,
                double *,double *,double *,double *,double **,double **);
void search(double,double *,int *,int *);
void var_interp(double *,double *,double *,double *,double *,double *,double *,
                double **,double *,double *,double *,double *,double *,double *,
                double *,double *,double **,double *,double *,double,
                double *,double,double *,double *,double *,double *,double **,double **,double *);
void new_coordinate(double,double *,double *,double *,double *,
                    double *,double *,double *);
void reduce_stepsize(double *,double *);
void stepsize_control(int,int,int,double *);
void coefficients(double *,double ***,double *,double **,double ***,
                  double **,double *,double *,double *,double *,double **,
                  double **,double *,double *,double *,double *,
                  double,double,double);/*,int/*NANNI);*/
void der1(double *,double *,double **,double **,double **);
void der2(double *,double *,double *,double *,double **,double **);
void stef_max(double **,double **,double *,double ***,double *,double *,
              double *,double **,double **,double,double,int,
              void (*)(double ,double *,double *,double *,double **,
                       double ,double ,double *,double *,double) );
void csi_expansion(double,double,double *,double *,double *);
void continuity(double *,double *);
void momentum(double *,double *,double *,double *,
              double *,double *,double *,double,double);
void energy(double *,double *,double *,double *,double **,double **,double *,
            double *,double *,double *,double *,double **,double **,double **,
            double **,double *,double,double,int,double,int/*NANNI*/);
void species(double **,double **,double **,double **,double *,double *,double **,
             double **,double ***,double *,double *,double **,double ***,
             double *,double *,double *,double *,double **,
             double,double,int,double,double *,double **,double *,int /*NANNI*/);

double *vector(int);
int *intvector(int);

double **matrix(int,int);
double ***trimatrix(int,int,int);
int **matrix_int(int,int);       /*NANNI*/
void ent2temp(double *,double *,double **,double,double);
double mixt_enthalpy(double *,double,double *);
double convergence(int,int,double *,double *,double *,double *,double **,
                   double **);
void inizialize_data(double *,double *,double *,double *,double *,double *,double *,
                     double *,double **,double **,double **,double **,
                     double *,double,double *,double *,double *,double *);
char *stringa(int);
void skin_and_heat(double,double,int,double *,double *,double *,double **,double **,
                   double *,double *,double *,double,double,double,int,
                   double *,double *,double */*NANNI*/);
void ableit(double *,double *,double,int);
void free_vector(double *);

void free_intvector(int *);

void set_1_(int*,int *,int *,int *,int *,int *,double *,double *,double *,double *,
            int *);
void set_2_(double *,double *,double *,int *);
void under_relax(int,int,double *,double *,double *,double *,double **,double **,double);
void finite_thickness_chi(double *,double,double *,double *);
/*void species_eq(double **,double **,double **,double **,double *,double *,double **,
             double **,double ***,double *,double *,double **,double ***,
             double *,double *,double *,double *,double **,
             double,double,int,double,double *,double **,double * /*NANNI,
	     double *, double *, int, double *); */

void elements_eq(double **, double **,double **,double **,double *,double *,
                   double **,double **,double ***,double *,double *,double **,
	           double ***,double *,double *,double *,double *,double **,
		   double ,double ,int ,double ,double *, double **, double *,
		   double *,int);/*, double *);*/


/*void species_eq_test(double **,double **,double **,double **,double *,double *,double **,
             double **,double ***,double *,double *,double **,double ***,
             double *,double *,double *,double *,double **,
             double,double,int,double,double *,double **,double * /*NANNI);*/
/*void find_ext_eq_conc(void);   /*NANNI*/
/*void read_nanni(void);   /*NANNI*/


#ifdef __GNUCO__

#define set_1_ set_1__
#define set_2_ set_2__
#endif

