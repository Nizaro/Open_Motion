/*
 * calibration.c
 *
 *  Created on: 30 May 2017
 *      Author: thomas
 */


#include "calibration.h"

/////////////////////////////////////////////
///////   Adaptive least square         /////
/////////////////////////////////////////////



void mat_to_vecs(struct omMatrix* A,struct omVector* vecs_A){

	int n = A->_columns;

	for(int i=0;i<n;++i){
		for(int j=0;j<n;++j){

			if( i <= j){
				double value = om_matrix_getValue(A, i, j);
				om_vector_setValue(vecs_A, ((j*(j+1))/2) + i, value);
			}

		}
	}

}


void vecs_to_mat(struct omVector* vecs_A,struct omMatrix* A){

	int nb = vecs_A->_length;
	int n = ((int)sqrt(1.0 + 8*nb) - 1)/2;

	for(int i=0;i<n;++i){
		for(int j=0;j<n;++j){

			if( i <= j){
				double value = om_vector_getValue(vecs_A, ((j*(j+1))/2) + i);
				om_matrix_setValue(A, i, j,value);
				om_matrix_setValue(A, j, i,value);
			}

		}
	}

}


double tensor_T(int k,int i,int l,omMatrix* X,double var){

	double x = om_matrix_getValue(X, i, l);
	double t_x = 0.0;

	switch(k){

		case 0:
			t_x = 1.0;
			break;

		case 1:
			t_x = x;
			break;

		case 2:
			t_x = (x*x) - var;
			break;

		case 3:
			t_x = (x*x*x) - (3.0*x*var);
			break;

		case 4:
			t_x = (x*x*x*x) - (6.0*x*x*var) + (3.0*var*var);
			break;

	}

	return t_x;
}


int tensor_R(int p,int q,int i,omMatrix* M){


	int res=0;

	res += (int)om_matrix_getValue(M, p, 0) == i ? 1 : 0;
	res += (int)om_matrix_getValue(M, p, 1) == i ? 1 : 0;
	res += (int)om_matrix_getValue(M, q, 0) == i ? 1 : 0;
	res += (int)om_matrix_getValue(M, q, 1) == i ? 1 : 0;

	return res;

}



void create_M(omMatrix* M,int n){

	int n_beta =  ((n*(n+1))/2) + n + 1 ;

	omMatrix A;
	omMatrix I;
	omMatrix A_t;
	omMatrix I_t;
	omMatrix A_tI;
	omMatrix I_tA;

	omVector vecs_I_tA;
	omVector vecs_A_tI;

	om_matrix_create(&A,1,n+1);
	om_matrix_create(&A_t,n+1,1);
	om_matrix_create(&I,1,n+1);
	om_matrix_create(&I_t,n+1,1);
	om_matrix_create(&A_tI,n+1,n+1);
	om_matrix_create(&I_tA,n+1,n+1);
	om_vector_create(&vecs_I_tA, n_beta);
	om_vector_create(&vecs_A_tI, n_beta);


	for(int i=0;i<n+1;i++){
		om_matrix_setValue(&A,0,i,1);
		om_matrix_setValue(&I,0,i,(double)i + 1.0);
	}
	om_matrix_setValue(&I,0,n,0.0);

	om_matrix_transpose(&A, &A_t);
	om_matrix_transpose(&I, &I_t);

	om_operator_matrix_mul(&A_t, &I, &A_tI);
	om_operator_matrix_mul(&I_t, &A, &I_tA);

	mat_to_vecs(&I_tA, &vecs_I_tA);
	mat_to_vecs(&A_tI, &vecs_A_tI);

	/*/
	printf("\n\nMatrix I_tA\n");om_matrix_display(&I_tA);
	printf("\nvecs_I_tA = ");om_vector_display(&vecs_I_tA);

	printf("\n\nMatrix A_tI\n");om_matrix_display(&A_tI);
	printf("\nvecs_A_tI = ");om_vector_display(&vecs_A_tI);
	//*/

	om_matrix_create(M,n_beta,2);

	for(int i=0;i<n_beta;i++){

		om_matrix_setValue(M, i, 0, om_vector_getValue(&vecs_I_tA, i));
		om_matrix_setValue(M, i, 1, om_vector_getValue(&vecs_A_tI, i));
	}

	om_matrix_free(&A);
	om_matrix_free(&A_t);
	om_matrix_free(&A_tI);
	om_matrix_free(&I);
	om_matrix_free(&I_t);
	om_matrix_free(&I_tA);
	om_vector_free(&vecs_I_tA);
	om_vector_free(&vecs_A_tI);

}



double function_n_als(int p,int q,omMatrix* M,omMatrix* X){

	int n = X->_rows;
	int m = X->_columns;

	double n_als = 0.0;

	for(int l=0;l<m;l++){

		double tmp = 1.0;

		for(int i=0;i<n;i++){

			//if(tensor_R(p, q, i+1, M) == 0)
			//	printf("R(%d,%d,%d) = %d\n",p,q,i,tensor_R(p, q, i+1, M));

			tmp *= tensor_T(tensor_R(p, q, i+1, M), i, l, X, 0.0075);
		}

		n_als += tmp;
	}

	return n_als;
}



int in_D(int integer){

	if(integer == 1 || integer == 3 || integer == 4 )
		return 1;
	else
		return 0;

}



void om_calibration_adjusted_least_square(omMatrix* X,omMatrix* Q,omVector* b,double* d,double Hm_square){


	int n=3;
	int n_beta =  ((n*(n+1))/2) + n + 1 ;
	omMatrix M;
	omMatrix Phi;

	create_M(&M,n);

	//printf("\n\nMatrix M\n");om_matrix_display(&M);

	om_matrix_create(&Phi, n_beta, n_beta);

	for(int p = 0;p<n_beta;p++){
		for(int q = 0;q<n_beta;q++){

			if( p <= q) {
				double value;

				if(in_D(p) == 1 && in_D(q) == 1){
					value = 4.0 * function_n_als(p, q,  &M, X);
				}else if(in_D(p) == 0 && in_D(q) == 0){
					value = function_n_als(p, q,  &M, X);
				}else{
					value = 2.0*function_n_als(p, q,  &M, X);
				}

				om_matrix_setValue(&Phi, p, q, value);
				om_matrix_setValue(&Phi, q, p, value);

			}

		}
	}

	omVector* eigen_vectors;
	double* eigen_values;
	om_matrix_getEingenValues(&Phi,&eigen_vectors,&eigen_values,50);


	int index_min;
	double min = 99999999999999999.9;

	printf("\n\n");
	for(int i=0;i<n_beta;i++){

		//printf("values %d = %f and vector ",i, eigen_values[i]);om_vector_display(&eigen_vectors[i]);

		if(min > eigen_values[i]){
			min = eigen_values[i];
			index_min = i;
		}

	}

	om_vector_normalize(&eigen_vectors[index_min]);
	//printf("\nbeta_als = ");om_vector_display(&eigen_vectors[index_min]);

	//om_operator_vector_scal_mul(&eigen_vectors[index_min], Hm_square, &eigen_vectors[index_min]);


	// beta_als to solution

	omMatrix M2;
	omMatrix M_inv;
	omVector vecs_beta_als_Q;
	omVector c;
	omVector Mc;

	om_matrix_create(&M2, n, n);
	om_matrix_create(&M_inv, n, n);
	om_vector_create(&vecs_beta_als_Q, (n*(n+1)/2));
	om_vector_create(&c, n);
	om_vector_create(&Mc, n);


	for(int i=0;i<(n*(n+1)/2);i++)
		om_vector_setValue(&vecs_beta_als_Q, i, om_vector_getValue(&eigen_vectors[index_min], i));

	vecs_to_mat(&vecs_beta_als_Q, &M2);

	om_matrix_inverse(&M2, &M_inv);


	for(int i=0;i<3;i++)
		om_vector_setValue(b, i, om_vector_getValue(&eigen_vectors[index_min],(n*(n+1)/2) + i));

	double beta = om_vector_getValue(&eigen_vectors[index_min],(n*(n+1)/2) + n);

	om_operator_matrix_vector_mul(&M_inv, b, &c);
	om_operator_vector_scal_div(&c, (-2.0), &c);

	om_operator_matrix_vector_mul(&M2, &c, &Mc);

	double theta = (om_vector_dotProduct(&c, &Mc) - beta);
	double gamma = Hm_square/ theta;

	om_operator_matrix_scal_mul(&M2, gamma, Q);
	om_operator_vector_scal_mul(b, gamma, b);

	(*d) = Hm_square * (  1.0 + (theta * (beta - 1.0) )) / theta;

	om_matrix_free(&M_inv);
	om_matrix_free(&M);
	om_matrix_free(&M2);
	om_matrix_free(&Phi);
	om_vector_free(&Mc);
	om_vector_free(&vecs_beta_als_Q);
	om_vector_free(&c);


}




/////////////////////////////////////////////
///////   Variance estimation           /////
/////////////////////////////////////////////




void init_fda(omVector** fda){


	//The idea here is to form a linear combination of successive elements
	//of the series. If the underlying form is locally nearly linear, then
	//a [1 -2 1] combination (for equally spaced data) will leave only
	//the noise remaining. Next, if we assume the measurement noise was
	//iid, N(0,s^2), then we can try to back out the noise variance.

	int nfda = 6;
	(*fda) = (omVector*)malloc(nfda * sizeof(omVector));


	om_vector_create(&(*fda)[0], 2,1.0,-1.0);
	om_vector_create(&(*fda)[1], 3,1.0,-2.0,1.0);
	om_vector_create(&(*fda)[2], 4,1.0,-3.0,3.0,-1.0);
	om_vector_create(&(*fda)[3], 5,1.0,-4.0,6.0,-4.0,1.0);
	om_vector_create(&(*fda)[4], 6,1.0,-5.0,10.0,-10.0,5.0,-1.0);
	om_vector_create(&(*fda)[5], 7,1.0,-6.0,15.0,-20.0,15.0,-6.0,1.0);

	for(int i = 0;i< nfda;i++){
		om_vector_normalize(&(*fda)[i]);

	}


}


void init_z(omVector* perc,omVector* z){

	//compute an interquantile range, like the distance between the 25%
	//and 75% points. This trims off the trash at each end, potentially
	//corrupted if there are discontinuities in the curve. It also deals
	//simply with a non-zero mean in this data. Actually do this for
	//several different interquantile ranges, then take a median.
	//NOTE: While I could have used other methods for the final variance
	//estimation, this method was chosen to avoid outlier issues when
	//the curve may have isolated discontinuities in function value or
	// a derivative.

	//The following points correspond to the central 90, 80, 75, 70, 65,
	//60, 55, 50, 45, 40, 35, 30, 25, and 20 percent ranges.

	om_vector_create(z,14);
	om_vector_create(perc,14);

	for(int i=0;i<14;i++){

		double perc_i = i == 0 ? 0.05 : 0.1 + ((double)(i-1)*0.025);
		double z_i = om_maths_erfinv( (1.0 - perc_i)*2.0 - 1.0 )*sqrt(2.0);

		om_vector_setValue(z, i, z_i);
		om_vector_setValue(perc, i, perc_i);

	}


}



void vec_change(omVector* m, int a, int b)

{

    double temp = om_vector_getValue(m,a);

    om_vector_setValue(m,a,om_vector_getValue(m, b));
    om_vector_setValue(m,b,temp);



}

void vec_quickSort(omVector* m, int debut, int fin){

    int gauche = debut-1;
    int droite = fin+1;

    const double pivot = om_vector_getValue(m,debut);//(*tableau)[debut];

    if(debut >= fin)
        return;


    while(1)
    {

        do droite--; while(/*(*tableau)[droite]*/ om_vector_getValue(m,droite)  > pivot);
        do gauche++; while(/*(*tableau)[gauche]*/ om_vector_getValue(m,gauche)  < pivot);


        if(gauche < droite)

        	vec_change(m, gauche, droite);

        else break;

    }


    /* Maintenant, tous les éléments inférieurs au pivot sont avant ceux
       supérieurs au pivot. On a donc deux groupes de cases à trier. On utilise
       pour cela... la méthode quickSort elle-même ! */
    vec_quickSort(m, debut, droite);
    vec_quickSort(m, droite+1, fin);

}

void mat_change(omMatrix* m, int a, int b){

    double temp = om_matrix_getValue(m,0,a);

    om_matrix_setValue(m,0,a,om_matrix_getValue(m, 0, b));
    om_matrix_setValue(m,0,b,temp);



}

void mat_quickSort(omMatrix* m, int debut, int fin){

    int gauche = debut-1;
    int droite = fin+1;

    const double pivot = om_matrix_getValue(m,0,debut);//(*tableau)[debut];

    if(debut >= fin)
        return;


    while(1)
    {

        do droite--; while(/*(*tableau)[droite]*/ om_matrix_getValue(m,0,droite)  > pivot);
        do gauche++; while(/*(*tableau)[gauche]*/ om_matrix_getValue(m,0,gauche)  < pivot);


        if(gauche < droite)

        	mat_change(m, gauche, droite);

        else break;

    }


    /* Maintenant, tous les éléments inférieurs au pivot sont avant ceux
       supérieurs au pivot. On a donc deux groupes de cases à trier. On utilise
       pour cela... la méthode quickSort elle-même ! */
    mat_quickSort(m, debut, droite);
    mat_quickSort(m, droite+1, fin);

}



void om_calibration_noise_estimation(omMatrix* X,omVector* v){

	omVector* fda;
	init_fda(&fda);

	omVector perc;
	omVector z;
	init_z(&perc, &z);



	for(int k=0;k<3;k++){

		omMatrix Xp;
		om_matrix_create(&Xp, 1, X->_columns);
		for(int i=0;i<X->_columns;i++)
			om_matrix_setValue(&Xp, 0, i, om_matrix_getValue(X, k, i));

		omVector sigmaest;
		om_vector_create(&sigmaest,6);
		for(int l=0;l<6;++l){

			omMatrix fda_l;
			om_matrix_create(&fda_l, 1, fda[l]._length);
			for(int i=0;i<fda[l]._length;i++)
				om_matrix_setValue(&fda_l, 0, i, om_vector_getValue(&fda[l], i));

			omMatrix noisedata;
			om_matrix_convolution2D_valid(&Xp, &fda_l, &noisedata);


			if(noisedata._columns > 2){

				mat_quickSort(&noisedata, 0, noisedata._columns - 1);

				omVector p;
				omVector noisedata2;

				p._length = noisedata._columns;
				noisedata2._length = noisedata._columns;

				p._values = (double*)malloc(noisedata._columns*sizeof(double));
				noisedata2._values = (double*)malloc(noisedata._columns*sizeof(double));

				for(int i=0;i<noisedata._columns;i++){
					om_vector_setValue(&p,i,  (0.5 + (double)(i+1))/(double)(noisedata._columns + 0.5));
					om_vector_setValue(&noisedata2,i,om_matrix_getValue(&noisedata,0,i));
				}

				omVector Q;
				om_vector_create(&Q,14);

				for(int j=0;j<14;j++){
					double q_j = om_vector_interpolation(&p, &noisedata2, 1.0 - om_vector_getValue(&perc, j)) - om_vector_interpolation(&p, &noisedata2, om_vector_getValue(&perc, j));
					q_j /= 2.0 * om_vector_getValue(&z, j);
					om_vector_setValue(&Q, j,q_j);
				}

				//printf("Q = ");om_vector_display(&Q);
				vec_quickSort(&Q,0,13);
				double median_q = (om_vector_getValue(&Q, 6) + om_vector_getValue(&Q, 7))/2.0;
				om_vector_setValue(&sigmaest, l, median_q);

				om_vector_free(&p);
				om_vector_free(&noisedata2);
				om_vector_free(&Q);

			}

			om_matrix_free(&noisedata);
			om_matrix_free(&fda_l);

		}

		//printf("sigmaest = ");om_vector_display(&sigmaest);
		vec_quickSort(&sigmaest,0,5);
		double var = (om_vector_getValue(&sigmaest, 2) + om_vector_getValue(&sigmaest, 3))/2.0;
		var *= var;
		var /= 1.0 + (15.0 * pow((double)(X->_columns) + 1.225,-1.245));

		om_vector_setValue(v, k, var);

		om_vector_free(&sigmaest);
		om_matrix_free(&Xp);

	}

}


