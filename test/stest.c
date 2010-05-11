/*
 * FEPC
 * Copyright (C) 2010 Stefan Handschuh (handschu@mis.mpg.de)
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
 
#include "fepc_easy.h"
#include "fepc_easy_helper.h"

fepc_real_t norm_x_sqr(vec_real_p x) {
    int n;
    
    fepc_real_t result = 0.;
    
    for (n = 0; n < x->dim; n++) {
        result += x->array[n]*x->array[n];
    }
    return result;
}

fepc_real_t const_2dim(vec_real_p x) {
    return x->dim*2.;
}

fepc_real_t two_xone(vec_real_p x) {
    return x->array[1]*2;
}

fepc_real_t nine_xone_sqr(vec_real_p x) {
    return 9*x->array[0]*x->array[0];
}

fepc_real_t function(vec_real_p x) {
    return 3*x->array[0]*x->array[0]*x->array[0] + x->array[1]*x->array[1]*x->array[1]*x->array[1];

}

fepc_real_t function_laplace(vec_real_p x) {
    return 18*x->array[0]+12*x->array[1]*x->array[1];
}
 


int main() {
    fepc_real_t stepping = 0.001;
    
    interval_p * intervals = (interval_p*) malloc(sizeof(interval_p)*1);
   
    int n, dimension = 2;
    
    intervals[0] = interval_new(dimension);
    
    for (n = 0; n < dimension; n++) {
        intervals[0]->start[n] = -.01;
        intervals[0]->end[n] = .01;
    }

    //func_p func_const_2dim =  setup_fepc_structure(const_2dim, intervals, 1, 0, stepping);
    
    
    //func_p func_norm_x_sqr =  setup_fepc_structure(norm_x_sqr, intervals, 1, 0, stepping);
    
    //func_p func_2x1 = setup_fepc_structure(two_xone, intervals, 1, 0, stepping);
    
    //func_p func_9x1_sqr = setup_fepc_structure(nine_xone_sqr, intervals, 1, 0, stepping);
    
    func_p func_function = setup_fepc_structure(function, intervals, 1, 0, stepping);
    
    printf("Funktion 1 fertig\n");
    func_p func_function_laplace = setup_fepc_structure(function_laplace, intervals, 1, 0, stepping);
    printf("Funktion 2 fertig\n");
    func_p func_function_laplace_calc = func_laplace(func_function, stepping);
    printf("Laplace fertig\n");
    
    //func_p func_laplace_norm_x_sqr = func_laplace(func_norm_x_sqr, stepping);
    
    //func_print(func_function, 3);
    //func_print(func_derive(func_function, 0, stepping), 3);
    //func_print(func_9x1_sqr, 3);
    
    //func_print(func_function_laplace_calc, 3);
    //func_print(func_function_laplace, 3);
    
    func_p temp = func_subtract(func_function_laplace_calc, func_function_laplace);
    
    
    fepc_real_t norm = func_norm_l2_sqr(temp, stepping);
    
    
    
    /*func_del(func_const_6);
    func_del(func_norm_x_sqr);
    func_del(func_laplace_norm_x_sqr);
    func_del(temp);
    */
    printf("Norm^2 = %1.15f\n", norm);

    return 0;
}
