/*
 * FEPC
 * Copyright (C) 2009 Peter Gerds (gerds@mis.mpg.de)
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

#include "faltung.h"

/**********************************************************
 *
 * lokale Funktionen (Def.)
 *
 **********************************************************/
void
func_plot_1d(func_p f, char dateiname1[20], fepc_real_t mesh);

void
func_plot_2d(func_p f, char dateiname1[20], fepc_real_t mesh);

/**********************************************************
 *
 * Hauptfunktion
 *
 **********************************************************/


int main(int argc, char* argv[]) {
	func_p  f, g, w, fepc, ref;
	fepc_real_t  norm;
	fepc_real_t  a,b,c, mesh, temp1, temp2;
	int  dim, i,j , l, pos;
	vec_p  grad, start, lang, r, s;
	vec_real_p  h;
	char dateiname1[20];



	dim = 2;
	a = 0.15;
	b = 4.;
	c = -3.5;
	mesh = 0.2;


	h = vec_real_new(5);
	for(l=0;l<5;l++) {
		h->array[l] = mesh*pow(2,-l);
	}





/*Die Funktion f erstellen*/
	f = func_new(2,dim);

	for(l=0;l<=f->maxlevel;l++) {
		grad = vec_new(dim);
		grad->array[0] = 0;
		grad->array[1] = 0;
		f->hierarchie[l] = folgen_vektor_new( grad );
	}

	start = vec_new(dim);
	lang = vec_new(dim);
	start->array[0] = 0;
	start->array[1] = 0;
	lang->array[0] = 5;
	lang->array[1] = 5;
	f->hierarchie[0]->vektor[0] = folge_new(start,lang);


	for(l=1;l<=f->maxlevel;l++) {
		start = vec_new(dim);
		lang = vec_new(dim);
		start->array[0] = 0;
		start->array[1] = 0;
		lang->array[0] = 4;
		lang->array[1] = pow(2,l)*5;
		f->hierarchie[l]->vektor[0] = folge_new(start,lang);
	}

	l=0;
	for(i=2;i<5;i++) {
		for(j=0;j<5;j++) {
			r = vec_new(dim);
			r->array[0] = i;
			r->array[1] = j;
			pos = entry_d2one(r,f->hierarchie[l]->vektor[0]->lang);
			vec_del(r);
			f->hierarchie[l]->vektor[0]->glied[pos] = log((i+1)*h->array[l]+a) - log(i*h->array[l]+a);
		}
	}


	for(l=1;l<f->maxlevel;l++) {
		for(i=2;i<4;i++) {
			for(j=0;j<(5*pow(2,l));j++) {
				r = vec_new(dim);
				r->array[0] = i;
				r->array[1] = j;
				pos = entry_d2one(r,f->hierarchie[l]->vektor[0]->lang);
				vec_del(r);
				f->hierarchie[l]->vektor[0]->glied[pos] = log((i+1)*h->array[l]+a) - log(i*h->array[l]+a);
			}
		}
	}


	l=f->maxlevel;
	for(i=0;i<4;i++) {
		for(j=0;j<(5*pow(2,l));j++) {
			r = vec_new(dim);
			r->array[0] = i;
			r->array[1] = j;
			pos = entry_d2one(r,f->hierarchie[l]->vektor[0]->lang);
			vec_del(r);
			f->hierarchie[l]->vektor[0]->glied[pos] = log((i+1)*h->array[l]+a) - log(i*h->array[l]+a);
		}
	}




	/*Die Funktion g erstellen*/

	g = func_new(0,dim);

	for(l=0;l<=g->maxlevel;l++) {
		grad = vec_new(dim);
		grad->array[0] = 0;
		grad->array[1] = 0;
		g->hierarchie[l] = folgen_vektor_new( grad );
	}

	for(l=0;l<=g->maxlevel;l++) {
		start = vec_new(dim);
		lang = vec_new(dim);
		start->array[0] = 0;
		lang->array[0] = 5*pow(2,l);
		start->array[1] = 0;
		lang->array[1] = 5*pow(2,l);
		g->hierarchie[l]->vektor[0] = folge_new(start,lang);
	}

	l=g->maxlevel;
	for(i=0;i<(5*pow(2,l));i++) {
		for(j=0;j<(5*pow(2,l));j++) {
			r = vec_new(dim);
			r->array[0] = i;
			r->array[1] = j;
			pos = entry_d2one(r,g->hierarchie[l]->vektor[0]->lang);
			vec_del(r);
			temp1 = b/2.*pow(h->array[l],2);
			temp2 = 2*i + 1;
			g->hierarchie[l]->vektor[0]->glied[pos] = temp1*temp2 + c*h->array[l];
		}
	}


	/*Die Funktion w erstellen*/
	w = func_new(2,dim);

	for(l=0;l<=w->maxlevel;l++) {
		grad = vec_new(dim);
		grad->array[0] = 0;
		grad->array[1] = 0;
		w->hierarchie[l] = folgen_vektor_new( grad );
	}

	start = vec_new(dim);
	lang = vec_new(dim);
	start->array[0] = 0;
	lang->array[0] = 5;
	start->array[1] = 0;
	lang->array[1] = 5;
	w->hierarchie[0]->vektor[0] = folge_new(start,lang);

	start = vec_new(dim);
	lang = vec_new(dim);
	start->array[0] = 2;
	lang->array[0] = 6;
	start->array[1] = 4;
	lang->array[1] = 6;
	w->hierarchie[1]->vektor[0] = folge_new(start,lang);

	start = vec_new(dim);
	lang = vec_new(dim);
	start->array[0] = 8;
	lang->array[0] = 4;
	start->array[1] = 16;
	lang->array[1] = 4;
	w->hierarchie[2]->vektor[0] = folge_new(start,lang);



	/*Durchführen der Berechnungen*/
	fepc = faltung_fepc(f,g,w,mesh);
	ref = faltung_ref(f,g,w,mesh);
	norm = faltung_hilfe_norm(fepc,ref);
	printf("\n Norm %lf \n",norm);


	sprintf(dateiname1,"Beispiel_2d/tab.xyz");
	func_plot_2d( fepc, dateiname1, mesh);

	return 0;
}

void
func_plot_1d(func_p f, char dateiname1[20], fepc_real_t mesh) {
	FILE  *tab;
	int  l,i, anfang, ende, a, b;
	vec_real_p  h;
	fepc_real_t  x, y, koef;

	h = vec_real_new(f->maxlevel+1);
	for(l=0;l<=f->maxlevel;l++) {
		h->array[l] = mesh*pow(2,-l);
	}

	tab = fopen(dateiname1,"w");
	for(l=0;l<=f->maxlevel;l++) {
		anfang = f->hierarchie[l]->vektor[0]->start->array[0];
		ende = anfang + f->hierarchie[l]->vektor[0]->lang->array[0]-1;
		for(i=anfang;i<=ende;i++) {
			x = (2.*i+1.) * h->array[l]/2.;
			koef = f->hierarchie[l]->vektor[0]->glied[i-anfang];
			y = koef/sqrt(h->array[l]);
			if(l<f->maxlevel) {
				a = f->hierarchie[l+1]->vektor[0]->start->array[0];
				b = a + f->hierarchie[l+1]->vektor[0]->lang->array[0]-1;
				if( (2*i<a)||(2*i>b) ){
					fprintf( tab, "%.5lf      %.5lf \n", x,y );
				}
			}
			else {
				fprintf( tab, "%.5lf      %.5lf \n", x,y );
			}
		}
	}
	fclose(tab);
}


void
func_plot_2d(func_p f, char dateiname1[20], fepc_real_t mesh) {
	FILE  *tab;
	int  l, i, j, anfang1, ende1, a1, b1, anfang0, ende0, a0, b0, dim, pos;
	vec_real_p  h;
	fepc_real_t  x, y, z, koef;
	vec_p  r, s;

	dim = 2;
	h = vec_real_new(f->maxlevel+1);
	for(l=0;l<=f->maxlevel;l++) {
		h->array[l] = mesh*pow(2,-l);
	}

	tab = fopen(dateiname1,"w");
	for(l=0;l<=f->maxlevel;l++) {

		anfang0 = f->hierarchie[l]->vektor[0]->start->array[0];
		ende0 = anfang0 + f->hierarchie[l]->vektor[0]->lang->array[0]-1;
		anfang1 = f->hierarchie[l]->vektor[0]->start->array[1];
		ende1 = anfang1 + f->hierarchie[l]->vektor[0]->lang->array[1]-1;

		for(i=anfang0;i<=ende0;i++) {
			for(j=anfang1;j<=ende1;j++) {
				x = (2.*i+1.) * h->array[l]/2.;
				y = (2.*j+1.) * h->array[l]/2.;
				r = vec_new(dim);
				r->array[0] = i;
				r->array[1] = j;
				s = vec_op(1,r,-1,f->hierarchie[l]->vektor[0]->start);
				pos = entry_d2one(s,f->hierarchie[l]->vektor[0]->lang);
				vec_del(r);
				vec_del(s);
				koef = f->hierarchie[l]->vektor[0]->glied[pos];
				z = koef/(h->array[l]);

				if(l<f->maxlevel) {
					a0 = f->hierarchie[l+1]->vektor[0]->start->array[0];
					b0 = a0 + f->hierarchie[l+1]->vektor[0]->lang->array[0]-1;
					a1 = f->hierarchie[l+1]->vektor[0]->start->array[1];
					b1 = a1 + f->hierarchie[l+1]->vektor[0]->lang->array[1]-1;
					if( (2*i<a0)||(2*i>b0)||(2*j<a1)||(2*j>b1) ){
						fprintf( tab, "%.5lf  %.5lf  %.5lf \n", x,y,z );
					}
				}
				else {
					fprintf( tab, "%.5lf  %.5lf  %.5lf \n", x,y,z );
				}
			}
		}
	}
	fclose(tab);
}
