/* For information on usage and redistribution, and for a DISCLAIMER OF ALL
* WARRANTIES, see the file, "LICENSE.txt," in this distribution.

iem_tab written by Thomas Musil, Copyright (c) IEM KUG Graz Austria 2000 - 2005 */

#ifdef NT
#pragma warning( disable : 4244 )
#pragma warning( disable : 4305 )
#endif


#include "m_pd.h"
#include "iemlib.h"

static t_class *iem_tab_class;

int iem_tab_check_arrays(t_symbol *obj_name, t_symbol *array_name, t_float **beg_mem, int *array_size, int max_index)
{
	int ok=1;
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(array_name, garray_class)))
	{
		error("%s: no such array", array_name->s_name);
		ok = 0;
	}
	else if(!garray_getfloatarray(a, array_size, beg_mem))
	{
		error("%s: bad template for %s", array_name->s_name, obj_name->s_name);
		ok = 0;
	}
	else if(*array_size < max_index)
	{
		error("%s: bad array-size: %d", array_name->s_name, *array_size);
		ok = 0;
	}
	return(ok);
}

static void *iem_tab_new(void)
{
	t_object *x = (t_object *)pd_new(iem_tab_class);
    
	return (x);
}

void tab_copy_setup(void);
void tab_reverse_setup(void);
void tab_min_max_setup(void);
void tab_min_index_setup(void);
void tab_max_index_setup(void);
void tab_find_peaks_setup(void);
void tab_abs_setup(void);
void tab_sqrt_setup(void);
void tab_sum_setup(void);
void tab_add_setup(void);
void tab_sub_setup(void);
void tab_mul_setup(void);
void tab_div_setup(void);
void tab_complex_mul_setup(void);
void tab_complex_inv_setup(void);
void tab_mul_scalar_setup(void);
void tab_add_scalar_setup(void);
void tab_const_setup(void);
void tab_fft_setup(void);
void tab_ifft_setup(void);
void tab_rfft_setup(void);
void tab_rifft_setup(void);
void tab_cross_corr_setup(void);
void tab_conv_setup(void);
void tab_gt_scalar_setup(void);
void tab_ge_scalar_setup(void);
void tab_lt_scalar_setup(void);
void tab_le_scalar_setup(void);
void tab_ne_scalar_setup(void);
void tab_eq_scalar_setup(void);
void tab_gt_setup(void);
void tab_ge_setup(void);
void tab_lt_setup(void);
void tab_le_setup(void);
void tab_ne_setup(void);
void tab_eq_setup(void);
void tab_counter_setup(void);
//void tab_mls_setup(void);

/* ------------------------ setup routine ------------------------- */

void iem_tab_setup(void)
{
	iem_tab_class = class_new(gensym("iem_tab"), iem_tab_new, 0,
    	sizeof(t_object), CLASS_NOINLET, 0);

		tab_copy_setup();
		tab_reverse_setup();
		tab_min_max_setup();
		tab_min_index_setup();
		tab_max_index_setup();
    tab_find_peaks_setup();
		tab_abs_setup();
		tab_sqrt_setup();
		tab_sum_setup();
		tab_add_setup();
		tab_sub_setup();
		tab_mul_setup();
		tab_div_setup();
		tab_complex_mul_setup();
		tab_complex_inv_setup();
		tab_mul_scalar_setup();
		tab_add_scalar_setup();
		tab_const_setup();
		tab_fft_setup();
		tab_ifft_setup();
		tab_rfft_setup();
		tab_rifft_setup();
    tab_cross_corr_setup();
    tab_conv_setup();
    tab_gt_scalar_setup();
    tab_ge_scalar_setup();
    tab_lt_scalar_setup();
    tab_le_scalar_setup();
    tab_ne_scalar_setup();
    tab_eq_scalar_setup();
    tab_gt_setup();
    tab_ge_setup();
    tab_lt_setup();
    tab_le_setup();
    tab_ne_setup();
    tab_eq_setup();
    tab_counter_setup();
//		tab_mls_setup();

    post("iem_tab (R-1.16) library loaded!   (c) Thomas Musil 05.2005");
	post("   musil%ciem.at iem KUG Graz Austria", '@');
}
