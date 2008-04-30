
/******************************************************
 *
 * parentposition - implementation file
 *
 * copyleft (c) IOhannes m zm�lnig
 *
 *   2007:forum::f�r::uml�ute:2007
 *
 *   institute of electronic music and acoustics (iem)
 *
 ******************************************************
 *
 * license: GNU General Public License v.2
 *
 ******************************************************/


/* 
 * this object provides a way to get the position of the containing abstraction
 * within the parent-patch
 */

#include "m_pd.h"

#include "g_canvas.h"
#include "m_imp.h"

/* ------------------------- parentposition ---------------------------- */

static t_class *parentposition_class;

typedef struct _parentposition
{
  t_object  x_obj;
  t_canvas  *x_canvas;

  t_outlet*xoutlet, *youtlet;
} t_parentposition;

static void parentposition_bang(t_parentposition *x)
{
  t_canvas*c=x->x_canvas;
  t_canvas*c0=0;

  int x1=0, y1=0, width=0, height=0;
  t_atom alist[2];

  if(!c) return;


  x1=c->gl_obj.te_xpix;
  y1=c->gl_obj.te_ypix;


  c0=c->gl_owner;
  if(c0!=0) {
    width= (int)(c0->gl_screenx2 - c0->gl_screenx1);
    height=(int)(c0->gl_screeny2 - c0->gl_screeny1);
  }

  SETFLOAT(alist, (t_float)width);
  SETFLOAT(alist+1, (t_float)height);
  outlet_list(x->youtlet, 0, 2, alist);

  //  outlet_float(x->youtlet, y1);
  SETFLOAT(alist, (t_float)x1);
  SETFLOAT(alist+1, (t_float)y1);
  outlet_list(x->xoutlet, 0, 2, alist);

}

static void parentposition_free(t_parentposition *x)
{
  outlet_free(x->xoutlet);
  outlet_free(x->youtlet);
}

static void *parentposition_new(void)
{
  t_parentposition *x = (t_parentposition *)pd_new(parentposition_class);
  t_glist *glist=(t_glist *)canvas_getcurrent();
  t_canvas *canvas=(t_canvas*)glist_getcanvas(glist);

  x->x_canvas = canvas;

  x->xoutlet=outlet_new(&x->x_obj, &s_list);
  x->youtlet=outlet_new(&x->x_obj, &s_list);

  return (x);
}

void parentposition_setup(void)
{
  parentposition_class = class_new(gensym("parentposition"), (t_newmethod)parentposition_new,
    (t_method)parentposition_free, sizeof(t_parentposition), 0, 0);
  class_addbang(parentposition_class, (t_method)parentposition_bang);
}
