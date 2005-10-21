// *********************(c)*2004*********************>
// -holzilib--holzilib--holzilib--holzilib--holzilib->
// ++++PD-External++by+Georg+Holzmann++grh@gmx.at++++>
//
// PDContainer: 
// this is a port of the containers from the C++ STL
// (Standard Template Library)
// for usage see the documentation and PD help files
// for license see readme.txt
//
// h_prioqueue.cpp


#include "include/HPrioQueue.h"


static t_class *h_priority_queue_class;
static t_class *proxy_class;

typedef struct _h_priority_queue 
{
  t_object  x_obj;
  t_canvas *x_canvas;
  t_outlet *out0, *out1, *out2;
  HPrioQueue *hpriority_queue;
  Element value;
  bool event_set;
} t_h_priority_queue;

typedef struct proxy
{
  t_object obj;
  t_int index;  // number of proxy inlet(s)
  t_h_priority_queue *x;	// we'll put the other struct in here
} t_proxy;

static void h_priority_queue_push(t_h_priority_queue *x, t_symbol *s, int argc, t_atom *argv)
{
  if(!x->event_set)
    {
      post("h_priority_queue, insert: you must first set a value at right inlet!");
      return;
    }

  float prio;
  if(argc && (argv[0].a_type == A_FLOAT))
    prio = static_cast<float>(argv[0].a_w.w_float);
  else
    {
      post("h_priority_queue, push: invalid priority!");
      return;
    }

  x->hpriority_queue->push(prio,x->value);
  x->event_set = false;
}

static void h_priority_queue_value(t_proxy *p, t_symbol *s, int argc, t_atom *argv)
{
  t_h_priority_queue *x = (t_h_priority_queue *)(p->x);

  // symbol without selector "symbol":
  if(argc == 0)
  {
    t_atom tmp;
    SETSYMBOL(&tmp, s);
    x->value.setAtoms(1, &tmp);
    x->event_set = true;
    return;
  }

  // input is a list without selector "list":
  if ( argc && (strcmp(s->s_name,"list")!=0) 
       && (strcmp(s->s_name,"float")!=0) 
       && (strcmp(s->s_name,"symbol")!=0)
       && (strcmp(s->s_name,"pointer")!=0) )
  {
    t_atom *atoms = (t_atom*)getbytes( (argc+1)*sizeof(t_atom) );

      // add the selector symbol to the list:
    SETSYMBOL(atoms, s);

    for(int i=0; i<argc; i++)
    {
      if(argv[i].a_type == A_FLOAT)
	SETFLOAT(&atoms[i+1],argv[i].a_w.w_float);
      if(argv[i].a_type == A_SYMBOL)
	SETSYMBOL(&atoms[i+1],argv[i].a_w.w_symbol);
      if(argv[i].a_type == A_POINTER)
	SETPOINTER(&atoms[i+1],argv[i].a_w.w_gpointer);
    }

    x->value.setAtoms(argc+1, atoms);

    x->event_set = true;
    freebytes(atoms, (argc+1)*sizeof(t_atom));
    return;
  }

  // "normal" input (list, float, pointer or symbol):
  if (argc)
  {
    x->value.setAtoms(argc, argv);
    x->event_set = true;
    return;
  }
}

static void h_priority_queue_top(t_h_priority_queue *x)
{
  if(x->hpriority_queue->getSize()==0)
    {
      // if there was no Element found, put out a bang at the right outlet
      outlet_bang(x->out2);
      return;
    }

  Element output = x->hpriority_queue->top();
 
  if(output.getLength() == 1) // symbol or float
    {
      if (output.getAtom()[0].a_type == A_FLOAT)
	outlet_float(x->out0, output.getAtom()[0].a_w.w_float);
      if (output.getAtom()[0].a_type == A_SYMBOL)
	outlet_symbol(x->out0, output.getAtom()[0].a_w.w_symbol);
      if (output.getAtom()[0].a_type == A_POINTER)
	outlet_pointer(x->out0, output.getAtom()[0].a_w.w_gpointer);
      return;
    }
  if(output.getLength() > 1) // list
    {
      outlet_list(x->out0,&s_list,output.getLength(),output.getAtom());
      return;
    }

  // if there was no Element found, put out a bang at the right outlet
  outlet_bang(x->out2);
}

static void h_priority_queue_pop(t_h_priority_queue *x)
{
  if(x->hpriority_queue->getSize()<=0)
    {
      post("h_priority_queue, pop: size is already 0 !");
      return;
    }

  x->hpriority_queue->pop();
}

static void h_priority_queue_getsize(t_h_priority_queue *x)
{
  outlet_float(x->out1,x->hpriority_queue->getSize());
}

static void h_priority_queue_help(t_h_priority_queue *x)
{
  x->hpriority_queue->help();
}

static void h_priority_queue_set_namespace(t_h_priority_queue *x, t_symbol *s)
{
  x->hpriority_queue->setNamespace(s->s_name);
}

static void h_priority_queue_get_namespace(t_h_priority_queue *x)
{
  post("h_priority_queue current namespace: %s",x->hpriority_queue->getNamespace().c_str());
}

static void h_priority_queue_clear(t_h_priority_queue *x)
{
  x->hpriority_queue->clearNamespace();
}

static void h_priority_queue_clear_all(t_h_priority_queue *x)
{
  x->hpriority_queue->clearAll();
}

static void *h_priority_queue_new(t_symbol *s, int argc, t_atom *argv) 
{
  t_h_priority_queue *x = (t_h_priority_queue *)pd_new(h_priority_queue_class);
  t_proxy *inlet = (t_proxy *)pd_new(proxy_class); // for the proxy inlet

  inlet->x = x;	// make x visible to the proxy inlets
  
  switch(argc)
    {
    default:
      post("h_priority_queue warning: only one argument for namespace is possible!");
   case 1:
      x->hpriority_queue = new HPrioQueue(atom_getsymbol(argv)->s_name);
      break;
    case 0:
      x->hpriority_queue = new HPrioQueue();
      break;
    }

  // we are going to create a proxy inlet no. 0
  // it belongs to the object x but the destination is t_proxy
  inlet->index = 0;
  inlet_new(&x->x_obj, &inlet->obj.ob_pd, 0,0);

  x->out0 = outlet_new(&x->x_obj, 0);
  x->out1 = outlet_new(&x->x_obj, &s_float);
  x->out2 = outlet_new(&x->x_obj, &s_bang);
  x->x_canvas = canvas_getcurrent();
    
  return (void *)x;
}

static void *h_priority_queue_free(t_h_priority_queue *x)
{
  delete x->hpriority_queue;
  return (void *)x;
}

void h_priority_queue_setup(void) 
{
  // the object class
  h_priority_queue_class = class_new(gensym("h_priority_queue"), (t_newmethod)h_priority_queue_new,
				(t_method)h_priority_queue_free, sizeof(t_h_priority_queue), 
				CLASS_DEFAULT, A_GIMME, 0);

  // a class for the proxy-inlet
  proxy_class = class_new(gensym("h_priority_queue_proxy"), NULL, NULL, sizeof(t_proxy),
			  CLASS_PD|CLASS_NOINLET, A_NULL);

  class_addmethod(h_priority_queue_class, (t_method)h_priority_queue_push, 
		  gensym("push"), A_GIMME, 0);
  class_addanything(proxy_class, (t_method)h_priority_queue_value); // the right inlet
  class_addmethod(h_priority_queue_class, (t_method)h_priority_queue_pop, 
		  gensym("pop"), A_DEFFLOAT, 0);
  class_addmethod(h_priority_queue_class, (t_method)h_priority_queue_top, 
		  gensym("top"), A_DEFFLOAT, 0);
  class_addmethod(h_priority_queue_class, (t_method)h_priority_queue_getsize, 
		  gensym("getsize"), A_DEFFLOAT , 0);
  class_addmethod(h_priority_queue_class, (t_method)h_priority_queue_set_namespace, 
		  gensym("namespace"), A_DEFSYMBOL , 0);
  class_addmethod(h_priority_queue_class, (t_method)h_priority_queue_get_namespace, 
		  gensym("getnamespace"), A_DEFFLOAT, 0);
  class_addmethod(h_priority_queue_class, (t_method)h_priority_queue_clear,  
		  gensym("clear"), A_DEFFLOAT, 0);
  class_addmethod(h_priority_queue_class, (t_method)h_priority_queue_clear_all,  
		  gensym("clearall"), A_DEFFLOAT, 0);

  // without an argument the following two methods wont work ??? why?? because of c++?
  class_addmethod(h_priority_queue_class, (t_method)h_priority_queue_help, gensym("help"),A_DEFFLOAT, 0);
}
