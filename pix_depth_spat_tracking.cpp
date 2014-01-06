////////////////////////////////////////////////////////
//
// GEM - Graphics Environment for Multimedia
//
// zmoelnig@iem.kug.ac.at
//
// Implementation file
//
//    Copyright (c) 1997-1998 Mark Danks.
//    Copyright (c) Günther Geiger.
//    Copyright (c) 2001-2011 IOhannes m zmölnig. forum::für::umläute. IEM. zmoelnig@iem.at
//    Copyright (c) 2002 James Tittle & Chris Clepper
//    For information on usage and redistribution, and for a DISCLAIMER OF ALL
//    WARRANTIES, see the file, "GEM.LICENSE.TERMS" in this distribution.
//
/////////////////////////////////////////////////////////

#include "pix_depth_spat_tracking.h"

CPPEXTERN_NEW_WITH_ONE_ARG(pix_depth_spat_tracking,t_floatarg, A_DEFFLOAT);

#define XtoZ 1.111466646194458
#define YtoZ 0.833599984645844

/////////////////////////////////////////////////////////
//
// pix_depth_spat_tracking
//
/////////////////////////////////////////////////////////
// Constructor
//
/////////////////////////////////////////////////////////
pix_depth_spat_tracking :: pix_depth_spat_tracking(t_floatarg f) : v1(3), v2(3), v3(3), v4(3),
                            T(3,3), debug_print(false), m_visualize(false)
{
  xsize = 640;
  ysize = 480;
  
  m_thresh = f;
  
  resolutionMess(4.0, 4.0, 4.0);
  
  // allocate space
  m_counter = (int*) malloc(m_res_x*m_res_y*m_res_z*sizeof(int));
  
  m_csize = 3;
  
  if (xsize < 0) xsize = 0;
  if (ysize < 0) ysize = 0;
  
  m_xsize = xsize;
  m_ysize = ysize;
  
  oldimagex = xsize;
  oldimagey = ysize;
  
  m_bufsize = m_xsize * m_ysize * m_csize;
  
  m_buffer = new t_atom[m_bufsize];

  
  // inlet threshold
  inlet_new(this->x_obj, &this->x_obj->ob_pd, gensym("float"), gensym("thresh"));
  //outlet
  m_outlet = outlet_new(this->x_obj, 0);
  
  m_outlet_temp = outlet_new(this->x_obj, 0);
}

/////////////////////////////////////////////////////////
// Destructor
//
/////////////////////////////////////////////////////////
pix_depth_spat_tracking :: ~pix_depth_spat_tracking()
{
	outlet_free(m_outlet);
  free(m_buffer);
  free(m_counter);
}


/////////////////////////////////////////////////////////
// processImage
//
/////////////////////////////////////////////////////////
void pix_depth_spat_tracking :: processImage(imageStruct &image)
{
    int x = m_xsize, y = m_ysize, c = m_csize;
    
    if (image.xsize != oldimagex) {
        oldimagex = image.xsize;
        m_xsize = ((!xsize) || (xsize > oldimagex))?oldimagex:xsize;
    }
    if (image.ysize != oldimagey) {
        oldimagey = image.ysize;
        m_ysize = ((!ysize) || (ysize > oldimagey))?oldimagey:ysize;
    }
    
    if (image.csize != m_csize) m_csize = image.csize;
    
    if ( (m_xsize != x) || (m_ysize != y) || (m_csize != c) ) {
        // resize the image buffer
        if(m_buffer)delete [] m_buffer;
        m_bufsize = m_xsize * m_ysize * m_csize;
        m_buffer = new t_atom[m_bufsize];
        
        m_xstep = m_csize * (static_cast<float>(image.xsize)/static_cast<float>(m_xsize));
        m_ystep = m_csize * (static_cast<float>(image.ysize)/static_cast<float>(m_ysize)) * image.xsize;
    }
  
  if (m_visualize) {
    int value = 0;
    
    int n = 0, m = 0;
    int i = 0;
    
    unsigned char *data, *line;
    
    Vector3d pos(3); // position vector
    
    Vector3d pos_transf(3);
    
    data = line = m_data;
    
    switch(m_csize){
      case 4:
        while (n < m_ysize) {
          while (m < m_xsize) {
            if (data[chAlpha]) // filter out zero alpha values
            {
              
              value = ((int)data[chRed] << 8) + (int)data[chGreen];
              bool set_zero = true;
              
              if (value) { // filter out zeros
                float real_x = ((m / 640.0) - 0.5) * (float)value * XtoZ; // real_x
                pos[0] = real_x - (float)v1[0]; // subtract origin of new vectorbase
                
                float real_y = (0.5 - (n / 480.0)) * (float)value * YtoZ; // real_y
                pos[1] = real_y - (float)v1[1]; // subtract origin of new vectorbase
                
                pos[2] = (float)value - (float)v1[2]; // subtract origin of new vectorbase
                
                pos_transf=T*pos; // transform with transformation matrix !!
                
                int pos_t_x = floor(pos_transf[0]);
                int pos_t_y = floor(pos_transf[1]);
                int pos_t_z = floor(pos_transf[2]);
                
                // look if inside parallelepiped
                if ((pos_t_x >= 0) && (pos_t_x < m_res_x))
                {
                  if ((pos_t_y >= 0) && (pos_t_y < m_res_y))
                  {
                    if ((pos_t_z >= 0) && (pos_t_z < m_res_z))
                    {
                      set_zero = false;
                    }
                  }
                }
              }
              
              if (set_zero) {
                data[chAlpha]=0;
              }

            }
            m++;
            
            data = line + static_cast<int>(m_xstep * static_cast<float>(m));
          }
          m = 0;
          n++;
          line = m_data + static_cast<int>(m_ystep*n);
          data = line;
        }        
        
        break;
        
      default:
        break;
    }
  }
    m_data = image.data;
}

/////////////////////////////////////////////////////////
// trigger
//
/////////////////////////////////////////////////////////
void pix_depth_spat_tracking :: trigger()
{
  if (!m_data) return;
  
  int value = 0;
  
  
  // set counter zero
  memset(m_counter,0,m_res_x*m_res_y*m_res_z*sizeof(int));
  
  int n = 0, m = 0;
  int i = 0;
  
  unsigned char *data, *line;
  
  Vector3d pos(3); // position vector
  
  Vector3d pos_transf(3);
  
  data = line = m_data;
  
  // outputMatrix(T);
  
  switch(m_csize){
    case 4:
      while (n < m_ysize) {
        while (m < m_xsize) {
          if ((int)data[chAlpha]) // filter out zero alpha values
          {
            value = ((int)data[chRed] << 8) + (int)data[chGreen];
            
            if (value) { // filter out zeros
              
              // calculate real x/y coordinates for trimming
              // x component
              //float FovH = 1.0144686707507438;
              //float XtoZ = tan(FovH / 2.0) * 2.0;
              float real_x = ((m / 640.0) - 0.5) * (float)value * XtoZ; // real_x
              
              pos[0] = real_x - (float)v1[0]; // subtract origin of new vectorbase
              
              // y component
              //float FovV=0.78980943449644714;
              //float YtoZ = tan(FovV / 2.0) * 2.0;
              float real_y = (0.5 - (n / 480.0)) * (float)value * YtoZ; // real_y
              
              pos[1] = real_y - (float)v1[1]; // subtract origin of new vectorbase
              
              pos[2] = (float)value - (float)v1[2]; // subtract origin of new vectorbase
              
              pos_transf=T*pos; // transform with transformation matrix !!
              
              // do the matrix vector mult by hand
              /*
              pos_transf[0]=T[0][0]*pos[0]+T[0][1]*pos[0]+T[0][2]*pos[0];
              pos_transf[1]=T[1][0]*pos[1]+T[1][1]*pos[1]+T[1][2]*pos[1];
              pos_transf[2]=T[2][0]*pos[2]+T[2][1]*pos[2]+T[2][2]*pos[2];
              */
              
              int pos_t_x = floor(pos_transf[0]);
              int pos_t_y = floor(pos_transf[1]);
              int pos_t_z = floor(pos_transf[2]);
              
              if (debug_print)
              {
                t_atom atom[3];
                SETFLOAT(&atom[0], pos_t_x);
                SETFLOAT(&atom[1], pos_t_y);
                SETFLOAT(&atom[2], pos_t_z);
                outlet_list(m_outlet_temp, gensym("list"), 3, atom);
              }
              
              // look if inside parallelepiped
              if ((pos_t_x >= 0) && (pos_t_x < m_res_x))
              {
                if ((pos_t_y >= 0) && (pos_t_y < m_res_y))
                {
                  if ((pos_t_z >= 0) && (pos_t_z < m_res_z))
                  {
                    // index: ((rowindex*col_size+colindex) * depth_size + depthindex)
                    //m_counter[pos_t_x*m_res_x+pos_t_y*m_res_y+pos_t_z]++; // increment counter
                    m_counter[(pos_t_x*m_res_y+pos_t_y)*m_res_z+pos_t_z]++; // increment counter
                  }
                }
              }
            }
          }
          m++;
          
          data = line + static_cast<int>(m_xstep * static_cast<float>(m));
        }
        m = 0;
        n++;
        line = m_data + static_cast<int>(m_ystep*n);
        data = line;
      }
      
      // output counter values
      for (int i = 0; i < m_res_x; ++i) {
        for (int j = 0; j < m_res_y; ++j) {
          for (int k = 0; k < m_res_z; ++k) {
            if (m_counter[(i*m_res_y+j)*m_res_z+k] >= m_thresh)
            {
              t_atom atom[4];
              SETFLOAT(&atom[0], i);
              SETFLOAT(&atom[1], j);
              SETFLOAT(&atom[2], k);
              SETFLOAT(&atom[3], m_counter[(i*m_res_y+j)*m_res_z+k]);
              outlet_list(m_outlet, gensym("list"), 4, atom);
              
            }
          }
        }
      }
      
      
      break;
      
      default:
      break;
  }
  
  debug_print = false;
   
}

// compute transformation matrix
// x1 y1 z1 x2 y2 z2 x3 y3 z3 x4 y4 z4
void pix_depth_spat_tracking :: verticesMess(int argc, t_atom *argv)
{
  if (argc < 12)
    return;

  // define 4 vectors
  v1[0] = atom_getfloat(&argv[0]);
  v1[1] = atom_getfloat(&argv[1]);
  v1[2] = atom_getfloat(&argv[2]);
  
  v2[0] = atom_getfloat(&argv[3]);
  v2[1] = atom_getfloat(&argv[4]);
  v2[2] = atom_getfloat(&argv[5]);
  
  v3[0] = atom_getfloat(&argv[6]);
  v3[1] = atom_getfloat(&argv[7]);
  v3[2] = atom_getfloat(&argv[8]);
  
  v4[0] = atom_getfloat(&argv[9]);
  v4[1] = atom_getfloat(&argv[10]);
  v4[2] = atom_getfloat(&argv[11]);
  
  
  computeInverse();
  
}

// compute transformation matrix
void pix_depth_spat_tracking :: computeInverse()
{
  // generate new vector space
  Matrix3d B(3,3);
  
  B.col(0)=(v2-v1)/m_res_x;
  B.col(1)=(v3-v1)/m_res_y;
  B.col(2)=(v4-v1)/m_res_z;
  
  // invert new vectorspace = transformation matrix
  
  bool invertible;
  B.computeInverseWithCheck(T,invertible);
  
  if (!invertible)
  {
    post("ERROR: Matrix inversion failed!");
  }
  
}

void pix_depth_spat_tracking :: outputMatrix(Matrix3d M)
{
  int rowno=M.rows();
  int colno=M.cols();
  
  // output matrix dimensions
  t_atom atom2[2];
  SETFLOAT(&atom2[0], rowno);
  SETFLOAT(&atom2[1], colno);
  outlet_list(m_outlet_temp, gensym("matrix"), 2, atom2);
  
  
  t_atom atom[colno];
  
  for (int i=0; i<rowno; i++) {
    int k=0;
    for (int j=0; j<colno; j++) {
      SETFLOAT(&atom[k], M(i,j));
      k++;
    }
    outlet_list(m_outlet_temp, gensym("matrix"), colno, atom);
  }
}

void pix_depth_spat_tracking :: threshMess(t_float thresh)
{
  m_thresh = thresh;
}

void pix_depth_spat_tracking :: resolutionMess(t_float res_x, t_float res_y, t_float res_z)
{
 // reallocate space
  m_counter = (int*) realloc(m_counter, (int)res_x*(int)res_y*(int)res_z*sizeof(int));
  
  m_res_x = (int)res_x;
  m_res_y = (int)res_y;
  m_res_z = (int)res_z;
  
  computeInverse();
}

void pix_depth_spat_tracking :: inverseMess()
{
  outputMatrix(T);
}

void pix_depth_spat_tracking :: visualizeMess(t_float value)
{
  if ((float)value < 0.5) {
    m_visualize = false;
  } else {
    m_visualize = true;
  }
}

/////////////////////////////////////////////////////////
// static member function
//
/////////////////////////////////////////////////////////
void pix_depth_spat_tracking :: obj_setupCallback(t_class *classPtr)
{
    class_addbang(classPtr, reinterpret_cast<t_method>(&pix_depth_spat_tracking::triggerMessCallback));
    class_addmethod(classPtr, reinterpret_cast<t_method>(&pix_depth_spat_tracking::threshMessCallback),
                    gensym("thresh"), A_FLOAT, A_NULL);
    class_addmethod(classPtr, reinterpret_cast<t_method>(&pix_depth_spat_tracking::resolutionMessCallback),
                    gensym("resolution"), A_FLOAT, A_FLOAT, A_FLOAT, A_NULL);
    class_addmethod(classPtr, reinterpret_cast<t_method>(&pix_depth_spat_tracking::verticesMessCallback),
                    gensym("vertices"), A_GIMME, A_NULL);
  class_addmethod(classPtr, reinterpret_cast<t_method>(&pix_depth_spat_tracking::printMessCallback),
                  gensym("print"), A_GIMME, A_NULL);
  class_addmethod(classPtr, reinterpret_cast<t_method>(&pix_depth_spat_tracking::inverseMessCallback),
                  gensym("inverse"), A_GIMME, A_NULL);
  class_addmethod(classPtr, reinterpret_cast<t_method>(&pix_depth_spat_tracking::visualizeMessCallback),
                  gensym("visualize"), A_FLOAT, A_NULL);
}


void pix_depth_spat_tracking :: triggerMessCallback(void *data)
{
  GetMyClass(data)->trigger();
}

void pix_depth_spat_tracking :: threshMessCallback(void *data, t_floatarg thresh)
{
    GetMyClass(data)->threshMess((float)thresh);
}
void pix_depth_spat_tracking :: resolutionMessCallback(void *data, t_floatarg res_x, t_floatarg res_y, t_floatarg res_z)
{
    GetMyClass(data)->resolutionMess((float)res_x, (float)res_y, (float)res_z);
}
void pix_depth_spat_tracking :: verticesMessCallback(void *data, t_symbol *, int argc, t_atom *argv)
{
    GetMyClass(data)->verticesMess(argc, argv);
}

void pix_depth_spat_tracking :: printMessCallback(void *data)
{
  GetMyClass(data)->debug_print=true;
}

void pix_depth_spat_tracking :: inverseMessCallback(void *data)
{
  GetMyClass(data)->inverseMess();
}

void pix_depth_spat_tracking :: visualizeMessCallback(void *data, t_floatarg value)
{
  GetMyClass(data)->visualizeMess((float)value);
}